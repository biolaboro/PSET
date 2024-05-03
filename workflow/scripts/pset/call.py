#!/usr/bin/env python3

import json
import re
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import Counter, defaultdict, namedtuple
from tempfile import NamedTemporaryFile

from Bio import SearchIO
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from taxa.taxa import ancestors, session_scope

from pset.assay import AlignType, Assay, align_type, decode_btop
from pset.util import fields_8CB, iter_hsps


def parse_argv(argv):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-assay", help="the input assay JSON file", required=True)
    parser.add_argument("-mapping", help="the input sequence JSON file", required=True)
    parser.add_argument("-taxa", help="the input taxa mapping TSV file", required=True)
    parser.add_argument("-alignment", help="the input alignment TSV file", required=True)
    parser.add_argument("-dburl", help="the path to the taxonomy SQLite database", required=True)
    parser.add_argument("-sim", help="the component similarity threshold", default=0.9, type=float)
    parser.add_argument("-dFR", help="the min/max distance range (inclusive) for F/R primers", default=(1, 1000), type=int, nargs=2)
    parser.add_argument("-dF3F2", help="the min/max distance range (inclusive) for F3/F2 primers", default=(20, 80), type=int, nargs=2)
    parser.add_argument("-dF2F1c", help="the min/max distance range (inclusive) for F2/F1c primers", default=(20, 80), type=int, nargs=2)
    parser.add_argument("-dF1cB1c", help="the min/max distance range (inclusive) for F1c/B1c primers", default=(1, 100), type=int, nargs=2)
    parser.add_argument("-xtaxa", help="the ancestral taxa to exclude", nargs="+", default=(81077,), type=int)
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def main(argv):
    args = parse_argv(argv[1:])

    xtaxa = set(map(int, args.xtaxa))

    # load assay
    assay = Assay.from_json(args.assay)
    # map duplicates
    with open(args.mapping) as file:
        mapped = {ele[0]: ele for ele in json.load(file).values()}
    # map accession to taxid
    with open(args.taxa) as file:
        taxa = {key: int(val) for key, val in (line.strip().split("\t") for line in file)}
    # map taxid to lineage
    engine = create_engine(args.dburl)
    with session_scope(sessionmaker(bind=engine)) as session:
        lineages = {int(val): set(ele.tax_id for ele in ancestors(session, val) if ele) for val in set(taxa.values())}
        nn_taxa = {ancestors(session, ele)[-1].parent_tax_id for ele in assay.targets}

    # process glsearch hits
    result = namedtuple("Result", ("qaln", "saln", "astr", "atypes", "qcov", "psim", "qsim"))
    sbj_results = defaultdict(list)
    hsp_results = {}
    # process BLAST+ results
    with NamedTemporaryFile() as temp, open(args.alignment, "rb") as file:
        temp.writelines((line for line in file if not line.startswith(b"#")))
        temp.flush()
        for hsp in iter_hsps(SearchIO.parse(temp.name, "blast-tab", comments=False, fields=fields_8CB)):
            # calculate query similarity percentage and mutations
            qry = next(assay.records(hsp.query_id)).seq[hsp.query_start: hsp.query_end]
            qry = str(qry if hsp.query_strand >= 0 else qry.reverse_complement())
            qaln = "".join(decode_btop(qry, hsp.btop, sbj=False))
            saln = "".join(decode_btop(qry, hsp.btop))
            qcov = (hsp.query_end - hsp.query_start) / len(qry)
            alignment = [align_type(*ele) for ele in zip(qaln, saln)]
            astr = "".join(str(ele.value) for ele in alignment)
            atypes = Counter(ele.name for ele in alignment)
            psim = sum(not (ele.islower() or ele == "-") for ele in saln) / hsp.aln_span
            qsim = qcov * psim
            sbj_id = re.sub(r":[0-9]+$", "", hsp.hit_id)  # acc.ver...:nth-alignment
            # aggregate HSPs by the subject accession)
            sbj_results[sbj_id].append(hsp)
            hsp_results[hsp.hit_id] = result(qaln, saln, astr, atypes, qcov, psim, qsim)

    # process hits
    near = set()
    obj = dict(
        assay=dict(zip(assay.key(), assay.val())),
        atype=[ele.name for ele in AlignType],
        hits=[],
        near=[],
    )
    dkwargs = dict(dFR=args.dFR, dF3F2=args.dF3F2, dF2F1c=args.dF2F1c, dF1cB1c=args.dF1cB1c)
    for key, val in sbj_results.items():  # iter subject id -> query result mapping
        # collect HSPs of queries
        hsps = sorted((ele for ele in val if ele.query_id in assay.components), key=lambda x: hsp_results[x.hit_id].qsim, reverse=True)
        hit = next(assay.hits(hsps, **dkwargs), ())
        if hit:
            # check if all components exceed threshold
            bind = {hsp.query_id: args.sim <= hsp_results[hsp.hit_id].qsim for hsp in hit}
            is_bind = len(bind) and all(bind.values())
            calls = defaultdict(list)
            add_n = "N" * any("n" in hsp_results[hsp.hit_id].saln for hsp in assay.components if hsp in hit)
            # make call for each duplicate sequence (with potentially different taxonomy)
            for acc in mapped[key]:
                tax = taxa[acc]
                is_target = bool(assay.targets & lineages[tax])
                is_exclude = bool(xtaxa & lineages[tax])
                call = "XX" if is_exclude else ((("TP" if is_target else "FP") if is_bind else ("FN" if is_target else "TN")) + add_n)
                calls[call].append(acc)
                if call[:2] in ("FP", "TN") and nn_taxa & lineages[tax]:
                    near.add(acc)
            # build evaluation object for all queries
            evals = {
                hsp.query_id: dict(
                    **dict(zip(result._fields, hsp_results[hsp.hit_id])),
                    qstr="-+"[hsp.query_strand >= 0],
                    bind=bind[hsp.query_id],
                )
                for hsp in hit
            }
            obj["hits"].append(dict(evals=evals, calls=calls))
    obj["near"] = list(near)

    # output
    with args.out as file:
        json.dump(obj, file, indent=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
