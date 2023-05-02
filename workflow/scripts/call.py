#!/usr/bin/env python3

import json
import re
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import Counter, defaultdict, namedtuple

from Bio import SearchIO
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from taxa.taxa import ancestors, session_scope

from pset.assay import Assay, align_type, decode_btop
from pset.util import fields_8CB, iter_hsps


def parse_argv(argv):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-assay", help="the input assay JSON file", required=True)
    parser.add_argument("-mapping", help="the input sequence JSON file", required=True)
    parser.add_argument("-taxa", help="the input taxa mapping TSV file", required=True)
    parser.add_argument("-alignment", help="the input alignment TSV file", required=True)
    parser.add_argument("-dburl", help="the path to the taxonomy SQLite database", required=True)
    parser.add_argument("-sim", help="the component similarity threshold", default=0.9, type=float)
    parser.add_argument("-dist", help="the maximum distance between outer primers", default=1000, type=int)
    parser.add_argument("-xtaxa", help="the ancestral taxa to exclude", nargs="+", default=(81077,), type=int)
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def main(argv):
    args = parse_argv(argv[1:])

    xtaxa = set(map(int, args.xtaxa))

    with open(args.assay) as file:
        db = json.load(file)["config"]["db"]
    # load assay
    assay = Assay.from_json(args.assay)
    # map duplicates
    with open(args.mapping) as file:
        mapped = {ele[0]: ele for ele in json.load(file).values()}
    # map accession to taxid
    with open(args.taxa) as file:
        taxa = dict(line.strip().split("\t") for line in file)
    # map taxid to lineage
    engine = create_engine(args.dburl)
    with session_scope(sessionmaker(bind=engine)) as session:
        lineages = {ele: set(ele[0] for ele in ancestors(session, ele)) for ele in set(taxa.values())}
        nn_taxa = {ancestors(session, ele)[-1][3] for ele in assay.targets}

    # process glsearch hits
    result = namedtuple("Result", ("qaln", "saln", "qcov", "atypes", "psim", "qsim"))
    results = defaultdict(list)
    # process BLAST+ results
    for hsp in iter_hsps(SearchIO.parse(args.alignment, "blast-tab", comments=True, fields=fields_8CB)):
        # calculate query similarity percentage and mutations
        qry = next(assay.records(hsp.query_id)).seq[hsp.query_start : hsp.query_end]
        qry = str(qry if hsp.query_strand >= 0 else qry.reverse_complement())
        qaln = "".join(decode_btop(qry, hsp.btop, sbj=False))
        saln = "".join(decode_btop(qry, hsp.btop))
        qcov = (hsp.query_end - hsp.query_start) / len(qry)
        atypes = Counter((align_type(*ele).name for ele in zip(qaln, saln)))
        psim = sum(not (ele.islower() or ele == "-") for ele in saln) / hsp.aln_span
        qsim = qcov * psim
        hit_id = re.sub(r":[0-9]+$", "", hsp.hit_id)  # acc.ver...:nth-alignment
        # tack the result onto the HSP object
        hsp.result = result(qaln, saln, qcov, atypes, psim, qsim)
        # aggregate HSPs by hit id (= the subject accession)
        results[hit_id].append(hsp)

    # process hits
    near = set()
    obj = dict(db=db, assay=dict(zip(assay.key(), assay.val())), hits=[], near=[])
    for key, val in results.items():  # iter subject id -> query result mapping
        # collect HSPs of primer queries
        hsps = sorted((ele for ele in val if ele.query_id in assay.primers), key=lambda x: x.result.qsim, reverse=True)
        hit = next(assay.hits(hsps, dist=args.dist), ())
        # check if all components exceed threshold
        bind = {hsp.query_id: args.sim <= hsp.result.qsim for hsp in hit}
        is_bind = len(bind) and all(bind.values())
        calls = defaultdict(list)
        add_n = "N" * any("n" in hsp.result.saln for hsp in assay.primers if hsp in hit)
        # make call for each sequence duplcate (with potentially different taxonomy)
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
                **dict(zip(result._fields, hsp.result)),
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
