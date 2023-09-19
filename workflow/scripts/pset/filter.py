#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from collections import defaultdict, namedtuple
from itertools import chain

from Bio import SearchIO

from pset.assay import Assay, decode_btop, is_similar
from pset.util import contextify, fields_8CB, iter_hsps, slice_aln


def parse_args(argv):
    parser = ArgumentParser(description="filter BLAST+", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("assay", help="the assay JSON file")
    parser.add_argument("blast", help="the BLAST+ results file")
    parser.add_argument("-tsv", help="the similarity output file", required=True)
    parser.add_argument("-txt", help="the entry batch output file", required=True)
    parser.add_argument("-sim", type=float, help="the similarity threshold", default=0.85)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with open(args.assay) as file:
        obj = json.load(file)

    ctx = obj["qry"]["ctx"]
    lcl = obj["qry"]["lcl"]
    assay = Assay.from_json(args.assay)
    rec = {key: str(next(assay.records(key, context=ctx[key])).seq) for key in lcl}

    # process BLAST+ results
    result = namedtuple("Result", ("hsp", "qcov", "psim", "qsim"))
    results = defaultdict(dict)
    for hsp in iter_hsps(SearchIO.parse(args.blast, "blast-tab", comments=True, fields=(*fields_8CB, "staxids"))):
        amb = rec[hsp.query_id]  # ambiguous query
        exp = lcl[hsp.query_id]  # expanded query
        # slice to remove context from alignment
        slc = slice(hsp.query_start, hsp.query_end)
        qaln, saln = "".join(decode_btop(amb[slc], hsp.btop)), "".join(decode_btop(exp[slc], hsp.btop))
        ctx5, ctx3 = ctx[hsp.query_id]
        aln = list(slice_aln(qaln, saln, hsp.query_start, ctx5, len(lcl[hsp.query_id]) - ctx3))
        # percentage query coverage
        qcov = sum("-" != ele[0] for ele in aln) / len(assay[hsp.query_id])
        # percentage alignment similarity
        psim = sum(is_similar(*ele) for ele in aln) / len(aln)
        qsim = qcov * psim
        # keep best hsp per subject/query
        if results[hsp.hit_id].get(hsp.query_id, (0,))[-1] < qsim:
            results[hsp.hit_id][hsp.query_id] = result(hsp, qcov, psim, qsim)

    # output similarity
    with open(args.tsv, "w") as file:
        print("hit_id", "query_id", *result._fields[1:], sep="\t", file=file)
        for key1, val1 in results.items():  # iter subject id -> query result mapping
            for key2, val2 in val1.items():  # iter query id -> result
                print(key1, key2, *val2[1:], sep="\t", file=file)

    # output entry batch for blastdbcmd
    with open(args.txt, "w") as file:
        for key, val in results.items():  # iter subject id -> query result mapping
            # check to see if all queries hit to the subject
            n = sum(args.sim <= ele.qsim for ele in val.values())
            if len(lcl) == n:
                hsps = (ele.hsp for ele in val.values())
                bat = [contextify(ele, len(lcl[ele.query_id])) for ele in hsps]
                start = min(chain.from_iterable((ele[1], ele[2]) for ele in bat))
                end = max(chain.from_iterable((ele[1], ele[2]) for ele in bat))
                print(key, f"{start}-{end}", bat[0][-1], file=file)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
