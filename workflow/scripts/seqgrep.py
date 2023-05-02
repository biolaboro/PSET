#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from csv import writer
from itertools import product
from multiprocessing import Pool

from Bio import SeqIO


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("qry")
    parser.add_argument("sbj")
    parser.add_argument("-proc", type=int, default=1)
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def func(qry, sbj):
    return int(qry.id), sbj.id, int(qry.seq in sbj or qry.reverse_complement().seq in sbj)


def main(argv):
    args = parse_args(argv[1:])

    with Pool(processes=args.proc) as pool:
        qry = SeqIO.parse(args.qry, "fasta")
        sbj = SeqIO.parse(args.sbj, "fasta")
        params = product(qry, sbj)
        results = pool.starmap(func, params)
        writer(args.out, delimiter="\t").writerows(results)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
