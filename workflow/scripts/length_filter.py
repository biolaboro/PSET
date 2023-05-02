#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType

from Bio import SeqIO


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", type=FileType())
    parser.add_argument("length", type=int)
    parser.add_argument("percent", type=float)
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    lower = args.length - args.length * args.percent
    upper = args.length + args.length * args.percent
    records = (record for record in SeqIO.parse(args.file, "fasta") if lower <= len(record) <= upper)
    SeqIO.write(records, args.out, "fasta")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
