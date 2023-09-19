#!/usr/bin/env python3

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import deque
from itertools import islice

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def sliding_window(iterable, n):
    # https://docs.python.org/3/library/itertools.html
    # sliding_window('ABCDEFG', 4) --> ABCD BCDE CDEF DEFG
    it = iter(iterable)
    window = deque(islice(it, n), maxlen=n)
    if len(window) == n:
        yield tuple(window)
    for x in it:
        window.append(x)
        yield tuple(window)


def main(argv):
    args = parse_args(argv[1:])
    record = SeqIO.read(args.file, "fasta")
    arr = [gc_fraction(ele, ambiguous=args.ambiguous) for ele in sliding_window(record.seq, args.window)]
    plt.title(f"GC%({record.id}, len={len(record)}, window={args.window}, ambiguous={args.ambiguous})")
    plt.plot(arr)
    plt.axhline(y=0.5, color="black")
    if args.show:
        plt.show()
    if args.out:
        plt.savefig(args.out)


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the input file", type=FileType(), default="-")
    parser.add_argument("-window", help="the window length", type=int, default=100)
    choices = ("remove", "ignore", "weighted")
    parser.add_argument("-ambiguous", help="the mode to deal with ambiguous letters", choices=choices, default=choices[-1])
    parser.add_argument("-show", action="store_true")
    parser.add_argument("-out", help="the figure output file (optional)")
    return parser.parse_args(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
