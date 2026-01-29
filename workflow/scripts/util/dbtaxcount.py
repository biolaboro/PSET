#!/usr/bin/env python3

import sys
from collections import Counter
from multiprocessing import Pool
from subprocess import PIPE, Popen


def count_taxa(volume):
    cmd = ("blastdbcmd", "-db", volume, "-entry", "all", "-outfmt", "%T")
    with Popen(cmd, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdout as file:
            return Counter(map(int, file))


def parse_volumes(lines):
    flag = 0
    for line in lines:
        if line.startswith("Volumes:"):
            flag = 1
        elif flag:
            yield line.strip()


def parse_args(argv):
    from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
    parser = ArgumentParser(
        description="Count the number of each taxon in the BLAST+ database...",
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("db", help="the BLAST+ database")
    parser.add_argument("-n", help="the number of processes", type=int, default=1)
    parser.add_argument("-f", help="the output fieldnames to use", nargs=2, default=("tax", "count"))
    parser.add_argument("-o", help="the output TSV file", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    cmd = ("blastdbcmd", "-info", "-db", args.db)
    with Popen(cmd, stdout=PIPE, universal_newlines=True) as proc:
        with proc.stdout as file:
            volumes = list(parse_volumes(file))

    with Pool(processes=args.n) as pool:
        result = sum(pool.map(count_taxa, volumes), start=Counter())
        with args.o as file:
            print(*args.f, sep="\t", file=file)
            for key, val in result.items():
                print(key, val, sep="\t", file=file)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
