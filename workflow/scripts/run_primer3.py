#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from functools import partial
from multiprocessing import Pool

import primer3
from Bio import SeqIO


def parse_args(argv):
    parser = ArgumentParser(description="FASTA file -> primers", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("seqs", help="the FASTA file", type=FileType())
    parser.add_argument("conf", help="the JSON file of Primer3 global input tags (https://primer3.org/manual.html#globalTags)")
    parser.add_argument("-proc", help="the number of processes", default=1, type=int)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with open(args.conf) as file:
        global_args = json.load(file)

    p_designPrimers = partial(primer3.bindings.designPrimers, global_args=global_args)
    with Pool(args.proc) as pool, args.seqs as file:
        records = list(SeqIO.parse(file, "fasta"))
        seq_args = (dict(SEQUENCE_TEMPLATE=str(record.seq)) for record in records)
        results = pool.map(p_designPrimers, seq_args)

        for record, result in zip(records, results):
            template = str(record.seq)
            result.update(
                SEQUENCE_IDENTIFIER=record.id,
                SEQUENCE_TEMPLATE=template,
            )

        json.dump(results, sys.stdout, indent=4)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
