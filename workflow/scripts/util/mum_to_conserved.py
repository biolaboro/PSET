#!/usr/bin/env python3

import json
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from itertools import product
from multiprocessing import Pool
from re import finditer

import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_ranges(file):
    for line in file:
        tokens = line.strip().split("\t")
        p, q = int(tokens[0]), int(tokens[1])
        yield p, q, *tokens[2:]


def parse_args(argv):
    parser = ArgumentParser(description="", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("ref")
    parser.add_argument("lib")
    parser.add_argument("tsv")
    parser.add_argument("-proc", type=int, default=1)
    parser.add_argument("-out", type=FileType("w"), default="-")
    return parser.parse_args(argv)


def func(qry, sbj):
    return int(qry.id), sbj.id, int(qry.seq in sbj or qry.reverse_complement().seq in sbj)


def main(argv):
    args = parse_args(argv[1:])
    ref = SeqIO.read(args.ref, "fasta").upper()

    with open(args.tsv) as file:
        # load accessions from the ranges file (start, finish, accession)
        accs = (ele.id for ele in SeqIO.parse(args.lib, "fasta"))
        # give each accession an index
        accs = {ele: idx for idx, ele in enumerate(accs)}
        # matrix of zerps n accessions (rows) by m reference positions (cols)
        cover = np.zeros(shape=(len(accs), len(ref) - 1))
        file.seek(0)
        for rng in parse_ranges(file):
            cover[accs[rng[2]], rng[0] - 1 : rng[1]] = 1

    # one way to get n accessions
    total = cover.shape[0]
    # get the column sums
    cover = cover.sum(axis=0)
    bits = "".join(map(str, (cover == total).astype(int)))

    seqs = []
    for idx, match in enumerate(finditer(rf"1+", bits)):
        start, finish = match.span()
        seq = ref[slice(*match.span())].seq
        seq = SeqRecord(seq, id=str(idx), name=ref.id, description=f"{start+1}-{finish}")
        seqs.append(seq.upper())
    else:
        idx = 0

    width = len(str(idx))
    for ele in seqs:
        ele.id = ele.id.zfill(width)

    # load the actual query sequences
    lib = SeqIO.parse(args.lib, "fasta")
    with Pool(processes=args.proc) as pool:
        params = product(seqs, lib)
        results = pool.starmap(func, params)

    # matrix: row = conserved sequence index, col = accession index, val = 0 or 1 (not) conserved
    hits = np.zeros(shape=(len(seqs), len(accs)))
    for key, acc, hit in results:
        hits[key][accs[acc]] = hit

    obj = dict(
        ref=ref.id,
        seqs=[str(ele.seq) for ele in seqs],
        coor=[tuple(map(int, ele.description.split("-"))) for ele in seqs],
        accs=list(accs),
        hits=hits.astype(int).tolist(),
    )
    json.dump(obj, args.out)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
