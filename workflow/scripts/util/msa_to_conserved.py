#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import Counter

import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def ref_pos(ref):
    n = 0
    for idx in range(len(ref)):
        # only increment non-gap positions
        n += ref[idx] != "-"
        yield n
    yield n + 1


def main(argv):
    args = parse_args(argv[1:])

    # load multiple alignment into a matrix/dataframe-like structure
    with args.file as file:
        msa = AlignIO.read(file, "fasta")

    # map ref aln pos to ref pos
    rpos = list(ref_pos(msa[0]))
    # calculate consensus values
    consensus = [Counter(msa[:, i]).most_common(1)[0] for i in range(msa.get_alignment_length())]
    keys = [ele for ele, _ in consensus]
    vals = [ele / len(msa) for _, ele in consensus]
    # process columns, output a 1 if is conserved and 0 otherwise
    bitstring = "".join(str(int(ele >= args.threshold)) for ele in vals)
    # process bitstring, locate positions of contiguous runs of 1s, and sort by length
    conserved = sorted(re.finditer(rf"1{{{args.min},}}", bitstring), key=lambda match: match.end() - match.start(), reverse=True)
    # extract conserved regions from multiple alignment and output FASTA
    width = len(str(msa.get_alignment_length()))
    context = tuple(map(int, args.context.split(",")))
    records = []
    pfx = args.pfx or msa[0].id
    for match in conserved:
        slc = slice(match.start() - context[0], match.end() + context[1])
        seq = "".join(keys[slc])
        val = vals[slc]
        if not args.gapless or (args.gapless and "-" not in seq):
            # inclusive range (b/c NCBI)
            record = SeqRecord(Seq(seq))
            record.id = f"{pfx}:{rpos[slc.start]:0{width}}-{rpos[slc.stop]-1:0{width}}"
            record.description = " ".join(
                (
                    f"len={len(seq)}",
                    f"min={np.min(val) * 100:06.2f}%",
                    f"avg={np.mean(val) * 100:06.2f}%",
                    f"max={np.max(val) * 100:06.2f}%",
                )
            )
            records.append(record)

    with args.out as file:
        SeqIO.write(records, file, "fasta")

    return 0


def parse_args(argv):
    parser = ArgumentParser(description="convert MSA -> conserved sequence FASTA file", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", help="the multiple sequence alignment as FASTA", type=FileType())
    parser.add_argument("-pfx", help="the id prefix for consensus records, defaults to accession of first MSA record")
    parser.add_argument("-min", help="the minimum acceptable length of conserved sequences", type=int, default=100)
    parser.add_argument("-threshold", help="the minimum conservation threshold per position", type=float, default=1)
    parser.add_argument("-context", help="the amount of 5'/3'-sequence context to include to the extracted sequences", default="0,0")
    parser.add_argument("-gapless", help="the flag to only output gapless sequences", action="store_true")
    parser.add_argument("-out", help="the output file", type=FileType("w"), default="-")
    return parser.parse_args(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
