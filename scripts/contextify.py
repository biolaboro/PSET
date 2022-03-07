#!/usr/bin/env python

import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from subprocess import PIPE, Popen

from Bio import SearchIO, SeqIO
from pset.assay import dna_val, parse_assays
from pset.util import fields_std, iter_hsps


def parse_argv(argv):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("assay", help="the assay file")
    parser.add_argument("reference", help="the reference FASTA file")
    parser.add_argument("n5", type=int, help="the amount of 5' nucleotides to prepend to the definition")
    parser.add_argument("n3", type=int, help="the amount of 3' nucleotides to append to the definition")
    parser.add_argument("-pident", default=95.0, type=float, help="the percent identity threshold [0, 100]")
    args = parser.parse_args(argv)
    return args


def main(argv):
    args = parse_argv(argv[1:])

    # load the reference
    lib = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    with open(args.assay) as file:
        print("id", "targets", "definition", "ref", "sstart", "send", "pident", sep="\t")
        # for each assay
        for assay in parse_assays(file):
            # load the assay amplicon record
            amp = next(assay.records("amplicon"))
            # search amplicon
            cmd = ("glsearch36", "-3", "-n", "-d", "1", "-m", "8C", "-", args.reference)
            with Popen(cmd, universal_newlines=True, stdin=PIPE, stdout=PIPE) as proc:
                with proc.stdin as file:
                    SeqIO.write(amp, file, "fasta")
                hsp = next(iter_hsps(SearchIO.parse(proc.stdout, "blast-tab", comments=True, fields=fields_std)))
                # guarantee full-query coverage above threshold
                assert args.pident <= hsp.ident_pct
                assert hsp.query_start == 0
                assert len(assay.amplicon()) == abs(hsp.query_end - hsp.query_start)
                span = (hsp.hit_start, hsp.hit_end)
            # add context
            ref = lib[hsp.hit_id]
            start, end = span
            definition = ref[start - args.n5:start].seq + assay.adefinition() + ref[end:end + args.n3].seq
            print(
                assay.id, ";".join(map(str, sorted(assay.targets))), definition,
                hsp.hit_id, hsp.hit_start, hsp.hit_end, hsp.ident_pct,
                sep="\t"
            )

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
