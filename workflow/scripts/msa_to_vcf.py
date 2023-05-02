#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, FileType
from collections import OrderedDict

from Bio import SeqIO


def ref_pos(ref):
    n = 0
    for idx in range(len(ref)):
        # only increment non-gap positions
        n += ref[idx] != "-"
        yield n
    yield n + 1


def msa_to_vcf(aln):
    # store the referece
    ref = next(aln).upper()  # ref aln seq
    rpos = list(ref_pos(ref))  # map ref aln pos to ref pos
    rseq = re.sub(r"[^ACGT]", "N", str(ref.seq).replace("-", ""))  # unamb/ungap ref seq
    # loop over the other sequences
    dna = set("ACGTN")
    for alt in aln:
        alt = alt.upper()  # alt aln seq
        apos = list(ref_pos(alt))  # map alt aln pos to alt pos
        aseq = re.sub(r"[^ACGT]", "N", str(alt.seq).replace("-", ""))  # unamb/ungap alt seq
        # SNPs
        for idx, pair in enumerate(zip(ref, alt), start=0):
            if pair[0] != pair[1] and pair[0] in dna and pair[1] in dna:
                yield alt.id, rpos[idx], *pair
        # insertions
        gaps = re.finditer(r"-+", "".join(".-"[x == "-" and y != "-"] for x, y in zip(ref, alt)))
        for gap in gaps:
            p1, p2 = gap.span()
            yield (
                (alt.id, 1, rseq[0], aseq[apos[p1] - 1 : apos[p2]])
                if p1 == 0
                else (alt.id, rpos[p1 - 1], rseq[rpos[p1] - 1], aseq[apos[p1] - 1 - 1 : apos[p2] - 1])
            )
        # deletions
        gaps = re.finditer(r"-+", "".join(".-"[x != "-" and y == "-"] for x, y in zip(ref, alt)))
        for gap in gaps:
            p1, p2 = gap.span()
            yield (
                (alt.id, 1, rseq[rpos[p1] - 1 : rpos[p2]], aseq[0])
                if p1 == 0
                else (alt.id, rpos[p1] - 1, rseq[rpos[p1] - 1 - 1 : rpos[p2] - 1], aseq[apos[p1] - 1])
            )


def parse_args(argv):
    parser = ArgumentParser(description="convert MSA -> VCF", formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", type=FileType(), help="the multiple sequence alignment as FASTA")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv[1:])

    with args.file as file:
        # read the reference
        ref = next(SeqIO.parse(file, "fasta"))
        file.seek(0)
        # read the alignment
        aln = SeqIO.parse(file, "fasta")
        # map the sequence id to its index
        keys = OrderedDict((ele.id, idx) for idx, ele in enumerate(aln))
        file.seek(0)
        # re-read
        aln = SeqIO.parse(file, "fasta")
        # output VCF header
        chrom = ref.id
        print("##fileformat=VCFv4.2")
        print(f"##contig=<ID={chrom},length={len(ref)}>")
        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *keys, sep="\t")
        # output variation
        for entry in msa_to_vcf(aln):
            key, pos, ref, alt = entry
            col = keys[key]
            print(chrom, pos, ".", ref, alt, ".", ".", ".", "GT", *(int(idx == col) for idx in range(len(keys))), sep="\t")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
