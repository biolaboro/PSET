#!/usr/bin/env python3

import sys

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

path = sys.argv[1]
path = sys.stdin if path == "-" else path

for record in SeqIO.parse(path, "fasta"):
    print(record.id, seguid(record), sep="\t")
