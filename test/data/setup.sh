#!/usr/bin/env bash

blastn -task blastn -query qry.fasta -subject lib.fasta -outfmt "7 qaccver saccver qstart qend qseq sseq length btop" > blastn.tsv
blastn -task blastn -query qry.fasta -subject lib.fasta > blastn.txt
fasta36 -m 8CB qry.fasta lib.fasta > fasta36.tsv
fasta36 -m 10 qry.fasta lib.fasta > fasta36.txt
glsearch36 -m 8CB qry.fasta lib.fasta > glsearch36.tsv
glsearch36 -m 10 qry.fasta lib.fasta > glsearch36.txt
