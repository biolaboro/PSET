#!/usr/bin/env bash

export PYTHONPATH="../.."

trap - SIGINT

cd "$(dirname "$0")" || exit

name="$(dirname "$0" | xargs basename)"
query="txid11266[Organism:exp] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP]"

esearch -db nuccore -query "$query" | \
  efetch -format uid | \
  head -n 500 | \
  esummary -db nuccore | \
  xtract -pattern DocumentSummary -element AccessionVersion TaxId > "$name.ssv"

cut -f 1 -d " " "$name.ssv" | \
    epost -db nuccore | \
    efetch -format fasta | \
    makeblastdb \
        -parse_seqids -hash_index \
        -blastdb_version 5 -dbtype nucl \
        -title "$name" -out "$name" -logfile "$name.log" -taxid_map "$name.ssv"
