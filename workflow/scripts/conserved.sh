#!/usr/bin/env bash

tax="$1"
len="$2"
met="${3:-0}"
dwn="${4:-0}"

min=20
pct=0.10
lab="$tax"
eva=10000
thr=64

# download and untar
[ "$dwn" -eq 1 ] && \
    datasets download genome taxon "$tax" --include genome --assembly-level complete --exclude-atypical --filename "$lab.zip" && \
    unzip -d "$lab" "$lab.zip"
# filter sequences
find "$lab" -type f -name "*.fna" | while read -r ele; do cat "$ele"; echo; done | \
    ../../../repo/pset/workflow/scripts/length_filter.py - "$len" "$pct" > "$lab.lib.fna"

tax=518
lab="$tax"
thr=8
min=20
i=0
find . -type f -name "dist.[0-9]*.tsv" | \
    sort | \
    while read -r ele; do
        # subset lib
        awk 'NR > 1 { print $2 }' "$ele" | xargs -L 8 samtools faidx "$lab.lib.fna" > "$i.lib.fna"
        # grab a random reference
        head -n 1 "$ele" | cut -f 1 -d ' ' | tr -d \> | xargs samtools faidx "$i.lib.fna" > "$i.ref.fna"
        # calculate and filter match ranges to the reference
        nucmer --maxmatch --minmatch "$min" --prefix "$i" --threads "$thr" "$i.ref.fna" "$i.lib.fna"
        # delta-filter -i 100 "$i.delta" > "$i.100.delta"
        # show-coords -T "$i.100.delta" | awk -v OFS=$"\t" 'NR > 4 { print $1, $2, $9; }' > "$i.100.tsv"
        # # calculate conserved ranges
        # ./workflow/scripts/mum_to_conserved.py "$i.ref.fna" "$i.lib.fna" "$i.100.tsv" -proc "$thr" > "$i.json"
        ((i++))
    done


# grab a random reference
head -n 1 "$lab.lib.fna" | cut -f 1 -d ' ' | tr -d \> | xargs samtools faidx "$lab.lib.fna" > "$lab.ref.fna"
# calculate and filter match ranges to the reference
nucmer --maxmatch --minmatch "$min" --prefix "$lab" --threads "$thr" "$lab.ref.fna" "$lab.lib.fna"
delta-filter -i 100 "$lab.delta" > "$lab.100.delta"
show-coords -T "$lab.100.delta" | awk -v OFS=$"\t" 'NR > 4 { print $1, $2, $9; }' > "$lab.100.tsv"
# calculate conserved ranges
./workflow/scripts/mum_to_conserved.py "$lab.ref.fna" "$lab.lib.fna" "$lab.100.tsv" -proc 64 > "$lab.json"
