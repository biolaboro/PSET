# PCR Signature Erosion Tool

## Abstract

The PCR Signature Erosion Tool (PSET) calculates the *in silico* detection capability of assays based on sequence alignment. The workflow generates a confusion matrix based on assay type, sequence alignment, and taxonomic lineage.

## Example

Create a [Conda](https://docs.conda.io/en/latest/) environment with [Mamba](https://github.com/mamba-org/mamba) and activate the **pset** environment.
```bash
mamba env create -f env/env.yml
conda activate pset
```

Download and install [taxa](https://github.com/biolaboro/taxa) and build the NCBI Taxonomy database.
```bash
URL="ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
target="$(basename "$URL")"
curl -O "$URL" && \
  python3 -m taxa.taxa -database taxa.db create "$target" && \
  rm -f "$target"
```

Download a data set of 500 *Ebolavirus* sequences and build a BLAST+ database.
```bash
bash -x data/ebola/ebola.sh
```

Set the **PYTHONPATH** environment variable.
```bash
export PYTHONPATH=.
```

Run the analysis.
```bash
snakemake \
  --cores 8 \
  --config \
    file=data/ebola/ebola.tsv \
    db=data/ebola/ebola \
    confb=task=blastn,num_threads=8,num_alignments=10000,subject_besthit \
    confg=E=10000,T=8 \
    confm=auto,thread=8 \
    dburl=sqlite:///taxa.db \
    out=ebola \
    nseq=100 \
    -- \
  all
```

Collect confusion matrix results.
```bash
find ebola \
  -name mat.tsv \
  -type f \
  -exec awk \
  -F $'\t' \
  -v OFS=$'\t' 'NR > 1 { split(FILENAME, t, "/"); m[t[2]"\t"$4]++; } END { for (k in m); print k, m[k]; }' {} \; | \
  sort -k 2,1
```

View results...
```
EBO1_2	TP	4
EBO3_4	FN	1
EBOGP	PT	451
EBO_GP	PT	1
ENZ	FNN	1
EboZNP	TN	4
Ebola_Bundibugyo_MGB	PT	2
Ebola_Bundibugyo_TM	PT	2
Ebola_Ivory_Coast_TM	TNN	2
Ebola_Sudan_MGB	TN	459
Ebola_Sudan_TM	TP	2
Filo_AB	FN	6
GAB_1	FNN	1
KGH	PT	447
Kulesh_MGB	PT	18
Kulesh_TM	FNN	3
NGDS_Primary_amplicon	TP	2
NGDS_Secondary_amplicon	PT	18
PanFiloL3_4	FN	9
PanFiloL_1_2	FN	1
Reston	TN	460
Sudan	TP	2
ZAI_NP	FN	5
ZebovGP	TP	455
pan_Ebola_Assay_MGB_EBOV	TP	426
pan_Ebola_Assay_MGB_RESTV	TN	433
pan_Ebola_Assay_MGB_SUDV	TP	2
```

## Assay

An assay shall consist of an **id**, **definition**, and **targets** list.

| Property | Description |
|-|-|
| id | the assay identifier and/or name |
| definition | the delimited DNA sequence that indicates component regions |
| targets | the list of NCBI Taxonomy identifiers that the assay targets |

### ID

The identifier shall be a unique name consisting of characters compatible with the file system.

### Definition

The assay definition implicitly defines the assay type. Each definition consists of IUPAC ambiguous DNA letter codes with delimiters that surround component regions. In other words, the definition is an amplicon sequence that indicates primer and probe regions. Definitions may also include additional sequence outside the amplicon region, termed context. There are four definition types.

| Assay | Format | Example | Description |
|-|-|-|-|
| amplicon | [] | G[ATTAC]A | amplicon only, no primers |
|primer | [][] | [GAT]T[AC]A | forward and reverse primer |
| probe | []\(\)[], [\(]\)[], []\([\)] | [G]A(TT[A)CA] | a primer assay with a probe that may overlap with the primers |
| nested | []{}{}[], [{]}{}[], []{}{[}], [{]}{[}] | [CA{T]GA}TT{ACA}[CAT] | a primer assay with another internal primer assay that may overlap with the outer primers |

### Targets

Each NCBI Taxonomy identifier (tax id) in the set of targets corresponds a node in the taxonomy tree. Any subject sequence with an assigned taxonomy in the set of assay targets is itself a target. Also, sequence with a taxonomic ancestor in the set of assay targets is itself a target. For example, if an assay targets set includes *Vibrio cholerae* (666), then any sequence with this ancestry is a target, such as *Vibrio cholerae O1 biovar El Tor* (686) with the ancestors "686 -> 127906 -> **666** -> 662 -> 641 -> 135623 -> 1236 -> 1224 -> 2 -> 131567 -> 1".

## Workflow

The PSET workflow consists of several sequence alignment phases.

### Local Alignment

The objective of this phase is to query the assay definition amplicon against a BLAST+ database to search for matching subjects. Any sequence context is also included in the search. The query sequence is expanded to remove any ambiguous DNA codes since they are incompatible with BLAST+. In other words, the query sequence represents the first permutation given the set of alternative letters represented at each ambiguous position.

|code|set|
|-|-|
|M|{A, C}|
|R|{A, G}|
|W|{A, T}|
|Y|{C, T}|
|S|{G, C}|
|K|{G, T}|
|V|{A, C, G}|
|H|{A, C, T}|
|D|{A, G, T}|
|B|{C, G, T}|
|N|{A, C, G, T}|

For example, "GA**W**TA**Y**A" has two ambiguous codes, representing four possible sequence expansions: "GA**A**TA**C**A", "GA**T**TA**C**A", "GA**A**TA**T**A", and "GA**T**TA**T**A". In such a case, the procedure only queries the first expansion.

This means that an additional step must re-evaluate the BLAST+ identity statistics using the original ambiguous query. Additionally, only the amplicon region is evaluated, meaning that any sequence context present in the alignment is ignored. Accordingly, this procedure scores only the amplicon region of the alignment for sequence similarity. Note however, that the unknown DNA character code **N** is always penalized if it is on the subject. The procedure will extend the alignment to cover the entire query range for similarity score.

### Paper

See publications for more details...
