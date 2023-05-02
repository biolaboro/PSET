# PCR Signature Erosion Tool

## Abstract

The PCR Signature Erosion Tool (PSET) calculates the *in silico* detection capability of assays based on assay type, sequence alignment, and taxonomic lineage of subject sequences.

## Setup

Create a [Conda](https://docs.conda.io/en/latest/) environment with [Mamba](https://github.com/mamba-org/mamba) and activate the **pset** environment.
```bash
mamba env create -f workflow/envs/env.yml
conda activate pset
```

Set the **PYTHONPATH** environment variable.
```bash
export PYTHONPATH=.
```

Build the NCBI Taxonomy database.
```bash
snakemake --cores 1 -s workflow/rules/setup.smk taxa
```

Download a data set of 500 *Ebolavirus* sequences and build a BLAST+ database.
```bash
snakemake --cores 1 -s workflow/rules/setup.smk data
```

## Example

Run the analysis.
```bash
snakemake \
  --printshellcmds \
  --cores 8 \
  --set-threads \
    local=4 \
    glocal=4 \
  --config \
    file=resources/assay/EBOV.tsv \
    db=resources/blast/EBOV/EBOV \
    out=results/EBOV \
    -- \
  target_tsv
```

View results (numbers may vary depending on database build date)...
```bash
cat results/EBOV/EBOV/con.tsv
```

```tsv
db	id	TP	FN	FP	TN
resources/blast/EBOV/EBOV	Ebola_Bundibugyo_MGB	3	0	0	0
resources/blast/EBOV/EBOV	Ebola_Bundibugyo_TM	3	0	0	0
resources/blast/EBOV/EBOV	EBO_GP	397	27	3	4
resources/blast/EBOV/EBOV	EBO1_2	421	4	0	12
resources/blast/EBOV/EBOV	EBO3_4	425	0	0	12
resources/blast/EBOV/EBOV	EBOGP	424	2	0	16
resources/blast/EBOV/EBOV	EboZNP	434	0	0	10
resources/blast/EBOV/EBOV	ENZ	425	4	0	0
resources/blast/EBOV/EBOV	GAB_1	408	3	0	27
resources/blast/EBOV/EBOV	KGH	423	0	0	0
resources/blast/EBOV/EBOV	Kulesh_MGB	425	1	0	8
resources/blast/EBOV/EBOV	Kulesh_TM	427	2	0	8
resources/blast/EBOV/EBOV	NGDS_Primary_amplicon	424	0	0	0
resources/blast/EBOV/EBOV	NGDS_Secondary_amplicon	425	1	0	0
resources/blast/EBOV/EBOV	pan_Ebola_Assay_MGB_EBOV	403	3	0	16
resources/blast/EBOV/EBOV	ZAI_NP	398	2	0	18
resources/blast/EBOV/EBOV	ZebovGP	424	4	0	12
resources/blast/EBOV/EBOV	Filo_AB	418	18	0	0
resources/blast/EBOV/EBOV	PanFiloL_1_2	430	6	0	0
resources/blast/EBOV/EBOV	PanFiloL3_4	437	3	0	0
resources/blast/EBOV/EBOV	Ebola_Reston_MGB	0	0	0	0
resources/blast/EBOV/EBOV	Ebola_Reston_TM	1	0	0	0
resources/blast/EBOV/EBOV	pan_Ebola_Assay_MGB_RESTV	0	0	0	409
resources/blast/EBOV/EBOV	Reston	1	0	0	441
resources/blast/EBOV/EBOV	Ebola_Sudan_MGB	5	0	0	437
resources/blast/EBOV/EBOV	Ebola_Sudan_TM	5	0	0	431
resources/blast/EBOV/EBOV	pan_Ebola_Assay_MGB_SUDV	5	0	0	410
resources/blast/EBOV/EBOV	Sudan	5	0	0	0
resources/blast/EBOV/EBOV	Ebola_Ivory_Coast_MGB	1	0	0	0
resources/blast/EBOV/EBOV	Ebola_Ivory_Coast_TM	1	0	0	442
```

Optionally, continue the analysis in report mode.
```bash
snakemake \
  --printshellcmds \
  --cores 8 \
  --set-threads \
    local=4 \
    glocal=4 \
  --config \
    file=resources/assay/EBOV.tsv \
    db=resources/blast/EBOV/EBOV \
    out=results/EBOV \
    -- \
  target_report
```

Now, generate the report with the `--report` flag.
```bash
snakemake \
  --report \
  --config \
    file=resources/assay/EBOV.tsv \
    db=resources/blast/EBOV/EBOV \
    out=results/EBOV \
    -- \
  target_report
```

Then open report.html.
