# Abstract

The PCR Signature Erosion Tool (PSET) calculates the *in silico* detection capability of PCR assays based on assay type, sequence alignment, and the taxonomic lineage of subject sequences.

# Methods

## Assay

An assay shall consist of an **id**, **definition**, and **targets**.

| Property | Description |
|-|-|
| id | the identifier |
| definition | the definition delimits all primer/probe regions, including optional 5'/3'-context |
| targets | the set of NCBI Taxonomy identifiers that the assay targets |

### id

The identifier shall be a unique name consisting of characters compatible with the file system.

### definition

The assay definition shall implicitly define the assay type based on its format. Each definition consists of IUPAC ambiguous DNA letter codes with delimiters surroundind assay component regions. The definition contains an amplicon sequence with primer and probe regions delimited. Any additional sequence outside the amplicon region is termed context. The general assay format consists of bracketed primer region(s).

Definition format:

```
5'-context-[primer1]-interprimer-(probe)-interprimer-[primer2]-context-3'
```

Definition terms:

- An assay **component** is a primer or probe sequence.
- A **delimiter** is a pair of open/close symbols to identify **component** regions, such as the "[]", "()", and "{}" bracket pairs, based on the corresponding assay format.
- The **amplicon** is the non-delimited amplicon sequence that includes all component/interprimer regions
- The **camplicon** (contextualized amplicon) is the non-delimited amplicon sequence that includes 5'/3'-context. This is equivalent to the non-delimited definition.

Definition types: amplicon, primer, probe, and nested.

| Assay | Examples | Description |
|-|-|-|
| amplicon | G[ATTAC]A | amplicon only, no primers |
| primer | [GAT]T[AC]A | forward and reverse primer |
| probe | [G]A(TT[A)CA] | a primer assay with a probe that may overlap with the primers |
| nested | [CA{T]GA}TT{ACA}[CAT] | a primer assay with another internal primer assay that may overlap with the outer primers |

### targets

Each NCBI Taxonomy identifier (tax id) in the set of targets corresponds to a node in the taxonomy tree. Any subject sequence bearing a tax id in the set is a positive hit. This is also true if the subject is a descendant of any of the tax id numbers in the set. For example, if an assay target set includes *Vibrio cholerae* (**666**), then any sequence with this ancestry is a target, such as *Vibrio cholerae O1 biovar El Tor* (**686**) having ancestors "**686** -> 127906 -> **666** -> 662 -> 641 -> 135623 -> 1236 -> 1224 -> 2 -> 131567 -> 1".

## Pipeline

The PSET workflow consists of several sequence alignment, scoring, and taxonomy evaluation phases. The objective is to determine whether the assay components alignmed with sufficient coverage and identity and with the correct arrangement and orientation to subjects bearing the targeted taxonomic identifier.

### Local Alignment

*Perform local alignment of the camplicon, which includes sequence context to promote local alignment near the 5'/3'-ends of the primers.*

The objective of this phase is to query the assay definition against a BLAST+ database to search for matching subjects. Any sequence context in the definition is also included in the search. The query sequence is expanded to remove any ambiguous DNA codes since they are incompatible with BLAST+. In other words, the query sequence represents the first permutation given the set of alternative letters represented at each ambiguous position.

|code|set|code|set|code|set|
|-|-|-|-|-|-|
|M|{A, C}|V|{A, C, G}|N|{A, C, G, T}|
|R|{A, G}|H|{A, C, T}|||||
|W|{A, T}|D|{A, G, T}|||||
|Y|{C, T}|B|{C, G, T}|||||
|S|{G, C}|||||||
|K|{G, T}|||||||

For example, "GA**W**TA**Y**A" has two ambiguous codes, representing four possible sequence expansions: "GA**A**TA**C**A", "GA**T**TA**C**A", "GA**A**TA**T**A", and "GA**T**TA**T**A". The procedure only queries the first expansion.

### Filter

*Score alignments within amplicon region and filter based on query coverage and similarity threshold.*

This step re-evaluates BLAST+ alignment statistics using the original, potentially ambiguous query. Only the amplicon region is evaluated for sequence similarity, meaning that any sequence context present in the alignment is ignored. Note however, that the unknown DNA character code **N** is always penalized if it is on the subject. Subject sequences meeting or exceeding the similarity threshold are then extracted such that it matches the original query length and coordinates.

```
       [       camplicon       ]
           [    amplicon   ]
qry   5'---[-----]---[-----]---3' ( id  ) x ( cov ) = ( sim )
sbj-1   101 11111 110 10111 111    9/10   x 10/10   =   90%   ✓
sbj-2       11110 111 10011 1      7/10   x 10/10   =   70%   x
sbj-3         101 111 111          5/6    x  6/10   =   50%   x
```

Only the unique set of sequences is extracted. Therefore, this step also outputs a mapping of all of the accessions bearing the extracted sequence.

### Global/Local Alignment

*Re-align components individually to the extracted subjects.*

Global/local (glocal) alignment guarantees complete alignent of each query component to the extracted subject sequences. This step re-aligns each primer/probe to the extracted set of subjects and are then scored similarly.

### Calls

*The final step outputs a confusion matrix call for each subject based on glocal alignment statistics and subject taxonomy.*

A call is made for each subject. The calculation is based on whether all primers aligned to a subject bearing the targeted taxonomy identifier in the set of assay targets (or is a descendant of one of them) with the correct arrange, orientation, and similarity at or above threshold.

|Alignment|Taxonomy|Call|Interpretation|
|-|-|-|-|
|`✓`|`✓`|TP|Good alignment, on-target.|
|`✓`|`x`|FP|Good alignment, off-target.
|`x`|`✓`|FN|Bad alignment, on-target.|
|`x`|`x`|TN|Bad alignment, off-target.|

A special call "XX" is reserved for alignments to synthetic constructs.

# Results

Result files are available on the sidebar and organized by under the individual "Assay" and aggregate "Summary" folders. The following figure is a heatmap of calls made for each assay against the unique set of extracted subject sequences. Color indicates the heat value, which is equal to the average similarity of each assay component to the subject.
