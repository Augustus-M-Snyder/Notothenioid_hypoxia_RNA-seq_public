# Notothenioid hypoxia RNA-seq

RNA-seq analysis of cardiac responses to hypoxia in five Antarctic and sub-Antarctic notothenioid fishes.

| Code | Species |
| --- | --- |
| ACE | *Chaenocephalus aceratus* |
| COR | *Notothenia coriiceps* |
| GEO | *Pseudochaenichthys georgianus* |
| GIB | *Gobionotothen gibberifrons* |
| MAC | *Eleginops maclovinus* |

The study compares normoxia with hypoxia and an additional low-oxygen treatment (ILOP, labeled `ILOS` in the analysis code). GEO includes the ILOP comparison only.

## Repository contents

| Directory | Contents |
| --- | --- |
| [`General_Pipeline/`](General_Pipeline/) | Read trimming, HISAT2 alignment and featureCounts scripts |
| [`EdgeR_Analysis/`](EdgeR_Analysis/) | Differential expression, WGCNA and comparative analyses in R |
| [`eggNOG/`](eggNOG/) | Orthology mapping and annotation results |
| [`Misc_Custom_Scripts/`](Misc_Custom_Scripts/) | Additional annotation and sequence-processing utilities |

## Analysis

The principal analysis is [`Ant_fish_RNA-seq_L2FC=0_revised.Rmd`](EdgeR_Analysis/Ant_fish_RNA-seq_L2FC=0_revised.Rmd), supported by [`rnaseq_helpers_revised.R`](EdgeR_Analysis/rnaseq_helpers_revised.R) and [`Map_to_zf.R`](EdgeR_Analysis/Map_to_zf.R).

The workflow includes read processing, gene-level quantification, edgeR differential expression, zebrafish ortholog mapping, functional enrichment and visualization.

Some analysis inputs and saved R objects are not included, and local file paths require adjustment before running the code. See [Reproducibility and data requirements](docs/REPRODUCIBILITY.md) for details.

## Citation

Analysis by Augustus M. Snyder. Portions of the R workflow were adapted from [Schweizer et al. (2023)](https://github.com/renaschweizer/pman_diaphragm_rnaseq).
