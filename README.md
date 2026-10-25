# Notothenioid hypoxia RNA-seq

Code and selected analysis files for a comparative RNA-seq study of hypoxia responses in Antarctic and sub-Antarctic notothenioid fishes. The repository covers sequencing-read processing, gene-level quantification, differential expression, gene annotation, cross-species orthology and downstream visualization.

**Reproducibility status:** This is a research analysis archive, not a self-contained, one-command pipeline. The main R Markdown notebook requires original sample metadata, count matrices, and local analysis files that are not all deposited here. See [Reproducibility and data requirements](docs/REPRODUCIBILITY.md).

## Study species

| Code | Species |
| --- | --- |
| ACE | *Chaenocephalus aceratus* |
| COR | *Notothenia coriiceps* |
| GEO | *Pseudochaenichthys georgianus* |
| GIB | *Gobionotothen gibberifrons* |
| MAC | *Eleginops maclovinus* |

The principal comparisons are **Hypoxia vs. Normoxia** and **ILOS vs. Normoxia**. GEO is analyzed for the ILOS contrast only. `ILOS` is the term used in many code objects; some later figures use **ILOP** for the same treatment.

## Repository contents

| Location | Description |
| --- | --- |
| [`General_Pipeline/`](General_Pipeline/) | Historical SLURM scripts for read trimming, HISAT2 alignment, index building and featureCounts |
| [`EdgeR_Analysis/`](EdgeR_Analysis/) | Main edgeR R Markdown analysis, helper functions and count-key mapping to zebrafish ortholog identifiers |
| [`eggNOG/`](eggNOG/) | eggNOG-mapper commands, annotation outputs and zebrafish orthology tables |
| [`Misc_Custom_Scripts/`](Misc_Custom_Scripts/) | Supplementary Bash and Python utilities for gene annotations, protein sequences and orthology |
| [`standardize_to_zf_one2one.R`](standardize_to_zf_one2one.R) | One-to-one zebrafish orthology table construction from eggNOG outputs |
| [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) | Input requirements, analysis order, software and known limitations |

Some directories contain exploratory or superseded scripts as well as later revisions; a filename alone should not be treated as proof that it produced a particular result.

## Analysis overview

```text
Paired-end RNA-seq FASTQs
  → adapter/quality trimming (BBDuk, Trimmomatic)
  → genome indexing and alignment (HISAT2)
  → gene-level fragment counting (featureCounts)
  → per-species differential expression (edgeR)
  → annotation and zebrafish orthology mapping (eggNOG)
  → within-/cross-species comparisons, enrichment and figures
```

The main expression-analysis files are:

- [`EdgeR_Analysis/Ant_fish_RNA-seq_L2FC=0_revised.Rmd`](EdgeR_Analysis/Ant_fish_RNA-seq_L2FC=0_revised.Rmd) — RNA-seq QC, species-specific differential expression and downstream analyses.
- [`EdgeR_Analysis/rnaseq_helpers_revised.R`](EdgeR_Analysis/rnaseq_helpers_revised.R) — functions called by the notebook.
- [`EdgeR_Analysis/Map_to_zf.R`](EdgeR_Analysis/Map_to_zf.R) — count-table gene keys to `DrOrtholog` mapping for cross-species comparisons.

In the version of the edgeR helpers currently deposited, the species-level model uses a treatment-only design (`~0 + treatment`), the control is Normoxia, and significance is evaluated at **FDR ≤ 0.05**, without a minimum absolute log2 fold-change threshold. The helper's initial count filter retains features with a mean count of at least 20. The notebook also contains additional and exploratory analyses; consult the code for the settings used in each section.

## Reuse

1. Obtain the corresponding FASTQs, genome assemblies/annotations, sample sheet, and species-level gene-count tables.
2. Install the versions of the command-line tools and R/Bioconductor packages required by the scripts.
3. Adapt the site-specific SLURM settings and absolute Windows/Chinook paths before execution.
4. Review [reproducibility requirements](docs/REPRODUCIBILITY.md) and run the analysis in stages, checking intermediate results.

**Important:** The R Markdown notebook is not expected to knit directly after cloning: it loads an external `.RData` workspace and references data files outside this repository. The committed eggNOG results support inspection of annotation and orthology outputs, but do not substitute for the original RNA-seq input data.

## Citation

Project code maintained by Augustus M. Snyder. Manuscript citation, raw-data accessions and verified software versions should be added when available.

The R Markdown analysis acknowledges code adapted from [Schweizer et al. (2023)](https://github.com/renaschweizer/pman_diaphragm_rnaseq).
