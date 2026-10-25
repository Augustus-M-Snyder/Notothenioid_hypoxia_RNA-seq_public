# Reproducibility and data requirements

This document describes the inputs and execution assumptions of the scripts currently committed to this repository. It is not a claim that the published or manuscript results can be regenerated from a fresh clone alone.

## Required inputs

| Input | Used by | Availability here |
| --- | --- | --- |
| Paired-end FASTQ files | `General_Pipeline/` | Not included |
| Species-specific reference genomes and GTF annotations | HISAT2 indexing, mapping, featureCounts | Not included |
| `Ant_fish_RNAseq_metadata.xlsx` | Main edgeR notebook | Not included |
| `Ant_fish_RNAseq_qc_metrics.xlsx` | Sequencing-QC analyses | Not included |
| `DEG_WGCNA_environment.RData` | Main notebook initialization | Not included |
| `ACE_counts_aggregated_by_GeneSymbol_corrected.txt` and `COR_counts.txt`, `GEO_counts.txt`, `GIB_counts.txt`, `MAC_counts.txt` | Main edgeR notebook | Expected at user-specific local paths; not supplied together as these named inputs |
| `ACE_protein_gene_mapping.csv` | `EdgeR_Analysis/Map_to_zf.R` | Expected at a user-specific path |
| `*_to_DR.ids.tsv` (one-to-one orthology mappings) | `Map_to_zf.R` | Related files are present in `eggNOG/emapper_zf_one2one/`; the script expects a different local directory layout |

The repository contains some historical featureCounts tables and eggNOG output files. Those are **not guaranteed to be identical to the five count-table inputs used by the current revised edgeR notebook**. Verify gene IDs, sample names, counting features, reference builds and strandedness before substituting any available counts.

## Analysis order and entry points

1. **Read processing:** `General_Pipeline/Step_1_modified.sh` (BBDuk), then `Step_2_Trimmomatic.sh`.
2. **Reference indexing and alignment:** `General_Pipeline/hisat2_builder.sh`, then `Step_3_HiSat2.sh`.
3. **Gene counts:** The `General_Pipeline/Step_4_featurecounts*.sh` files represent different counting approaches (grouped/separate samples, exon/gene features). Confirm the actual count-generation settings and files before choosing a version.
4. **Differential expression:** `EdgeR_Analysis/Ant_fish_RNA-seq_L2FC=0_revised.Rmd` sources `rnaseq_helpers_revised.R`. Initial input, QC and edgeR model sections precede figures and enrichment.
5. **Zebrafish orthology:** `eggNOG/`, `standardize_to_zf_one2one.R`, and `EdgeR_Analysis/Map_to_zf.R` contain distinct mapping stages. The revised notebook uses the latter's `DrOrtholog` mapping tables for cross-species DEG overlaps.

These are *logical* steps, not a validated makefile or an automatically executable workflow. Some preprocessing scripts were written for earlier HPC environments, and some orthology outputs were generated through separate pipelines.

## Differential-expression assumptions

The deposited edgeR helpers implement the following key choices:

- Count-table filtering: feature mean count **≥20** before model fitting.
- Normalization: `edgeR::calcNormFactors()`.
- Design: `model.matrix(~0 + treatment)`; Normoxia is the control. Sex is recorded in metadata but **is not a covariate in this design**.
- Fitting/testing: `glmQLFit()` and `glmQLFTest()`, with FDR ≤0.05 for significant DEGs.
- Hypoxia and ILOS contrasts for ACE, COR, GIB and MAC; ILOS only for GEO.
- The `L2FC=0` analysis does **not** impose a minimum magnitude of absolute log2 fold change.

Other sections of the notebook perform WGCNA-related work, orthology overlaps and functional enrichment; their input objects and analysis-specific thresholds must be checked separately. The script may depend on objects loaded from the saved `.RData` workspace.

## Orthology and annotation cautions

Cross-species DEG comparisons depend on how species-specific featureCounts identifiers are mapped to shared zebrafish-related keys.

- `EdgeR_Analysis/Map_to_zf.R` explicitly labels its cross-species matching key **`DrOrtholog`** and cautions that this is not necessarily a verified zebrafish Ensembl gene ID.
- One of the earlier normalization scripts (`standardize_to_zf_one2one.R`) attempts protein/transcript-to-gene ID conversion by replacing Ensembl identifier prefixes. This approach should **not** be assumed to establish a genuine protein/transcript-to-gene relationship without independent annotation lookup.
- One-to-one orthology filters, ambiguous mappings, species-specific identifiers and mapping losses should be recorded before calculating shared DEG sets or performing enrichment.
- GO/KEGG analysis requires an explicitly verified identifier type and an appropriate gene universe; a `DrOrtholog` label alone is not evidence of a valid Ensembl/Entrez mapping.

## Software and environment

The scripts reference SLURM, Bash, BBDuk, Trimmomatic, HISAT2, samtools, Subread/featureCounts, Python, eggNOG-mapper, and R with edgeR and other CRAN/Bioconductor packages. The notebook loads or calls packages including `edgeR`, `limma`, `DESeq2`, `WGCNA`, `tidyverse`, `lme4`, `lmerTest`, `gprofiler2`, `UpSetR`, `VennDiagram`, `ggplot2`, `ggpubr`, and various plotting/enrichment packages.

No complete version-locked environment accompanies the main analysis. Before a reproducibility release, record `sessionInfo()` from the working R session; record `hisat2 --version`, `featureCounts -v`, `emapper.py --version`, and reference assembly/annotation accessions and checksums. External annotation services and databases should also be versioned where possible.

## Path setup and safe execution

The main notebook sets a hard-coded Windows working directory, then calls `source("rnaseq_helpers_revised.R")` and `load("DEG_WGCNA_environment.RData")`. Other scripts contain local absolute paths and HPC-specific modules, resource requests and email settings. These paths and settings must be changed to fit the user's system.

Recommended sequence for reproducing a result:

1. Begin from a clean copy of the input data, preserving original FASTQs and metadata.
2. Check that count-table sample columns match the metadata table in the same order.
3. Confirm which reference assemblies, count features and counting parameters produced the chosen matrices.
4. Run QC and edgeR sections before any analysis dependent on `sig.*`, `DrOrtholog` or WGCNA objects.
5. Save the Git commit hash, package/tool versions, input file checksums and exact commands for each completed run.
6. Write outputs to a separate analysis directory instead of overwriting historical files.

**Current status:** Code and selected derived results are available, but the archive lacks the complete input dataset, local workspaces and environment lockfile required for verified end-to-end reruns.
