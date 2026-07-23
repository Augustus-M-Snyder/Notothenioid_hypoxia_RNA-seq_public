#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ----------- CLI -----------
args <- commandArgs(trailingOnly = TRUE)
argval <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0 || i == length(args)) return(default)
  args[i + 1]
}
indir    <- argval("--in",  "emapper_zf")
outdir   <- argval("--out", "emapper_zf")
species  <- strsplit(argval("--species", "ACE,COR,GEO,GIB,MAC"), ",")[[1]] %>% trimws()
map_tsv  <- argval("--map_tsv", NULL)  # optional offline DR mapping TSV

if (!dir.exists(indir)) stop("Input dir not found: ", indir)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cat("[biomaRt] Attempting to fetch DR mappings from Ensembl OR load provided TSV...\n")

# ----------- Build zebrafish ID→gene mapping -----------
get_map_df <- function() {
  if (!is.null(map_tsv) && file.exists(map_tsv)) {
    m <- readr::read_tsv(map_tsv, show_col_types = FALSE)
    # Expect columns: ensembl_gene_id, ensembl_transcript_id, ensembl_peptide_id, external_gene_name
  } else {
    suppressPackageStartupMessages({
      library(biomaRt)
    })
    mart <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "drerio_gene_ensembl")
    m <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                     "ensembl_peptide_id", "external_gene_name"),
      mart = mart
    )
  }

  # Keep a stable copy of gene ID, then pivot to long so ANY id can map to gene id
  m <- m %>% dplyr::mutate(gene_id = .data$ensembl_gene_id)
  m_long <- m %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(c("ensembl_gene_id","ensembl_transcript_id",
                             "ensembl_peptide_id","external_gene_name")),
      names_to = "id_type", values_to = "id"
    ) %>%
    dplyr::filter(!is.na(.data$id), .data$id != "") %>%
    dplyr::select(.data$id, .data$gene_id) %>%
    dplyr::distinct()

  # return a named vector id -> gene_id
  stats <- m_long %>% dplyr::summarise(n_ids = dplyr::n_distinct(.data$id),
                                       n_genes = dplyr::n_distinct(.data$gene_id))
  cat(sprintf("[INFO] DR map: %s ids -> %s genes\n", stats$n_ids, stats$n_genes))
  setNames(m_long$gene_id, m_long$id)
}

id2gene <- get_map_df()

# ----------- Helpers -----------
read_ortholog_any_ids <- function(sp) {
  f <- file.path(indir, sprintf("%s_to_DR.emapper.orthologs", sp))
  if (!file.exists(f) || file.size(f) == 0) return(character())

  # Read lines, keep Danio rerio(7955) lines, grab last column (target id)
  ln <- readLines(f, warn = FALSE)
  ln <- ln[!startsWith(ln, "#")]
  ln_dr <- ln[grepl("Danio rerio\\(7955\\)", ln)]
  if (!length(ln_dr)) return(character())

  last_field <- vapply(strsplit(ln_dr, "\t"), function(x) x[length(x)], character(1))
  last_field <- sub("^\\*", "", last_field)

  # Extract any Ensembl DR tokens: gene (ENSDARG), protein (ENSDARP), transcript (ENSDART), plus fallback gene symbols
  ens_any <- stringr::str_extract(last_field, "(ENSDARG[0-9]+|ENSDARP[0-9]+|ENSDART[0-9]+)")
  ens_any <- ens_any[!is.na(ens_any)]

  # If no Ensembl tokens found on a line, also try a conservative symbol token (letters+digits+._-), but not needed if Ensembl present
  if (!length(ens_any)) {
    syms <- stringr::str_extract(last_field, "[A-Za-z0-9._-]+")
    ens_any <- syms[!is.na(syms)]
  }

  unique(ens_any)
}

count_fasta_headers <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  con <- file(path, "r"); on.exit(close(con))
  n <- 0L
  repeat {
    chunk <- readLines(con, n = 100000, warn = FALSE)
    if (!length(chunk)) break
    n <- n + sum(substr(chunk, 1, 1) == ">")
  }
  n
}

find_input_fasta <- function(sp) {
  cands <- c(
    sprintf("%s_cds_longest_per_gene.prot.fa", sp),
    sprintf("%s_cds_longest_per_gene.faa", sp),
    sprintf("%s_cds_longest_per_gene.fa", sp)
  )
  existing <- cands[file.exists(cands)]
  if (length(existing)) existing[[1]] else NA_character_
}

# ----------- Main loop over species -----------
summary_rows <- list()
per_species_genes <- list()

for (sp in species) {
  any_ids <- read_ortholog_any_ids(sp)

  # Map ANY → ENSEMBL GENE using id2gene (gene, transcript, protein, or symbol if present in map)
  gene_ids <- unname(id2gene[any_ids])
  gene_ids <- sort(unique(na.omit(gene_ids)))

  # Write outputs
  readr::write_lines(gene_ids, file.path(outdir, sprintf("%s_zfish_keys.txt", sp)))
  readr::write_lines(any_ids,  file.path(outdir, sprintf("%s_any_ids.txt", sp)))

  # Count input sequences
  fasta_path <- find_input_fasta(sp)
  fasta_n <- if (is.na(fasta_path)) NA_integer_ else count_fasta_headers(fasta_path)

  summary_rows[[sp]] <- tibble::tibble(
    species = sp,
    input_fasta = ifelse(is.na(fasta_path), "", fasta_path),
    input_seqs  = fasta_n,
    mapped_any_ids   = length(any_ids),
    mapped_gene_ids  = length(gene_ids)
  )

  per_species_genes[[sp]] <- gene_ids
  cat(sprintf("[OK] %s: input=%s; ANY=%d; GENE=%d\n",
              sp, ifelse(is.na(fasta_path), "NA", fasta_path),
              length(any_ids), length(gene_ids)))
}

summary_df <- dplyr::bind_rows(summary_rows)
readr::write_tsv(summary_df, file.path(outdir, "mapping_summary.tsv"))

# ----------- Presence/absence matrix -----------
all_genes <- sort(unique(unlist(per_species_genes)))
if (length(all_genes)) {
  mat <- dplyr::tibble(zfish_key = all_genes)
  for (sp in species) {
    have <- per_species_genes[[sp]]
    mat[[sp]] <- as.integer(all_genes %in% have)
  }
  readr::write_tsv(mat, file.path(outdir, "presence_matrix.tsv"))
  cat("[OK] Wrote presence_matrix.tsv and mapping_summary.tsv to: ", outdir, "\n", sep = "")
} else {
  cat("[WARN] No gene-level DR IDs produced. Check ortholog files and mapping.\n")
}

