#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "$"), args)
  if (length(hit) == 0 || hit == length(args)) return(default)
  args[hit + 1]
}
arg_flag <- function(flag) any(args == flag)

IN_DIR  <- arg_val("--in",  stop("Missing --in"))
OUT_DIR <- arg_val("--out", stop("Missing --out"))
SPECIES <- strsplit(arg_val("--species", stop("Missing --species")), ",")[[1]] |> trimws()

if (!dir.exists(IN_DIR))  stop("Input dir not found: ", IN_DIR)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

message("[CONFIG] IN_DIR=", IN_DIR, "  OUT_DIR=", OUT_DIR, "  SPECIES={", paste(SPECIES, collapse=","), "}")

# ---------- helpers ----------
read_orthologs <- function(path) {
  # emapper *.orthologs are 4+ tab cols; we only need first 4
  if (!file.exists(path)) return(NULL)
  df <- suppressMessages(readr::read_tsv(
    path, comment = "#", col_names = FALSE,
    col_types = readr::cols(.default = "c"), progress = FALSE
  ))
  if (nrow(df) == 0) return(df)
  # pad or trim to 4 columns
  maxc <- max(4, ncol(df))
  if (ncol(df) < 4) {
    for (k in (ncol(df)+1):4) df[[paste0("X", k)]] <- NA_character_
  }
  df <- df[, 1:4, drop = FALSE]
  names(df) <- c("query_name","relation","target","target_id")
  df
}

count_fasta_seqs <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  # count header lines '>'
  as.integer(system2("bash", c("-lc", sprintf("grep -c '^>' %s || echo 0", shQuote(path))), stdout = TRUE))
}

# Convert Ensembl transcript/protein → gene by stable-id prefix swap
to_gene_id <- function(x) {
  x <- str_remove(x, "^\\*") # emapper marks chosen target with '*'
  case_when(
    str_detect(x, "^ENSDARG\\d+$") ~ x,
    str_detect(x, "^ENSDARP\\d+$") ~ str_replace(x, "^ENSDARP", "ENSDARG"),
    str_detect(x, "^ENSDART\\d+$") ~ str_replace(x, "^ENSDART", "ENSDARG"),
    TRUE                           ~ NA_character_
  )
}

# choose input FASTA (for summary)
pick_input_fasta <- function(sp) {
  cands <- c(
    sprintf("%s_cds_longest_per_gene.prot.fa", sp),
    sprintf("%s_cds_longest_per_gene.faa", sp),
    sprintf("%s_cds_longest_per_gene.fa",  sp)
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

# ---------- main per-species loop ----------
per_sp_keys <- list()
summary_rows <- list()

for (sp in SPECIES) {
  ortho_path <- file.path(IN_DIR, sprintf("%s_to_DR.emapper.orthologs", sp))
  ann <- read_orthologs(ortho_path)

  if (is.null(ann) || nrow(ann) == 0) {
    warning("[WARN] No orthologs for ", sp, " at ", ortho_path)
    per_sp_keys[[sp]] <- character(0)
    summary_rows[[sp]] <- tibble(
      species = sp,
      input_fasta   = pick_input_fasta(sp),
      input_seqs    = count_fasta_seqs(pick_input_fasta(sp)),
      mapped_any_ids  = 0L,
      mapped_gene_ids = 0L
    )
    next
  }

  # Only zebrafish & one-to-one
  ann_dr <- ann %>%
    filter(relation == "one2one",
           str_detect(target, "Danio rerio\\(7955\\)"))

  # ANY DR IDs (protein/transcript/gene)
  ids_any <- ann_dr$target_id %>% to_gene_id() # normalize to gene where possible by prefix swap
  ids_gene <- ids_any[!is.na(ids_any)] %>% unique()

  # write per-species results
  write_lines(ids_gene, file.path(OUT_DIR, sprintf("%s_zfish_keys.txt", sp)))

  # also write pairs: zfish_key \t species_gene
  pairs <- tibble(
    zfish_key    = ids_gene,
    species_gene = ann_dr$query_name[!is.na(ids_any)][match(ids_gene, ids_any[!is.na(ids_any)])]
  ) %>% distinct()
  readr::write_tsv(pairs, file.path(OUT_DIR, sprintf("%s_to_DR.ids.tsv", sp)))

  per_sp_keys[[sp]] <- ids_gene

  # summary
  infa <- pick_input_fasta(sp)
  summary_rows[[sp]] <- tibble(
    species          = sp,
    input_fasta      = infa,
    input_seqs       = count_fasta_seqs(infa),
    mapped_any_ids   = length(unique(ann_dr$target_id)), # before prefix-swap NA filtering
    mapped_gene_ids  = length(ids_gene)                  # after gene-level standardization
  )

  message(sprintf("[OK] %s: one2one zebrafish genes = %d", sp, length(ids_gene)))
}

# ---------- presence/absence matrix ----------
all_keys <- sort(unique(unlist(per_sp_keys)))
mat <- tibble(zfish_key = all_keys)
for (sp in SPECIES) {
  present <- all_keys %in% per_sp_keys[[sp]]
  mat[[sp]] <- as.integer(present)
}
readr::write_tsv(mat, file.path(OUT_DIR, "presence_matrix.tsv"))
message("[OK] Wrote matrix: ", file.path(OUT_DIR, "presence_matrix.tsv"),
        " (", nrow(mat), " rows)")

# ---------- summary ----------
summary_df <- bind_rows(summary_rows) %>%
  arrange(species)
readr::write_tsv(summary_df, file.path(OUT_DIR, "mapping_summary.tsv"))
message("[OK] Wrote summary: ", file.path(OUT_DIR, "mapping_summary.tsv"))

# ---------- quick console peek ----------
print(summary_df, n = nrow(summary_df))

