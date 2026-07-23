# ================================================================
# Map_to_zf.R — SIMPLE / REPRODUCIBLE VERSION
#
# Purpose:
#   Map species-specific count-table keys to shared zebrafish/EggNOG
#   ortholog IDs for cross-species overlap analyses.
#
# What this script DOES:
#   - Reads species count table keys
#   - Reads *_to_DR.ids.tsv files
#   - Extracts gene_id and gene_name/gene_symbol from species_gene
#   - Maps count keys to shared DrOrtholog IDs
#   - Writes clean map, unmapped keys, ambiguous keys, and summary table
#
# What this script DOES NOT do:
#   - Does not convert DrOrtholog IDs to Entrez
#   - Does not run GO/KEGG/pathway analysis
#   - Does not assume DrOrtholog IDs are true ENSDARG gene IDs
#
# Downstream rule:
#   Use DrOrtholog only for overlap / shared-gene-set analyses.
#   Use separate, explicit annotation steps for GO/KEGG/Entrez.
# ================================================================


# ================================================================
# 0. Config
# ================================================================

base_dir    <- "C:/Users/augus/Desktop/UAF/Antarctic Fish Project/RNAseq"
emapper_dir <- file.path(base_dir, "emapper_zf_one2one")
counts_dir  <- base_dir
out_dir     <- file.path(base_dir, "mapped_counts")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Do not write legacy *_key_to_ENSDARG_minimal.tsv by default.
# We want downstream code to be forced to use DrOrtholog explicitly.
write_legacy_compatibility_files <- FALSE


# ================================================================
# 1. Packages
# ================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
})


# ================================================================
# 2. Helpers
# ================================================================

read_auto <- function(path) {
  delim <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  
  readr::read_delim(
    path,
    delim = delim,
    show_col_types = FALSE,
    progress = FALSE
  )
}

read_counts_keys <- function(path) {
  df <- readr::read_delim(
    path,
    delim = "\t",
    comment = "#",
    col_select = 1,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  key_header <- names(df)[1]
  
  out <- df %>%
    transmute(key = as.character(.data[[key_header]])) %>%
    filter(!is.na(key), nzchar(key)) %>%
    distinct()
  
  attr(out, "key_header") <- key_header
  
  out
}

norm_key <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_ace_aux <- function(df) {
  df <- df[, !grepl("^\\.\\.\\.", names(df)), drop = FALSE]
  names(df) <- tolower(names(df))
  
  if ("geneid" %in% names(df)) {
    names(df)[names(df) == "geneid"] <- "gene_id"
  }
  
  if ("genesymbol" %in% names(df)) {
    names(df)[names(df) == "genesymbol"] <- "gene_symbol"
  }
  
  if ("proteinaccession" %in% names(df)) {
    names(df)[names(df) == "proteinaccession"] <- "protein_accession"
  }
  
  first_col <- names(df)[1]
  
  df %>%
    filter(!(tolower(.data[[first_col]]) %in% c("gene_id", "geneid")))
}

parse_species_gene_field <- function(x) {
  x <- as.character(x)
  
  tibble(
    species_gene = x,
    
    gene_id = stringr::str_extract(
      x,
      "(?i)(?<=\\bgene[:=])[A-Za-z0-9._:-]+"
    ),
    
    gene_symbol = stringr::str_extract(
      x,
      "(?i)(?<=\\bgene_name[:=])[^:|\\s]+"
    )
  ) %>%
    mutate(
      gene_id = na_if(trimws(gene_id), ""),
      gene_symbol = na_if(trimws(gene_symbol), ""),
      gene_symbol = ifelse(
        gene_symbol %in% c("-", ".", "::"),
        NA_character_,
        gene_symbol
      )
    )
}


# ================================================================
# 3. Build candidate map for one species
# ================================================================

make_candidates <- function(species, ids_tsv, ace_aux_csv = NA_character_) {
  ids <- read_auto(ids_tsv)
  
  if (!all(c("zfish_key", "species_gene") %in% names(ids))) {
    stop(
      "[", species, "] ",
      basename(ids_tsv),
      " must contain columns: zfish_key and species_gene"
    )
  }
  
  parsed <- parse_species_gene_field(ids$species_gene) %>%
    bind_cols(DrOrtholog = as.character(ids$zfish_key)) %>%
    filter(!is.na(DrOrtholog), nzchar(DrOrtholog)) %>%
    distinct()
  
  # Route 1: species gene_id -> DrOrtholog
  by_gene_id <- parsed %>%
    filter(!is.na(gene_id), nzchar(gene_id)) %>%
    transmute(
      match_norm = norm_key(gene_id),
      matched_value = gene_id,
      DrOrtholog,
      picked_from = "gene_id",
      priority = 1L,
      species_gene,
      species_gene_id = gene_id,
      species_symbol = gene_symbol
    ) %>%
    distinct()
  
  # Route 2: species gene_name/symbol -> DrOrtholog
  by_symbol <- parsed %>%
    filter(!is.na(gene_symbol), nzchar(gene_symbol)) %>%
    transmute(
      match_norm = norm_key(gene_symbol),
      matched_value = gene_symbol,
      DrOrtholog,
      picked_from = "symbol",
      priority = 3L,
      species_gene,
      species_gene_id = gene_id,
      species_symbol = gene_symbol
    ) %>%
    distinct()
  
  # Route 3: ACE-specific symbol -> gene_id -> DrOrtholog bridge
  by_ace_aux <- tibble()
  
  if (
    species == "ACE" &&
    !is.na(ace_aux_csv) &&
    file.exists(ace_aux_csv)
  ) {
    aux <- read_auto(ace_aux_csv) %>%
      clean_ace_aux()
    
    if (all(c("gene_id", "gene_symbol") %in% names(aux))) {
      aux_clean <- aux %>%
        transmute(
          species_gene_id = as.character(gene_id),
          aux_symbol = as.character(gene_symbol)
        ) %>%
        filter(
          !is.na(species_gene_id), nzchar(species_gene_id),
          !is.na(aux_symbol), nzchar(aux_symbol),
          aux_symbol != "-"
        ) %>%
        distinct()
      
      by_ace_aux <- aux_clean %>%
        mutate(gene_id_norm = norm_key(species_gene_id)) %>%
        inner_join(
          by_gene_id %>%
            mutate(gene_id_norm = norm_key(species_gene_id)),
          by = "gene_id_norm"
        ) %>%
        transmute(
          match_norm = norm_key(aux_symbol),
          matched_value = aux_symbol,
          DrOrtholog,
          picked_from = "ace_aux_symbol_to_gene_id",
          priority = 2L,
          species_gene,
          species_gene_id = species_gene_id.x,
          species_symbol = aux_symbol
        ) %>%
        distinct()
    }
  }
  
  bind_rows(
    by_gene_id,
    by_ace_aux,
    by_symbol
  ) %>%
    distinct()
}


# ================================================================
# 4. Choose one ortholog per count key
# ================================================================

choose_one_mapping <- function(df) {
  if (!nrow(df)) {
    return(tibble())
  }
  
  priorities <- sort(unique(df$priority))
  
  for (p in priorities) {
    sub <- df %>%
      filter(priority == p)
    
    orthologs <- unique(sub$DrOrtholog[!is.na(sub$DrOrtholog) & nzchar(sub$DrOrtholog)])
    
    if (length(orthologs) == 1) {
      return(
        sub %>%
          filter(DrOrtholog == orthologs[1]) %>%
          slice(1)
      )
    }
    
    # If the best-priority route gives multiple orthologs, flag as ambiguous.
    # Do not fall through to lower-priority routes.
    if (length(orthologs) > 1) {
      return(tibble())
    }
  }
  
  tibble()
}


# ================================================================
# 5. Make map for one species
# ================================================================

make_species_map <- function(species, counts_path, ids_tsv, ace_aux_csv = NA_character_) {
  message("\n== ", species, " ==")
  
  if (!file.exists(counts_path)) {
    stop("[", species, "] Missing counts file: ", counts_path)
  }
  
  if (!file.exists(ids_tsv)) {
    stop("[", species, "] Missing ids file: ", ids_tsv)
  }
  
  keys <- read_counts_keys(counts_path)
  key_header <- attr(keys, "key_header")
  
  message(
    "[", species, "] key column = ", key_header,
    "; first keys = ", paste(head(keys$key, 3), collapse = ", ")
  )
  
  candidates <- make_candidates(
    species = species,
    ids_tsv = ids_tsv,
    ace_aux_csv = ace_aux_csv
  )
  
  joined <- keys %>%
    mutate(match_norm = norm_key(key)) %>%
    left_join(candidates, by = "match_norm") %>%
    filter(!is.na(DrOrtholog), nzchar(DrOrtholog))
  
  chosen <- bind_rows(
    lapply(
      split(joined, joined$key),
      choose_one_mapping
    )
  )
  
  if (!nrow(chosen)) {
    chosen <- joined[0, ]
  }
  
  ambiguous <- joined %>%
    group_by(key) %>%
    summarise(
      n_candidate_orthologs = n_distinct(DrOrtholog),
      candidate_orthologs = paste(sort(unique(DrOrtholog)), collapse = ";"),
      candidate_routes = paste(sort(unique(picked_from)), collapse = ";"),
      .groups = "drop"
    ) %>%
    filter(
      !key %in% chosen$key,
      n_candidate_orthologs > 1
    )
  
  unmapped <- keys %>%
    filter(!key %in% chosen$key) %>%
    left_join(
      ambiguous %>%
        select(key, n_candidate_orthologs),
      by = "key"
    ) %>%
    mutate(
      unmapped_reason = ifelse(
        is.na(n_candidate_orthologs),
        "no_unique_ortholog_mapping",
        "ambiguous_multiple_orthologs"
      )
    ) %>%
    select(key, unmapped_reason)
  
  full_map <- chosen %>%
    transmute(
      species,
      key,
      DrOrtholog,
      picked_from,
      matched_value,
      species_gene,
      species_gene_id,
      species_symbol
    ) %>%
    distinct() %>%
    arrange(key)
  
  minimal_map <- full_map %>%
    select(key, DrOrtholog, picked_from)
  
  readr::write_tsv(
    full_map,
    file.path(out_dir, paste0(species, "_key_to_DrOrtholog_map.tsv"))
  )
  
  readr::write_tsv(
    minimal_map,
    file.path(out_dir, paste0(species, "_key_to_DrOrtholog_minimal.tsv"))
  )
  
  readr::write_tsv(
    unmapped,
    file.path(out_dir, paste0(species, "_unmapped_keys.tsv"))
  )
  
  readr::write_tsv(
    ambiguous,
    file.path(out_dir, paste0(species, "_ambiguous_keys.tsv"))
  )
  
  if (write_legacy_compatibility_files) {
    legacy_map <- minimal_map %>%
      rename(ENSDARG = DrOrtholog)
    
    readr::write_tsv(
      legacy_map,
      file.path(out_dir, paste0(species, "_key_to_ENSDARG_minimal.tsv"))
    )
  }
  
  tibble(
    species = species,
    key_header = key_header,
    total_keys = nrow(keys),
    mapped_keys = nrow(full_map),
    unmapped_keys = nrow(unmapped),
    ambiguous_keys = nrow(ambiguous),
    unique_DrOrtholog = n_distinct(full_map$DrOrtholog),
    mapped_by_gene_id = sum(full_map$picked_from == "gene_id", na.rm = TRUE),
    mapped_by_symbol = sum(full_map$picked_from == "symbol", na.rm = TRUE),
    mapped_by_ace_aux = sum(full_map$picked_from == "ace_aux_symbol_to_gene_id", na.rm = TRUE)
  )
}


# ================================================================
# 6. Species configuration
# ================================================================

cfg <- tibble::tribble(
  ~species, ~counts, ~ids_tsv, ~ace_aux,
  
  "ACE",
  file.path(counts_dir, "ACE_counts_aggregated_by_GeneSymbol_corrected.txt"),
  file.path(emapper_dir, "ACE_to_DR.ids.tsv"),
  file.path(base_dir, "ACE_protein_gene_mapping.csv"),
  
  "COR",
  file.path(counts_dir, "COR_counts.txt"),
  file.path(emapper_dir, "COR_to_DR.ids.tsv"),
  NA_character_,
  
  "GEO",
  file.path(counts_dir, "GEO_counts.txt"),
  file.path(emapper_dir, "GEO_to_DR.ids.tsv"),
  NA_character_,
  
  "GIB",
  file.path(counts_dir, "GIB_counts.txt"),
  file.path(emapper_dir, "GIB_to_DR.ids.tsv"),
  NA_character_,
  
  "MAC",
  file.path(counts_dir, "MAC_counts.txt"),
  file.path(emapper_dir, "MAC_to_DR.ids.tsv"),
  NA_character_
)

missing_files <- c(
  cfg$counts[!file.exists(cfg$counts)],
  cfg$ids_tsv[!file.exists(cfg$ids_tsv)],
  cfg$ace_aux[!is.na(cfg$ace_aux) & !file.exists(cfg$ace_aux)]
)

if (length(missing_files)) {
  stop(
    "Missing required files:\n",
    paste0(" - ", missing_files, collapse = "\n")
  )
}


# ================================================================
# 7. Run
# ================================================================

summary_df <- purrr::pmap_dfr(
  cfg,
  function(species, counts, ids_tsv, ace_aux) {
    make_species_map(
      species = species,
      counts_path = counts,
      ids_tsv = ids_tsv,
      ace_aux_csv = ace_aux
    )
  }
)

readr::write_tsv(
  summary_df,
  file.path(out_dir, "mapping_counts_summary.tsv")
)

cat("\n=== Mapping summary ===\n")
print(summary_df, n = Inf)

cat(
  "\nDone.\n",
  "Maps written to: ", out_dir, "\n\n",
  "Important outputs:\n",
  "  <SPEC>_key_to_DrOrtholog_map.tsv\n",
  "  <SPEC>_key_to_DrOrtholog_minimal.tsv\n",
  "  <SPEC>_unmapped_keys.tsv\n",
  "  <SPEC>_ambiguous_keys.tsv\n",
  "  mapping_counts_summary.tsv\n\n",
  sep = ""
)