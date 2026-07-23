read_featurecounts_table <- function(path, rename_first_col = FALSE) {
  df <- read.delim(path, comment.char = "#")
  cn <- colnames(df)
  cn[7:length(cn)] <- sub("bam\\.([A-Z]+\\d+)_.*", "\\1", cn[7:length(cn)])
  if (rename_first_col) cn[1] <- "Geneid"
  colnames(df) <- cn
  df
}

featurecounts_to_matrix <- function(df) {
  out <- df[, -c(2:6), drop = FALSE]
  rownames(out) <- df$Geneid
  out$Geneid <- NULL
  out
}

align_metadata_to_counts <- function(id, counts_list) {
  id <- id %>%
    dplyr::mutate(Species = sub("([A-Z]+)\\d+", "\\1", SampleID)) %>%
    dplyr::relocate(Species, .after = SampleID)

  species_levels <- intersect(names(counts_list), unique(id$Species))
  aligned_blocks <- lapply(species_levels, function(sp) {
    block <- dplyr::filter(id, Species == sp)
    block[match(colnames(counts_list[[sp]]), block$SampleID), , drop = FALSE]
  })

  dplyr::bind_rows(aligned_blocks)
}

filter_counts_by_mean <- function(counts_list, min_mean = 20) {
  lapply(counts_list, function(x) {
    keep <- rowMeans(x) >= min_mean
    x[keep, , drop = FALSE]
  })
}

make_species_design <- function(species, keep_counts, id, reference = "Normoxia") {
  stopifnot(identical(colnames(keep_counts), dplyr::filter(id, Species == species)$SampleID))

  treatment <- dplyr::filter(id, Species == species) %>% dplyr::pull(Treatment)
  sex <- dplyr::filter(id, Species == species) %>% dplyr::pull(Sex)

  sample_table <- data.frame(
    treatment = factor(treatment),
    sex = factor(sex)
  )
  sample_table$treatment <- stats::relevel(sample_table$treatment, ref = reference)
  group <- factor(sample_table$treatment)

  design <- stats::model.matrix(~0 + treatment, data = sample_table)
  colnames(design) <- sub("^treatment", "", colnames(design))
  colnames(design)[colnames(design) == reference] <- "Control"

  y <- edgeR::DGEList(counts = keep_counts, group = group)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateGLMCommonDisp(y, design, verbose = TRUE)
  y <- edgeR::estimateGLMTrendedDisp(y, design)
  y <- edgeR::estimateGLMTagwiseDisp(y, design)

  list(
    treatment = treatment,
    sex = factor(sex),
    table = sample_table,
    group = group,
    design = design,
    y = y
  )
}

make_mds_plot <- function(y, group, sex, species, run_date,
                          manual_colors = NULL,
                          manual_shapes = NULL) {
  mds_out <- limma::plotMDS(y, plot = FALSE, top = nrow(y$counts))
  mds_df <- data.frame(
    Dim1 = mds_out$x,
    Dim2 = mds_out$y,
    SampleID = colnames(y),
    Treatment = group,
    Sex = sex
  )

  p <- ggplot2::ggplot(
    mds_df,
    ggplot2::aes(x = Dim1, y = Dim2, color = Treatment, shape = Treatment)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::stat_ellipse(ggplot2::aes(group = Treatment), level = 0.95, linetype = "solid", linewidth = 1) +
    ggplot2::labs(
      x = paste0("Dimension 1 \n (", round(mds_out$var.explained[1], 2), "% of leading logFC distance captured)"),
      y = paste0("Dimension 2 \n (", round(mds_out$var.explained[2], 2), "% of leading logFC distance captured)"),
      title = paste0(species, " MDS Plot with 95% Confidence Ellipses")
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "right")

  if (!is.null(manual_colors)) {
    p <- p + ggplot2::scale_color_manual(values = manual_colors)
  }
  if (!is.null(manual_shapes)) {
    p <- p + ggplot2::scale_shape_manual(values = manual_shapes)
  }

  out_dir <- file.path("edgeR_output", species)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(file = file.path(out_dir, paste0("MDS_plot_with_ellipses_", run_date, ".pdf")), width = 8, height = 6)
  print(p)
  grDevices::dev.off()
  print(p)

  list(mds_out = mds_out, mds_df = mds_df, p_mds = p)
}

plot_full_table_scatter <- function(full_table, species, run_date) {
  out_dir <- file.path("edgeR_output", species, "FDR=0.05_L2FC=0")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  myCol <- c(
    rgb(68/255, 1/255, 84/255, alpha = 0.5),
    rgb(31/255, 120/255, 180/255, alpha = 0.5),
    rgb(178/255, 223/255, 138/255, alpha = 0.7)
  )

  draw_panel <- function() {
    graphics::par(mfrow = c(1, 1), mar = c(5, 5, 4, 2) + 0.1)
    graphics::plot(
      full_table$logFC_Hypoxia,
      full_table$logFC_ILOS,
      pch = ".",
      col = "gray",
      main = paste(species, "Significant DEGs"),
      xlab = "logFC Hypoxia",
      ylab = "logFC ILOS",
      xlim = c(-6, 6),
      ylim = c(-6, 6),
      cex.lab = 2
    )
    graphics::abline(coef = c(0, 1))
    graphics::abline(v = 0, lty = 2)
    graphics::abline(h = 0, lty = 2)

    for (i in seq_len(nrow(full_table))) {
      tmp_row <- full_table[i, ]
      if (tmp_row$FDR_ILOS <= 0.05 & tmp_row$FDR_Hypoxia > 0.05) {
        graphics::points(tmp_row$logFC_Hypoxia, tmp_row$logFC_ILOS, pch = 15, col = myCol[1], cex = 0.5)
      }
      if (tmp_row$FDR_Hypoxia <= 0.05 & tmp_row$FDR_ILOS > 0.05) {
        graphics::points(tmp_row$logFC_Hypoxia, tmp_row$logFC_ILOS, pch = 17, col = myCol[2], cex = 0.5)
      }
      if (tmp_row$FDR_Hypoxia <= 0.05 & tmp_row$FDR_ILOS <= 0.05) {
        graphics::points(tmp_row$logFC_Hypoxia, tmp_row$logFC_ILOS, pch = 21, col = rgb(0, 0.5, 0, alpha = 0.8), cex = 0.8)
      }
    }

    graphics::legend(
      "bottomright",
      legend = c("Hypoxia/ILOS Intersection", "Hypoxia DE", "ILOS DE"),
      col = c(rgb(0, 0.5, 0, alpha = 0.8), myCol[2], myCol[1]),
      pch = c(21, 17, 15),
      cex = 0.8
    )
  }

  grDevices::pdf(
    file = file.path(out_dir, paste0("Figure_X_logFC_pairwiseAll_", species, "_", run_date, ".pdf")),
    height = 8,
    width = 16
  )
  draw_panel()
  grDevices::dev.off()
  draw_panel()
}

run_species_deg_pipeline <- function(species,
                                     keep_counts,
                                     id,
                                     run_date,
                                     contrast_names,
                                     manual_colors = NULL,
                                     manual_shapes = NULL) {
  setup <- make_species_design(species, keep_counts, id)
  mds <- make_mds_plot(
    y = setup$y,
    group = setup$group,
    sex = setup$sex,
    species = species,
    run_date = run_date,
    manual_colors = manual_colors,
    manual_shapes = manual_shapes
  )

  fit <- edgeR::glmQLFit(setup$y, setup$design)
  contrast_defs <- stats::setNames(paste0(contrast_names, "-Control"), contrast_names)
  my.contrasts <- do.call(limma::makeContrasts, c(as.list(contrast_defs), list(levels = setup$design)))

  out <- list(
    mds[["p_mds"]],
    fit,
    my.contrasts
  )
  names(out) <- c(
    paste0("p_mds_", species),
    paste0("fit_", species),
    paste0("my.contrasts_", species)
  )

  out[[paste0("treatment_", species)]] <- setup$treatment
  out[[paste0("sex_", species)]] <- setup$sex
  out[[paste0("table_", species)]] <- setup$table
  out[[paste0("group_", species)]] <- setup$group
  out[[paste0("design_", species)]] <- setup$design
  out[[paste0("y_", species)]] <- setup$y

  out_dir <- file.path("edgeR_output", species, "FDR=0.05_L2FC=0")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (cond in contrast_names) {
    qlf <- edgeR::glmQLFTest(fit, contrast = my.contrasts[, cond])
    de <- limma::decideTests(qlf, adjust.method = "fdr", p.value = 0.05, lfc = 0)
    detags <- rownames(setup$y)[as.logical(de)]

    edgeR::plotSmear(qlf, de.tags = detags, main = paste(cond, "vs. Control"))
    graphics::abline(h = c(-1, 1), col = "blue")

    grDevices::pdf(file.path(out_dir, paste0("smearPlot_", cond, "_", species, "_de_", run_date, ".pdf")), w = 7, h = 6)
    edgeR::plotSmear(qlf, de.tags = detags, main = paste(cond, "vs. Control"))
    graphics::abline(h = c(-1, 1), col = "blue")
    grDevices::dev.off()

    tbl <- data.frame(edgeR::topTags(qlf, n = nrow(keep_counts)))
    sig_p <- subset(tbl, PValue <= 0.05)
    sig_fdr <- subset(tbl, FDR <= 0.05)
    sig_fdr_l2fc <- dplyr::filter(sig_fdr, logFC < 0 | logFC > 0)
    sig_fdr_up <- rownames(sig_fdr[sig_fdr$logFC > 0, , drop = FALSE])
    sig_fdr_down <- rownames(sig_fdr[sig_fdr$logFC < 0, , drop = FALSE])

    out[[paste0("qlf.", cond, "_", species)]] <- qlf
    out[[paste0("de", match(cond, contrast_names), "_", species)]] <- de
    out[[paste0("detags", match(cond, contrast_names), "_", species)]] <- detags
    out[[paste0(cond, "_table_", species)]] <- tbl
    out[[paste0("sig.", cond, "_", species)]] <- sig_p
    out[[paste0("sig.", cond, "_", species, ".fdr")]] <- sig_fdr
    out[[paste0("sig.", cond, "_", species, ".fdr.l2fc")]] <- sig_fdr_l2fc
    out[[paste0("sig.", cond, "_", species, ".fdr_up")]] <- sig_fdr_up
    out[[paste0("sig.", cond, "_", species, ".fdr_down")]] <- sig_fdr_down

    utils::write.csv(sig_fdr, file = file.path(out_dir, paste0(species, "_", cond, "_DEGs_", run_date, ".csv")))
  }

  if (length(contrast_names) == 2 && all(c("Hypoxia", "ILOS") %in% contrast_names)) {
    hypoxia_df <- out[[paste0("sig.Hypoxia_", species, ".fdr")]]
    ilos_df <- out[[paste0("sig.ILOS_", species, ".fdr")]]

    listInput <- list(
      stats::setNames(list(out[[paste0("sig.Hypoxia_", species, ".fdr_up")]]), paste0("Hypoxia_", species, "_up"))[[1]],
      stats::setNames(list(out[[paste0("sig.Hypoxia_", species, ".fdr_down")]]), paste0("Hypoxia_", species, "_down"))[[1]],
      stats::setNames(list(out[[paste0("sig.ILOS_", species, ".fdr_up")]]), paste0("ILOS_", species, "_up"))[[1]],
      stats::setNames(list(out[[paste0("sig.ILOS_", species, ".fdr_down")]]), paste0("ILOS_", species, "_down"))[[1]]
    )
    names(listInput) <- c(
      paste0("Hypoxia_", species, "_up"),
      paste0("Hypoxia_", species, "_down"),
      paste0("ILOS_", species, "_up"),
      paste0("ILOS_", species, "_down")
    )
    out[[paste0("listInput_", species)]] <- listInput

    UpSetR::upset(
      UpSetR::fromList(listInput),
      order.by = "freq",
      mainbar.y.label = "DE Genes",
      sets.x.label = "DE Genes Per Stage Comparison",
      sets.bar.color = "#56B4E9"
    )
    grDevices::pdf(file.path(out_dir, paste0("Intersection_of_DE_Genes_", species, "_", run_date, ".pdf")), w = 8, h = 6)
    UpSetR::upset(
      UpSetR::fromList(listInput),
      order.by = "freq",
      mainbar.y.label = "DE Genes",
      sets.x.label = "DE Genes Per Stage Comparison",
      sets.bar.color = "#56B4E9"
    )
    grDevices::dev.off()

    venn_grob <- VennDiagram::venn.diagram(
      x = listInput,
      category.names = c("Hypoxia \n Up", "Hypoxia \n Down", "ILOS Up", "ILOS Down"),
      filename = NULL,
      lwd = 0,
      lty = "blank",
      fill = c("#9b59b6", "#CD69C9", "#1f78b4", "#56B4E9"),
      cex = 1.5,
      main = paste(species, "Significant DEGs"),
      main.pos = c(0.5, 1.05),
      main.fontfamily = "sans",
      main.fontface = "bold",
      main.cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontface = "plain",
      cat.default.pos = "outer",
      cat.fontfamily = "sans",
      cat.dist = c(0.2, 0.2, 0.1, 0.1),
      cat.pos = c(-10, 15, 0, 0)
    )
    grid::grid.newpage()
    grid::grid.draw(venn_grob)
    VennDiagram::venn.diagram(
      x = listInput,
      category.names = c("Hypoxia \n Up", "Hypoxia \n Down", "ILOS Up", "ILOS Down"),
      filename = file.path(out_dir, paste0("venn_diagramm_of_DE_genes_", species, "_", run_date, ".png")),
      output = TRUE,
      imagetype = "png",
      height = 1000,
      width = 1600,
      resolution = 300,
      compression = "lzw",
      lwd = 2,
      lty = "blank",
      fill = c("#9b59b6", "#CD69C9", "#1f78b4", "#56B4E9"),
      cex = 1.5,
      main = paste(species, "Significant DEGs"),
      main.pos = c(0.5, 1.05),
      main.fontfamily = "sans",
      main.fontface = "bold",
      main.cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontface = "plain",
      cat.default.pos = "outer",
      cat.fontfamily = "sans",
      cat.dist = c(0.2, 0.2, 0.1, 0.1),
      cat.pos = c(-10, 15, 0, 0)
    )

    hypoxia_df$Gene <- rownames(hypoxia_df)
    ilos_df$Gene <- rownames(ilos_df)
    sig_intersect <- data.frame(Gene = intersect(rownames(hypoxia_df), rownames(ilos_df)))
    tmp1 <- merge(sig_intersect, hypoxia_df[, c("Gene", "logFC")], by = "Gene")
    colnames(tmp1) <- c("Gene", "logFC.Hypoxia")
    tmp2 <- merge(tmp1, ilos_df[, c("Gene", "logFC")], by = "Gene")
    colnames(tmp2) <- c("Gene", "logFC.Hypoxia", "logFC.ILOS")
    intersect_FC <- tmp2

    out[[paste0("sig.intersect.fdr_", species)]] <- sig_intersect
    out[[paste0("intersect_FC_", species)]] <- intersect_FC

    utils::write.table(
      intersect_FC,
      file = file.path(out_dir, paste0("intersectGenes_FC_", species, "_", run_date, ".txt")),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )

    intersect_FC_plot <- intersect_FC[stats::complete.cases(intersect_FC[, c("logFC.Hypoxia", "logFC.ILOS")]), , drop = FALSE]

    grDevices::pdf(file = file.path(out_dir, paste0("intersectGenes_FC_Hypoxia-vs-ILOS_", run_date, ".pdf")))
    if (nrow(intersect_FC_plot) > 0) {
      graphics::plot(intersect_FC_plot$logFC.Hypoxia,
                     intersect_FC_plot$logFC.ILOS,
                     main = paste("Hypoxia vs. ILOS (", nrow(intersect_FC_plot), ")", sep = ""),
                     xlab = "Hypoxia Log2FC",
                     ylab = "ILOS Log2FC")
      graphics::abline(coef = c(0, 1))
      graphics::abline(v = 0, lty = 2)
      graphics::abline(h = 0, lty = 2)
    } else {
      graphics::plot.new()
      graphics::title(main = paste("Hypoxia vs. ILOS (", species, ")", sep = ""))
      graphics::text(0.5, 0.5, "No overlapping DEGs with finite logFC values to plot")
    }
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out_dir, paste0("intersectGenes_boxplot_Hypoxia-vs-ILOS_", species, "_", run_date, ".pdf")))
    if (nrow(intersect_FC_plot) > 0) {
      boxplot_color <- rgb(0.5, 0.5, 1, alpha = 0.5)
      yvals <- unlist(intersect_FC_plot[, 2:3], use.names = FALSE)
      yvals <- yvals[is.finite(yvals)]
      if (!length(yvals)) yvals <- c(-1, 1)
      graphics::plot(NULL,
                     xlim = c(0.5, 2.5),
                     ylim = range(yvals),
                     ylab = "logFC relative to control",
                     xlab = "",
                     xaxt = "n")
      graphics::points(jitter(rep(1:2, each = nrow(intersect_FC_plot))),
                       c(intersect_FC_plot[, 2], intersect_FC_plot[, 3]),
                       col = "black",
                       pch = 21)
      graphics::boxplot(intersect_FC_plot[, 2:3], add = TRUE, col = boxplot_color, border = "darkblue", outline = FALSE, xaxt = "n")
      graphics::axis(1, at = 1:2, labels = colnames(intersect_FC_plot[, 2:3]))
    } else {
      graphics::plot.new()
      graphics::text(0.5, 0.5, "No overlapping DEGs with finite logFC values to summarize")
    }
    grDevices::dev.off()

    grDevices::pdf(file = file.path(out_dir, paste0("allGenes_FC_Hypoxia-vs-ILOS_", species, "_", run_date, ".pdf")))
    graphics::plot(out[[paste0("qlf.Hypoxia_", species)]]$table$logFC,
                   out[[paste0("qlf.ILOS_", species)]]$table$logFC,
                   main = paste("Hypoxia vs. ILOS (n=", nrow(out[[paste0("qlf.Hypoxia_", species)]]$table), ")", sep = ""),
                   xlab = "Hypoxia Log2FC",
                   ylab = "ILOS Log2FC")
    graphics::abline(coef = c(0, 1))
    graphics::abline(v = 0, lty = 2)
    graphics::abline(h = 0, lty = 2)
    grDevices::dev.off()

    hypoxia_table <- out[[paste0("Hypoxia_table_", species)]]
    ilos_table <- out[[paste0("ILOS_table_", species)]]
    hypoxia_table$Gene <- rownames(hypoxia_table)
    ilos_table$Gene <- rownames(ilos_table)
    full_table <- merge(hypoxia_table, ilos_table, by = "Gene", suffixes = c("_Hypoxia", "_ILOS"))
    full_table <- full_table[, c(1, 2, 6, 7, 11)]
    out[[paste0("Full_table_", species)]] <- full_table
    plot_full_table_scatter(full_table, species, run_date)

    out[[paste0("sig.intersect_Hypoxia_ILOS_", species)]] <- intersect(hypoxia_df$Gene, ilos_df$Gene)
    out[[paste0("only_Hypoxia_", species)]] <- setdiff(hypoxia_df$Gene, ilos_df$Gene)
    out[[paste0("only_ILOS_", species)]] <- setdiff(ilos_df$Gene, hypoxia_df$Gene)
  }

  out
}


flatten_list_columns <- function(df) {
  if (is.null(df)) return(df)
  as.data.frame(lapply(df, function(x) {
    if (is.list(x)) vapply(x, toString, character(1)) else x
  }), stringsAsFactors = FALSE)
}

map_symbol_keys_to_gprofiler_ids <- function(keys, label = "query") {
  keys <- clean_set(keys)
  
  if (!length(keys)) {
    return(tibble(
      input_key = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      ENSDARG = character()
    ))
  }
  
  candidate_symbols <- clean_set(c(
    keys,
    tolower(keys),
    toupper(keys)
  ))
  
  valid_symbols <- AnnotationDbi::keys(org.Dr.eg.db, keytype = "SYMBOL")
  query_symbols <- intersect(candidate_symbols, valid_symbols)
  
  if (!length(query_symbols)) {
    message("[", label, "] no valid zebrafish symbols found.")
    return(tibble(
      input_key = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      ENSDARG = character()
    ))
  }
  
  sym_lookup <- suppressMessages(
    AnnotationDbi::select(
      org.Dr.eg.db,
      keys = query_symbols,
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENTREZID", "ENSEMBL")
    )
  ) |>
    transmute(
      SYMBOL = as.character(SYMBOL),
      ENTREZID = as.character(ENTREZID),
      ENSDARG = as.character(ENSEMBL),
      symbol_norm = tolower(SYMBOL)
    ) |>
    filter(
      !is.na(ENSDARG), nzchar(ENSDARG)
    ) |>
    distinct(symbol_norm, ENSDARG, .keep_all = TRUE)
  
  out <- tibble(
    input_key = keys,
    symbol_norm = tolower(keys)
  ) |>
    left_join(sym_lookup, by = "symbol_norm", relationship = "many-to-many") |>
    filter(!is.na(ENSDARG), nzchar(ENSDARG)) |>
    distinct(input_key, ENSDARG, .keep_all = TRUE) |>
    select(input_key, SYMBOL, ENTREZID, ENSDARG)
  
  message(
    "[", label, "] ",
    length(keys), " input keys -> ",
    dplyr::n_distinct(out$ENSDARG), " zebrafish Ensembl gene IDs"
  )
  
  out
}

get_bp_terms_with_fallback <- function(query_entrez, bg_entrez,
                                       fdr = 0.10, topN = NA,
                                       sim_method = c("Wang", "Rel"),
                                       significant = TRUE) {
  sim_method <- match.arg(sim_method)
  query_entrez <- unique(as.character(query_entrez))
  bg_entrez <- unique(c(as.character(bg_entrez), query_entrez))

  if (!length(query_entrez) || !length(bg_entrez)) {
    empty_bp <- data.frame(term_id = character(), stringsAsFactors = FALSE)
    return(list(bp = empty_bp, scores = numeric(), source = "empty_input"))
  }

  gp <- try(gprofiler2::gost(
    query = query_entrez,
    organism = "drerio",
    sources = "GO",
    correction_method = "fdr",
    significant = significant,
    user_threshold = fdr,
    custom_bg = bg_entrez,
    domain_scope = "custom"
  ), silent = TRUE)

  if (inherits(gp, "try-error") || is.null(gp$result) || !nrow(gp$result)) {
    alt_sig <- !significant
    gp <- try(gprofiler2::gost(
      query = query_entrez,
      organism = "drerio",
      sources = "GO",
      correction_method = "fdr",
      significant = alt_sig,
      user_threshold = fdr,
      custom_bg = bg_entrez,
      domain_scope = "custom"
    ), silent = TRUE)
  }

  if (!inherits(gp, "try-error") && !is.null(gp$result) && nrow(gp$result)) {
    bp <- dplyr::filter(gp$result, source == "GO:BP") |>
      dplyr::distinct(term_id, .keep_all = TRUE) |>
      dplyr::arrange(p_value)

    if (!is.na(topN)) {
      bp <- dplyr::slice_head(bp, n = min(topN, nrow(bp)))
    }

    scores <- setNames(-log10(pmax(bp$p_value, .Machine$double.eps)), bp$term_id)
    return(list(bp = bp, scores = scores, source = "gprofiler"))
  }

  ann <- AnnotationDbi::select(
    org.Dr.eg.db::org.Dr.eg.db,
    keys = unique(query_entrez),
    keytype = "ENTREZID",
    columns = c("GO", "ONTOLOGY")
  ) |>
    dplyr::filter(ONTOLOGY == "BP", !is.na(GO)) |>
    dplyr::distinct(ENTREZID, GO)

  bp <- ann |>
    dplyr::count(GO, name = "k") |>
    dplyr::arrange(dplyr::desc(k))
  names(bp)[1] <- "term_id"

  if (!is.na(topN)) {
    bp <- dplyr::slice_head(bp, n = min(topN, nrow(bp)))
  }

  scores <- setNames(as.numeric(bp$k), bp$term_id)
  list(bp = bp, scores = scores, source = "offline_coverage")
}

make_rrvgo_treemap <- function(bp_tbl, scores, title_txt,
                               sim_method = c("Wang", "Rel"),
                               threshold = 0.7,
                               out_stub = "rrvgo",
                               return = c("list", "plot", "grob"),
                               write_tsv_file = TRUE,
                               save_plot_rds = FALSE,
                               out_dir = ".") {
  sim_method <- match.arg(sim_method)
  return <- match.arg(return)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(bp_tbl) || !nrow(bp_tbl)) {
    return(NULL)
  }

  term_ids <- bp_tbl$term_id
  sim <- rrvgo::calculateSimMatrix(
    term_ids,
    orgdb = "org.Dr.eg.db",
    ont = "BP",
    method = sim_method
  )
  red <- rrvgo::reduceSimMatrix(
    sim,
    scores,
    threshold = threshold,
    orgdb = "org.Dr.eg.db"
  )

  rt <- as.data.frame(red)
  rt$cluster <- as.factor(rt$cluster)
  p <- rrvgo::treemapPlot(reducedTerms = rt, size = "score", title = title_txt)

  if (write_tsv_file) {
    readr::write_tsv(rt, file.path(out_dir, paste0(out_stub, "_reduced.tsv")))
  }
  if (save_plot_rds) {
    saveRDS(p, file.path(out_dir, paste0(out_stub, "_plot.rds")))
  }

  if (return == "plot") return(p)
  if (return == "grob") return(ggplotify::as.grob(p))

  list(plot = p, grob = ggplotify::as.grob(p), data = rt, sim = sim, reduced = red)
}

only_ensdarg <- function(x) {
  unique(grep("^ENSDARG", as.character(x), value = TRUE))
}

ensdarg_to_symbol <- function(overlap_ensdarg,
                              map_df,
                              symbol_col = "Symbol",
                              key_col = "key") {
  overlap_ensdarg <- only_ensdarg(overlap_ensdarg)
  if (!length(overlap_ensdarg)) return(character())

  use_symbol_col <- symbol_col %in% names(map_df) &&
    any(!is.na(map_df[[symbol_col]]) & nzchar(as.character(map_df[[symbol_col]])))

  hits <- map_df[map_df$ENSDARG %in% overlap_ensdarg, , drop = FALSE]
  vals <- if (use_symbol_col) hits[[symbol_col]] else hits[[key_col]]
  vals <- as.character(vals)
  unique(vals[!is.na(vals) & nzchar(vals)])
}

run_overlap_go_bp_analysis <- function(overlap_ensdarg,
                                       map_df,
                                       background_symbols,
                                       title_txt,
                                       out_stub,
                                       out_dir = file.path("edgeR_output", "overlap_GO"),
                                       topN = 100,
                                       fdr = 0.10,
                                       sim_method = "Rel",
                                       significant = FALSE) {
  overlap_ensdarg <- only_ensdarg(overlap_ensdarg)
  overlap_symbols <- ensdarg_to_symbol(overlap_ensdarg, map_df)

  query_map <- map_symbol_keys_to_gprofiler_ids(
    mac_key_table$MAC_key,
    label = paste0(comparison, " overlap MAC keys")
  )
  
  background_map <- map_symbol_keys_to_gprofiler_ids(
    background_keys,
    label = paste0(comparison, " MAC DrOrtholog-tested background")
  )
  
  query_gprof <- clean_set(query_map$ENSDARG)
  bg_gprof <- clean_set(background_map$ENSDARG)
  
  # Keep query in the universe because query genes came from tested genes.
  bg_gprof <- clean_set(c(bg_gprof, query_gprof))

  plot <- NULL
  if (!is.null(bp_res$bp) && nrow(bp_res$bp)) {
    plot <- make_rrvgo_treemap(
      bp_tbl = bp_res$bp,
      scores = bp_res$scores,
      title_txt = title_txt,
      sim_method = sim_method,
      out_stub = out_stub,
      return = "plot",
      out_dir = out_dir
    )
  }

  list(
    overlap_ensdarg = overlap_ensdarg,
    overlap_symbols = overlap_symbols,
    query_map = query_map,
    bg_map = bg_map,
    query_entrez = query_entrez,
    bg_entrez = bg_entrez,
    bp_res = bp_res,
    plot = plot
  )
}

build_go_bubble_plot <- function(top_terms,
                                 title_txt,
                                 palette = "vik",
                                 base_size = 18) {
  if (is.null(top_terms) || !nrow(top_terms)) return(NULL)

  ggplot2::ggplot(
    top_terms,
    ggplot2::aes(x = percent_DEGs, y = reorder(term_name, percent_DEGs))
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = intersection_size, fill = -log10(p_value)),
      shape = 21,
      color = "black",
      stroke = 0.6,
      alpha = 0.9
    ) +
    scico::scale_fill_scico(palette = palette, name = expression(-log[10](p))) +
    ggplot2::scale_size_area(max_size = 12, name = "DEG count") +
    ggplot2::labs(
      x = "% of DEGs in GO term",
      y = "GO term",
      title = title_txt
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "bold", size = 14, lineheight = 0.95),
      axis.text.x = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(face = "bold", size = 18),
      plot.title = ggplot2::element_text(face = "bold", size = 20, margin = ggplot2::margin(b = 8)),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 12),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ================================================================
# GO enrichment helper using mapped Entrez IDs
# ================================================================

run_species_go_enrichment_mapped <- function(species,
                                             condition,
                                             sig_mapped_df,
                                             background_mapped_df,
                                             run_date,
                                             out_root = "edgeR_output",
                                             correction_method = "g_SCS",
                                             user_threshold = 0.05,
                                             sources = c("GO:BP", "GO:MF", "GO:CC"),
                                             palette = "vik",
                                             top_n = 20,
                                             min_query_genes = 5,
                                             return_all_terms = TRUE) {
  
  query_genes <- clean_set(sig_mapped_df$ENTREZID)
  background_genes <- clean_set(background_mapped_df$ENTREZID)
  
  out_dir <- file.path(out_root, species, "GO_enrichment")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  empty <- list(
    species = species,
    condition = condition,
    query_genes = query_genes,
    background_genes = background_genes,
    gost = NULL,
    result_table = data.frame(),
    top_terms = data.frame(),
    plot = NULL,
    out_dir = out_dir,
    status = "empty_or_failed"
  )
  
  if (length(query_genes) < min_query_genes) {
    message(
      species, " ", condition,
      ": only ", length(query_genes),
      " mapped query genes; skipping enrichment."
    )
    empty$status <- "too_few_query_genes"
    return(empty)
  }
  
  if (!length(background_genes)) {
    message(species, " ", condition, ": no mapped background genes; skipping enrichment.")
    empty$status <- "no_background"
    return(empty)
  }
  
  missing_from_bg <- setdiff(query_genes, background_genes)
  
  if (length(missing_from_bg)) {
    message(
      species, " ", condition,
      ": ", length(missing_from_bg),
      " query Entrez IDs were absent from the mapped background. ",
      "Adding them to background because query genes came from tested genes."
    )
    
    background_genes <- unique(c(background_genes, query_genes))
  }
  
  # Important:
  # significant = FALSE returns the full table and lets us count significant terms afterward.
  # This prevents g:Profiler's "No results to show" message from hiding useful diagnostics.
  gostres <- try(
    gprofiler2::gost(
      query = query_genes,
      organism = "drerio",
      sources = sources,
      correction_method = correction_method,
      significant = !return_all_terms,
      user_threshold = user_threshold,
      custom_bg = background_genes,
      domain_scope = "custom"
    ),
    silent = TRUE
  )
  
  if (inherits(gostres, "try-error") || is.null(gostres$result)) {
    message(species, " ", condition, ": g:Profiler returned no result table.")
    return(empty)
  }
  
  res_df <- as.data.frame(gostres$result)
  
  if (!nrow(res_df)) {
    message(species, " ", condition, ": g:Profiler result table had 0 rows.")
    empty$gost <- gostres
    empty$status <- "no_terms_returned"
    return(empty)
  }
  
  res_df$percent_DEGs <- 100 * res_df$intersection_size / res_df$query_size
  res_df$GeneRatio <- res_df$intersection_size / res_df$term_size
  
  top_terms <- res_df
  
  if ("significant" %in% names(top_terms)) {
    sig_terms <- top_terms[top_terms$significant %in% TRUE, , drop = FALSE]
    
    if (nrow(sig_terms)) {
      top_terms <- sig_terms
    }
  }
  
  top_terms <- top_terms[order(-top_terms$percent_DEGs, top_terms$p_value), , drop = FALSE]
  top_terms <- utils::head(top_terms, top_n)
  
  plot <- build_go_bubble_plot(
    top_terms = top_terms,
    title_txt = paste(condition, "GO Enrichment (", species, ")", sep = ""),
    palette = palette
  )
  
  res_out <- flatten_list_columns(res_df)
  
  utils::write.csv(
    res_out,
    file = file.path(out_dir, paste0(species, "_", condition, "_full_GO_enrichment_", run_date, ".csv")),
    row.names = FALSE
  )
  
  if (!is.null(plot)) {
    grDevices::pdf(
      file = file.path(out_dir, paste0(condition, "_", species, "_GO_bubbleplot_", run_date, ".pdf")),
      width = 10,
      height = 8
    )
    print(plot)
    grDevices::dev.off()
  }
  
  list(
    species = species,
    condition = condition,
    query_genes = query_genes,
    background_genes = background_genes,
    gost = gostres,
    result_table = res_df,
    top_terms = top_terms,
    plot = plot,
    out_dir = out_dir,
    status = "ok"
  )
}