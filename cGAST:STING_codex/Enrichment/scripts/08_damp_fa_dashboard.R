run_damp_fa_dashboard_module <- function(
  tbl,
  padj_cutoff,
  logFC_cutoff,
  ego,
  ego_MF,
  ekegg,
  ereact,
  gsea_go,
  gsea_go_MF,
  gsea_kegg,
  damp_pattern,
  lipid_pattern,
  sterile_pattern,
  dirs,
  save_plot,
  register_if_exists
) {
  damp_fa_panel <- if (exists("damp_fa_panel", inherits = TRUE)) {
    get("damp_fa_panel", inherits = TRUE)
  } else {
    list(
      DAMP_alarmins = c("HMGB1", "S100A8", "S100A9", "HSPA1A", "HSP90AA1", "IL33"),
      PRR_sensing = c("TLR2", "TLR4", "NOD1", "NOD2", "NLRP3", "AIM2", "TMEM173", "MB21D1"),
      Inflammasome_IL1 = c("CASP1", "IL1B", "PYCARD", "NLRC4", "GSDMD"),
      Lipid_FA = c("PPARA", "PPARG", "CPT1A", "ACSL1", "ACOX1", "FASN", "SCD", "ELOVL6"),
      Eicosanoid = c("PLA2G4A", "PTGS1", "PTGS2", "ALOX5", "ALOX12", "ALOX15")
    )
  }

  if (!"SYMBOL" %in% names(tbl)) tbl$SYMBOL <- NA_character_
  if (!"gene_label" %in% names(tbl)) tbl$gene_label <- tbl$gene

  classify_reg <- function(logfc, padj) {
    dplyr::case_when(
      !is.na(logfc) & !is.na(padj) & logfc >= logFC_cutoff & padj < padj_cutoff ~ "Up",
      !is.na(logfc) & !is.na(padj) & logfc <= -logFC_cutoff & padj < padj_cutoff ~ "Down",
      TRUE ~ "NoChange"
    )
  }

  panel_rows <- lapply(names(damp_fa_panel), function(bucket) {
    genes <- unique(as.character(damp_fa_panel[[bucket]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    lapply(genes, function(g) {
      g_upper <- toupper(g)
      idx <- which(
        toupper(as.character(tbl$gene_label)) == g_upper |
        toupper(as.character(tbl$SYMBOL)) == g_upper |
        toupper(as.character(tbl$gene)) == g_upper
      )
      if (length(idx) == 0) {
        data.frame(
          bucket = bucket,
          gene_label = g,
          ENTREZID = NA_character_,
          logFC = NA_real_,
          padj = NA_real_,
          Regulation = "NoChange",
          Present_in_input = FALSE,
          stringsAsFactors = FALSE
        )
      } else {
        hit <- tbl[idx, , drop = FALSE]
        data.frame(
          bucket = bucket,
          gene_label = as.character(hit$gene_label),
          ENTREZID = as.character(hit$ENTREZID),
          logFC = as.numeric(hit$logFC),
          padj = as.numeric(hit$padj),
          Regulation = classify_reg(as.numeric(hit$logFC), as.numeric(hit$padj)),
          Present_in_input = TRUE,
          stringsAsFactors = FALSE
        )
      }
    })
  })

  panel_hits <- dplyr::bind_rows(unlist(panel_rows, recursive = FALSE))
  panel_hits <- panel_hits %>% dplyr::distinct(bucket, gene_label, ENTREZID, logFC, padj, Regulation, Present_in_input, .keep_all = TRUE)

  panel_csv <- file.path(dirs$tables_sig, "damp_fa_gene_panel_hits.csv")
  panel_xlsx <- file.path(dirs$tables_sig, "damp_fa_gene_panel_hits.xlsx")
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(panel_hits, panel_csv)
  } else {
    utils::write.csv(panel_hits, panel_csv, row.names = FALSE)
  }
  openxlsx::write.xlsx(panel_hits, panel_xlsx, overwrite = TRUE)
  register_if_exists("table", "damp_fa_gene_panel_hits_csv", panel_csv, "Curated DAMP/FA panel hits.")
  register_if_exists("table", "damp_fa_gene_panel_hits_xlsx", panel_xlsx, "Curated DAMP/FA panel hits.")

  panel_labels <- unique(panel_hits$gene_label[panel_hits$Present_in_input %in% TRUE])
  panel_labels <- panel_labels[!is.na(panel_labels) & nzchar(panel_labels)]
  if (length(panel_labels) > 0) {
    volc_panel <- EnhancedVolcano::EnhancedVolcano(
      tbl,
      lab = tbl$gene_label,
      x = "logFC",
      y = "padj",
      pCutoff = padj_cutoff,
      FCcutoff = logFC_cutoff,
      selectLab = panel_labels,
      title = "Volcano: DAMP/FA panel genes"
    )
    save_plot(volc_panel, file.path(dirs$dashboard, "volcano_damp_fa_panel.pdf"), 9, 8)
  }

  panel_hits_plot <- panel_hits %>%
    dplyr::filter(Present_in_input %in% TRUE, !is.na(logFC), !is.na(gene_label), nzchar(gene_label)) %>%
    dplyr::arrange(logFC)

  if (nrow(panel_hits_plot) > 0) {
    panel_hits_plot <- panel_hits_plot %>%
      dplyr::mutate(gene_bucket = paste(gene_label, bucket, sep = " | "))
    p_bar <- ggplot2::ggplot(panel_hits_plot, ggplot2::aes(x = stats::reorder(gene_bucket, logFC), y = logFC, fill = bucket)) +
      ggplot2::geom_col(width = 0.8) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "DAMP/FA panel genes in input",
        x = "Gene (bucket)",
        y = "log2 fold change"
      ) +
      ggplot2::theme_bw(base_size = 10)
    save_plot(p_bar, file.path(dirs$dashboard, "damp_fa_panel_logFC_bar.pdf"), 9, 7)
  }

  as_df_safe <- function(x) {
    if (is.null(x)) return(data.frame())
    out <- tryCatch(as.data.frame(x), error = function(e) data.frame())
    if (is.null(out)) return(data.frame())
    out
  }

  clean_subset <- function(df, source_name, tag_name) {
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(source = source_name, subset = tag_name, Note = "No results", stringsAsFactors = FALSE))
    }
    dplyr::bind_cols(data.frame(source = source_name, subset = tag_name, stringsAsFactors = FALSE), as.data.frame(df))
  }

  subset_plot <- function(df, source_name, tag_name) {
    if (is.null(df) || nrow(df) < 3 || !"Description" %in% names(df)) return(invisible(NULL))

    d <- as.data.frame(df)
    d <- d[seq_len(min(nrow(d), 20)), , drop = FALSE]
    d$Description <- stringr::str_wrap(d$Description, width = 50)

    if ("GeneRatio" %in% names(d)) {
      ratio_num <- suppressWarnings(vapply(strsplit(as.character(d$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1)))
      d$xval <- ratio_num
      xlab <- "Gene ratio"
      size_col <- if ("Count" %in% names(d)) "Count" else NULL
    } else if ("NES" %in% names(d)) {
      d$xval <- d$NES
      xlab <- "NES"
      size_col <- if ("setSize" %in% names(d)) "setSize" else NULL
    } else if ("Count" %in% names(d)) {
      d$xval <- d$Count
      xlab <- "Count"
      size_col <- "Count"
    } else {
      d$xval <- seq_len(nrow(d))
      xlab <- "Rank"
      size_col <- NULL
    }

    color_col <- if ("p.adjust" %in% names(d)) "p.adjust" else if ("padj" %in% names(d)) "padj" else NULL
    p <- ggplot2::ggplot(d, ggplot2::aes(x = xval, y = stats::reorder(Description, xval)))
    if (!is.null(size_col) && !is.null(color_col)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]], color = .data[[color_col]]))
    } else if (!is.null(size_col)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(size = .data[[size_col]]), color = "steelblue")
    } else if (!is.null(color_col)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color_col]]), size = 3)
    } else {
      p <- p + ggplot2::geom_point(size = 3, color = "steelblue")
    }

    p <- p +
      ggplot2::labs(title = paste0("Term subset: ", source_name, " (", tag_name, ")"), x = xlab, y = NULL) +
      ggplot2::theme_bw(base_size = 10)

    out_pdf <- file.path(dirs$dashboard, paste0("term_subset_", tolower(gsub("[^A-Za-z0-9]+", "_", source_name)), "_", tolower(tag_name), ".pdf"))
    save_plot(p, out_pdf, 9, 6)
  }

  source_tables <- list(
    ORA_GO_BP = as_df_safe(ego),
    ORA_GO_MF = as_df_safe(ego_MF),
    ORA_KEGG = as_df_safe(ekegg),
    ORA_Reactome = as_df_safe(ereact),
    GSEA_GO_BP = as_df_safe(gsea_go),
    GSEA_GO_MF = as_df_safe(gsea_go_MF),
    GSEA_KEGG = as_df_safe(gsea_kegg)
  )

  subset_tables <- list()
  subset_stacked <- list()

  for (src in names(source_tables)) {
    df <- source_tables[[src]]
    if (!"Description" %in% names(df)) df$Description <- rep(NA_character_, nrow(df))

    damp_df <- df[!is.na(df$Description) & grepl(damp_pattern, df$Description, ignore.case = TRUE), , drop = FALSE]
    lipid_df <- df[!is.na(df$Description) & grepl(lipid_pattern, df$Description, ignore.case = TRUE), , drop = FALSE]
    sterile_df <- df[!is.na(df$Description) & grepl(sterile_pattern, df$Description, ignore.case = TRUE), , drop = FALSE]

    out_xlsx <- file.path(dirs$tables_sig, paste0("term_subsets_", tolower(gsub("[^A-Za-z0-9]+", "_", src)), ".xlsx"))
    openxlsx::write.xlsx(
      list(
        DAMP = if (nrow(damp_df) > 0) damp_df else data.frame(Note = "No results", stringsAsFactors = FALSE),
        LIPID = if (nrow(lipid_df) > 0) lipid_df else data.frame(Note = "No results", stringsAsFactors = FALSE),
        STERILE = if (nrow(sterile_df) > 0) sterile_df else data.frame(Note = "No results", stringsAsFactors = FALSE)
      ),
      file = out_xlsx,
      overwrite = TRUE
    )
    register_if_exists("table", paste0("term_subsets_", tolower(src)), out_xlsx, "DAMP/lipid/sterile term subset workbook.")

    subset_tables[[src]] <- list(DAMP = damp_df, LIPID = lipid_df, STERILE = sterile_df)

    subset_stacked[[paste0(src, "_DAMP")]] <- clean_subset(damp_df, src, "DAMP")
    subset_stacked[[paste0(src, "_LIPID")]] <- clean_subset(lipid_df, src, "LIPID")
    subset_stacked[[paste0(src, "_STERILE")]] <- clean_subset(sterile_df, src, "STERILE")

    subset_plot(damp_df, src, "DAMP")
    subset_plot(lipid_df, src, "LIPID")
    subset_plot(sterile_df, src, "STERILE")
  }

  list(
    panel_hits = panel_hits,
    subset_tables = subset_tables,
    subset_stacked = dplyr::bind_rows(subset_stacked),
    panel_xlsx = panel_xlsx
  )
}
