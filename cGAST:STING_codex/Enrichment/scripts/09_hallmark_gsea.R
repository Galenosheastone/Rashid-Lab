run_hallmark_gsea_module <- function(
  ranks,
  hallmark_species,
  hallmark_selected,
  minGSSize,
  maxGSSize,
  dirs,
  save_plot,
  register_if_exists
) {
  empty_out <- list(
    fgsea = data.frame(),
    selected = data.frame(),
    selected_summary = data.frame(),
    selected_plot = NULL,
    all_xlsx = NULL,
    selected_xlsx = NULL
  )

  if (!requireNamespace("msigdbr", quietly = TRUE) || !requireNamespace("fgsea", quietly = TRUE)) {
    note_file <- file.path(dirs$tables_sig, "hallmark_gsea_note.txt")
    writeLines("Hallmark GSEA skipped because msigdbr and/or fgsea is not installed.", note_file)
    register_if_exists("table", "hallmark_gsea_note", note_file, "Hallmark module skip note.")
    return(empty_out)
  }

  rank_names <- names(ranks)
  ranks <- suppressWarnings(as.numeric(ranks))
  names(ranks) <- as.character(rank_names)
  keep <- !is.na(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))
  ranks <- ranks[keep]
  ranks <- sort(ranks, decreasing = TRUE)

  if (length(ranks) < 10) {
    note_file <- file.path(dirs$tables_sig, "hallmark_gsea_note.txt")
    writeLines("Hallmark GSEA skipped because ranks vector is too small after filtering.", note_file)
    register_if_exists("table", "hallmark_gsea_note", note_file, "Hallmark module skip note.")
    return(empty_out)
  }

  msig_df <- msigdbr::msigdbr(species = hallmark_species, category = "H")
  if (nrow(msig_df) == 0) {
    note_file <- file.path(dirs$tables_sig, "hallmark_gsea_note.txt")
    writeLines(paste0("No Hallmark gene sets were returned for species: ", hallmark_species), note_file)
    register_if_exists("table", "hallmark_gsea_note", note_file, "Hallmark module skip note.")
    return(empty_out)
  }

  gene_col <- if ("entrez_gene" %in% names(msig_df)) {
    "entrez_gene"
  } else if ("ncbi_gene" %in% names(msig_df)) {
    "ncbi_gene"
  } else {
    NULL
  }
  if (is.null(gene_col)) {
    stop("msigdbr output lacks an entrez gene column (expected 'entrez_gene' or 'ncbi_gene').")
  }

  hallmark_term2gene <- msig_df %>%
    dplyr::transmute(gs_name = .data$gs_name, entrez_gene = as.character(.data[[gene_col]])) %>%
    dplyr::filter(!is.na(entrez_gene), nzchar(entrez_gene)) %>%
    dplyr::distinct(gs_name, entrez_gene)

  hallmark_list <- split(hallmark_term2gene$entrez_gene, hallmark_term2gene$gs_name)
  hallmark_list <- hallmark_list[lengths(hallmark_list) > 0]

  fgsea_res <- fgsea::fgsea(
    pathways = hallmark_list,
    stats = ranks,
    minSize = minGSSize,
    maxSize = maxGSSize
  )
  fgsea_df <- as.data.frame(fgsea_res)
  if (nrow(fgsea_df) == 0) {
    fgsea_df <- data.frame(Note = "No results", stringsAsFactors = FALSE)
  }

  fgsea_export <- fgsea_df
  if ("leadingEdge" %in% names(fgsea_export)) {
    fgsea_export$leadingEdge_length <- vapply(fgsea_export$leadingEdge, length, integer(1))
    fgsea_export$leadingEdge_genes <- vapply(fgsea_export$leadingEdge, function(x) paste(x, collapse = "/"), character(1))
    fgsea_export$leadingEdge <- NULL
  }
  if ("pathway" %in% names(fgsea_export) && "padj" %in% names(fgsea_export)) {
    fgsea_export <- fgsea_export %>% dplyr::arrange(.data$padj)
  }

  all_xlsx <- file.path(dirs$tables_sig, "hallmark_fgsea_all.xlsx")
  openxlsx::write.xlsx(fgsea_export, all_xlsx, overwrite = TRUE)
  register_if_exists("table", "hallmark_fgsea_all", all_xlsx, "Hallmark fgsea results.")

  selected_df <- if ("pathway" %in% names(fgsea_export)) {
    fgsea_export %>% dplyr::filter(.data$pathway %in% hallmark_selected)
  } else {
    data.frame()
  }
  if (nrow(selected_df) > 0) {
    selected_summary <- selected_df
    if (!"leadingEdge_length" %in% names(selected_summary)) {
      selected_summary$leadingEdge_length <- NA_real_
    } else {
      selected_summary$leadingEdge_length <- suppressWarnings(as.numeric(selected_summary$leadingEdge_length))
    }
    selected_summary <- selected_summary %>%
      dplyr::select(dplyr::any_of(c("pathway", "NES", "padj", "leadingEdge_length")))
  } else {
    selected_summary <- data.frame()
  }

  selected_xlsx <- file.path(dirs$tables_sig, "hallmark_selected.xlsx")
  openxlsx::write.xlsx(
    list(
      selected = if (nrow(selected_df) > 0) selected_df else data.frame(Note = "No results", stringsAsFactors = FALSE),
      summary = if (nrow(selected_summary) > 0) selected_summary else data.frame(Note = "No results", stringsAsFactors = FALSE)
    ),
    selected_xlsx,
    overwrite = TRUE
  )
  register_if_exists("table", "hallmark_selected", selected_xlsx, "Selected Hallmark summary.")

  selected_plot <- NULL
  if (nrow(selected_summary) > 0 && "NES" %in% names(selected_summary)) {
    sel_plot_df <- selected_summary %>% dplyr::arrange(.data$NES)
    p <- ggplot2::ggplot(sel_plot_df, ggplot2::aes(x = stats::reorder(pathway, NES), y = NES, fill = NES > 0)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "royalblue"), guide = "none") +
      ggplot2::labs(title = "Selected Hallmark GSEA NES", x = NULL, y = "NES") +
      ggplot2::theme_bw(base_size = 11)
    selected_plot <- file.path(dirs$signatures, "hallmark_selected_NES_bar.pdf")
    save_plot(p, selected_plot, 9, 6)
  }

  if (nrow(selected_df) > 0 && "padj" %in% names(selected_df) && "pathway" %in% names(selected_df)) {
    top_sel <- selected_df %>%
      dplyr::arrange(.data$padj) %>%
      dplyr::filter(!is.na(.data$padj), .data$padj < 0.1) %>%
      dplyr::slice_head(n = 3)
    if (nrow(top_sel) > 0) {
      running_pdf <- file.path(dirs$signatures, "hallmark_selected_running_scores.pdf")
      grDevices::pdf(running_pdf, width = 8, height = 6)
      for (pw in top_sel$pathway) {
        gene_set <- hallmark_list[[pw]]
        if (is.null(gene_set)) next
        p <- tryCatch(
          fgsea::plotEnrichment(gene_set, ranks) + ggplot2::ggtitle(pw),
          error = function(e) NULL
        )
        if (!is.null(p)) print(p)
      }
      grDevices::dev.off()
      register_if_exists("plot", "hallmark_selected_running_scores", running_pdf, "Top selected hallmark running-score plots.")
    }
  }

  list(
    fgsea = fgsea_export,
    selected = selected_df,
    selected_summary = selected_summary,
    selected_plot = selected_plot,
    all_xlsx = all_xlsx,
    selected_xlsx = selected_xlsx
  )
}
