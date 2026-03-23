labelize_contrast <- function(group1, group0) {
  gsub("[^A-Za-z0-9._-]+", "_", paste0(group1, "_vs_", group0))
}

find_coef_name <- function(dds_obj, grp1, grp0) {
  result_names <- DESeq2::resultsNames(dds_obj)
  exact_match <- paste0("condition_", make.names(grp1), "_vs_", make.names(grp0))

  if (exact_match %in% result_names) {
    return(exact_match)
  }

  grp1_key <- gsub("[^A-Za-z0-9.]+", ".", grp1)
  grp0_key <- gsub("[^A-Za-z0-9.]+", ".", grp0)
  pattern <- paste0("^condition[_:]*", grp1_key, ".*_vs_.*", grp0_key, "$")
  hits <- grep(pattern, result_names, value = TRUE)

  if (length(hits) >= 1L) {
    return(hits[1L])
  }

  stop(
    "Could not find a DESeq2 coefficient for ", grp1, " vs ", grp0,
    " in resultsNames(dds).",
    call. = FALSE
  )
}

compute_contrast_results <- function(dds_obj, group1, group0, alpha = 0.05) {
  DESeq2::results(
    dds_obj,
    contrast = c("condition", group1, group0),
    alpha = alpha,
    parallel = TRUE
  )
}

shrink_contrast_results <- function(dds_obj, group1, group0, shrink_method = "apeglm") {
  method <- tolower(shrink_method)
  contrast_vector <- c("condition", group1, group0)

  if (identical(method, "ashr")) {
    shrunk <- DESeq2::lfcShrink(
      dds_obj,
      contrast = contrast_vector,
      type = "ashr",
      parallel = TRUE
    )

    return(list(
      result = shrunk,
      method = method,
      coef_name = NA_character_,
      used_refit = FALSE,
      refit_reference = NA_character_,
      error = NA_character_
    ))
  }

  coef_name <- tryCatch(
    find_coef_name(dds_obj, group1, group0),
    error = function(e) NULL
  )

  if (!is.null(coef_name)) {
    shrunk <- DESeq2::lfcShrink(
      dds_obj,
      coef = coef_name,
      type = "apeglm",
      parallel = TRUE
    )

    return(list(
      result = shrunk,
      method = method,
      coef_name = coef_name,
      used_refit = FALSE,
      refit_reference = NA_character_,
      error = NA_character_
    ))
  }

  message(
    "apeglm shrinkage for ", group1, " vs ", group0,
    " requires an explicit refit with reference level '", group0, "'."
  )

  dds_refit <- dds_obj
  dds_refit$condition <- stats::relevel(dds_refit$condition, ref = group0)
  dds_refit <- DESeq2::DESeq(dds_refit, parallel = TRUE)
  coef_name <- find_coef_name(dds_refit, group1, group0)
  shrunk <- DESeq2::lfcShrink(
    dds_refit,
    coef = coef_name,
    type = "apeglm",
    parallel = TRUE
  )

  list(
    result = shrunk,
    method = method,
    coef_name = coef_name,
    used_refit = TRUE,
    refit_reference = group0,
    error = NA_character_
  )
}

build_contrast_export_table <- function(raw_results, shrunk_results = NULL) {
  if (!is.null(shrunk_results)) {
    data.frame(
      baseMean = raw_results$baseMean,
      log2FC = raw_results$log2FoldChange,
      log2FC_shrunken = shrunk_results$log2FoldChange,
      lfcSE_shrunken = shrunk_results$lfcSE,
      stat = raw_results$stat,
      pvalue = raw_results$pvalue,
      padj = raw_results$padj,
      row.names = rownames(raw_results),
      check.names = FALSE
    )
  } else {
    data.frame(
      baseMean = raw_results$baseMean,
      log2FC = raw_results$log2FoldChange,
      stat = raw_results$stat,
      pvalue = raw_results$pvalue,
      padj = raw_results$padj,
      row.names = rownames(raw_results),
      check.names = FALSE
    )
  }
}

make_labeled_export_tables <- function(export_table, gene_annotation) {
  gene_subset <- gene_annotation[rownames(export_table), , drop = FALSE]

  if (!"SYMBOL" %in% colnames(gene_subset)) {
    gene_subset$SYMBOL <- NA_character_
  }

  if (!"Ensembl_GeneID" %in% colnames(gene_subset)) {
    gene_subset$Ensembl_GeneID <- NA_character_
  }

  ensembl_labels <- trimws(as.character(gene_subset$Ensembl_GeneID))
  ensembl_labels[is.na(ensembl_labels) | !nzchar(ensembl_labels)] <- rownames(export_table)
  ensembl_labels <- make.unique(ensembl_labels)

  symbol_labels <- trimws(as.character(gene_subset$SYMBOL))
  symbol_labels[is.na(symbol_labels) | !nzchar(symbol_labels)] <- rownames(export_table)
  symbol_labels <- make.unique(symbol_labels)

  export_by_ensembl <- export_table
  rownames(export_by_ensembl) <- ensembl_labels

  export_by_symbol <- export_table
  rownames(export_by_symbol) <- symbol_labels

  list(
    by_ensembl = export_by_ensembl,
    by_symbol = export_by_symbol
  )
}

summarize_sig_gene_counts <- function(raw_results, alpha = 0.05) {
  is_sig <- !is.na(raw_results$padj) & raw_results$padj < alpha

  data.frame(
    up = sum(is_sig & raw_results$log2FoldChange > 0, na.rm = TRUE),
    down = sum(is_sig & raw_results$log2FoldChange < 0, na.rm = TRUE)
  )
}

make_de_call <- function(log2fc, padj, group1, group0, alpha = 0.05) {
  out <- rep("not_sig", length(log2fc))
  significant <- !is.na(padj) & padj < alpha
  out[significant & log2fc > 0] <- paste0(group1, "_up")
  out[significant & log2fc < 0] <- paste0(group0, "_up")
  out
}

write_sig_gene_counts <- function(sig_gene_counts, file_path) {
  utils::write.table(
    sig_gene_counts,
    file = file_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

make_ma_plot <- function(raw_results, file_path) {
  grDevices::pdf(file_path)
  DESeq2::plotMA(raw_results, ylim = c(-5, 5))
  grDevices::dev.off()
}

make_volcano_plot <- function(raw_results_df, group1, group0, alpha, file_path) {
  plot_df <- transform(
    raw_results_df,
    neglog10padj = -log10(pmin(padj, 1)),
    log2FC = log2FoldChange
  )
  plot_df <- plot_df[is.finite(plot_df$neglog10padj) & is.finite(plot_df$log2FC), , drop = FALSE]

  grDevices::pdf(file_path)
  graphics::plot(
    plot_df$log2FC,
    plot_df$neglog10padj,
    pch = 20,
    xlab = "log2 fold-change",
    ylab = "-log10(FDR)",
    main = paste("Volcano -", group1, "vs", group0)
  )
  graphics::abline(h = -log10(alpha), lty = 2)
  grDevices::dev.off()
}

run_requested_contrasts <- function(fit, cleaned, config) {
  contrast_plan <- cleaned$contrast_plan
  contrast_results <- vector("list", nrow(contrast_plan))
  names(contrast_results) <- vapply(
    seq_len(nrow(contrast_plan)),
    function(i) labelize_contrast(contrast_plan$group1[i], contrast_plan$group0[i]),
    character(1)
  )

  for (i in seq_len(nrow(contrast_plan))) {
    group1 <- contrast_plan$group1[i]
    group0 <- contrast_plan$group0[i]
    label <- labelize_contrast(group1, group0)

    message("Running contrast: ", group1, " vs ", group0)

    raw_results <- compute_contrast_results(
      fit$dds_wald,
      group1 = group1,
      group0 = group0,
      alpha = config$results_alpha
    )

    shrinkage <- tryCatch(
      shrink_contrast_results(
        fit$dds_wald,
        group1 = group1,
        group0 = group0,
        shrink_method = config$shrink_method
      ),
      error = function(e) {
        message("Shrinkage failed for ", label, ": ", conditionMessage(e))
        list(
          result = NULL,
          method = config$shrink_method,
          coef_name = NA_character_,
          used_refit = FALSE,
          refit_reference = NA_character_,
          error = conditionMessage(e)
        )
      }
    )

    export_table <- build_contrast_export_table(raw_results, shrinkage$result)
    labeled_tables <- make_labeled_export_tables(export_table, cleaned$gene_annotation)
    gsea_result <- run_contrast_gsea(
      raw_results = raw_results,
      gene_annotation = cleaned$gene_annotation[rownames(raw_results), , drop = FALSE],
      config = config,
      label = label
    )

    contrast_results[[label]] <- list(
      label = label,
      group1 = group1,
      group0 = group0,
      contrast = c("condition", group1, group0),
      raw_results = raw_results,
      raw_results_df = as.data.frame(raw_results),
      shrunken_results = shrinkage$result,
      shrinkage = shrinkage[setdiff(names(shrinkage), "result")],
      export_table = export_table,
      export_table_by_ensembl = labeled_tables$by_ensembl,
      export_table_by_symbol = labeled_tables$by_symbol,
      sig_gene_counts = summarize_sig_gene_counts(raw_results, config$results_alpha),
      gsea = gsea_result
    )
  }

  contrast_results
}
