build_combined_contrast_export <- function(fit, cleaned, contrast_results, config) {
  gene_ids <- rownames(fit$dds_wald)
  gene_annotation <- cleaned$gene_annotation[gene_ids, , drop = FALSE]
  annotation_block <- gene_annotation[, setdiff(colnames(gene_annotation), "Geneid_clean"), drop = FALSE]

  raw_counts <- as.data.frame(DESeq2::counts(fit$dds_wald, normalized = FALSE))
  raw_counts <- raw_counts[gene_ids, , drop = FALSE]
  colnames(raw_counts) <- paste0("raw_", colnames(raw_counts))

  normalized_counts <- as.data.frame(DESeq2::counts(fit$dds_wald, normalized = TRUE))
  normalized_counts <- normalized_counts[gene_ids, , drop = FALSE]
  colnames(normalized_counts) <- paste0("norm_", colnames(normalized_counts))

  norm_matrix <- DESeq2::counts(fit$dds_wald, normalized = TRUE)
  mean_norm_by_condition <- lapply(levels(fit$dds_wald$condition), function(level_name) {
    idx <- fit$dds_wald$condition == level_name
    values <- rowMeans(norm_matrix[, idx, drop = FALSE])
    out <- data.frame(values, row.names = gene_ids, check.names = FALSE)
    colnames(out) <- paste0("mean_norm_", level_name)
    out
  })
  mean_norm_by_condition <- do.call(cbind, mean_norm_by_condition)

  first_label <- names(contrast_results)[1]
  if (is.null(first_label) || is.null(contrast_results[[first_label]])) {
    stop("No contrast results are available for the combined export.", call. = FALSE)
  }

  combined <- data.frame(
    Geneid_clean = gene_ids,
    annotation_block,
    baseMean = contrast_results[[first_label]]$export_table[gene_ids, "baseMean", drop = TRUE],
    raw_counts,
    normalized_counts,
    mean_norm_by_condition,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (label in names(contrast_results)) {
    contrast_result <- contrast_results[[label]]
    export_table <- contrast_result$export_table[gene_ids, , drop = FALSE]
    export_columns <- intersect(
      c("baseMean", "log2FC", "log2FC_shrunken", "lfcSE_shrunken", "stat", "pvalue", "padj"),
      colnames(export_table)
    )
    export_block <- export_table[, export_columns, drop = FALSE]
    colnames(export_block) <- paste0(label, "__", colnames(export_block))
    export_block[[paste0(label, "__DE_call")]] <- make_de_call(
      log2fc = export_table$log2FC,
      padj = export_table$padj,
      group1 = contrast_result$group1,
      group0 = contrast_result$group0,
      alpha = config$results_alpha
    )

    combined <- cbind(combined, export_block)
  }

  rownames(combined) <- gene_ids
  combined
}

write_repair_log <- function(repair_log, file_path) {
  utils::write.csv(repair_log, file_path, row.names = FALSE)
}

write_contrast_outputs <- function(contrast_result, config) {
  label <- contrast_result$label
  output_dir <- config$output_dir

  utils::write.csv(
    contrast_result$export_table,
    file.path(output_dir, paste0("DESeq2_", label, ".csv"))
  )
  utils::write.csv(
    contrast_result$export_table_by_ensembl,
    file.path(output_dir, paste0("DESeq2_", label, "_by_ENS.csv"))
  )
  utils::write.csv(
    contrast_result$export_table_by_symbol,
    file.path(output_dir, paste0("DESeq2_", label, "_by_SYMBOL.csv"))
  )

  write_sig_gene_counts(
    contrast_result$sig_gene_counts,
    file.path(output_dir, paste0("sig_gene_counts_", label, ".txt"))
  )
  make_ma_plot(
    contrast_result$raw_results,
    file.path(output_dir, paste0("MA_", label, ".pdf"))
  )
  make_volcano_plot(
    contrast_result$raw_results_df,
    contrast_result$group1,
    contrast_result$group0,
    config$results_alpha,
    file.path(output_dir, paste0("volcano_", label, ".pdf"))
  )
  write_gsea_outputs(contrast_result$gsea, output_dir, label)
}

write_all_outputs <- function(fit, cleaned, contrast_results, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  utils::write.csv(
    DESeq2::counts(fit$dds_wald, normalized = TRUE),
    file.path(config$output_dir, "normalized_counts.csv")
  )
  utils::write.csv(
    as.data.frame(SummarizedExperiment::colData(fit$dds_wald)),
    file.path(config$output_dir, "coldata_used.csv")
  )
  utils::write.csv(
    as.data.frame(cleaned$gene_annotation[rownames(fit$dds_wald), , drop = FALSE]),
    file.path(config$output_dir, "gene_annotation_used.csv")
  )
  utils::write.csv(
    cleaned$contrast_plan,
    file.path(config$output_dir, "contrast_plan_used.csv"),
    row.names = FALSE
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(config$output_dir, "sessionInfo.txt")
  )
  write_repair_log(
    cleaned$repair_log,
    file.path(config$output_dir, config$repair_log_file)
  )

  combined_export <- build_combined_contrast_export(fit, cleaned, contrast_results, config)
  utils::write.csv(
    combined_export,
    file.path(config$output_dir, config$combined_export_file),
    row.names = FALSE
  )

  if (!is.null(fit$lrt_results)) {
    utils::write.csv(
      as.data.frame(fit$lrt_results),
      file.path(config$output_dir, "LRT_condition_vs_null.csv")
    )
  }

  for (contrast_result in contrast_results) {
    write_contrast_outputs(contrast_result, config)
  }

  invisible(list(combined_export = combined_export))
}

print_run_summary <- function(fit, cleaned, contrast_results, config) {
  cat("\n=== Run summary ===\n")
  cat("Samples:", ncol(fit$dds_wald), "| Genes (after filter):", nrow(fit$dds_wald), "\n")
  cat("Condition levels:", paste(levels(fit$dds_wald$condition), collapse = ", "), "\n")
  cat("Reference level:", fit$ref_level, "\n")
  cat("Abundance filter:", describe_filter(config), "\n")
  cat("Results alpha:", config$results_alpha, "\n")
  cat(
    "Contrasts:",
    paste(
      paste(cleaned$contrast_plan$group1, "vs", cleaned$contrast_plan$group0),
      collapse = "; "
    ),
    "\n"
  )
  cat("Combined export:", file.path(config$output_dir, config$combined_export_file), "\n")
  cat("Repair log:", file.path(config$output_dir, config$repair_log_file), "\n")
  cat("Outputs in:", normalizePath(config$output_dir), "\n")
  invisible(contrast_results)
}
