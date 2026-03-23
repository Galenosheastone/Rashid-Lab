save_size_factor_plot <- function(dds, file_path) {
  grDevices::pdf(file_path)
  graphics::barplot(
    DESeq2::sizeFactors(dds),
    las = 2,
    ylab = "sizeFactor",
    main = "Size factors"
  )
  grDevices::dev.off()
}

save_dispersion_plot <- function(dds, file_path) {
  grDevices::pdf(file_path)
  DESeq2::plotDispEsts(dds)
  grDevices::dev.off()
}

save_sample_distance_heatmap <- function(vsd, file_path) {
  sample_distance <- as.matrix(stats::dist(t(SummarizedExperiment::assay(vsd))))
  grDevices::pdf(file_path)
  pheatmap::pheatmap(
    sample_distance,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    main = "Sample-sample distance (VST/rlog)"
  )
  grDevices::dev.off()
}

save_cooks_distance_plot <- function(dds, file_path) {
  cook_values <- SummarizedExperiment::assays(dds)[["cooks"]]
  if (is.null(cook_values)) {
    return(invisible(NULL))
  }

  grDevices::pdf(file_path)
  graphics::boxplot(
    log10(cook_values + 1e-8),
    range = 0,
    outline = FALSE,
    las = 2,
    ylab = "log10 Cook's distance",
    main = "Influential observations"
  )
  grDevices::dev.off()
}

save_mean_sd_plot <- function(vsd, file_path) {
  grDevices::pdf(file_path)
  vsn::meanSdPlot(SummarizedExperiment::assay(vsd))
  grDevices::dev.off()
}

save_pvalue_histogram <- function(dds, coef_names, alpha, file_path) {
  first_coef <- grep("^condition", coef_names, value = TRUE)[1]
  if (is.na(first_coef)) {
    return(invisible(NULL))
  }

  pvalue_results <- DESeq2::results(dds, name = first_coef, alpha = alpha, parallel = TRUE)

  grDevices::pdf(file_path)
  graphics::hist(
    pvalue_results$pvalue,
    breaks = 50,
    col = "grey",
    main = paste("P-value distribution -", first_coef),
    xlab = "p-value"
  )
  grDevices::dev.off()
}

make_qc_outputs <- function(fit, config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  pca_plot <- DESeq2::plotPCA(fit$vsd, intgroup = "condition") +
    ggplot2::ggtitle("PCA - variance-stabilized counts")

  ggplot2::ggsave(
    filename = file.path(config$output_dir, "QC_PCA_vst_or_rlog.pdf"),
    plot = pca_plot,
    width = 6.5,
    height = 5
  )

  save_size_factor_plot(fit$dds_wald, file.path(config$output_dir, "QC_size_factors.pdf"))
  save_dispersion_plot(fit$dds_wald, file.path(config$output_dir, "QC_dispersion_trend.pdf"))
  save_sample_distance_heatmap(
    fit$vsd,
    file.path(config$output_dir, "QC_sample_distance_heatmap.pdf")
  )
  save_cooks_distance_plot(
    fit$dds_wald,
    file.path(config$output_dir, "QC_cooks_distance_boxplot.pdf")
  )
  save_mean_sd_plot(fit$vsd, file.path(config$output_dir, "QC_meanSD_vst_or_rlog.pdf"))
  save_pvalue_histogram(
    fit$dds_wald,
    fit$coef_names,
    config$results_alpha,
    file.path(config$output_dir, "QC_pvalue_histogram.pdf")
  )

  invisible(TRUE)
}
