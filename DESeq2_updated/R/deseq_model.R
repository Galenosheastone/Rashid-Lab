configure_parallel <- function(n_workers) {
  workers <- as.integer(n_workers)
  if (is.na(workers) || workers < 1L) {
    workers <- 1L
  }

  param <- tryCatch(
    {
      if (workers == 1L) {
        candidate <- BiocParallel::SerialParam()
      } else if (.Platform$OS.type == "windows") {
        candidate <- BiocParallel::SnowParam(workers = workers)
      } else {
        candidate <- BiocParallel::MulticoreParam(workers = workers)
      }

      started_param <- BiocParallel::bpstart(candidate)
      BiocParallel::bpstop(started_param)
      candidate
    },
    error = function(e) {
      warning(
        "Falling back to SerialParam() because parallel worker setup failed: ",
        conditionMessage(e),
        call. = FALSE
      )
      BiocParallel::SerialParam()
    }
  )

  BiocParallel::register(param, default = TRUE)
  invisible(param)
}

keep_by_group_abundance <- function(count_matrix, condition, min_count = 10L, min_reps = 3L) {
  stopifnot(min_count >= 1L, min_reps >= 1L)

  condition <- factor(condition)
  keep_by_level <- lapply(levels(condition), function(level_name) {
    idx <- condition == level_name
    rowSums(count_matrix[, idx, drop = FALSE] >= min_count) >= min_reps
  })

  Reduce("|", keep_by_level)
}

fit_deseq_pipeline <- function(cleaned, config) {
  configure_parallel(config$n_workers)

  count_matrix <- round(cleaned$counts_matrix)
  storage.mode(count_matrix) <- "integer"
  coldata <- cleaned$metadata

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = ~ condition
  )

  if (identical(config$filter_mode, "legacy")) {
    keep <- rowSums(DESeq2::counts(dds) > 0) > 0
  } else {
    keep <- keep_by_group_abundance(
      count_matrix = DESeq2::counts(dds),
      condition = dds$condition,
      min_count = config$min_count,
      min_reps = config$min_reps_within_group
    )
  }

  if (!any(keep)) {
    stop("No genes passed the configured abundance filter.", call. = FALSE)
  }

  dds <- dds[keep, ]
  dds_wald <- DESeq2::DESeq(dds, parallel = TRUE)
  coef_names <- DESeq2::resultsNames(dds_wald)

  dds_lrt <- NULL
  lrt_results <- NULL
  if (isTRUE(config$run_lrt)) {
    dds_lrt <- DESeq2::DESeq(dds, test = "LRT", reduced = ~ 1, parallel = TRUE)
    lrt_results <- DESeq2::results(dds_lrt, alpha = config$results_alpha, parallel = TRUE)
  }

  transformed <- if (isTRUE(config$use_rlog)) {
    DESeq2::rlog(dds_wald, blind = FALSE)
  } else {
    DESeq2::vst(dds_wald, blind = FALSE)
  }

  list(
    dds_wald = dds_wald,
    dds_lrt = dds_lrt,
    lrt_results = lrt_results,
    vsd = transformed,
    coef_names = coef_names,
    ref_level = cleaned$ref_level,
    n_genes_before_filter = nrow(cleaned$counts_matrix),
    n_genes_after_filter = nrow(dds_wald)
  )
}
