#!/usr/bin/env Rscript

# pathway_scoring.R
# -----------------
# Command-line pathway scoring for bulk RNA-seq expression tables that include
# DESeq2-style result columns plus normalized expression columns.
#
# Example usage with CSV/TSV gene sets:
#   Rscript pathway_scoring.R \
#     --input DESeq2_full_combined_ssRNAseq.csv \
#     --genesets curated_gene_sets.csv \
#     --gene-id-col SYMBOL \
#     --expr-prefix norm_ \
#     --group-pattern '^norm_([A-Za-z]+)_(\\d+)$' \
#     --outdir pathway_scoring_results \
#     --methods gsva,ssgsea,singscore \
#     --min-set-size 5 \
#     --max-set-size 500 \
#     --low-expr-threshold 1 \
#     --low-expr-min-samples 2 \
#     --stats-test kruskal \
#     --disable-pairwise-correction
#
# Example usage with GMT gene sets:
#   Rscript pathway_scoring.R \
#     --input DESeq2_full_combined_ssRNAseq.csv \
#     --genesets curated_gene_sets.gmt \
#     --outdir pathway_scoring_results
#
# Notes:
# - Expression input is restricted to columns whose names start with --expr-prefix.
# - Gene IDs default to SYMBOL with fallback to Geneid_clean and Ensembl_GeneID.
# - Duplicate gene IDs are collapsed by keeping the row with the highest mean
#   normalized expression across samples.
# - Gene sets are read from an external CSV/TSV or GMT file; no pathways are
#   hard-coded in this script.
# - Pairwise p-value correction is enabled by default and can be disabled with
#   --disable-pairwise-correction.

options(stringsAsFactors = FALSE, warn = 1)

timestamp <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

msg <- function(fmt, ...) {
  cat(sprintf("[%s] %s\n", timestamp(), sprintf(fmt, ...)))
}

warnf <- function(fmt, ...) {
  warning(sprintf(fmt, ...), call. = FALSE, immediate. = TRUE)
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

trim_to_char <- function(x) {
  out <- trimws(as.character(x))
  out[is.na(out)] <- ""
  out
}

is_non_empty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

unique_keep_order <- function(x) {
  x <- as.character(x)
  x[!duplicated(x)]
}

find_ci_col <- function(df_names, target_name) {
  hits <- df_names[tolower(df_names) == tolower(target_name)]
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[[1]]
}

safe_first_non_na <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

check_required_packages <- function(methods, geneset_file, plot_method_correlation) {
  required <- c("readr", "tibble", "dplyr", "tidyr", "stringr", "purrr", "ggplot2", "pheatmap")
  if (tolower(tools::file_ext(geneset_file)) == "gmt") {
    required <- c(required, "GSEABase")
  }
  if (any(c("gsva", "ssgsea") %in% methods)) {
    required <- c(required, "GSVA")
  }
  if ("singscore" %in% methods) {
    required <- c(required, "singscore")
  }
  if (plot_method_correlation) {
    required <- c(required, "pheatmap")
  }
  missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stopf(
      "Missing required package(s): %s\nInstall them before running the script.",
      paste(unique(missing), collapse = ", ")
    )
  }
}

parse_methods <- function(method_string) {
  methods <- unlist(strsplit(method_string, ",", fixed = TRUE))
  methods <- trimws(tolower(methods))
  methods <- methods[nzchar(methods)]
  valid <- c("gsva", "ssgsea", "singscore")
  unknown <- setdiff(methods, valid)
  if (length(methods) == 0) {
    stopf("No methods were requested. Use --methods with one or more of: %s", paste(valid, collapse = ", "))
  }
  if (length(unknown) > 0) {
    stopf("Unsupported method(s): %s. Valid choices are: %s", paste(unknown, collapse = ", "), paste(valid, collapse = ", "))
  }
  unique(methods)
}

pathway_scoring_defaults <- function() {
  list(
    input = NULL,
    genesets = NULL,
    gene_id_col = "SYMBOL",
    expr_prefix = "norm_",
    group_pattern = "^norm_([A-Za-z]+)_([0-9]+)$",
    outdir = "pathway_scoring_results",
    methods = c("gsva", "ssgsea", "singscore"),
    min_set_size = 5L,
    max_set_size = 500L,
    low_expr_threshold = 1,
    low_expr_min_samples = 2L,
    disable_low_expr_filter = FALSE,
    stats_test = "kruskal",
    pairwise_test = "auto",
    pairwise_correction = TRUE,
    adjust_method = "BH",
    plot_method_correlation = FALSE,
    boxplot_significance = TRUE
  )
}

normalize_run_options <- function(opts) {
  if (!is.list(opts)) {
    stopf("Options must be provided as a named list.")
  }

  defaults <- pathway_scoring_defaults()
  unknown <- setdiff(names(opts), names(defaults))
  if (length(unknown) > 0) {
    stopf("Unknown option(s): %s", paste(unknown, collapse = ", "))
  }

  for (name in names(defaults)) {
    if (is.null(opts[[name]])) {
      opts[[name]] <- defaults[[name]]
    }
  }

  opts$input <- if (length(opts$input) == 0 || is.null(opts$input)) "" else trimws(as.character(opts$input[[1]]))
  opts$genesets <- if (length(opts$genesets) == 0 || is.null(opts$genesets)) "" else trimws(as.character(opts$genesets[[1]]))
  opts$gene_id_col <- if (length(opts$gene_id_col) == 0 || is.null(opts$gene_id_col)) "SYMBOL" else trimws(as.character(opts$gene_id_col[[1]]))
  opts$expr_prefix <- if (length(opts$expr_prefix) == 0 || is.null(opts$expr_prefix)) "norm_" else as.character(opts$expr_prefix[[1]])
  opts$group_pattern <- if (length(opts$group_pattern) == 0 || is.null(opts$group_pattern)) {
    pathway_scoring_defaults()$group_pattern
  } else {
    as.character(opts$group_pattern[[1]])
  }
  opts$outdir <- if (length(opts$outdir) == 0 || is.null(opts$outdir)) {
    pathway_scoring_defaults()$outdir
  } else {
    trimws(as.character(opts$outdir[[1]]))
  }

  if (length(opts$methods) == 0 || is.null(opts$methods)) {
    opts$methods <- pathway_scoring_defaults()$methods
  }
  if (length(opts$methods) == 1) {
    opts$methods <- parse_methods(as.character(opts$methods[[1]]))
  } else {
    opts$methods <- parse_methods(paste(as.character(opts$methods), collapse = ","))
  }

  opts$min_set_size <- as.integer(opts$min_set_size[[1]])
  opts$max_set_size <- as.integer(opts$max_set_size[[1]])
  opts$low_expr_threshold <- as.numeric(opts$low_expr_threshold[[1]])
  opts$low_expr_min_samples <- as.integer(opts$low_expr_min_samples[[1]])
  opts$disable_low_expr_filter <- isTRUE(as.logical(opts$disable_low_expr_filter[[1]]))
  opts$stats_test <- trimws(tolower(as.character(opts$stats_test[[1]])))
  opts$pairwise_test <- trimws(tolower(as.character(opts$pairwise_test[[1]])))
  opts$pairwise_correction <- isTRUE(as.logical(opts$pairwise_correction[[1]]))
  opts$adjust_method <- as.character(opts$adjust_method[[1]])
  opts$plot_method_correlation <- isTRUE(as.logical(opts$plot_method_correlation[[1]]))
  opts$boxplot_significance <- isTRUE(as.logical(opts$boxplot_significance[[1]]))

  if (!nzchar(opts$input)) {
    stopf("Missing required argument: input")
  }
  if (!nzchar(opts$genesets)) {
    stopf("Missing required argument: genesets")
  }
  if (!file.exists(opts$input)) {
    stopf("Input file does not exist: %s", opts$input)
  }
  if (!file.exists(opts$genesets)) {
    stopf("Gene-set file does not exist: %s", opts$genesets)
  }
  if (is.na(opts$min_set_size) || opts$min_set_size < 1) {
    stopf("min_set_size must be at least 1.")
  }
  if (is.na(opts$max_set_size) || opts$max_set_size < opts$min_set_size) {
    stopf("max_set_size (%d) must be >= min_set_size (%d).", opts$max_set_size, opts$min_set_size)
  }
  if (is.na(opts$low_expr_threshold) || opts$low_expr_threshold < 0) {
    stopf("low_expr_threshold must be >= 0.")
  }
  if (is.na(opts$low_expr_min_samples) || opts$low_expr_min_samples < 0) {
    stopf("low_expr_min_samples must be >= 0.")
  }

  valid_stats <- c("kruskal", "anova", "none")
  valid_pairwise <- c("auto", "wilcox", "t", "none")
  if (!opts$stats_test %in% valid_stats) {
    stopf("stats_test must be one of: %s", paste(valid_stats, collapse = ", "))
  }
  if (!opts$pairwise_test %in% valid_pairwise) {
    stopf("pairwise_test must be one of: %s", paste(valid_pairwise, collapse = ", "))
  }
  if (identical(opts$pairwise_test, "auto")) {
    opts$pairwise_test <- switch(
      opts$stats_test,
      kruskal = "wilcox",
      anova = "t",
      none = "none",
      "wilcox"
    )
  }

  tryCatch(
    stats::p.adjust(c(0.01, 0.02), method = opts$adjust_method),
    error = function(e) stopf("Invalid adjust_method '%s': %s", opts$adjust_method, conditionMessage(e))
  )

  opts
}

build_pathway_scoring_options <- function(input, genesets, ...) {
  opts <- pathway_scoring_defaults()
  opts$input <- input
  opts$genesets <- genesets
  overrides <- list(...)
  unknown <- setdiff(names(overrides), names(opts))
  if (length(unknown) > 0) {
    stopf("Unknown option override(s): %s", paste(unknown, collapse = ", "))
  }
  for (name in names(overrides)) {
    opts[[name]] <- overrides[[name]]
  }
  normalize_run_options(opts)
}

parse_options <- function() {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required for command-line use. Install it with install.packages('optparse').", call. = FALSE)
  }

  option_list <- list(
    optparse::make_option("--input", type = "character", help = "Input DESeq2-style CSV file."),
    optparse::make_option("--genesets", type = "character", help = "Gene-set file in GMT, CSV, or TSV format."),
    optparse::make_option("--gene-id-col", type = "character", default = "SYMBOL", dest = "gene_id_col",
                help = "Primary gene ID column to prefer before fallback columns [default: %default]."),
    optparse::make_option("--expr-prefix", type = "character", default = "norm_", dest = "expr_prefix",
                help = "Prefix used to detect normalized expression columns [default: %default]."),
    optparse::make_option("--group-pattern", type = "character", default = "^norm_([A-Za-z]+)_([0-9]+)$", dest = "group_pattern",
                help = "Regex with capture groups for sample group and replicate [default: %default]."),
    optparse::make_option("--outdir", type = "character", default = "pathway_scoring_results",
                help = "Output directory [default: %default]."),
    optparse::make_option("--methods", type = "character", default = "gsva,ssgsea,singscore",
                help = "Comma-separated methods to run [default: %default]."),
    optparse::make_option("--min-set-size", type = "integer", default = 5L, dest = "min_set_size",
                help = "Minimum matched genes per gene set after filtering [default: %default]."),
    optparse::make_option("--max-set-size", type = "integer", default = 500L, dest = "max_set_size",
                help = "Maximum matched genes per gene set after filtering [default: %default]."),
    optparse::make_option("--low-expr-threshold", type = "double", default = 1, dest = "low_expr_threshold",
                help = "Minimum normalized expression threshold used by the low-expression gene filter [default: %default]."),
    optparse::make_option("--low-expr-min-samples", type = "integer", default = 2L, dest = "low_expr_min_samples",
                help = "Minimum number of samples that must meet --low-expr-threshold [default: %default]."),
    optparse::make_option("--disable-low-expr-filter", action = "store_true", default = FALSE, dest = "disable_low_expr_filter",
                help = "Disable low-expression gene filtering before scoring."),
    optparse::make_option("--stats-test", type = "character", default = "kruskal", dest = "stats_test",
                help = "Omnibus test: kruskal, anova, or none [default: %default]."),
    optparse::make_option("--pairwise-test", type = "character", default = "auto", dest = "pairwise_test",
                help = "Pairwise test: auto, wilcox, t, or none [default: %default]."),
    optparse::make_option("--disable-pairwise-correction", action = "store_true", default = FALSE, dest = "disable_pairwise_correction",
                help = "Turn off multiple-testing correction for pairwise comparisons."),
    optparse::make_option("--adjust-method", type = "character", default = "BH", dest = "adjust_method",
                help = "P-value adjustment method passed to p.adjust() [default: %default]."),
    optparse::make_option("--plot-method-correlation", action = "store_true", default = FALSE, dest = "plot_method_correlation",
                help = "Also create a method-correlation heatmap when multiple methods are run."),
    optparse::make_option("--disable-boxplot-significance", action = "store_true", default = FALSE, dest = "disable_boxplot_significance",
                help = "Turn off significance brackets and labels on boxplots.")
  )

  parser <- optparse::OptionParser(
    usage = "Rscript %prog --input expression.csv --genesets curated_sets.csv --outdir results_dir [options]",
    option_list = option_list,
    description = "Single-sample pathway scoring for bulk RNA-seq DESeq2-style tables using GSVA, ssGSEA, and singscore."
  )

  opts <- optparse::parse_args(parser)
  if (is.null(opts$input) || !nzchar(opts$input) || is.null(opts$genesets) || !nzchar(opts$genesets)) {
    optparse::print_help(parser)
  }
  opts$pairwise_correction <- !isTRUE(opts$disable_pairwise_correction)
  opts$boxplot_significance <- !isTRUE(opts$disable_boxplot_significance)
  opts$disable_pairwise_correction <- NULL
  opts$disable_boxplot_significance <- NULL
  opts$help <- NULL
  normalize_run_options(opts)
}

create_output_dirs <- function(outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  plots_dir <- file.path(outdir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  }
  list(
    outdir = normalizePath(outdir, winslash = "/", mustWork = TRUE),
    plots_dir = normalizePath(plots_dir, winslash = "/", mustWork = TRUE)
  )
}

read_expression_data <- function(input_file, expr_prefix) {
  msg("Reading expression table: %s", input_file)
  df <- readr::read_csv(
    input_file,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )

  if (nrow(df) == 0) {
    stopf("Input table contains zero rows: %s", input_file)
  }

  expr_cols <- names(df)[startsWith(names(df), expr_prefix)]
  if (length(expr_cols) == 0) {
    stopf("No normalized expression columns found with prefix '%s'.", expr_prefix)
  }

  raw_expr <- as.data.frame(df[, expr_cols, drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  raw_expr_char <- as.matrix(data.frame(lapply(raw_expr, as.character), check.names = FALSE, stringsAsFactors = FALSE))
  expr_numeric <- as.data.frame(
    lapply(raw_expr, function(col) suppressWarnings(as.numeric(col))),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expr_numeric_mat <- as.matrix(expr_numeric)
  bad_mask <- is.na(expr_numeric_mat) & !is.na(raw_expr_char) & nzchar(trimws(raw_expr_char))
  if (sum(bad_mask) > 0) {
    bad_cols <- unique(col(expr_numeric_mat)[bad_mask])
    bad_col_names <- names(expr_numeric)[bad_cols]
    stopf(
      "Expression columns contain non-numeric values after coercion. Example affected columns: %s",
      paste(utils::head(bad_col_names, 10), collapse = ", ")
    )
  }

  na_count <- sum(is.na(as.matrix(expr_numeric)))
  if (na_count > 0) {
    stopf("Expression matrix contains %d missing value(s). Please fill or remove missing normalized values before scoring.", na_count)
  }

  df[, expr_cols] <- tibble::as_tibble(expr_numeric)
  msg("Detected %d normalized sample column(s).", length(expr_cols))

  list(
    data = tibble::as_tibble(df),
    expr_cols = expr_cols
  )
}

resolve_gene_ids <- function(df, expr_cols, primary_gene_id_col = "SYMBOL") {
  fallback_cols <- c("SYMBOL", "Geneid_clean", "Ensembl_GeneID", "Geneid")
  if (!identical(primary_gene_id_col, "auto") && !primary_gene_id_col %in% names(df)) {
    stopf("Requested --gene-id-col '%s' was not found in the input table.", primary_gene_id_col)
  }

  id_priority <- unique(c(primary_gene_id_col, fallback_cols))
  id_priority <- id_priority[id_priority %in% names(df)]

  if (length(id_priority) == 0) {
    stopf("No usable gene ID columns were found. Looked for: %s", paste(unique(c(primary_gene_id_col, fallback_cols)), collapse = ", "))
  }

  msg("Gene ID priority: %s", paste(id_priority, collapse = " -> "))

  resolved_id <- rep(NA_character_, nrow(df))
  resolved_source <- rep(NA_character_, nrow(df))
  for (col_name in id_priority) {
    values <- trim_to_char(df[[col_name]])
    take <- is.na(resolved_id) & is_non_empty(values)
    resolved_id[take] <- values[take]
    resolved_source[take] <- col_name
  }

  df$.resolved_gene_id <- resolved_id
  df$.resolved_gene_id_source <- resolved_source
  df$.row_index <- seq_len(nrow(df))
  df$.mean_normalized_expr <- rowMeans(as.matrix(df[, expr_cols, drop = FALSE]), na.rm = TRUE)
  df$.mean_normalized_expr[!is.finite(df$.mean_normalized_expr)] <- -Inf

  missing_gene_id_rows <- sum(!is_non_empty(df$.resolved_gene_id))
  with_ids <- df[is_non_empty(df$.resolved_gene_id), , drop = FALSE]
  duplicate_rows_before_collapse <- sum(duplicated(with_ids$.resolved_gene_id))

  with_ids <- with_ids[order(-with_ids$.mean_normalized_expr, with_ids$.row_index), , drop = FALSE]
  collapsed <- with_ids[!duplicated(with_ids$.resolved_gene_id), , drop = FALSE]
  collapsed <- collapsed[order(collapsed$.row_index), , drop = FALSE]

  summary_overall <- tibble::tibble(
    summary_type = "overall",
    item = c(
      "input_rows",
      "rows_missing_resolved_gene_id",
      "rows_after_gene_id_filter",
      "duplicate_rows_before_collapse",
      "rows_removed_by_duplicate_collapse",
      "unique_gene_ids_after_collapse"
    ),
    value = c(
      nrow(df),
      missing_gene_id_rows,
      nrow(with_ids),
      duplicate_rows_before_collapse,
      nrow(with_ids) - nrow(collapsed),
      nrow(collapsed)
    )
  )

  source_counts <- tibble::tibble(
    summary_type = "gene_id_source",
    item = id_priority,
    value = vapply(
      id_priority,
      function(col_name) sum(collapsed$.resolved_gene_id_source == col_name, na.rm = TRUE),
      FUN.VALUE = numeric(1)
    )
  )

  list(
    data = tibble::as_tibble(collapsed),
    summary = dplyr::bind_rows(summary_overall, source_counts),
    id_priority = id_priority
  )
}

filter_expression_matrix <- function(df, expr_cols, low_expr_threshold, low_expr_min_samples, disable_low_expr_filter) {
  expr_mat <- as.matrix(df[, expr_cols, drop = FALSE])
  storage.mode(expr_mat) <- "double"
  rownames(expr_mat) <- df$.resolved_gene_id

  if (any(expr_mat <= -1, na.rm = TRUE)) {
    stopf("Expression matrix contains value(s) <= -1, so log2(x + 1) cannot be computed safely.")
  }

  keep_low_expr <- rep(TRUE, nrow(expr_mat))
  if (!disable_low_expr_filter && low_expr_min_samples > 0) {
    keep_low_expr <- rowSums(expr_mat >= low_expr_threshold) >= low_expr_min_samples
    msg(
      "Low-expression filter kept %d of %d genes (threshold >= %.4f in at least %d sample(s)).",
      sum(keep_low_expr), nrow(expr_mat), low_expr_threshold, low_expr_min_samples
    )
  } else {
    msg("Low-expression filter disabled.")
  }

  expr_mat <- expr_mat[keep_low_expr, , drop = FALSE]
  if (nrow(expr_mat) == 0) {
    stopf("No genes remain after low-expression filtering.")
  }

  expr_log <- log2(expr_mat + 1)
  row_variance <- apply(expr_log, 1, stats::var)
  keep_variance <- is.finite(row_variance) & row_variance > 0

  removed_zero_variance <- sum(!keep_variance)
  expr_mat <- expr_mat[keep_variance, , drop = FALSE]
  expr_log <- expr_log[keep_variance, , drop = FALSE]

  if (nrow(expr_log) == 0) {
    stopf("No genes remain after zero-variance filtering.")
  }

  filter_summary <- tibble::tibble(
    summary_type = "expression_filtering",
    item = c(
      "genes_before_expression_filtering",
      "genes_after_low_expression_filter",
      "genes_removed_by_low_expression_filter",
      "genes_removed_by_zero_variance_filter",
      "genes_after_all_filters"
    ),
    value = c(
      length(keep_low_expr),
      sum(keep_low_expr),
      length(keep_low_expr) - sum(keep_low_expr),
      removed_zero_variance,
      nrow(expr_log)
    )
  )

  list(
    expression_matrix = expr_mat,
    log_expression_matrix = expr_log,
    summary = filter_summary
  )
}

parse_sample_metadata <- function(sample_names, group_pattern) {
  if (anyDuplicated(sample_names) > 0) {
    dupes <- unique(sample_names[duplicated(sample_names)])
    stopf("Duplicate normalized sample names found: %s", paste(dupes, collapse = ", "))
  }

  matches <- stringr::str_match(sample_names, group_pattern)
  if (ncol(matches) < 2) {
    stopf("--group-pattern must contain at least one capture group for sample group.")
  }

  unmatched <- is.na(matches[, 1])
  if (any(unmatched)) {
    stopf(
      "The following sample names did not match --group-pattern '%s': %s",
      group_pattern,
      paste(sample_names[unmatched], collapse = ", ")
    )
  }

  group <- matches[, 2]
  replicate <- rep(NA_integer_, length(sample_names))
  if (ncol(matches) >= 3) {
    replicate_raw <- matches[, 3]
    replicate_num <- suppressWarnings(as.integer(replicate_raw))
    bad_reps <- is.na(replicate_num) & !is.na(replicate_raw)
    if (any(bad_reps)) {
      stopf("Replicate capture group from --group-pattern must be numeric. Offending sample(s): %s", paste(sample_names[bad_reps], collapse = ", "))
    }
    replicate <- replicate_num
  }

  tibble::tibble(
    sample = sample_names,
    group = factor(group, levels = unique(group)),
    replicate = replicate,
    sample_order = seq_along(sample_names)
  )
}

read_gene_sets <- function(geneset_file) {
  ext <- tolower(tools::file_ext(geneset_file))
  msg("Reading gene sets: %s", geneset_file)

  if (ext == "gmt") {
    gsc <- GSEABase::getGmt(geneset_file)
    if (length(gsc) == 0) {
      stopf("GMT file contains zero gene sets: %s", geneset_file)
    }
    rows <- purrr::imap_dfr(gsc, function(gs, idx) {
      gene_ids <- GSEABase::geneIds(gs)
      tibble::tibble(
        set_name = GSEABase::setName(gs),
        gene_symbol = gene_ids,
        direction = NA_character_,
        description = GSEABase::shortDescription(gs),
        pathway_order = idx,
        row_order = seq_along(gene_ids)
      )
    })
    source_format <- "gmt"
  } else {
    delim <- if (grepl("\\.(tsv|tab|txt)$", tolower(geneset_file))) "\t" else ","
    reader <- if (identical(delim, "\t")) readr::read_tsv else readr::read_csv
    raw_tbl <- reader(
      geneset_file,
      show_col_types = FALSE,
      progress = FALSE,
      name_repair = "minimal"
    )
    if (nrow(raw_tbl) == 0) {
      stopf("Gene-set file contains zero rows: %s", geneset_file)
    }

    set_name_col <- find_ci_col(names(raw_tbl), "set_name")
    gene_symbol_col <- find_ci_col(names(raw_tbl), "gene_symbol")
    direction_col <- find_ci_col(names(raw_tbl), "direction")
    description_col <- find_ci_col(names(raw_tbl), "description")

    if (is.na(set_name_col) || is.na(gene_symbol_col)) {
      stopf("CSV/TSV gene-set file must contain columns named set_name and gene_symbol.")
    }

    rows <- tibble::tibble(
      set_name = raw_tbl[[set_name_col]],
      gene_symbol = raw_tbl[[gene_symbol_col]],
      direction = if (!is.na(direction_col)) raw_tbl[[direction_col]] else NA_character_,
      description = if (!is.na(description_col)) raw_tbl[[description_col]] else NA_character_,
      pathway_order = match(raw_tbl[[set_name_col]], unique(raw_tbl[[set_name_col]])),
      row_order = seq_len(nrow(raw_tbl))
    )
    source_format <- if (identical(delim, "\t")) "tsv" else "csv"
  }

  rows$set_name <- trim_to_char(rows$set_name)
  rows$gene_symbol <- trim_to_char(rows$gene_symbol)
  rows$description <- trim_to_char(rows$description)
  rows$description[!nzchar(rows$description)] <- NA_character_
  rows$direction <- tolower(trim_to_char(rows$direction))
  rows$direction[!nzchar(rows$direction)] <- NA_character_

  if (any(!is_non_empty(rows$set_name))) {
    stopf("Gene-set file contains blank set_name entries.")
  }

  blank_gene_rows <- sum(!is_non_empty(rows$gene_symbol))
  if (blank_gene_rows > 0) {
    warnf("Removing %d gene-set row(s) with blank gene symbols.", blank_gene_rows)
    rows <- rows[is_non_empty(rows$gene_symbol), , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    stopf("No valid gene-set rows remain after removing blank gene symbols.")
  }

  invalid_directions <- setdiff(stats::na.omit(unique(rows$direction)), c("up", "down"))
  if (length(invalid_directions) > 0) {
    stopf(
      "Invalid values found in the direction column: %s. Allowed values are 'up' and 'down'.",
      paste(invalid_directions, collapse = ", ")
    )
  }

  pathway_order <- unique(rows$set_name)
  descriptions <- rows |>
    dplyr::group_by(set_name) |>
    dplyr::summarise(description = safe_first_non_na(description), .groups = "drop")
  description_map <- stats::setNames(descriptions$description, descriptions$set_name)

  aggregated_sets <- stats::setNames(vector("list", length(pathway_order)), pathway_order)
  singscore_input <- stats::setNames(vector("list", length(pathway_order)), pathway_order)

  for (set_name in pathway_order) {
    subset_tbl <- rows[rows$set_name == set_name, , drop = FALSE]
    aggregated_sets[[set_name]] <- unique_keep_order(subset_tbl$gene_symbol)

    has_direction <- any(!is.na(subset_tbl$direction))
    if (!has_direction) {
      singscore_input[[set_name]] <- list(
        set_name = set_name,
        input_mode = "unsigned",
        up = unique_keep_order(subset_tbl$gene_symbol),
        down = character(0)
      )
    } else {
      singscore_input[[set_name]] <- list(
        set_name = set_name,
        input_mode = "directional",
        up = unique_keep_order(subset_tbl$gene_symbol[subset_tbl$direction == "up"]),
        down = unique_keep_order(subset_tbl$gene_symbol[subset_tbl$direction == "down"])
      )
    }
  }

  msg("Loaded %d pathway(s) from %s gene-set input.", length(pathway_order), toupper(source_format))

  list(
    source_format = source_format,
    rows = tibble::as_tibble(rows),
    pathway_order = pathway_order,
    descriptions = description_map,
    aggregated_sets = aggregated_sets,
    singscore_input = singscore_input
  )
}

prepare_gene_sets_for_scoring <- function(gene_set_data, expr_gene_ids, min_set_size, max_set_size) {
  expr_gene_ids <- unique(expr_gene_ids)

  coverage_rows <- vector("list", length(gene_set_data$pathway_order))
  gsva_sets <- list()
  singscore_specs <- list()

  for (i in seq_along(gene_set_data$pathway_order)) {
    set_name <- gene_set_data$pathway_order[[i]]
    all_genes <- unique_keep_order(gene_set_data$aggregated_sets[[set_name]])
    matched_union <- all_genes[all_genes %in% expr_gene_ids]
    unmatched_union <- all_genes[!all_genes %in% expr_gene_ids]

    gsva_reason <- "included"
    if (length(matched_union) < min_set_size) {
      gsva_reason <- sprintf("excluded: matched genes below min-set-size (%d < %d)", length(matched_union), min_set_size)
    } else if (length(matched_union) > max_set_size) {
      gsva_reason <- sprintf("excluded: matched genes above max-set-size (%d > %d)", length(matched_union), max_set_size)
    } else {
      gsva_sets[[set_name]] <- matched_union
    }

    score_input <- gene_set_data$singscore_input[[set_name]]
    up_input <- unique_keep_order(score_input$up)
    down_input <- unique_keep_order(score_input$down)
    up_matched <- up_input[up_input %in% expr_gene_ids]
    down_matched <- down_input[down_input %in% expr_gene_ids]

    up_pass <- length(up_matched) >= min_set_size && length(up_matched) <= max_set_size
    down_pass <- length(down_matched) >= min_set_size && length(down_matched) <= max_set_size

    singscore_mode <- "excluded"
    if (identical(score_input$input_mode, "unsigned")) {
      if (up_pass) {
        singscore_specs[[set_name]] <- list(
          set_name = set_name,
          mode = "unsigned",
          up = up_matched,
          down = character(0)
        )
        singscore_mode <- "unsigned"
      }
    } else {
      if (up_pass && down_pass) {
        singscore_specs[[set_name]] <- list(
          set_name = set_name,
          mode = "paired_up_down",
          up = up_matched,
          down = down_matched
        )
        singscore_mode <- "paired_up_down"
      } else if (up_pass) {
        singscore_specs[[set_name]] <- list(
          set_name = set_name,
          mode = "up_only",
          up = up_matched,
          down = character(0)
        )
        singscore_mode <- if (length(down_input) > 0) "up_only_after_filter" else "up_only"
      } else if (down_pass) {
        singscore_specs[[set_name]] <- list(
          set_name = set_name,
          mode = "down_only",
          up = character(0),
          down = down_matched
        )
        singscore_mode <- if (length(up_input) > 0) "down_only_after_filter" else "down_only"
      }
    }

    coverage_rows[[i]] <- tibble::tibble(
      set_name = set_name,
      description = gene_set_data$descriptions[[set_name]],
      input_mode = score_input$input_mode,
      input_gene_set_size = length(all_genes),
      matched_gene_count = length(matched_union),
      unmatched_gene_count = length(unmatched_union),
      percent_coverage = if (length(all_genes) == 0) 0 else round((length(matched_union) / length(all_genes)) * 100, 2),
      matched_genes = paste(matched_union, collapse = ";"),
      unmatched_genes = paste(unmatched_union, collapse = ";"),
      gsva_ssgsea_status = gsva_reason,
      singscore_mode = singscore_mode,
      up_input_size = length(up_input),
      up_matched_gene_count = length(up_matched),
      down_input_size = length(down_input),
      down_matched_gene_count = length(down_matched)
    )
  }

  coverage_tbl <- dplyr::bind_rows(coverage_rows)

  list(
    coverage = coverage_tbl,
    gsva_sets = gsva_sets,
    singscore_specs = singscore_specs
  )
}

score_with_gsva <- function(expr_log, gene_sets, method_name, min_set_size, max_set_size) {
  if (!method_name %in% c("gsva", "ssgsea")) {
    stopf("Unsupported GSVA-family method requested: %s", method_name)
  }
  if (length(gene_sets) == 0) {
    stopf("No valid gene sets are available for %s scoring after filtering.", toupper(method_name))
  }

  msg("Running %s on %d pathway(s) across %d sample(s).", toupper(method_name), length(gene_sets), ncol(expr_log))

  if (identical(method_name, "gsva")) {
    param <- GSVA::gsvaParam(
      exprData = expr_log,
      geneSets = gene_sets,
      minSize = min_set_size,
      maxSize = max_set_size,
      kcdf = "Gaussian",
      sparse = FALSE,
      checkNA = "no",
      filterRows = TRUE,
      ondisk = "no",
      verbose = FALSE
    )
  } else {
    param <- GSVA::ssgseaParam(
      exprData = expr_log,
      geneSets = gene_sets,
      minSize = min_set_size,
      maxSize = max_set_size,
      alpha = 0.25,
      normalize = TRUE,
      checkNA = "no",
      verbose = FALSE
    )
  }

  scores <- GSVA::gsva(param)
  scores <- as.matrix(scores)

  expected_rows <- names(gene_sets)
  missing_rows <- setdiff(expected_rows, rownames(scores))
  if (length(missing_rows) > 0) {
    stopf("GSVA output is missing expected pathway row(s): %s", paste(utils::head(missing_rows, 10), collapse = ", "))
  }

  scores <- scores[expected_rows, colnames(expr_log), drop = FALSE]
  scores
}

detect_singscore_score_column <- function(score_df) {
  preferred <- c("TotalScore", "Total.Score", "score", "Score")
  present <- preferred[preferred %in% names(score_df)]
  if (length(present) > 0) {
    return(present[[1]])
  }

  numeric_cols <- names(score_df)[vapply(score_df, is.numeric, logical(1))]
  numeric_score_cols <- numeric_cols[grepl("score", numeric_cols, ignore.case = TRUE) & !grepl("dispersion", numeric_cols, ignore.case = TRUE)]
  if (length(numeric_score_cols) > 0) {
    return(numeric_score_cols[[1]])
  }

  if (length(numeric_cols) > 0) {
    return(numeric_cols[[1]])
  }

  stopf("Could not identify a numeric singscore output column from columns: %s", paste(names(score_df), collapse = ", "))
}

extract_singscore_vector <- function(score_obj, sample_names) {
  score_df <- as.data.frame(score_obj, check.names = FALSE, stringsAsFactors = FALSE)

  if (!all(sample_names %in% rownames(score_df))) {
    sample_col <- names(score_df)[tolower(names(score_df)) %in% c("sample", "samples")]
    if (length(sample_col) > 0) {
      rownames(score_df) <- as.character(score_df[[sample_col[[1]]]])
    }
  }

  if (!all(sample_names %in% rownames(score_df))) {
    stopf("Could not align singscore output back to the input sample names.")
  }

  score_col <- detect_singscore_score_column(score_df)
  as.numeric(score_df[sample_names, score_col])
}

score_with_singscore <- function(expr_log, singscore_specs) {
  if (length(singscore_specs) == 0) {
    stopf("No valid gene sets are available for singscore after filtering.")
  }

  msg("Running singscore on %d pathway(s) across %d sample(s).", length(singscore_specs), ncol(expr_log))
  rank_data <- singscore::rankGenes(expr_log)
  sample_names <- colnames(expr_log)
  score_mat <- matrix(NA_real_, nrow = length(singscore_specs), ncol = length(sample_names))
  rownames(score_mat) <- names(singscore_specs)
  colnames(score_mat) <- sample_names

  for (set_name in names(singscore_specs)) {
    spec <- singscore_specs[[set_name]]
    if (identical(spec$mode, "paired_up_down")) {
      score_obj <- singscore::simpleScore(rankData = rank_data, upSet = spec$up, downSet = spec$down)
    } else if (identical(spec$mode, "unsigned") || identical(spec$mode, "up_only")) {
      score_obj <- singscore::simpleScore(rankData = rank_data, upSet = spec$up)
    } else if (identical(spec$mode, "down_only")) {
      score_obj <- tryCatch(
        singscore::simpleScore(rankData = rank_data, downSet = spec$down),
        error = function(e) {
          warnf(
            "Pathway '%s': singscore down-only scoring failed (%s). Falling back to negated one-sided scoring of the down gene set.",
            set_name, conditionMessage(e)
          )
          fallback <- singscore::simpleScore(rankData = rank_data, upSet = spec$down)
          fallback_df <- as.data.frame(fallback, check.names = FALSE, stringsAsFactors = FALSE)
          score_col <- detect_singscore_score_column(fallback_df)
          fallback_df[[score_col]] <- -1 * fallback_df[[score_col]]
          fallback_df
        }
      )
    } else {
      stopf("Unsupported singscore mode for pathway '%s': %s", set_name, spec$mode)
    }

    score_mat[set_name, ] <- extract_singscore_vector(score_obj, sample_names)
  }

  score_mat
}

score_matrix_to_long <- function(score_mat, sample_meta, method_name) {
  score_df <- as.data.frame(score_mat, check.names = FALSE, stringsAsFactors = FALSE)
  score_df$pathway <- rownames(score_mat)
  score_df <- score_df[, c("pathway", colnames(score_mat)), drop = FALSE]
  long_tbl <- tidyr::pivot_longer(
    tibble::as_tibble(score_df),
    cols = dplyr::all_of(setdiff(colnames(score_df), "pathway")),
    names_to = "sample",
    values_to = "score"
  )
  long_tbl <- dplyr::left_join(long_tbl, sample_meta, by = "sample")
  long_tbl$method <- method_name
  long_tbl$pathway <- factor(long_tbl$pathway, levels = rownames(score_mat))
  long_tbl$group <- factor(long_tbl$group, levels = levels(sample_meta$group))
  long_tbl
}

empty_stats_table <- function() {
  tibble::tibble(
    pathway = character(),
    test_scope = character(),
    test_method = character(),
    group1 = character(),
    group2 = character(),
    n_group1 = integer(),
    n_group2 = integer(),
    statistic = double(),
    parameter = double(),
    estimate = double(),
    effect_type = character(),
    p_value = double(),
    p_adj = double(),
    adjust_method = character(),
    note = character()
  )
}

adjust_pvalues <- function(tbl, p_col = "p_value", new_col = "p_adj", method = "BH") {
  if (nrow(tbl) == 0) {
    tbl[[new_col]] <- numeric(0)
    return(tbl)
  }
  pvals <- tbl[[p_col]]
  adjusted <- rep(NA_real_, length(pvals))
  keep <- !is.na(pvals)
  if (any(keep)) {
    adjusted[keep] <- stats::p.adjust(pvals[keep], method = method)
  }
  tbl[[new_col]] <- adjusted
  tbl
}

run_statistics <- function(score_long, stats_test = "kruskal", pairwise_test = "wilcox", pairwise_correction = TRUE, adjust_method = "BH") {
  if (identical(stats_test, "none")) {
    return(empty_stats_table())
  }

  pathways <- levels(score_long$pathway)
  group_levels <- levels(score_long$group)
  omnibus_rows <- vector("list", length(pathways))
  pairwise_rows <- list()

  for (i in seq_along(pathways)) {
    pathway_name <- pathways[[i]]
    pathway_df <- score_long[score_long$pathway == pathway_name, , drop = FALSE]
    pathway_df <- pathway_df[!is.na(pathway_df$score) & !is.na(pathway_df$group), , drop = FALSE]
    present_groups <- unique(as.character(pathway_df$group))

    if (length(present_groups) < 2) {
      omnibus_rows[[i]] <- tibble::tibble(
        pathway = pathway_name,
        test_scope = "omnibus",
        test_method = stats_test,
        group1 = NA_character_,
        group2 = NA_character_,
        n_group1 = NA_integer_,
        n_group2 = NA_integer_,
        statistic = NA_real_,
        parameter = NA_real_,
        estimate = NA_real_,
        effect_type = NA_character_,
        p_value = NA_real_,
        adjust_method = adjust_method,
        note = "Fewer than two groups contained non-missing scores."
      )
      next
    }

    if (identical(stats_test, "kruskal")) {
      fit <- stats::kruskal.test(score ~ group, data = pathway_df)
      omnibus_rows[[i]] <- tibble::tibble(
        pathway = pathway_name,
        test_scope = "omnibus",
        test_method = "kruskal",
        group1 = NA_character_,
        group2 = NA_character_,
        n_group1 = NA_integer_,
        n_group2 = NA_integer_,
        statistic = unname(fit$statistic),
        parameter = unname(fit$parameter),
        estimate = NA_real_,
        effect_type = NA_character_,
        p_value = fit$p.value,
        adjust_method = adjust_method,
        note = NA_character_
      )
    } else {
      fit <- stats::aov(score ~ group, data = pathway_df)
      anova_tbl <- summary(fit)[[1]]
      omnibus_rows[[i]] <- tibble::tibble(
        pathway = pathway_name,
        test_scope = "omnibus",
        test_method = "anova",
        group1 = NA_character_,
        group2 = NA_character_,
        n_group1 = NA_integer_,
        n_group2 = NA_integer_,
        statistic = unname(anova_tbl[1, "F value"]),
        parameter = unname(anova_tbl[1, "Df"]),
        estimate = NA_real_,
        effect_type = NA_character_,
        p_value = unname(anova_tbl[1, "Pr(>F)"]),
        adjust_method = adjust_method,
        note = NA_character_
      )
    }

    if (!identical(pairwise_test, "none")) {
      valid_groups <- group_levels[group_levels %in% present_groups]
      if (length(valid_groups) >= 2) {
        pairs <- utils::combn(valid_groups, 2, simplify = FALSE)
        for (pair in pairs) {
          g1 <- pair[[1]]
          g2 <- pair[[2]]
          sub_df <- pathway_df[pathway_df$group %in% c(g1, g2), , drop = FALSE]
          sub_df$group <- droplevels(sub_df$group)
          scores_g1 <- sub_df$score[sub_df$group == g1]
          scores_g2 <- sub_df$score[sub_df$group == g2]

          pair_result <- tryCatch(
            {
              if (identical(pairwise_test, "wilcox")) {
                fit_pw <- suppressWarnings(stats::wilcox.test(score ~ group, data = sub_df))
                tibble::tibble(
                  pathway = pathway_name,
                  test_scope = "pairwise",
                  test_method = "wilcox",
                  group1 = g1,
                  group2 = g2,
                  n_group1 = length(scores_g1),
                  n_group2 = length(scores_g2),
                  statistic = unname(fit_pw$statistic),
                  parameter = NA_real_,
                  estimate = stats::median(scores_g1) - stats::median(scores_g2),
                  effect_type = "median_difference",
                  p_value = fit_pw$p.value,
                  adjust_method = if (isTRUE(pairwise_correction)) adjust_method else "none",
                  note = NA_character_
                )
              } else {
                fit_pw <- stats::t.test(score ~ group, data = sub_df)
                tibble::tibble(
                  pathway = pathway_name,
                  test_scope = "pairwise",
                  test_method = "t",
                  group1 = g1,
                  group2 = g2,
                  n_group1 = length(scores_g1),
                  n_group2 = length(scores_g2),
                  statistic = unname(fit_pw$statistic),
                  parameter = unname(fit_pw$parameter),
                  estimate = mean(scores_g1) - mean(scores_g2),
                  effect_type = "mean_difference",
                  p_value = fit_pw$p.value,
                  adjust_method = if (isTRUE(pairwise_correction)) adjust_method else "none",
                  note = NA_character_
                )
              }
            },
            error = function(e) {
              tibble::tibble(
                pathway = pathway_name,
                test_scope = "pairwise",
                test_method = pairwise_test,
                group1 = g1,
                group2 = g2,
                n_group1 = length(scores_g1),
                n_group2 = length(scores_g2),
                statistic = NA_real_,
                parameter = NA_real_,
                estimate = NA_real_,
                effect_type = if (identical(pairwise_test, "wilcox")) "median_difference" else "mean_difference",
                p_value = NA_real_,
                adjust_method = if (isTRUE(pairwise_correction)) adjust_method else "none",
                note = conditionMessage(e)
              )
            }
          )
          pairwise_rows[[length(pairwise_rows) + 1]] <- pair_result
        }
      }
    }
  }

  omnibus_tbl <- dplyr::bind_rows(omnibus_rows)
  omnibus_tbl <- adjust_pvalues(omnibus_tbl, method = adjust_method)
  pairwise_tbl <- if (length(pairwise_rows) > 0) dplyr::bind_rows(pairwise_rows) else empty_stats_table()
  if (nrow(pairwise_tbl) > 0) {
    if (isTRUE(pairwise_correction)) {
      pairwise_tbl <- adjust_pvalues(pairwise_tbl, method = adjust_method)
    } else {
      pairwise_tbl$p_adj <- pairwise_tbl$p_value
      pairwise_tbl$adjust_method <- "none"
    }
  }

  dplyr::bind_rows(omnibus_tbl, pairwise_tbl)
}

write_csv_safe <- function(tbl, output_file) {
  readr::write_csv(tbl, output_file, na = "")
}

write_score_matrix <- function(score_mat, output_file) {
  score_df <- as.data.frame(score_mat, check.names = FALSE, stringsAsFactors = FALSE)
  score_df <- tibble::rownames_to_column(score_df, var = "pathway")
  readr::write_csv(score_df, output_file, na = "")
}

build_group_palette <- function(groups) {
  groups <- unique(as.character(groups))
  palette_vals <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
  if (length(groups) > length(palette_vals)) {
    palette_vals <- grDevices::colorRampPalette(palette_vals)(length(groups))
  } else {
    palette_vals <- palette_vals[seq_along(groups)]
  }
  stats::setNames(palette_vals, groups)
}

p_value_to_significance_label <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  if (p_value <= 1e-4) {
    return("****")
  }
  if (p_value <= 1e-3) {
    return("***")
  }
  if (p_value <= 1e-2) {
    return("**")
  }
  if (p_value <= 5e-2) {
    return("*")
  }
  "ns"
}

build_boxplot_annotations <- function(score_long, stats_tbl) {
  if (is.null(stats_tbl) || nrow(stats_tbl) == 0) {
    return(NULL)
  }

  pairwise_tbl <- stats_tbl[
    stats_tbl$test_scope == "pairwise" &
      !is.na(stats_tbl$group1) &
      !is.na(stats_tbl$group2),
    ,
    drop = FALSE
  ]
  if (nrow(pairwise_tbl) == 0) {
    return(NULL)
  }

  group_levels <- levels(score_long$group)
  x_map <- stats::setNames(seq_along(group_levels), group_levels)
  pathway_levels <- levels(score_long$pathway)
  annotations <- vector("list", length(pathway_levels))

  for (i in seq_along(pathway_levels)) {
    pathway_name <- pathway_levels[[i]]
    pathway_scores <- score_long[score_long$pathway == pathway_name, , drop = FALSE]
    pathway_pairs <- pairwise_tbl[pairwise_tbl$pathway == pathway_name, , drop = FALSE]
    if (nrow(pathway_scores) == 0 || nrow(pathway_pairs) == 0) {
      next
    }

    score_min <- min(pathway_scores$score, na.rm = TRUE)
    score_max <- max(pathway_scores$score, na.rm = TRUE)
    score_range <- score_max - score_min
    if (!is.finite(score_range) || score_range <= 0) {
      score_range <- max(abs(score_max), 1)
    }

    pathway_pairs$xmin <- x_map[as.character(pathway_pairs$group1)]
    pathway_pairs$xmax <- x_map[as.character(pathway_pairs$group2)]
    pathway_pairs$xleft <- pmin(pathway_pairs$xmin, pathway_pairs$xmax)
    pathway_pairs$xright <- pmax(pathway_pairs$xmin, pathway_pairs$xmax)
    pathway_pairs$comparison_width <- pathway_pairs$xright - pathway_pairs$xleft
    pathway_pairs <- pathway_pairs[order(pathway_pairs$comparison_width, pathway_pairs$xleft), , drop = FALSE]

    base_offset <- score_range * 0.12
    step_offset <- score_range * 0.16
    tip_length <- score_range * 0.05
    text_offset <- score_range * 0.03

    pathway_pairs$pathway <- factor(pathway_name, levels = pathway_levels)
    pathway_pairs$y_position <- score_max + base_offset + (seq_len(nrow(pathway_pairs)) - 1) * step_offset
    pathway_pairs$y_tip <- pathway_pairs$y_position - tip_length
    pathway_pairs$label_y <- pathway_pairs$y_position + text_offset
    pathway_pairs$label <- vapply(
      ifelse(is.na(pathway_pairs$p_adj), pathway_pairs$p_value, pathway_pairs$p_adj),
      p_value_to_significance_label,
      FUN.VALUE = character(1)
    )

    annotations[[i]] <- pathway_pairs
  }

  annot_tbl <- dplyr::bind_rows(annotations)
  if (nrow(annot_tbl) == 0) {
    return(NULL)
  }

  annot_tbl
}

plot_heatmaps <- function(score_mat, sample_meta, method_name, output_file) {
  annotation_col <- as.data.frame(sample_meta[, c("group"), drop = FALSE], stringsAsFactors = FALSE)
  rownames(annotation_col) <- sample_meta$sample
  score_mat <- score_mat[, sample_meta$sample, drop = FALSE]

  scaled_mat <- t(scale(t(score_mat)))
  scaled_mat[!is.finite(scaled_mat)] <- 0

  heat_colors <- grDevices::colorRampPalette(c("#053061", "#F7F7F7", "#67001F"))(100)

  grDevices::pdf(
    output_file,
    width = max(8, ncol(score_mat) * 0.7),
    height = max(6, nrow(score_mat) * 0.35)
  )

  ph_raw <- pheatmap::pheatmap(
    score_mat,
    color = heat_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = sprintf("%s pathway scores", toupper(method_name)),
    silent = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(ph_raw$gtable)

  ph_scaled <- pheatmap::pheatmap(
    scaled_mat,
    color = heat_colors,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    border_color = NA,
    main = sprintf("%s pathway scores (row-scaled)", toupper(method_name)),
    silent = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(ph_scaled$gtable)

  grDevices::dev.off()
}

plot_boxplots <- function(score_long, method_name, output_file, stats_tbl = NULL, show_significance = TRUE) {
  n_pathways <- dplyr::n_distinct(score_long$pathway)
  ncol_wrap <- min(4, max(1, ceiling(sqrt(n_pathways))))
  nrow_wrap <- ceiling(n_pathways / ncol_wrap)
  group_palette <- build_group_palette(levels(score_long$group))
  annot_tbl <- if (isTRUE(show_significance)) build_boxplot_annotations(score_long, stats_tbl) else NULL

  p <- ggplot2::ggplot(score_long, ggplot2::aes(x = group, y = score, fill = group, color = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.65) +
    ggplot2::geom_jitter(width = 0.12, size = 2.3, alpha = 0.9) +
    ggplot2::facet_wrap(stats::as.formula("~ pathway"), scales = "free_y", ncol = ncol_wrap) +
    ggplot2::scale_fill_manual(values = group_palette, drop = FALSE) +
    ggplot2::scale_color_manual(values = group_palette, drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s pathway scores by group", toupper(method_name)),
      x = "Group",
      y = "Pathway score"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      strip.background = ggplot2::element_rect(fill = "grey95"),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (!is.null(annot_tbl) && nrow(annot_tbl) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = annot_tbl,
        ggplot2::aes(x = xleft, xend = xright, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = annot_tbl,
        ggplot2::aes(x = xleft, xend = xleft, y = y_tip, yend = y_position),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_segment(
        data = annot_tbl,
        ggplot2::aes(x = xright, xend = xright, y = y_tip, yend = y_position),
        inherit.aes = FALSE,
        color = "black",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        data = annot_tbl,
        ggplot2::aes(x = (xleft + xright) / 2, y = label_y, label = label),
        inherit.aes = FALSE,
        color = "black",
        size = 3.3
      )
  }

  ggplot2::ggsave(
    filename = output_file,
    plot = p,
    width = max(9, 4 * ncol_wrap),
    height = max(6, 3.4 * nrow_wrap),
    units = "in",
    limitsize = FALSE
  )
}

plot_pca <- function(score_mat, sample_meta, method_name, output_file) {
  group_palette <- build_group_palette(levels(sample_meta$group))
  score_mat <- score_mat[, sample_meta$sample, drop = FALSE]

  if (ncol(score_mat) < 2 || nrow(score_mat) < 2) {
    note_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "PCA requires at least 2 samples and 2 pathways after filtering.",
        size = 4.2
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(sprintf("%s PCA", toupper(method_name)))
    ggplot2::ggsave(output_file, note_plot, width = 7, height = 5, units = "in")
    return(invisible(NULL))
  }

  pca <- stats::prcomp(t(score_mat), center = TRUE, scale. = TRUE)
  if (ncol(pca$x) < 2) {
    note_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text",
        x = 0,
        y = 0,
        label = "PCA returned fewer than 2 principal components.",
        size = 4.2
      ) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(sprintf("%s PCA", toupper(method_name)))
    ggplot2::ggsave(output_file, note_plot, width = 7, height = 5, units = "in")
    return(invisible(NULL))
  }

  variance_explained <- summary(pca)$importance[2, 1:2] * 100
  pca_df <- tibble::tibble(
    sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  pca_df <- dplyr::left_join(pca_df, sample_meta, by = "sample")

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::geom_text(size = 3, vjust = -0.9, show.legend = FALSE, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = group_palette, drop = FALSE) +
    ggplot2::labs(
      title = sprintf("%s PCA", toupper(method_name)),
      x = sprintf("PC1 (%.1f%%)", variance_explained[[1]]),
      y = sprintf("PC2 (%.1f%%)", variance_explained[[2]]),
      color = "Group"
    ) +
    ggplot2::theme_bw(base_size = 11)

  ggplot2::ggsave(output_file, p, width = 7.5, height = 5.5, units = "in")
}

plot_method_correlation <- function(score_mats, output_file) {
  if (length(score_mats) < 2) {
    return(invisible(NULL))
  }

  common_pathways <- Reduce(intersect, lapply(score_mats, rownames))
  common_samples <- Reduce(intersect, lapply(score_mats, colnames))
  if (length(common_pathways) == 0 || length(common_samples) == 0) {
    warnf("Skipping method-correlation heatmap because the methods do not share any common pathways and samples.")
    return(invisible(NULL))
  }

  flattened <- lapply(score_mats, function(mat) {
    as.vector(mat[common_pathways, common_samples, drop = FALSE])
  })
  corr_mat <- stats::cor(do.call(cbind, flattened), use = "pairwise.complete.obs", method = "pearson")
  colnames(corr_mat) <- names(score_mats)
  rownames(corr_mat) <- names(score_mats)

  grDevices::pdf(output_file, width = 6.5, height = 5.5)
  ph <- pheatmap::pheatmap(
    corr_mat,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    border_color = NA,
    main = "Correlation between pathway-scoring methods",
    silent = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off()
}

write_run_parameters <- function(output_file, opts, output_dirs, input_info, clean_info, filter_info, gene_set_info, prepared_sets) {
  package_list <- unique(c("optparse", "readr", "dplyr", "tidyr", "tibble", "stringr", "purrr", "ggplot2", "pheatmap", "GSVA", "GSEABase", "singscore"))
  installed_pkgs <- package_list[vapply(package_list, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  package_lines <- vapply(
    installed_pkgs,
    function(pkg) sprintf("  - %s: %s", pkg, as.character(utils::packageVersion(pkg))),
    FUN.VALUE = character(1)
  )

  param_lines <- c(
    sprintf("Run timestamp: %s", timestamp()),
    sprintf("Working directory: %s", getwd()),
    sprintf("Command line: %s", paste(commandArgs(trailingOnly = TRUE), collapse = " ")),
    sprintf("R version: %s", R.version.string),
    "",
    sprintf("Input file: %s", normalizePath(opts$input, winslash = "/", mustWork = TRUE)),
    sprintf("Gene-set file: %s", normalizePath(opts$genesets, winslash = "/", mustWork = TRUE)),
    sprintf("Output directory: %s", output_dirs$outdir),
    sprintf("Methods: %s", paste(opts$methods, collapse = ", ")),
    sprintf("Expression prefix: %s", opts$expr_prefix),
    sprintf("Group pattern: %s", opts$group_pattern),
    sprintf("Primary gene ID column: %s", opts$gene_id_col),
    sprintf("Gene ID priority used: %s", paste(clean_info$id_priority, collapse = " -> ")),
    sprintf("Min set size: %d", opts$min_set_size),
    sprintf("Max set size: %d", opts$max_set_size),
    sprintf("Low-expression threshold: %.4f", opts$low_expr_threshold),
    sprintf("Low-expression min samples: %d", opts$low_expr_min_samples),
    sprintf("Low-expression filter disabled: %s", opts$disable_low_expr_filter),
    sprintf("Statistics test: %s", opts$stats_test),
    sprintf("Pairwise test: %s", opts$pairwise_test),
    sprintf("Pairwise correction applied: %s", opts$pairwise_correction),
    sprintf("P-value adjustment: %s", opts$adjust_method),
    sprintf("Boxplot significance annotations: %s", opts$boxplot_significance),
    "",
    sprintf("Detected normalized sample columns: %d", length(input_info$expr_cols)),
    sprintf("Samples: %s", paste(input_info$expr_cols, collapse = ", ")),
    sprintf("Rows after gene ID collapse: %d", nrow(clean_info$data)),
    sprintf("Genes after all expression filters: %d", nrow(filter_info$log_expression_matrix)),
    sprintf("Pathways loaded from input file: %d", length(gene_set_info$pathway_order)),
    sprintf("Pathways retained for GSVA/ssGSEA: %d", length(prepared_sets$gsva_sets)),
    sprintf("Pathways retained for singscore: %d", length(prepared_sets$singscore_specs)),
    "",
    "Package versions:",
    package_lines
  )

  writeLines(param_lines, con = output_file)
}

output_file_map <- function(outdir, plots_dir) {
  list(
    run_parameters = file.path(outdir, "00_run_parameters.txt"),
    sample_metadata = file.path(outdir, "01_sample_metadata.csv"),
    mapping_summary = file.path(outdir, "02_gene_id_mapping_summary.csv"),
    coverage = file.path(outdir, "03_gene_set_coverage.csv"),
    methods = list(
      gsva = list(
        score = file.path(outdir, "04_scores_gsva.csv"),
        long = file.path(outdir, "07_scores_long_gsva.csv"),
        stats = file.path(outdir, "10_stats_gsva.csv"),
        heatmap = file.path(plots_dir, "heatmap_gsva.pdf"),
        boxplots = file.path(plots_dir, "boxplots_gsva.pdf"),
        pca = file.path(plots_dir, "pca_gsva.pdf")
      ),
      ssgsea = list(
        score = file.path(outdir, "05_scores_ssgsea.csv"),
        long = file.path(outdir, "08_scores_long_ssgsea.csv"),
        stats = file.path(outdir, "11_stats_ssgsea.csv"),
        heatmap = file.path(plots_dir, "heatmap_ssgsea.pdf"),
        boxplots = file.path(plots_dir, "boxplots_ssgsea.pdf"),
        pca = file.path(plots_dir, "pca_ssgsea.pdf")
      ),
      singscore = list(
        score = file.path(outdir, "06_scores_singscore.csv"),
        long = file.path(outdir, "09_scores_long_singscore.csv"),
        stats = file.path(outdir, "12_stats_singscore.csv"),
        heatmap = file.path(plots_dir, "heatmap_singscore.pdf"),
        boxplots = file.path(plots_dir, "boxplots_singscore.pdf"),
        pca = file.path(plots_dir, "pca_singscore.pdf")
      )
    ),
    method_correlation = file.path(plots_dir, "method_correlation_heatmap.pdf")
  )
}

run_pathway_scoring <- function(opts) {
  opts <- normalize_run_options(opts)
  check_required_packages(opts$methods, opts$genesets, opts$plot_method_correlation)
  output_dirs <- create_output_dirs(opts$outdir)
  output_paths <- output_file_map(output_dirs$outdir, output_dirs$plots_dir)

  input_info <- read_expression_data(opts$input, opts$expr_prefix)
  clean_info <- resolve_gene_ids(input_info$data, input_info$expr_cols, opts$gene_id_col)
  filter_info <- filter_expression_matrix(
    clean_info$data,
    input_info$expr_cols,
    opts$low_expr_threshold,
    opts$low_expr_min_samples,
    opts$disable_low_expr_filter
  )
  sample_meta <- parse_sample_metadata(colnames(filter_info$log_expression_matrix), opts$group_pattern)
  gene_set_info <- read_gene_sets(opts$genesets)
  prepared_sets <- prepare_gene_sets_for_scoring(
    gene_set_info,
    expr_gene_ids = rownames(filter_info$log_expression_matrix),
    min_set_size = opts$min_set_size,
    max_set_size = opts$max_set_size
  )

  mapping_summary <- dplyr::bind_rows(clean_info$summary, filter_info$summary)
  write_csv_safe(sample_meta, output_paths$sample_metadata)
  write_csv_safe(mapping_summary, output_paths$mapping_summary)
  write_csv_safe(prepared_sets$coverage, output_paths$coverage)

  if (any(c("gsva", "ssgsea") %in% opts$methods) && length(prepared_sets$gsva_sets) == 0) {
    stopf("No gene sets passed overlap and size filters for GSVA/ssGSEA. Review %s for details.", output_paths$coverage)
  }
  if ("singscore" %in% opts$methods && length(prepared_sets$singscore_specs) == 0) {
    stopf("No gene sets passed overlap and size filters for singscore. Review %s for details.", output_paths$coverage)
  }

  score_mats <- list()
  long_tables <- list()
  stats_tables <- list()
  if ("gsva" %in% opts$methods) {
    score_mats$gsva <- score_with_gsva(
      expr_log = filter_info$log_expression_matrix,
      gene_sets = prepared_sets$gsva_sets,
      method_name = "gsva",
      min_set_size = opts$min_set_size,
      max_set_size = opts$max_set_size
    )
  }
  if ("ssgsea" %in% opts$methods) {
    score_mats$ssgsea <- score_with_gsva(
      expr_log = filter_info$log_expression_matrix,
      gene_sets = prepared_sets$gsva_sets,
      method_name = "ssgsea",
      min_set_size = opts$min_set_size,
      max_set_size = opts$max_set_size
    )
  }
  if ("singscore" %in% opts$methods) {
    score_mats$singscore <- score_with_singscore(
      expr_log = filter_info$log_expression_matrix,
      singscore_specs = prepared_sets$singscore_specs
    )
  }

  for (method_name in names(score_mats)) {
    msg("Writing outputs for %s.", method_name)
    score_mat <- score_mats[[method_name]]
    score_mat <- score_mat[, sample_meta$sample, drop = FALSE]

    long_tbl <- score_matrix_to_long(score_mat, sample_meta, method_name)
        stats_tbl <- run_statistics(long_tbl, opts$stats_test, opts$pairwise_test, opts$pairwise_correction, opts$adjust_method)
    long_tables[[method_name]] <- long_tbl
    stats_tables[[method_name]] <- stats_tbl

    write_score_matrix(score_mat, output_paths$methods[[method_name]]$score)
    write_csv_safe(long_tbl, output_paths$methods[[method_name]]$long)
    write_csv_safe(stats_tbl, output_paths$methods[[method_name]]$stats)

    plot_heatmaps(
      score_mat = score_mat,
      sample_meta = sample_meta,
      method_name = method_name,
      output_file = output_paths$methods[[method_name]]$heatmap
    )
    plot_boxplots(
      score_long = long_tbl,
      method_name = method_name,
      output_file = output_paths$methods[[method_name]]$boxplots,
      stats_tbl = stats_tbl,
      show_significance = opts$boxplot_significance
    )
    plot_pca(
      score_mat = score_mat,
      sample_meta = sample_meta,
      method_name = method_name,
      output_file = output_paths$methods[[method_name]]$pca
    )
  }

  if (opts$plot_method_correlation && length(score_mats) > 1) {
    plot_method_correlation(
      score_mats = score_mats,
      output_file = output_paths$method_correlation
    )
  }

  write_run_parameters(
    output_file = output_paths$run_parameters,
    opts = opts,
    output_dirs = output_dirs,
    input_info = input_info,
    clean_info = clean_info,
    filter_info = filter_info,
    gene_set_info = gene_set_info,
    prepared_sets = prepared_sets
  )

  msg("Pathway scoring completed successfully.")
  msg("Results written to: %s", output_dirs$outdir)

  invisible(list(
    options = opts,
    output_dirs = output_dirs,
    output_paths = output_paths,
    input_info = input_info,
    sample_metadata = sample_meta,
    mapping_summary = mapping_summary,
    gene_set_info = gene_set_info,
    coverage = prepared_sets$coverage,
    score_mats = score_mats,
    score_long = long_tables,
    stats = stats_tables
  ))
}

main <- function() {
  opts <- parse_options()
  invisible(run_pathway_scoring(opts))
}

if (sys.nframe() == 0) {
  tryCatch(
    main(),
    error = function(e) {
      cat(sprintf("[%s] ERROR: %s\n", timestamp(), conditionMessage(e)), file = stderr())
      quit(save = "no", status = 1)
    }
  )
}
