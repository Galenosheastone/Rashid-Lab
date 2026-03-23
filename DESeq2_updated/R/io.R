`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    y
  } else {
    x
  }
}

normalize_optional_scalar <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(NULL)
  }

  if (is.character(x) && length(x) == 1L && trimws(x) %in% c("", "NA", "NULL", "null")) {
    return(NULL)
  }

  x
}

default_analysis_config <- function() {
  list(
    output_dir = "outputs",
    counts_file = NULL,
    metadata_file = NULL,
    ref_level = NULL,
    results_alpha = 0.05,
    filter_mode = "strict",
    min_count = 10L,
    min_reps_within_group = 3L,
    use_rlog = FALSE,
    run_lrt = TRUE,
    run_gsea = FALSE,
    OrgDb_pkg = NULL,
    kegg_org = NULL,
    shrink_method = "apeglm",
    n_workers = 1L,
    infer_metadata_from_sample_names = FALSE,
    allow_condition_aliases = FALSE,
    contrast_plan = list()
  )
}

coerce_contrast_plan <- function(plan) {
  if (is.null(plan) || length(plan) == 0L) {
    return(data.frame(group1 = character(), group0 = character(), stringsAsFactors = FALSE))
  }

  if (is.data.frame(plan)) {
    out <- as.data.frame(plan, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.list(plan)) {
    rows <- lapply(plan, function(entry) {
      as.data.frame(as.list(entry), stringsAsFactors = FALSE, check.names = FALSE)
    })
    out <- do.call(rbind, rows)
  } else {
    stop("contrast_plan must be a list of group1/group0 entries or a data frame.", call. = FALSE)
  }

  required_cols <- c("group1", "group0")
  if (!all(required_cols %in% colnames(out))) {
    stop("contrast_plan must contain 'group1' and 'group0' columns.", call. = FALSE)
  }

  out$group1 <- as.character(out$group1)
  out$group0 <- as.character(out$group0)
  out
}

load_config <- function(path = "config/analysis_config.yaml") {
  if (!file.exists(path)) {
    stop("Config file not found: ", path, call. = FALSE)
  }

  user_config <- yaml::read_yaml(path)
  if (is.null(user_config)) {
    user_config <- list()
  }

  contrast_plan_raw <- if ("contrast_plan" %in% names(user_config)) {
    user_config$contrast_plan
  } else {
    default_analysis_config()$contrast_plan
  }

  config <- utils::modifyList(
    default_analysis_config(),
    user_config[setdiff(names(user_config), "contrast_plan")]
  )
  config$output_dir <- normalize_optional_scalar(config$output_dir) %||% "outputs"
  config$counts_file <- normalize_optional_scalar(config$counts_file)
  config$metadata_file <- normalize_optional_scalar(config$metadata_file)
  config$ref_level <- normalize_optional_scalar(config$ref_level)
  config$OrgDb_pkg <- normalize_optional_scalar(config$OrgDb_pkg)
  config$kegg_org <- normalize_optional_scalar(config$kegg_org)
  config$filter_mode <- tolower(as.character(config$filter_mode %||% "strict"))
  config$shrink_method <- tolower(as.character(config$shrink_method %||% "apeglm"))
  config$results_alpha <- as.numeric(config$results_alpha %||% 0.05)
  config$min_count <- as.integer(config$min_count %||% 10L)
  config$min_reps_within_group <- as.integer(config$min_reps_within_group %||% 3L)
  config$n_workers <- as.integer(config$n_workers %||% 1L)
  config$contrast_plan <- coerce_contrast_plan(contrast_plan_raw)
  config$repair_log_file <- "input_repair_log.csv"
  config$combined_export_file <- "DESeq2_combined_contrast_export.csv"
  config$config_path <- path
  config$OrgDb <- NULL

  if (is.null(config$counts_file)) {
    stop("config$counts_file must be set.", call. = FALSE)
  }

  if (!config$filter_mode %in% c("legacy", "strict")) {
    stop("config$filter_mode must be either 'legacy' or 'strict'.", call. = FALSE)
  }

  if (!config$shrink_method %in% c("apeglm", "ashr")) {
    stop("config$shrink_method must be either 'apeglm' or 'ashr'.", call. = FALSE)
  }

  if (is.na(config$results_alpha) || config$results_alpha <= 0 || config$results_alpha >= 1) {
    stop("config$results_alpha must be a number between 0 and 1.", call. = FALSE)
  }

  if (is.na(config$min_count) || config$min_count < 1L) {
    stop("config$min_count must be >= 1.", call. = FALSE)
  }

  if (is.na(config$min_reps_within_group) || config$min_reps_within_group < 1L) {
    stop("config$min_reps_within_group must be >= 1.", call. = FALSE)
  }

  if (is.na(config$n_workers) || config$n_workers < 1L) {
    config$n_workers <- 1L
  }

  if (nrow(config$contrast_plan) == 0L) {
    stop("config$contrast_plan must contain at least one contrast.", call. = FALSE)
  }

  if (isTRUE(config$run_gsea) && is.null(config$OrgDb_pkg)) {
    stop("config$OrgDb_pkg must be set when run_gsea is TRUE.", call. = FALSE)
  }

  if (isTRUE(config$run_gsea) && is.null(config$kegg_org)) {
    stop("config$kegg_org must be set when run_gsea is TRUE.", call. = FALSE)
  }

  config
}

required_packages_for_config <- function(config) {
  pkgs <- c("yaml", "DESeq2", "BiocParallel", "ggplot2", "pheatmap", "vsn")

  if (identical(config$shrink_method, "apeglm")) {
    pkgs <- c(pkgs, "apeglm")
  }

  if (identical(config$shrink_method, "ashr")) {
    pkgs <- c(pkgs, "ashr")
  }

  if (isTRUE(config$run_gsea)) {
    pkgs <- c(pkgs, "AnnotationDbi", "clusterProfiler", "enrichplot", "pathview", config$OrgDb_pkg)
  }

  unique(pkgs)
}

resolve_orgdb_object <- function(pkg_name) {
  tryCatch(
    get(pkg_name, envir = asNamespace(pkg_name)),
    error = function(e) {
      tryCatch(get(pkg_name, envir = .GlobalEnv), error = function(e2) NULL)
    }
  )
}

load_required_packages <- function(config) {
  required_pkgs <- required_packages_for_config(config)
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]

  if (length(missing_pkgs) > 0L) {
    stop(
      "Missing required package(s): ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running the pipeline.",
      call. = FALSE
    )
  }

  for (pkg in required_pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }

  if (isTRUE(config$run_gsea)) {
    config$OrgDb <- resolve_orgdb_object(config$OrgDb_pkg)
    if (is.null(config$OrgDb)) {
      warning(
        "Could not retrieve the OrgDb object from package '", config$OrgDb_pkg,
        "'. GSEA outputs will be skipped.",
        call. = FALSE
      )
    }
  }

  config
}

strip_bom <- function(x) {
  sub("^\ufeff", "", x)
}

read_delimited_table <- function(path, row_names = FALSE) {
  is_csv <- grepl("\\.csv(?:\\.gz)?$", path, ignore.case = TRUE)

  if (is_csv) {
    utils::read.csv(
      path,
      row.names = if (isTRUE(row_names)) 1 else NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
  } else {
    utils::read.delim(
      path,
      row.names = if (isTRUE(row_names)) 1 else NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      fileEncoding = "UTF-8-BOM"
    )
  }
}

read_metadata_table <- function(path) {
  metadata_table <- read_delimited_table(path, row_names = TRUE)
  colnames(metadata_table) <- strip_bom(colnames(metadata_table))

  if ("condition" %in% colnames(metadata_table)) {
    return(metadata_table)
  }

  raw_metadata <- read_delimited_table(path, row_names = FALSE)
  colnames(raw_metadata) <- strip_bom(colnames(raw_metadata))
  sample_col_candidates <- c("sample_id", "sample", "sample_name", "SampleID", "Sample", "Sample_Name")
  sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(raw_metadata)][1]

  if (!is.na(sample_col) && "condition" %in% colnames(raw_metadata)) {
    rownames(raw_metadata) <- as.character(raw_metadata[[sample_col]])
    raw_metadata[[sample_col]] <- NULL
    return(raw_metadata)
  }

  stop(
    "Metadata file must either use the first column as sample IDs or include a ",
    "'sample_id' column together with a 'condition' column: ", path,
    call. = FALSE
  )
}

load_inputs <- function(config) {
  if (!file.exists(config$counts_file)) {
    stop("Counts file not found: ", config$counts_file, call. = FALSE)
  }

  counts_table <- read_delimited_table(config$counts_file, row_names = FALSE)
  colnames(counts_table) <- strip_bom(colnames(counts_table))

  if (!"Geneid" %in% colnames(counts_table)) {
    stop("Counts file must contain a 'Geneid' column: ", config$counts_file, call. = FALSE)
  }

  metadata_table <- NULL
  metadata_source <- "missing"

  if (!is.null(config$metadata_file) && file.exists(config$metadata_file)) {
    metadata_table <- read_metadata_table(config$metadata_file)
    metadata_source <- "file"
  }

  list(
    counts_table = counts_table,
    metadata_table = metadata_table,
    metadata_source = metadata_source,
    counts_file = config$counts_file,
    metadata_file = config$metadata_file
  )
}

describe_filter <- function(config) {
  if (identical(config$filter_mode, "legacy")) {
    "legacy (>0 reads in at least 1 sample)"
  } else {
    paste0(
      "strict (>=", config$min_count,
      " reads in >=", config$min_reps_within_group,
      " samples within at least 1 condition)"
    )
  }
}
