is_integer_like <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  abs(x - round(x)) <= tolerance
}

format_example_values <- function(x, n = 5L) {
  paste(utils::head(unique(as.character(x)), n), collapse = ", ")
}

validate_inputs <- function(cleaned, config) {
  counts_matrix <- cleaned$counts_matrix
  coldata <- cleaned$metadata
  contrast_plan <- cleaned$contrast_plan

  if (!is.matrix(counts_matrix) || !is.numeric(counts_matrix)) {
    stop("Counts matrix must be numeric.", call. = FALSE)
  }

  if (nrow(counts_matrix) == 0L || ncol(counts_matrix) == 0L) {
    stop("Counts matrix is empty after input cleaning.", call. = FALSE)
  }

  if (anyNA(counts_matrix)) {
    stop("Counts matrix contains NA values.", call. = FALSE)
  }

  if (!all(is.finite(counts_matrix))) {
    stop("Counts matrix contains non-finite values.", call. = FALSE)
  }

  if (any(counts_matrix < 0)) {
    stop("Counts matrix contains negative values.", call. = FALSE)
  }

  if (anyDuplicated(rownames(counts_matrix))) {
    stop("Duplicate gene IDs remain after cleaning.", call. = FALSE)
  }

  if (anyDuplicated(colnames(counts_matrix))) {
    dupes <- unique(colnames(counts_matrix)[duplicated(colnames(counts_matrix))])
    stop(
      "Duplicate sample IDs remain in the counts matrix after cleaning: ",
      paste(dupes, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(rownames(coldata)) || any(!nzchar(rownames(coldata)))) {
    stop("Metadata row names must contain non-empty sample IDs.", call. = FALSE)
  }

  if (anyDuplicated(rownames(coldata))) {
    dupes <- unique(rownames(coldata)[duplicated(rownames(coldata))])
    stop(
      "Duplicate sample IDs remain in metadata after cleaning: ",
      paste(dupes, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"condition" %in% colnames(coldata)) {
    stop("Metadata is missing the required 'condition' column.", call. = FALSE)
  }

  if (!identical(colnames(counts_matrix), rownames(coldata))) {
    if (setequal(colnames(counts_matrix), rownames(coldata))) {
      stop(
        "Counts columns and metadata row names contain the same samples but are not aligned.",
        call. = FALSE
      )
    }

    stop(
      "Counts columns and metadata row names do not match after cleaning.\n",
      "Not in metadata: ", paste(setdiff(colnames(counts_matrix), rownames(coldata)), collapse = ", "), "\n",
      "Not in counts: ", paste(setdiff(rownames(coldata), colnames(counts_matrix)), collapse = ", "),
      call. = FALSE
    )
  }

  bad_conditions <- is.na(coldata$condition) | !nzchar(as.character(coldata$condition))
  if (any(bad_conditions)) {
    bad_labels <- unique(coldata$condition_original[bad_conditions])
    stop(
      "Some metadata condition labels could not be canonicalized safely: ",
      paste(bad_labels, collapse = ", "),
      ". Fix the metadata or set allow_condition_aliases = TRUE if those labels are expected.",
      call. = FALSE
    )
  }

  if (anyNA(contrast_plan$group1) || anyNA(contrast_plan$group0)) {
    bad_groups <- unique(
      c(
        contrast_plan$group1_original[is.na(contrast_plan$group1)],
        contrast_plan$group0_original[is.na(contrast_plan$group0)]
      )
    )
    stop(
      "Some contrast plan labels could not be canonicalized safely: ",
      paste(bad_groups, collapse = ", "),
      call. = FALSE
    )
  }

  self_comparisons <- contrast_plan$group1 == contrast_plan$group0
  if (any(self_comparisons, na.rm = TRUE)) {
    bad_pairs <- paste(
      contrast_plan$group1[self_comparisons],
      "vs",
      contrast_plan$group0[self_comparisons]
    )
    stop(
      "Contrast plan contains self-comparisons: ",
      paste(bad_pairs, collapse = "; "),
      call. = FALSE
    )
  }

  missing_groups <- setdiff(
    unique(c(contrast_plan$group1, contrast_plan$group0)),
    levels(coldata$condition)
  )
  if (length(missing_groups) > 0L) {
    stop(
      "Requested contrast groups are missing from metadata condition levels: ",
      paste(missing_groups, collapse = ", "),
      ". Available levels: ",
      paste(levels(coldata$condition), collapse = ", "),
      call. = FALSE
    )
  }

  integer_like_mask <- is_integer_like(counts_matrix)
  if (!all(integer_like_mask, na.rm = TRUE)) {
    non_integer_values <- counts_matrix[!integer_like_mask]
    stop(
      "Non-integer count values detected before coercion. Example value(s): ",
      format_example_values(non_integer_values),
      ". DESeq2 expects raw counts.",
      call. = FALSE
    )
  }

  replicate_counts <- table(coldata$condition)
  low_replicate_groups <- replicate_counts[replicate_counts < 2L]
  if (length(low_replicate_groups) > 0L) {
    warning(
      "Some condition groups have fewer than 2 replicates: ",
      paste(names(low_replicate_groups), low_replicate_groups, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }

  if (identical(config$filter_mode, "strict")) {
    requested_counts <- replicate_counts[
      names(replicate_counts) %in% unique(c(contrast_plan$group1, contrast_plan$group0))
    ]
    strict_warning_groups <- requested_counts[requested_counts < config$min_reps_within_group]
    if (length(strict_warning_groups) > 0L) {
      warning(
        "Strict filtering is configured for >=", config$min_reps_within_group,
        " samples within a condition. These requested groups have fewer replicates: ",
        paste(names(strict_warning_groups), strict_warning_groups, sep = "=", collapse = ", "),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
