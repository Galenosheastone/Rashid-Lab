new_repair_log <- function() {
  data.frame(
    category = character(),
    item = character(),
    original = character(),
    cleaned = character(),
    value = character(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}

add_repair_entry <- function(
  repair_log,
  category,
  item,
  original = NA_character_,
  cleaned = NA_character_,
  value = NA_character_,
  detail = NA_character_
) {
  rbind(
    repair_log,
    data.frame(
      category = as.character(category),
      item = as.character(item),
      original = as.character(original),
      cleaned = as.character(cleaned),
      value = as.character(value),
      detail = as.character(detail),
      stringsAsFactors = FALSE
    )
  )
}

add_repair_metric <- function(repair_log, item, value, detail = NA_character_) {
  add_repair_entry(
    repair_log = repair_log,
    category = "summary",
    item = item,
    value = as.character(value),
    detail = detail
  )
}

log_value_changes <- function(repair_log, category, item, original, cleaned, detail = NA_character_) {
  original_chr <- as.character(original)
  cleaned_chr <- as.character(cleaned)

  changed <- (is.na(cleaned) & !is.na(original) & nzchar(trimws(original_chr))) |
    (!is.na(original) & !is.na(cleaned) & original_chr != cleaned_chr)

  if (!any(changed)) {
    return(repair_log)
  }

  rows <- data.frame(
    category = rep(as.character(category), sum(changed)),
    item = rep(as.character(item), sum(changed)),
    original = original_chr[changed],
    cleaned = cleaned_chr[changed],
    value = rep(NA_character_, sum(changed)),
    detail = rep(as.character(detail), sum(changed)),
    stringsAsFactors = FALSE
  )

  rbind(repair_log, rows)
}

clean_names <- function(x) {
  x <- trimws(as.character(x))
  x <- basename(x)
  x <- sub("\\.(bam|sam|counts(?:\\.txt)?|txt|fq|fastq|gz)$", "", x, ignore.case = TRUE)
  x <- sub("\\.sorted$", "", x, ignore.case = TRUE)
  x <- sub("(_S\\d+)?(_R[12])?(_00[12])?$", "", x, ignore.case = TRUE)
  x <- sub("(_S\\d+)?(_R[12])?(_001)?$", "", x, ignore.case = TRUE)
  x <- sub("(\\.aligned)?$", "", x, ignore.case = TRUE)
  x
}

normalize_condition_key <- function(x) {
  gsub("[^a-z0-9]+", "", tolower(trimws(as.character(x))))
}

canonicalize_condition <- function(x, allow_aliases = FALSE) {
  original <- trimws(as.character(x))
  original[is.na(original)] <- ""
  key <- normalize_condition_key(original)
  out <- rep(NA_character_, length(original))

  exact_map <- c(sacral = "Sacral", free = "Free", pygo = "Pygo")
  exact_hits <- key %in% names(exact_map)
  out[exact_hits] <- unname(exact_map[key[exact_hits]])

  if (isTRUE(allow_aliases)) {
    out[is.na(out) & grepl("^sacral[0-9]+$", key)] <- "Sacral"
    out[is.na(out) & grepl("^free[0-9]+$", key)] <- "Free"
    out[is.na(out) & grepl("^pygo[0-9]+$", key)] <- "Pygo"
  }

  out[key == ""] <- NA_character_
  out
}

infer_condition_from_sample <- function(x, allow_aliases = FALSE) {
  prefixes <- sub("([._-].*)$", "", trimws(as.character(x)))
  canonicalize_condition(prefixes, allow_aliases = allow_aliases)
}

first_token <- function(x) {
  sub(";.*$", "", x)
}

strip_ens_version <- function(x) {
  sub("\\.[0-9]+$", "", x)
}

annotation_columns <- function(expression_table) {
  intersect(
    colnames(expression_table),
    c(
      "SYMBOL", "NCBI_GeneID", "Ensembl_GeneID", "gene_biotype",
      "Chr", "Start", "End", "Strand", "Length"
    )
  )
}

maybe_fill_symbols_from_orgdb <- function(gene_annotation, config, repair_log) {
  if (is.null(config$OrgDb) || !"Ensembl_GeneID" %in% colnames(gene_annotation)) {
    return(list(gene_annotation = gene_annotation, repair_log = repair_log))
  }

  symbol_values <- trimws(as.character(gene_annotation$SYMBOL))
  symbol_values[is.na(symbol_values)] <- ""
  if (any(symbol_values != "")) {
    return(list(gene_annotation = gene_annotation, repair_log = repair_log))
  }

  ensembl_ids <- trimws(as.character(gene_annotation$Ensembl_GeneID))
  ensembl_ids[is.na(ensembl_ids)] <- ""
  ensembl_ids <- strip_ens_version(ensembl_ids)
  valid_ids <- ensembl_ids != ""

  if (!any(valid_ids)) {
    return(list(gene_annotation = gene_annotation, repair_log = repair_log))
  }

  symbol_map <- tryCatch(
    AnnotationDbi::mapIds(
      config$OrgDb,
      keys = unique(ensembl_ids[valid_ids]),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    ),
    error = function(e) NULL
  )

  if (is.null(symbol_map)) {
    repair_log <- add_repair_metric(
      repair_log,
      "symbol_mapped_from_orgdb_count",
      0L,
      "OrgDb mapping to SYMBOL was attempted but returned no results."
    )
    return(list(gene_annotation = gene_annotation, repair_log = repair_log))
  }

  mapped_symbols <- unname(symbol_map[ensembl_ids[valid_ids]])
  gene_annotation$SYMBOL[valid_ids] <- mapped_symbols
  mapped_count <- sum(!is.na(mapped_symbols) & nzchar(mapped_symbols))

  repair_log <- add_repair_metric(
    repair_log,
    "symbol_mapped_from_orgdb_count",
    mapped_count,
    "Filled missing SYMBOL values from OrgDb using Ensembl_GeneID."
  )

  list(gene_annotation = gene_annotation, repair_log = repair_log)
}

extract_numeric_count_matrix <- function(expression_table, annotation_cols, repair_log) {
  candidate_cols <- setdiff(colnames(expression_table), c("Geneid", annotation_cols))
  if (length(candidate_cols) == 0L) {
    stop("No sample columns were found after excluding annotation columns.", call. = FALSE)
  }

  candidate_df <- expression_table[, candidate_cols, drop = FALSE]
  numeric_mask <- vapply(
    candidate_df,
    function(col) is.numeric(col) || is.integer(col),
    logical(1)
  )

  if (!all(numeric_mask)) {
    dropped_cols <- candidate_cols[!numeric_mask]
    warning(
      "Dropping non-numeric columns from the count matrix: ",
      paste(dropped_cols, collapse = ", "),
      call. = FALSE
    )

    for (col_name in dropped_cols) {
      repair_log <- add_repair_entry(
        repair_log,
        category = "counts_matrix",
        item = "dropped_non_numeric_column",
        original = col_name,
        cleaned = NA_character_,
        value = "dropped",
        detail = "Excluded from count matrix because the column is not numeric."
      )
    }
  }

  sample_cols <- candidate_cols[numeric_mask]
  if (length(sample_cols) == 0L) {
    stop("No numeric sample columns remain after dropping non-numeric columns.", call. = FALSE)
  }

  count_matrix <- as.matrix(candidate_df[, sample_cols, drop = FALSE])
  storage.mode(count_matrix) <- "double"

  list(
    count_matrix = count_matrix,
    repair_log = repair_log
  )
}

clean_inputs <- function(inputs, config) {
  repair_log <- new_repair_log()
  expression <- inputs$counts_table
  colnames(expression) <- strip_bom(colnames(expression))

  geneid_original <- as.character(expression$Geneid)
  geneid_clean <- trimws(geneid_original)
  geneid_clean[is.na(geneid_clean)] <- ""

  blank_geneid_count <- sum(geneid_clean == "")
  repair_log <- add_repair_metric(repair_log, "blank_Geneid_count", blank_geneid_count)

  ensembl_ids <- rep("", nrow(expression))
  if ("Ensembl_GeneID" %in% colnames(expression)) {
    ensembl_ids <- trimws(as.character(expression$Ensembl_GeneID))
    ensembl_ids[is.na(ensembl_ids)] <- ""
    ensembl_ids <- first_token(ensembl_ids)
  }

  use_ensembl_fallback <- geneid_clean == "" & ensembl_ids != ""
  geneid_clean[use_ensembl_fallback] <- ensembl_ids[use_ensembl_fallback]
  repair_log <- add_repair_metric(repair_log, "ensembl_fallback_used", any(use_ensembl_fallback))
  repair_log <- add_repair_metric(repair_log, "ensembl_fallback_count", sum(use_ensembl_fallback))

  still_missing_geneid <- geneid_clean == ""
  if (any(still_missing_geneid)) {
    geneid_clean[still_missing_geneid] <- sprintf(
      "UNLABELED_%06d",
      seq_len(sum(still_missing_geneid))
    )
  }

  repair_log <- add_repair_metric(
    repair_log,
    "generated_UNLABELED_count",
    sum(still_missing_geneid)
  )

  repair_log <- add_repair_metric(
    repair_log,
    "duplicate_Geneid_count_before_make.unique",
    sum(duplicated(geneid_clean))
  )

  geneid_unique <- make.unique(geneid_clean)
  repair_log <- add_repair_metric(
    repair_log,
    "make.unique_adjustment_count",
    sum(geneid_unique != geneid_clean)
  )

  geneid_clean <- geneid_unique
  rownames(expression) <- geneid_clean
  anno_cols <- annotation_columns(expression)

  gene_annotation <- expression[, unique(c("Geneid", anno_cols)), drop = FALSE]
  gene_annotation$Geneid_original <- geneid_original
  gene_annotation$Geneid_clean <- geneid_clean
  rownames(gene_annotation) <- geneid_clean

  placeholder_cols <- c("SYMBOL", "Ensembl_GeneID", "NCBI_GeneID")
  for (col_name in placeholder_cols) {
    if (!col_name %in% colnames(gene_annotation)) {
      gene_annotation[[col_name]] <- NA_character_
      repair_log <- add_repair_entry(
        repair_log,
        category = "annotation",
        item = "added_placeholder_column",
        original = NA_character_,
        cleaned = col_name,
        value = "NA",
        detail = "Added missing annotation column as an NA placeholder."
      )
    }
  }

  symbol_fill <- maybe_fill_symbols_from_orgdb(gene_annotation, config, repair_log)
  gene_annotation <- symbol_fill$gene_annotation
  repair_log <- symbol_fill$repair_log

  count_info <- extract_numeric_count_matrix(expression, anno_cols, repair_log)
  counts_matrix <- count_info$count_matrix
  repair_log <- count_info$repair_log

  original_count_names <- colnames(counts_matrix)
  cleaned_count_names <- clean_names(original_count_names)
  colnames(counts_matrix) <- cleaned_count_names
  repair_log <- log_value_changes(
    repair_log,
    category = "sample_name",
    item = "counts_column_cleaned",
    original = original_count_names,
    cleaned = cleaned_count_names,
    detail = "Counts column name cleaned with clean_names()."
  )

  metadata_inferred <- FALSE
  if (!is.null(inputs$metadata_table)) {
    coldata <- inputs$metadata_table
    if (is.null(rownames(coldata))) {
      stop("Metadata table must have sample IDs as row names.", call. = FALSE)
    }
    original_metadata_names <- rownames(coldata)
    cleaned_metadata_names <- clean_names(original_metadata_names)
    rownames(coldata) <- cleaned_metadata_names
    repair_log <- log_value_changes(
      repair_log,
      category = "sample_name",
      item = "metadata_rowname_cleaned",
      original = original_metadata_names,
      cleaned = cleaned_metadata_names,
      detail = "Metadata row name cleaned with clean_names()."
    )
  } else {
    if (!isTRUE(config$infer_metadata_from_sample_names)) {
      stop(
        "Metadata file was not loaded and infer_metadata_from_sample_names is FALSE.",
        call. = FALSE
      )
    }

    metadata_inferred <- TRUE
    coldata <- data.frame(
      condition = infer_condition_from_sample(
        colnames(counts_matrix),
        allow_aliases = config$allow_condition_aliases
      ),
      row.names = colnames(counts_matrix),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (!"condition" %in% colnames(coldata)) {
    if (!isTRUE(config$infer_metadata_from_sample_names)) {
      stop(
        "Metadata must contain a 'condition' column unless inference is explicitly enabled.",
        call. = FALSE
      )
    }

    coldata$condition <- infer_condition_from_sample(
      rownames(coldata),
      allow_aliases = config$allow_condition_aliases
    )
    metadata_inferred <- TRUE
  }

  repair_log <- add_repair_metric(repair_log, "metadata_inferred", metadata_inferred)
  repair_log <- add_repair_metric(repair_log, "metadata_source", inputs$metadata_source)

  coldata$condition_original <- as.character(coldata$condition)
  canonical_condition <- canonicalize_condition(
    coldata$condition_original,
    allow_aliases = config$allow_condition_aliases
  )

  repair_log <- log_value_changes(
    repair_log,
    category = "condition",
    item = "metadata_condition_canonicalized",
    original = coldata$condition_original,
    cleaned = canonical_condition,
    detail = "Condition label canonicalized with canonicalize_condition()."
  )

  sanitized_condition <- ifelse(is.na(canonical_condition), NA_character_, make.names(canonical_condition))
  coldata$condition <- factor(sanitized_condition)

  if (setequal(colnames(counts_matrix), rownames(coldata))) {
    metadata_reordered <- !identical(colnames(counts_matrix), rownames(coldata))
    coldata <- coldata[colnames(counts_matrix), , drop = FALSE]
  } else {
    metadata_reordered <- FALSE
  }

  repair_log <- add_repair_metric(
    repair_log,
    "metadata_reordered_to_match_counts",
    metadata_reordered
  )

  contrast_plan <- config$contrast_plan
  contrast_plan$group1_original <- contrast_plan$group1
  contrast_plan$group0_original <- contrast_plan$group0

  contrast_group1 <- canonicalize_condition(
    contrast_plan$group1_original,
    allow_aliases = config$allow_condition_aliases
  )
  contrast_group0 <- canonicalize_condition(
    contrast_plan$group0_original,
    allow_aliases = config$allow_condition_aliases
  )

  repair_log <- log_value_changes(
    repair_log,
    category = "condition",
    item = "contrast_group1_canonicalized",
    original = contrast_plan$group1_original,
    cleaned = contrast_group1,
    detail = "Contrast group1 canonicalized with canonicalize_condition()."
  )
  repair_log <- log_value_changes(
    repair_log,
    category = "condition",
    item = "contrast_group0_canonicalized",
    original = contrast_plan$group0_original,
    cleaned = contrast_group0,
    detail = "Contrast group0 canonicalized with canonicalize_condition()."
  )

  contrast_plan$group1 <- ifelse(is.na(contrast_group1), NA_character_, make.names(contrast_group1))
  contrast_plan$group0 <- ifelse(is.na(contrast_group0), NA_character_, make.names(contrast_group0))
  contrast_plan <- unique(contrast_plan)

  ref_level <- NULL
  if (nlevels(coldata$condition) > 0L) {
    if (is.null(config$ref_level)) {
      ref_level <- if ("Free" %in% levels(coldata$condition)) "Free" else levels(coldata$condition)[1]
    } else {
      ref_candidate <- canonicalize_condition(
        config$ref_level,
        allow_aliases = config$allow_condition_aliases
      )

      if (is.na(ref_candidate)) {
        stop(
          "Configured ref_level could not be canonicalized safely: ",
          config$ref_level,
          call. = FALSE
        )
      }

      ref_level <- make.names(ref_candidate)
      if (!ref_level %in% levels(coldata$condition)) {
        stop(
          "Requested ref_level '", ref_level,
          "' is not present in metadata condition levels: ",
          paste(levels(coldata$condition), collapse = ", "),
          call. = FALSE
        )
      }
    }

    coldata$condition <- stats::relevel(coldata$condition, ref = ref_level)
  }

  repair_log <- add_repair_metric(
    repair_log,
    "reference_level",
    if (is.null(ref_level)) "" else ref_level
  )

  list(
    counts_matrix = counts_matrix,
    metadata = coldata,
    gene_annotation = gene_annotation,
    contrast_plan = contrast_plan,
    repair_log = repair_log,
    ref_level = ref_level
  )
}
