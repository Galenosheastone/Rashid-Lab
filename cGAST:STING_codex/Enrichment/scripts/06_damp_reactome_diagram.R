source(file.path("R", "damp_sbgnview.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

.merge_lists <- function(defaults, overrides) {
  if (is.null(overrides) || length(overrides) == 0) {
    return(defaults)
  }
  out <- defaults
  for (nm in names(overrides)) {
    ov <- overrides[[nm]]
    if (nm %in% names(out) && is.list(out[[nm]]) && is.list(ov)) {
      out[[nm]] <- .merge_lists(out[[nm]], ov)
    } else {
      out[[nm]] <- ov
    }
  }
  out
}

.is_abs_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

.sanitize_filename <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "contrast", x)
}

.detect_gene_col <- function(df) {
  candidates <- c(
    "gene_id", "gene", "ENSEMBL", "ensembl", "gene_symbol",
    "symbol", "SYMBOL", "Gene", "GeneID"
  )
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) {
    return(hit[1])
  }
  if (!is.null(rownames(df))) {
    rn <- rownames(df)
    is_default <- identical(rn, as.character(seq_len(nrow(df))))
    if (!is_default && any(nzchar(rn))) {
      return(".rownames")
    }
  }
  NULL
}

.detect_lfc_col <- function(df) {
  candidates <- c("log2FC", "logFC", "log2FoldChange", "LFC", "stat")
  hit <- intersect(candidates, names(df))
  if (length(hit) > 0) {
    return(hit[1])
  }
  NULL
}

.detect_input_type <- function(ids, configured) {
  configured <- toupper(as.character(configured %||% "ENSEMBL"))
  if (configured %in% c("ENSEMBL", "SYMBOL")) {
    return(configured)
  }
  ids <- stats::na.omit(as.character(ids))
  if (length(ids) == 0) {
    return("ENSEMBL")
  }
  if (all(grepl("^ENS", ids))) {
    return("ENSEMBL")
  }
  "SYMBOL"
}

.detect_organism_from_ids <- function(de_list) {
  ids <- character(0)
  for (de_tbl in de_list) {
    de_tbl <- as.data.frame(de_tbl)
    gene_col <- .detect_gene_col(de_tbl)
    if (identical(gene_col, ".rownames")) {
      ids <- c(ids, head(rownames(de_tbl), 4000))
    } else if (!is.null(gene_col) && gene_col %in% names(de_tbl)) {
      ids <- c(ids, head(as.character(de_tbl[[gene_col]]), 4000))
    }
  }
  ids <- unique(stats::na.omit(ids))
  if (length(ids) == 0) {
    return(NULL)
  }
  if (any(grepl("^ENSGALG", ids, ignore.case = TRUE))) {
    return("chicken")
  }
  if (any(grepl("^ENSMUSG", ids, ignore.case = TRUE))) {
    return("mouse")
  }
  if (any(grepl("^ENSG", ids, ignore.case = TRUE))) {
    return("human")
  }
  NULL
}

.load_damp_config <- function(config_path = file.path("config", "pathways.yml")) {
  defaults <- list(
    sbgnview = list(
      enabled = TRUE,
      organism = "mouse",
      reactome_id_human = "R-HSA-5686938",
      reactome_id_mouse = "R-MMU-5686938",
      gene_id_input = "ENSEMBL",
      gene_id_for_sbgnview = "ENTREZ",
      output_formats = c("png", "svg"),
      out_dir = "figures/pathways/damp_reactome"
    )
  )

  cfg <- defaults
  if (file.exists(config_path) && requireNamespace("yaml", quietly = TRUE)) {
    loaded <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
    cfg <- .merge_lists(defaults, loaded)
  }
  cfg$sbgnview
}

run_damp_reactome_diagrams <- function(
  de_list,
  config_path = file.path("config", "pathways.yml"),
  output_root = "output"
) {
  cfg <- .load_damp_config(config_path = config_path)
  if (!isTRUE(cfg$enabled %||% TRUE)) {
    message("DAMP Reactome SBGNview module disabled in config.")
    return(invisible(NULL))
  }

  if (is.null(de_list) || length(de_list) == 0) {
    warning("DAMP Reactome module: de_list is missing; skipping.")
    return(invisible(NULL))
  }

  if (is.null(names(de_list)) || any(names(de_list) == "")) {
    names(de_list) <- paste0("contrast_", seq_along(de_list))
  }

  inferred_org <- .detect_organism_from_ids(de_list)
  organism <- as.character(cfg$organism %||% inferred_org %||% "mouse")
  if (tolower(organism) == "auto" && !is.null(inferred_org)) {
    organism <- inferred_org
  }
  if (!is.null(inferred_org) && !identical(tolower(organism), tolower(inferred_org))) {
    message("DAMP module: using organism inferred from IDs (", inferred_org, ") instead of config (", organism, ").")
    organism <- inferred_org
  }

  org_code <- infer_org_code(organism)
  if (is.na(org_code)) {
    warning("DAMP Reactome module skipped due to unsupported organism='", organism, "'.")
    return(invisible(NULL))
  }

  reactome_stable_id <- if (identical(org_code, "mmu")) {
    as.character(cfg$reactome_id_mouse %||% "R-MMU-5686938")
  } else {
    as.character(cfg$reactome_id_human %||% "R-HSA-5686938")
  }

  pathway_id <- resolve_sbgn_pathway_id(reactome_stable_id = reactome_stable_id, org_code = org_code)
  if (is.na(pathway_id) || !nzchar(pathway_id)) {
    warning(
      "DAMP Reactome module could not resolve pathway ID and will skip. ",
      "Attempted stable_id=", reactome_stable_id, ", organism=", organism, " (", org_code, ")."
    )
    return(invisible(NULL))
  }

  out_formats <- tolower(as.character(unlist(cfg$output_formats %||% c("png", "svg"))))
  out_formats <- unique(out_formats[out_formats %in% c("png", "svg", "pdf")])
  if (length(out_formats) == 0) {
    out_formats <- c("png", "svg")
  }

  fig_rel <- as.character(cfg$out_dir %||% "figures/pathways/damp_reactome")
  fig_dir <- if (.is_abs_path(fig_rel)) fig_rel else file.path(output_root, fig_rel)
  results_dir <- file.path(output_root, "results", "pathways", "damp_reactome")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  summary_rows <- list()
  for (contrast in names(de_list)) {
    de_tbl <- as.data.frame(de_list[[contrast]])
    gene_col <- .detect_gene_col(de_tbl)
    lfc_col <- .detect_lfc_col(de_tbl)

    if (is.null(lfc_col) || (is.null(gene_col) && is.null(rownames(de_tbl)))) {
      warning("DAMP module: missing gene/log2FC columns for contrast '", contrast, "'. Skipping.")
      next
    }

    if (is.null(gene_col)) {
      warning("DAMP module: no gene ID column detected for contrast '", contrast, "'. Skipping.")
      next
    }
    if (identical(gene_col, ".rownames")) {
      gene_ids <- rownames(de_tbl)
    } else {
      gene_ids <- as.character(de_tbl[[gene_col]])
    }
    log2fc <- suppressWarnings(as.numeric(de_tbl[[lfc_col]]))
    input_type <- .detect_input_type(gene_ids, cfg$gene_id_input)
    entrez <- map_ids_to_entrez(gene_ids = gene_ids, organism = organism, input_type = input_type)

    mapping_tbl <- data.frame(
      input_id = as.character(gene_ids),
      entrez_id = as.character(entrez),
      log2FC = log2fc,
      stringsAsFactors = FALSE
    )
    mapping_tbl$mapping_status <- ifelse(
      is.na(mapping_tbl$entrez_id) | mapping_tbl$entrez_id == "",
      "unmapped",
      "mapped"
    )

    mapped_tbl <- mapping_tbl[
      mapping_tbl$mapping_status == "mapped" & !is.na(mapping_tbl$log2FC),
      c("entrez_id", "log2FC"),
      drop = FALSE
    ]

    if (nrow(mapped_tbl) > 0) {
      fc_tbl <- stats::aggregate(log2FC ~ entrez_id, data = mapped_tbl, FUN = mean)
      fc_named_entrez <- fc_tbl$log2FC
      names(fc_named_entrez) <- fc_tbl$entrez_id
    } else {
      fc_named_entrez <- numeric(0)
    }

    contrast_file <- .sanitize_filename(contrast)
    map_tsv <- file.path(results_dir, paste0(contrast_file, "_entrez_mapping.tsv"))
    utils::write.table(mapping_tbl, file = map_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

    out_prefix <- file.path(fig_dir, paste0(contrast_file, "_DAMP_TLR"))
    rendered <- render_damp_tlr_diagram(
      fc_named_entrez = fc_named_entrez,
      org_code = org_code,
      pathway_id = pathway_id,
      out_prefix = out_prefix,
      formats = out_formats
    )

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      contrast = contrast,
      organism = organism,
      org_code = org_code,
      reactome_stable_id = reactome_stable_id,
      sbgn_pathway_id = pathway_id,
      n_input = nrow(mapping_tbl),
      n_mapped = sum(mapping_tbl$mapping_status == "mapped"),
      n_unmapped = sum(mapping_tbl$mapping_status == "unmapped"),
      n_unique_entrez = length(unique(stats::na.omit(mapping_tbl$entrez_id))),
      rendered = isTRUE(rendered),
      output_prefix = out_prefix,
      mapping_tsv = map_tsv,
      stringsAsFactors = FALSE
    )
  }

  if (length(summary_rows) > 0) {
    summary_tbl <- do.call(rbind, summary_rows)
    utils::write.table(
      summary_tbl,
      file = file.path(results_dir, "damp_reactome_render_summary.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  invisible(NULL)
}

if (exists("de_list", inherits = TRUE)) {
  run_damp_reactome_diagrams(
    de_list = get("de_list", inherits = TRUE),
    config_path = file.path("config", "pathways.yml"),
    output_root = if (exists("outdir", inherits = TRUE)) get("outdir", inherits = TRUE) else "output"
  )
} else {
  warning("DAMP Reactome script expected object 'de_list' in parent environment. Nothing was run.")
}
