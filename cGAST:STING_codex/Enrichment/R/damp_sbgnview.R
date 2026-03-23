`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

infer_org_code <- function(organism) {
  org <- tolower(trimws(as.character(organism %||% "")))
  if (org %in% c("mouse", "mus musculus", "mmu", "m.musculus", "mm")) {
    return("mmu")
  }
  if (org %in% c("human", "homo sapiens", "hsa", "h.sapiens", "hs")) {
    return("hsa")
  }
  warning(
    "DAMP Reactome SBGNview currently supports only mouse/human. ",
    "Received organism='", organism, "'."
  )
  NA_character_
}

.resolve_orgdb <- function(org_code) {
  pkg <- switch(
    org_code,
    mmu = "org.Mm.eg.db",
    hsa = "org.Hs.eg.db",
    NULL
  )
  if (is.null(pkg)) {
    return(NULL)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning(pkg, " is required for ID mapping but is not installed.")
    return(NULL)
  }
  tryCatch(
    getExportedValue(pkg, pkg),
    error = function(e) {
      warning("Failed to load OrgDb object from ", pkg, ": ", conditionMessage(e))
      NULL
    }
  )
}

map_ids_to_entrez <- function(gene_ids, organism, input_type = c("ENSEMBL", "SYMBOL")) {
  input_type <- match.arg(input_type)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    warning("AnnotationDbi is required for ID mapping but is not installed.")
    return(rep(NA_character_, length(gene_ids)))
  }

  gene_ids <- as.character(gene_ids)
  gene_ids_clean <- gene_ids
  if (identical(input_type, "ENSEMBL")) {
    gene_ids_clean <- sub("\\.\\d+$", "", gene_ids_clean)
  }

  org_code <- infer_org_code(organism)
  orgdb <- .resolve_orgdb(org_code)
  if (is.na(org_code) || is.null(orgdb)) {
    out <- rep(NA_character_, length(gene_ids_clean))
    names(out) <- gene_ids
    return(out)
  }

  mapped <- tryCatch(
    AnnotationDbi::mapIds(
      x = orgdb,
      keys = unique(gene_ids_clean),
      column = "ENTREZID",
      keytype = input_type,
      multiVals = "first"
    ),
    error = function(e) {
      warning("Failed to map IDs to ENTREZID: ", conditionMessage(e))
      NULL
    }
  )

  out <- rep(NA_character_, length(gene_ids_clean))
  names(out) <- gene_ids
  if (is.null(mapped) || length(mapped) == 0) {
    return(out)
  }

  mapped_chr <- as.character(mapped)
  names(mapped_chr) <- names(mapped)
  out <- unname(mapped_chr[gene_ids_clean])
  names(out) <- gene_ids
  out
}

.extract_pathway_table <- function(pathway_obj) {
  if (is.null(pathway_obj)) {
    return(data.frame())
  }
  if (is.data.frame(pathway_obj)) {
    return(pathway_obj)
  }
  if (is.list(pathway_obj)) {
    for (x in pathway_obj) {
      if (is.data.frame(x)) {
        return(x)
      }
    }
  }
  tryCatch(as.data.frame(pathway_obj), error = function(e) data.frame())
}

.first_matching_col <- function(tbl, candidates) {
  hits <- intersect(candidates, names(tbl))
  if (length(hits) == 0) {
    return(NULL)
  }
  hits[1]
}

.pick_pathway_row <- function(tbl, reactome_stable_id) {
  if (nrow(tbl) == 0) {
    return(tbl)
  }

  id_cols <- intersect(
    c("pathway.id", "pathway_id", "pathwayId", "reactome.id", "reactome_id", "id"),
    names(tbl)
  )

  if (length(id_cols) > 0) {
    stable <- tolower(reactome_stable_id)
    for (col in id_cols) {
      vals <- tolower(as.character(tbl[[col]]))
      idx <- which(vals == stable)
      if (length(idx) > 0) {
        return(tbl[idx[1], , drop = FALSE])
      }
    }
  }

  score_col <- .first_matching_col(tbl, c("relevance", "score", "Rank", "rank"))
  if (!is.null(score_col)) {
    suppressWarnings({
      score_vals <- as.numeric(tbl[[score_col]])
    })
    if (any(!is.na(score_vals))) {
      idx <- which.max(score_vals)
      return(tbl[idx, , drop = FALSE])
    }
  }

  tbl[1, , drop = FALSE]
}

resolve_sbgn_pathway_id <- function(reactome_stable_id, org_code) {
  if (!requireNamespace("SBGNview", quietly = TRUE)) {
    warning("SBGNview is not installed; cannot resolve Reactome pathway.")
    return(NA_character_)
  }

  res_by_id <- tryCatch(
    SBGNview::findPathways(
      keywords = reactome_stable_id,
      keyword.type = "pathway.id",
      org = org_code
    ),
    error = function(e) NULL
  )
  tbl <- .extract_pathway_table(res_by_id)

  if (nrow(tbl) == 0) {
    res_by_name <- tryCatch(
      SBGNview::findPathways(
        keywords = "Regulation of TLR by endogenous ligand",
        keyword.type = "pathway.name",
        org = org_code
      ),
      error = function(e) NULL
    )
    tbl <- .extract_pathway_table(res_by_name)
  }

  if (nrow(tbl) == 0) {
    warning(
      "Could not resolve SBGN pathway ID for Reactome stable ID ",
      reactome_stable_id,
      " (org=", org_code, ")."
    )
    return(NA_character_)
  }

  picked <- .pick_pathway_row(tbl, reactome_stable_id)
  id_col <- .first_matching_col(picked, c("pathway.id", "pathway_id", "pathwayId", "id"))
  if (is.null(id_col)) {
    warning(
      "Resolved pathway table is missing a usable pathway ID column for ",
      reactome_stable_id, " (org=", org_code, ")."
    )
    return(NA_character_)
  }

  as.character(picked[[id_col]][1])
}

render_damp_tlr_diagram <- function(fc_named_entrez, org_code, pathway_id, out_prefix, formats = c("png", "svg")) {
  if (!requireNamespace("SBGNview", quietly = TRUE)) {
    warning("SBGNview is not installed; skipping DAMP/TLR diagram rendering.")
    return(invisible(FALSE))
  }

  fc_vals <- suppressWarnings(as.numeric(fc_named_entrez))
  names(fc_vals) <- names(fc_named_entrez)
  valid <- !is.na(fc_vals) & is.finite(fc_vals) & !is.na(names(fc_vals)) & names(fc_vals) != ""
  fc_named_entrez <- fc_vals[valid]

  if (length(fc_named_entrez) == 0) {
    warning("No mapped ENTREZ log2FC values available for SBGNview rendering.")
    return(invisible(FALSE))
  }

  out_dir <- dirname(out_prefix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sbgn_cache <- file.path(out_dir, "sbgnml_cache")
  dir.create(sbgn_cache, recursive = TRUE, showWarnings = FALSE)

  sbgn_obj <- tryCatch(
    SBGNview::SBGNview(
      gene.data = fc_named_entrez,
      input.sbgn = pathway_id,
      gene.id.type = "entrez",
      org = org_code,
      output.file = out_prefix,
      output.formats = formats,
      sbgn.dir = sbgn_cache
    ),
    error = function(e) {
      warning("SBGNview rendering failed for ", out_prefix, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(sbgn_obj)) {
    return(invisible(FALSE))
  }

  print(sbgn_obj)
  invisible(TRUE)
}
