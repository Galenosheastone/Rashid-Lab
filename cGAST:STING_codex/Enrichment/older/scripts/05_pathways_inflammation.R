source(file.path("R", "pathways_inflammation.R"))

.matrix_to_tsv_df <- function(mat, row_col = "feature") {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(data.frame())
  }
  out <- data.frame(row_id = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
  names(out)[1] <- row_col
  out
}

.default_de_list_from_tbl <- function(env = parent.frame()) {
  if (!exists("tbl", envir = env, inherits = TRUE)) {
    return(NULL)
  }

  de_tbl <- as.data.frame(get("tbl", envir = env, inherits = TRUE))
  if (!"gene_id" %in% names(de_tbl)) {
    if ("SYMBOL" %in% names(de_tbl)) {
      de_tbl$gene_id <- ifelse(is.na(de_tbl$SYMBOL) | de_tbl$SYMBOL == "", de_tbl$gene, de_tbl$SYMBOL)
    } else if ("gene" %in% names(de_tbl)) {
      de_tbl$gene_id <- as.character(de_tbl$gene)
    } else if (!is.null(rownames(de_tbl))) {
      de_tbl$gene_id <- rownames(de_tbl)
    }
  }

  if (!"log2FC" %in% names(de_tbl) && "logFC" %in% names(de_tbl)) {
    de_tbl$log2FC <- de_tbl$logFC
  }
  if (!"pvalue" %in% names(de_tbl) && "padj" %in% names(de_tbl)) {
    de_tbl$pvalue <- de_tbl$padj
  }

  contrast_id <- "contrast_1"
  if (exists("input_file", envir = env, inherits = TRUE)) {
    contrast_id <- tools::file_path_sans_ext(basename(get("input_file", envir = env, inherits = TRUE)))
  }

  out <- list()
  out[[contrast_id]] <- de_tbl
  out
}

.prepare_meta <- function(meta, expr_norm) {
  if (is.null(meta)) {
    return(NULL)
  }
  meta <- as.data.frame(meta)
  if (is.null(rownames(meta)) && !is.null(expr_norm) && nrow(meta) == ncol(expr_norm)) {
    rownames(meta) <- colnames(expr_norm)
  }
  meta
}

.resolve_objects <- function(expr_norm, meta, de_list, env = parent.frame()) {
  if (is.null(expr_norm) && exists("expr_norm", envir = env, inherits = TRUE)) {
    expr_norm <- get("expr_norm", envir = env, inherits = TRUE)
  }
  if (is.null(meta) && exists("meta", envir = env, inherits = TRUE)) {
    meta <- get("meta", envir = env, inherits = TRUE)
  }
  if (is.null(de_list) && exists("de_list", envir = env, inherits = TRUE)) {
    de_list <- get("de_list", envir = env, inherits = TRUE)
  }
  if (is.null(de_list)) {
    de_list <- .default_de_list_from_tbl(env = env)
  }

  list(
    expr_norm = expr_norm,
    meta = meta,
    de_list = de_list
  )
}

.guess_organism <- function(cfg, expr_norm, de_list, config_exists) {
  ids <- character(0)
  if (!is.null(expr_norm) && !is.null(rownames(expr_norm))) {
    ids <- c(ids, head(rownames(expr_norm), 5000))
  }
  if (!is.null(de_list) && length(de_list) > 0) {
    first_tbl <- as.data.frame(de_list[[1]])
    gene_col <- intersect(c("gene_symbol", "gene_id", "gene", "symbol", "SYMBOL", "ENSEMBL"), names(first_tbl))
    if (length(gene_col) > 0) {
      ids <- c(ids, head(as.character(first_tbl[[gene_col[1]]]), 5000))
    } else if (!is.null(rownames(first_tbl))) {
      ids <- c(ids, head(rownames(first_tbl), 5000))
    }
  }

  guessed <- .infer_organism_from_ids(ids)
  cfg_org <- .readable_organism(cfg$organism)

  if (!is.null(guessed) && !identical(cfg_org, guessed)) {
    if (isTRUE(config_exists)) {
      warning(
        "Config organism (", cfg_org, ") does not match detected IDs (", guessed, "). ",
        "Using detected organism for this run."
      )
    }
    return(guessed)
  }

  if (isTRUE(config_exists) && !is.null(cfg_org) && nzchar(cfg_org)) {
    return(cfg_org)
  }

  guessed %||% cfg_org
}

run_pathway_module <- function(
  expr_norm = NULL,
  meta = NULL,
  de_list = NULL,
  config_path = file.path("config", "pathways.yml"),
  output_root = "output"
) {
  cfg <- read_pathway_config(config_path = config_path)
  if (!isTRUE(cfg$run_module %||% TRUE)) {
    message("Pathway module disabled by config (run_module: FALSE).")
    return(invisible(NULL))
  }

  resolved <- .resolve_objects(expr_norm = expr_norm, meta = meta, de_list = de_list, env = parent.frame())
  expr_norm <- resolved$expr_norm
  meta <- resolved$meta
  de_list <- resolved$de_list

  if (is.null(de_list) || length(de_list) == 0) {
    warning("Pathway module: de_list is missing; skipping module run.")
    return(invisible(NULL))
  }

  if (is.null(names(de_list)) || any(names(de_list) == "")) {
    names(de_list) <- paste0("contrast_", seq_along(de_list))
  }

  meta <- .prepare_meta(meta = meta, expr_norm = expr_norm)

  cfg$organism <- .guess_organism(
    cfg = cfg,
    expr_norm = expr_norm,
    de_list = de_list,
    config_exists = file.exists(config_path)
  )

  message("Pathway module organism: ", cfg$organism)

  genesets <- load_genesets_msigdb(
    organism = cfg$organism,
    categories = cfg$msigdb_categories,
    include_gene_sets = cfg$include_gene_sets
  )
  genesets <- add_custom_sets(genesets, cfg$custom_sets)

  if (length(genesets) == 0) {
    warning("No pathway gene sets available after loading/filtering; stopping pathway module.")
    return(invisible(NULL))
  }

  pathway_root <- file.path(output_root, "pathways")
  results_dir <- file.path(pathway_root, "results")
  figures_dir <- file.path(pathway_root, "figures")
  dir.create(file.path(results_dir, "gsea"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(results_dir, "ssgsea"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(results_dir, "decoupleR"), recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  expr_sym <- NULL
  if (!is.null(expr_norm)) {
    harmonized <- ensure_symbols(expr_norm = expr_norm, de_table = de_list[[1]], organism = cfg$organism)
    expr_sym <- harmonized$expr_sym
  }

  fgsea_out <- list()
  if (isTRUE(cfg$run_fgsea %||% TRUE)) {
    for (contrast in names(de_list)) {
      de_tbl <- as.data.frame(de_list[[contrast]])
      de_sym <- ensure_symbols(expr_norm = NULL, de_table = de_tbl, organism = cfg$organism)$de_sym
      ranks <- make_ranks(de_sym)
      fgsea_res <- run_fgsea_for_contrast(ranks = ranks, genesets = genesets, fgsea_params = cfg$fgsea)

      contrast_file <- sanitize_filename(contrast)
      fgsea_tbl <- as.data.frame(fgsea_res$results_tbl)
      leading_tbl <- as.data.frame(fgsea_res$leading_edge_long)
      leading_plot_tbl <- leading_tbl

      if (nrow(leading_plot_tbl) > 0 && "pathway" %in% names(leading_plot_tbl)) {
        max_le <- as.integer(cfg$plots$leading_edge_genes %||% 30)
        if (is.finite(max_le) && max_le > 0) {
          split_tbl <- split(leading_plot_tbl, leading_plot_tbl$pathway)
          split_tbl <- lapply(split_tbl, function(x) utils::head(x, n = max_le))
          leading_plot_tbl <- do.call(rbind, split_tbl)
        }
      }

      safe_write_tsv(
        if (nrow(fgsea_tbl) > 0) fgsea_tbl else data.frame(),
        file.path(results_dir, "gsea", paste0(contrast_file, "_fgsea.tsv"))
      )

      if (nrow(leading_tbl) > 0) {
        safe_write_tsv(
          leading_tbl,
          file.path(results_dir, "gsea", paste0(contrast_file, "_leading_edge.tsv"))
        )
      }

      plot_fgsea_dotplot(
        fgsea_tbl = fgsea_tbl,
        out_png = file.path(figures_dir, paste0(contrast_file, "_fgsea_dotplot.png")),
        top_n = as.integer(cfg$plots$top_pathways %||% 20)
      )

      if (!is.null(expr_sym) && nrow(leading_plot_tbl) > 0) {
        plot_leading_edge_heatmap(
          expr_sym = expr_sym,
          meta = meta,
          leading_edge_long = leading_plot_tbl,
          out_png = file.path(figures_dir, paste0(contrast_file, "_leading_edge_heatmap.png"))
        )
      }

      fgsea_out[[contrast]] <- fgsea_tbl
    }
  }

  ssgsea_scores <- matrix(numeric(0), nrow = 0, ncol = 0)
  if (isTRUE(cfg$run_ssgsea %||% TRUE)) {
    if (is.null(expr_sym) || nrow(expr_sym) == 0 || ncol(expr_sym) == 0) {
      warning("ssGSEA requested but expr_norm is unavailable or empty; skipping ssGSEA.")
    } else {
      ssgsea_scores <- run_ssgsea(expr_sym = expr_sym, genesets = genesets, ssgsea_params = cfg$ssgsea)
      if (nrow(ssgsea_scores) > 0) {
        safe_write_tsv(
          .matrix_to_tsv_df(ssgsea_scores, row_col = "pathway"),
          file.path(results_dir, "ssgsea", "ssgsea_scores.tsv")
        )

        plot_score_heatmap(
          score_mat = ssgsea_scores,
          meta = meta,
          out_png = file.path(figures_dir, "ssgsea_heatmap.png"),
          title = "ssGSEA pathway scores"
        )

        if (!is.null(meta) && "time" %in% names(meta)) {
          default_paths <- c(
            "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE",
            "HALLMARK_FATTY_ACID_METABOLISM",
            "DAMP_ALARMS",
            "DAMP_SENSING_TLR_RAGE",
            "INFLAMMASOME_CORE",
            "PYROPTOSIS_CORE",
            "NECROPTOSIS_CORE",
            "FERROPTOSIS_CORE",
            "CGAS_STING_CORE"
          )
          selected_paths <- intersect(unique(c(cfg$include_gene_sets, default_paths)), rownames(ssgsea_scores))
          selected_paths <- utils::head(selected_paths, n = 12)

          for (pw in selected_paths) {
            plot_trajectory(
              score_mat = ssgsea_scores,
              meta = meta,
              pathway = pw,
              out_png = file.path(figures_dir, paste0("trajectory_", sanitize_filename(pw), ".png")),
              time_col = "time",
              group_col = if ("group" %in% names(meta)) "group" else "time"
            )
          }
        }
      }
    }
  }

  progeny_activities <- matrix(numeric(0), nrow = 0, ncol = 0)
  tf_activities <- matrix(numeric(0), nrow = 0, ncol = 0)
  if (isTRUE(cfg$run_decoupleR %||% TRUE)) {
    if (is.null(expr_sym) || nrow(expr_sym) == 0 || ncol(expr_sym) == 0) {
      warning("decoupleR requested but expr_norm is unavailable or empty; skipping activity inference.")
    } else {
      progeny_activities <- run_progeny(expr_sym = expr_sym, organism = cfg$organism)
      if (nrow(progeny_activities) > 0) {
        safe_write_tsv(
          .matrix_to_tsv_df(progeny_activities, row_col = "pathway"),
          file.path(results_dir, "decoupleR", "progeny_activities.tsv")
        )
        plot_score_heatmap(
          score_mat = progeny_activities,
          meta = meta,
          out_png = file.path(figures_dir, "progeny_activity_heatmap.png"),
          title = "PROGENy pathway activity"
        )
      }

      tf_activities <- run_tf_activity(expr_sym = expr_sym, organism = cfg$organism, conf_levels = c("A", "B", "C"))
      if (nrow(tf_activities) > 0) {
        safe_write_tsv(
          .matrix_to_tsv_df(tf_activities, row_col = "tf"),
          file.path(results_dir, "decoupleR", "tf_activities.tsv")
        )
        plot_score_heatmap(
          score_mat = tf_activities,
          meta = meta,
          out_png = file.path(figures_dir, "tf_activity_heatmap.png"),
          title = "DoRothEA TF activity"
        )
      }
    }
  }

  manifest_files <- list.files(pathway_root, recursive = TRUE, full.names = TRUE)
  manifest_tbl <- data.frame(file = manifest_files, stringsAsFactors = FALSE)
  safe_write_tsv(manifest_tbl, file.path(pathway_root, "pathways_manifest.tsv"))

  invisible(list(
    config = cfg,
    genesets = genesets,
    fgsea = fgsea_out,
    ssgsea = ssgsea_scores,
    progeny = progeny_activities,
    tf = tf_activities,
    output_root = pathway_root
  ))
}

if (identical(sys.nframe(), 0L)) {
  run_pathway_module()
}
