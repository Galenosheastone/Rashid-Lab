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

.infer_organism_from_ids <- function(ids) {
  ids <- unique(stats::na.omit(as.character(ids)))
  if (length(ids) == 0) {
    return(NULL)
  }
  if (any(grepl("^ENSMUSG", ids, ignore.case = TRUE))) {
    return("Mus musculus")
  }
  if (any(grepl("^ENSGALG", ids, ignore.case = TRUE))) {
    return("Gallus gallus")
  }
  if (any(grepl("^ENSG", ids, ignore.case = TRUE))) {
    return("Homo sapiens")
  }
  NULL
}

.readable_organism <- function(organism) {
  organism <- tolower(trimws(as.character(organism %||% "")))
  if (organism %in% c("mus musculus", "mouse", "mm", "m.musculus")) {
    return("Mus musculus")
  }
  if (organism %in% c("homo sapiens", "human", "hs", "h.sapiens")) {
    return("Homo sapiens")
  }
  if (organism %in% c("gallus gallus", "chicken", "gga")) {
    return("Gallus gallus")
  }
  if (nchar(organism) > 0) {
    return(tools::toTitleCase(organism))
  }
  "Mus musculus"
}

.to_decoupler_organism <- function(organism) {
  org <- .readable_organism(organism)
  if (identical(org, "Mus musculus")) {
    return("mouse")
  }
  if (identical(org, "Homo sapiens")) {
    return("human")
  }
  NULL
}

.to_orgdb_package <- function(organism) {
  org <- .readable_organism(organism)
  if (identical(org, "Mus musculus")) {
    return("org.Mm.eg.db")
  }
  if (identical(org, "Homo sapiens")) {
    return("org.Hs.eg.db")
  }
  if (identical(org, "Gallus gallus")) {
    return("org.Gg.eg.db")
  }
  NULL
}

read_pathway_config <- function(config_path = NULL) {
  defaults <- list(
    run_module = TRUE,
    organism = "Mus musculus",
    msigdb_categories = c("H"),
    include_gene_sets = c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_FATTY_ACID_METABOLISM"
    ),
    custom_sets = list(
      DAMP_ALARMS = c("HMGB1", "S100A8", "S100A9", "IL33", "HSP90AA1", "HSPA1A", "HSPA1B"),
      DAMP_SENSING_TLR_RAGE = c("TLR2", "TLR4", "AGER", "MYD88", "TICAM1", "IRAK1", "IRAK4", "TRAF6"),
      INFLAMMASOME_CORE = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD"),
      PYROPTOSIS_CORE = c("CASP1", "CASP4", "CASP11", "GSDMD", "IL1B"),
      NECROPTOSIS_CORE = c("RIPK1", "RIPK3", "MLKL", "TNF", "FADD", "CASP8"),
      FERROPTOSIS_CORE = c("ACSL4", "SLC7A11", "GPX4", "FTH1", "TFRC"),
      CGAS_STING_CORE = c("MB21D1", "TMEM173", "TBK1", "IRF3", "IFNB1")
    ),
    fgsea = list(
      nperm = 10000,
      minSize = 10,
      maxSize = 500,
      seed = 666
    ),
    ssgsea = list(
      normalize = TRUE
    ),
    plots = list(
      top_pathways = 20,
      leading_edge_genes = 30
    ),
    run_decoupleR = TRUE,
    run_ssgsea = TRUE,
    run_fgsea = TRUE
  )

  cfg_path <- config_path %||% file.path("config", "pathways.yml")
  cfg <- defaults

  if (file.exists(cfg_path)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      warning("yaml is not installed; using default pathway config.")
    } else {
      loaded <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) {
        warning("Failed to read config at ", cfg_path, ": ", conditionMessage(e))
        NULL
      })
      cfg <- .merge_lists(defaults, loaded)
    }
  }

  cfg$organism <- .readable_organism(cfg$organism)
  cfg$msigdb_categories <- as.character(unlist(cfg$msigdb_categories %||% c("H")))
  cfg$include_gene_sets <- as.character(unlist(cfg$include_gene_sets %||% character(0)))
  cfg
}

.load_orgdb <- function(organism) {
  pkg <- .to_orgdb_package(organism)
  if (is.null(pkg)) {
    return(NULL)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    warning("AnnotationDbi is not installed; cannot map IDs to SYMBOL.")
    return(NULL)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning("", pkg, " is not installed; proceeding without ID mapping.")
    return(NULL)
  }
  tryCatch(getExportedValue(pkg, pkg), error = function(e) {
    warning("Could not load OrgDb object from ", pkg, ": ", conditionMessage(e))
    NULL
  })
}

.detect_keytype <- function(ids) {
  ids <- stats::na.omit(as.character(ids))
  if (length(ids) == 0) {
    return(NULL)
  }
  if (all(grepl("^ENS[A-Z]*G", ids, ignore.case = TRUE))) {
    return("ENSEMBL")
  }
  NULL
}

map_ids_to_symbols <- function(ids, organism) {
  ids <- as.character(ids)
  ids <- sub("\\.\\d+$", "", ids)
  if (length(ids) == 0) {
    return(ids)
  }

  keytype <- .detect_keytype(ids)
  if (is.null(keytype)) {
    return(ids)
  }

  orgdb <- .load_orgdb(organism)
  if (is.null(orgdb)) {
    return(ids)
  }

  mapped <- tryCatch(
    AnnotationDbi::select(
      orgdb,
      keys = unique(ids),
      columns = c("SYMBOL"),
      keytype = keytype
    ),
    error = function(e) {
      warning("ID mapping failed; using original IDs: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(mapped) || nrow(mapped) == 0 || !("SYMBOL" %in% names(mapped))) {
    return(ids)
  }

  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", , drop = FALSE]
  if (nrow(mapped) == 0) {
    return(ids)
  }

  idx <- match(ids, mapped[[keytype]])
  symbols <- mapped$SYMBOL[idx]
  symbols[is.na(symbols) | symbols == ""] <- ids[is.na(symbols) | symbols == ""]
  symbols
}

.aggregate_rows_mean <- function(mat, groups) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  groups <- as.character(groups)

  mat0 <- mat
  mat0[is.na(mat0)] <- 0
  sums <- rowsum(mat0, group = groups, reorder = FALSE)

  counts <- rowsum(!is.na(mat), group = groups, reorder = FALSE)
  counts <- pmax(counts, 1)

  sums / counts
}

.detect_de_gene_col <- function(de_table) {
  if (is.null(de_table)) {
    return(NULL)
  }
  nms <- names(de_table)
  prefs <- c("gene_symbol", "gene_id", "gene", "symbol", "SYMBOL", "ensembl", "ENSEMBL")
  hit <- prefs[prefs %in% nms]
  if (length(hit) > 0) {
    return(hit[1])
  }
  NULL
}

ensure_symbols <- function(expr_norm, de_table, organism) {
  expr_sym <- expr_norm
  de_sym <- de_table

  if (!is.null(expr_norm)) {
    expr_sym <- as.matrix(expr_norm)
    if (is.null(rownames(expr_sym))) {
      warning("expr_norm has no rownames; cannot map to SYMBOL.")
    } else {
      mapped <- map_ids_to_symbols(rownames(expr_sym), organism)
      keep <- !is.na(mapped) & mapped != ""
      expr_sym <- expr_sym[keep, , drop = FALSE]
      mapped <- mapped[keep]
      expr_sym <- .aggregate_rows_mean(expr_sym, mapped)
      expr_sym <- expr_sym[!is.na(rownames(expr_sym)) & rownames(expr_sym) != "", , drop = FALSE]
    }
  }

  if (!is.null(de_table)) {
    gene_col <- .detect_de_gene_col(de_table)
    if (is.null(gene_col)) {
      if (!is.null(rownames(de_table))) {
        de_sym$gene_id <- rownames(de_table)
        gene_col <- "gene_id"
      }
    }
    if (!is.null(gene_col)) {
      gene_ids <- map_ids_to_symbols(de_sym[[gene_col]], organism)
      de_sym$gene_symbol <- gene_ids
    } else {
      warning("Could not detect a gene ID column in DE table; rankings may fail.")
    }
  }

  list(expr_sym = expr_sym, de_sym = de_sym)
}

load_genesets_msigdb <- function(organism, categories, include_gene_sets = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    warning("msigdbr is not installed; no MSigDB pathways loaded.")
    return(list())
  }

  org <- .readable_organism(organism)
  cats <- unique(as.character(categories %||% "H"))
  msig_formals <- names(formals(msigdbr::msigdbr))
  use_collection <- "collection" %in% msig_formals

  msig_list <- lapply(cats, function(cat) {
    tryCatch(
      if (isTRUE(use_collection)) {
        msigdbr::msigdbr(species = org, collection = cat)
      } else {
        msigdbr::msigdbr(species = org, category = cat)
      },
      error = function(e) {
        warning("msigdbr failed for category ", cat, ": ", conditionMessage(e))
        NULL
      }
    )
  })

  msig_df <- do.call(rbind, msig_list)
  if (is.null(msig_df) || nrow(msig_df) == 0) {
    warning("No MSigDB sets loaded for organism ", org, ".")
    return(list())
  }

  if (!("gs_name" %in% names(msig_df)) || !("gene_symbol" %in% names(msig_df))) {
    warning("msigdbr output is missing required columns gs_name/gene_symbol.")
    return(list())
  }

  msig_df <- msig_df[!is.na(msig_df$gs_name) & !is.na(msig_df$gene_symbol), , drop = FALSE]

  sets <- split(as.character(msig_df$gene_symbol), as.character(msig_df$gs_name))
  sets <- lapply(sets, function(x) unique(stats::na.omit(as.character(x))))
  sets <- sets[lengths(sets) > 0]

  include <- as.character(unlist(include_gene_sets %||% character(0)))
  if (length(include) > 0) {
    missing_sets <- setdiff(include, names(sets))
    if (length(missing_sets) > 0) {
      warning(
        "Requested include_gene_sets were not found for organism ", org, ": ",
        paste(missing_sets, collapse = ", ")
      )
    }
    sets <- sets[intersect(include, names(sets))]
  }

  sets
}

add_custom_sets <- function(genesets, custom_sets) {
  out <- genesets %||% list()
  if (is.null(custom_sets) || length(custom_sets) == 0) {
    return(out)
  }

  for (nm in names(custom_sets)) {
    vals <- unique(stats::na.omit(as.character(custom_sets[[nm]])))
    if (length(vals) == 0) {
      next
    }
    out[[nm]] <- vals
  }
  out
}

.detect_col <- function(x, candidates) {
  nms <- names(x)
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) {
    return(NULL)
  }
  hit[1]
}

make_ranks <- function(de_tbl) {
  if (is.null(de_tbl) || nrow(de_tbl) == 0) {
    return(setNames(numeric(0), character(0)))
  }

  gene_col <- .detect_col(de_tbl, c("gene_symbol", "gene_id", "gene", "symbol", "SYMBOL", "ENSEMBL"))
  if (is.null(gene_col) && !is.null(rownames(de_tbl))) {
    de_tbl$gene_id <- rownames(de_tbl)
    gene_col <- "gene_id"
  }
  if (is.null(gene_col)) {
    warning("make_ranks: No gene column found in DE table.")
    return(setNames(numeric(0), character(0)))
  }

  stat_col <- .detect_col(de_tbl, c("stat", "t", "wald", "Wald", "signed_stat"))
  lfc_col <- .detect_col(de_tbl, c("log2FC", "logFC", "log2FoldChange", "lfc", "log2fc", "log2FC_shrunken", "log2fc_shrunken"))
  p_col <- .detect_col(de_tbl, c("pvalue", "pval", "padj", "adj.P.Val", "FDR"))

  ranks <- NULL
  if (!is.null(stat_col) && any(is.finite(as.numeric(de_tbl[[stat_col]])))) {
    ranks <- suppressWarnings(as.numeric(de_tbl[[stat_col]]))
  } else {
    if (is.null(lfc_col)) {
      warning("make_ranks: Missing stat and log2FC columns; cannot create ranks.")
      return(setNames(numeric(0), character(0)))
    }
    lfc <- suppressWarnings(as.numeric(de_tbl[[lfc_col]]))
    if (!is.null(p_col)) {
      pvals <- suppressWarnings(as.numeric(de_tbl[[p_col]]))
      pvals[is.na(pvals)] <- 1
      ranks <- sign(lfc) * abs(lfc) * (-log10(pmax(pvals, 1e-300)))
    } else {
      warning("make_ranks: p-value column missing; using signed log2FC fallback.")
      ranks <- lfc
    }
  }

  genes <- as.character(de_tbl[[gene_col]])
  genes <- sub("\\.\\d+$", "", genes)

  keep <- !is.na(genes) & genes != "" & is.finite(ranks)
  if (!any(keep)) {
    return(setNames(numeric(0), character(0)))
  }

  rank_df <- data.frame(
    gene = genes[keep],
    rank = as.numeric(ranks[keep]),
    stringsAsFactors = FALSE
  )

  rank_df <- rank_df[order(abs(rank_df$rank), decreasing = TRUE), , drop = FALSE]
  rank_df <- rank_df[!duplicated(rank_df$gene), , drop = FALSE]
  rank_df <- rank_df[order(rank_df$rank, decreasing = TRUE), , drop = FALSE]

  out <- rank_df$rank
  names(out) <- rank_df$gene
  out
}

.filter_genesets_by_ranks <- function(genesets, rank_names, min_size = 10, max_size = 500) {
  gs <- lapply(genesets, function(x) intersect(unique(as.character(x)), rank_names))
  sz <- lengths(gs)
  gs[sz >= min_size & sz <= max_size]
}

run_fgsea_for_contrast <- function(ranks, genesets, fgsea_params) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    warning("fgsea is not installed; skipping contrast-level fgsea.")
    return(list(
      results_tbl = data.frame(),
      leading_edge_long = data.frame()
    ))
  }

  min_size <- as.integer(fgsea_params$minSize %||% 10)
  max_size <- as.integer(fgsea_params$maxSize %||% 500)
  nperm <- as.integer(fgsea_params$nperm %||% 10000)
  seed <- as.integer(fgsea_params$seed %||% 666)

  ranks <- ranks[is.finite(ranks)]
  if (length(ranks) == 0) {
    warning("run_fgsea_for_contrast: Empty rank vector.")
    return(list(results_tbl = data.frame(), leading_edge_long = data.frame()))
  }

  gs_use <- .filter_genesets_by_ranks(
    genesets = genesets,
    rank_names = names(ranks),
    min_size = min_size,
    max_size = max_size
  )

  if (length(gs_use) == 0) {
    warning("No gene sets passed overlap/minSize/maxSize filters for fgsea.")
    return(list(results_tbl = data.frame(), leading_edge_long = data.frame()))
  }

  fgsea_fun <- fgsea::fgsea
  fgsea_args <- list(
    pathways = gs_use,
    stats = ranks,
    minSize = min_size,
    maxSize = max_size
  )
  if ("nperm" %in% names(formals(fgsea_fun))) {
    fgsea_args$nperm <- nperm
  }

  set.seed(seed)
  fgsea_res <- tryCatch(
    do.call(fgsea_fun, fgsea_args),
    error = function(e) {
      warning("fgsea failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(fgsea_res)) {
    return(list(results_tbl = data.frame(), leading_edge_long = data.frame()))
  }

  fgsea_tbl <- as.data.frame(fgsea_res)
  if (nrow(fgsea_tbl) == 0) {
    return(list(results_tbl = data.frame(), leading_edge_long = data.frame()))
  }

  fgsea_tbl <- fgsea_tbl[order(fgsea_tbl$padj, -abs(fgsea_tbl$NES)), , drop = FALSE]

  leading_edge_long <- data.frame()
  if ("leadingEdge" %in% names(fgsea_tbl)) {
    lead_list <- lapply(seq_len(nrow(fgsea_tbl)), function(i) {
      genes <- as.character(unlist(fgsea_tbl$leadingEdge[[i]]))
      if (length(genes) == 0) {
        return(NULL)
      }
      data.frame(
        pathway = fgsea_tbl$pathway[i],
        gene_symbol = genes,
        NES = fgsea_tbl$NES[i],
        padj = fgsea_tbl$padj[i],
        stringsAsFactors = FALSE
      )
    })
    lead_list <- lead_list[!vapply(lead_list, is.null, logical(1))]
    if (length(lead_list) > 0) {
      leading_edge_long <- do.call(rbind, lead_list)
    }
  }

  list(results_tbl = fgsea_tbl, leading_edge_long = leading_edge_long)
}

run_ssgsea <- function(expr_sym, genesets, ssgsea_params) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    warning("GSVA is not installed; skipping ssGSEA.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  expr_sym <- as.matrix(expr_sym)
  if (nrow(expr_sym) == 0 || ncol(expr_sym) == 0) {
    warning("run_ssgsea: empty expression matrix.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  gs_use <- lapply(genesets, function(x) intersect(unique(as.character(x)), rownames(expr_sym)))
  gs_use <- gs_use[lengths(gs_use) > 0]
  if (length(gs_use) == 0) {
    warning("run_ssgsea: no gene-set overlaps with expression matrix.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  normalize <- isTRUE(ssgsea_params$normalize %||% TRUE)

  score_mat <- tryCatch({
    if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
      param <- GSVA::ssgseaParam(exprData = expr_sym, geneSets = gs_use, normalize = normalize)
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(
        expr = expr_sym,
        gset.idx.list = gs_use,
        method = "ssgsea",
        ssgsea.norm = normalize,
        verbose = FALSE
      )
    }
  }, error = function(e) {
    warning("ssGSEA failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(score_mat)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  as.matrix(score_mat)
}

.harmonize_network <- function(network) {
  network <- as.data.frame(network)
  nms <- names(network)

  if (!("source" %in% nms)) {
    src <- intersect(c("tf", "pathway", "regulator", "source"), nms)
    if (length(src) > 0) {
      network$source <- network[[src[1]]]
    }
  }

  if (!("target" %in% nms)) {
    tgt <- intersect(c("target", "gene", "genesymbol", "symbol"), nms)
    if (length(tgt) > 0) {
      network$target <- network[[tgt[1]]]
    }
  }

  if (!("weight" %in% nms)) {
    w <- intersect(c("weight", "mor", "likelihood"), nms)
    if (length(w) > 0) {
      network$weight <- suppressWarnings(as.numeric(network[[w[1]]]))
    } else {
      network$weight <- 1
    }
  }

  if (!("source" %in% names(network)) || !("target" %in% names(network))) {
    return(data.frame())
  }

  network <- network[, c("source", "target", "weight"), drop = FALSE]
  network$source <- as.character(network$source)
  network$target <- as.character(network$target)
  network$weight <- suppressWarnings(as.numeric(network$weight))
  network$weight[is.na(network$weight)] <- 1

  network <- network[
    !is.na(network$source) & network$source != "" &
      !is.na(network$target) & network$target != "",
    , drop = FALSE
  ]
  unique(network)
}

.decouple_to_matrix <- function(res, sample_order) {
  if (is.null(res)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  if (is.matrix(res)) {
    mat <- as.matrix(res)
    keep <- intersect(sample_order, colnames(mat))
    if (length(keep) > 0) {
      mat <- mat[, keep, drop = FALSE]
    }
    return(mat)
  }

  if (!is.data.frame(res)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  x <- as.data.frame(res)
  names(x) <- sub("^\\.", "", names(x))

  source_col <- intersect(c("source", "pathway", "tf", "regulator"), names(x))
  cond_col <- intersect(c("condition", "sample"), names(x))
  value_col <- intersect(c("score", "statistic", "estimate", "activity", "value"), names(x))

  if (length(source_col) == 0 || length(cond_col) == 0 || length(value_col) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  source_col <- source_col[1]
  cond_col <- cond_col[1]
  value_col <- value_col[1]

  x <- x[
    !is.na(x[[source_col]]) & !is.na(x[[cond_col]]) & is.finite(as.numeric(x[[value_col]])),
    , drop = FALSE
  ]
  if (nrow(x) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  vals <- tapply(
    as.numeric(x[[value_col]]),
    list(as.character(x[[source_col]]), as.character(x[[cond_col]])),
    mean,
    na.rm = TRUE
  )

  mat <- as.matrix(vals)
  if (length(sample_order) > 0 && !is.null(colnames(mat))) {
    keep <- intersect(sample_order, colnames(mat))
    if (length(keep) > 0) {
      mat <- mat[, keep, drop = FALSE]
    }
  }
  mat
}

.run_decouple_method <- function(expr_sym, network, method_name) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    warning("decoupleR is not installed.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  method_fun <- tryCatch(getExportedValue("decoupleR", method_name), error = function(e) NULL)
  if (is.null(method_fun)) {
    warning("decoupleR method not available: ", method_name)
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  fml <- names(formals(method_fun))
  args <- list()

  if ("mat" %in% fml) args$mat <- as.matrix(expr_sym)
  if ("x" %in% fml) args$x <- as.matrix(expr_sym)
  if ("network" %in% fml) args$network <- network
  if ("net" %in% fml) args$net <- network
  if (".source" %in% fml) args$.source <- "source"
  if ("source" %in% fml) args$source <- "source"
  if (".target" %in% fml) args$.target <- "target"
  if ("target" %in% fml) args$target <- "target"
  if (".mor" %in% fml) args$.mor <- "weight"
  if ("mor" %in% fml) args$mor <- "weight"
  if ("minsize" %in% fml) args$minsize <- 5L
  if ("times" %in% fml) args$times <- 0L
  if ("verbose" %in% fml) args$verbose <- FALSE

  res <- tryCatch(do.call(method_fun, args), error = function(e) {
    warning("decoupleR::", method_name, " failed: ", conditionMessage(e))
    NULL
  })

  .decouple_to_matrix(res, sample_order = colnames(expr_sym))
}

run_progeny <- function(expr_sym, organism) {
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    warning("decoupleR is not installed; skipping PROGENy activity inference.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  dec_org <- .to_decoupler_organism(organism)
  if (is.null(dec_org)) {
    warning("PROGENy currently supported here for mouse/human only; organism was ", organism)
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  net <- tryCatch(decoupleR::get_progeny(organism = dec_org, top = 500), error = function(e) {
    warning("Could not load PROGENy network: ", conditionMessage(e))
    NULL
  })

  if (is.null(net)) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  net <- .harmonize_network(net)
  if (nrow(net) == 0) {
    warning("PROGENy network is empty after harmonization.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  methods <- getNamespaceExports("decoupleR")
  method_name <- if ("run_wmean" %in% methods) "run_wmean" else if ("run_mlm" %in% methods) "run_mlm" else NULL
  if (is.null(method_name)) {
    warning("No supported decoupleR method found for PROGENy.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  .run_decouple_method(expr_sym, net, method_name)
}

run_tf_activity <- function(expr_sym, organism, conf_levels = c("A", "B", "C")) {
  if (!requireNamespace("decoupleR", quietly = TRUE) || !requireNamespace("dorothea", quietly = TRUE)) {
    warning("decoupleR and/or dorothea is missing; skipping TF activity inference.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  org <- .readable_organism(organism)
  regulons <- NULL
  if (identical(org, "Mus musculus")) {
    regulons <- tryCatch(dorothea::dorothea_mm, error = function(e) NULL)
  }
  if (identical(org, "Homo sapiens")) {
    regulons <- tryCatch(dorothea::dorothea_hs, error = function(e) NULL)
  }

  if (is.null(regulons)) {
    warning("DoRothEA regulons are only configured here for mouse/human; organism was ", org)
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  reg <- as.data.frame(regulons)
  if ("confidence" %in% names(reg)) {
    reg <- reg[reg$confidence %in% conf_levels, , drop = FALSE]
  }

  reg <- .harmonize_network(reg)
  if (nrow(reg) == 0) {
    warning("DoRothEA regulon network is empty after filtering.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  methods <- getNamespaceExports("decoupleR")
  method_name <- if ("run_viper" %in% methods) "run_viper" else if ("run_wmean" %in% methods) "run_wmean" else NULL
  if (is.null(method_name)) {
    warning("No supported decoupleR method found for TF activity.")
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  tf_mat <- .run_decouple_method(expr_sym, reg, method_name)
  if (nrow(tf_mat) == 0 && method_name != "run_wmean" && "run_wmean" %in% methods) {
    tf_mat <- .run_decouple_method(expr_sym, reg, "run_wmean")
  }
  tf_mat
}

plot_fgsea_dotplot <- function(fgsea_tbl, out_png, top_n = 20) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping fgsea dotplot.")
    return(invisible(NULL))
  }

  df <- as.data.frame(fgsea_tbl)
  needed <- c("pathway", "NES", "padj")
  if (nrow(df) == 0 || !all(needed %in% names(df))) {
    return(invisible(NULL))
  }

  df <- df[is.finite(df$NES) & is.finite(df$padj), , drop = FALSE]
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }

  df <- df[order(df$padj, -abs(df$NES)), , drop = FALSE]
  df <- utils::head(df, n = min(top_n, nrow(df)))
  df$pathway <- factor(df$pathway, levels = rev(df$pathway))
  df$direction <- ifelse(df$NES >= 0, "Activated", "Suppressed")
  df$neglog10_padj <- -log10(pmax(df$padj, 1e-300))

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = NES, y = pathway, size = neglog10_padj, color = direction)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_manual(values = c(Activated = "#B2182B", Suppressed = "#2166AC")) +
    ggplot2::labs(
      x = "NES",
      y = NULL,
      size = "-log10(padj)",
      color = "Direction",
      title = "fgsea inflammatory/immunometabolic pathways"
    ) +
    ggplot2::theme_bw(base_size = 11)

  ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  invisible(out_png)
}

.build_col_annotation <- function(meta, sample_ids) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(meta) || nrow(meta) == 0 || is.null(rownames(meta))) {
    return(NULL)
  }

  keep <- intersect(sample_ids, rownames(meta))
  if (length(keep) == 0) {
    return(NULL)
  }

  ann_cols <- intersect(c("group", "time", "batch"), names(meta))
  if (length(ann_cols) == 0) {
    return(NULL)
  }

  ann_df <- meta[keep, ann_cols, drop = FALSE]
  if (nrow(ann_df) == 0 || ncol(ann_df) == 0) {
    return(NULL)
  }

  ComplexHeatmap::HeatmapAnnotation(df = ann_df)
}

.zscore_rows <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(mat)
  }
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  z
}

plot_leading_edge_heatmap <- function(expr_sym, meta, leading_edge_long, out_png) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    warning("ComplexHeatmap/circlize missing; skipping leading-edge heatmap.")
    return(invisible(NULL))
  }

  expr_sym <- as.matrix(expr_sym)
  if (nrow(expr_sym) == 0 || ncol(expr_sym) == 0) {
    return(invisible(NULL))
  }

  le <- as.data.frame(leading_edge_long)
  if (nrow(le) == 0 || !("pathway" %in% names(le)) || !("gene_symbol" %in% names(le))) {
    return(invisible(NULL))
  }

  if ("padj" %in% names(le)) {
    pathway_order <- unique(le$pathway[order(le$padj)])
  } else {
    pathway_order <- unique(le$pathway)
  }
  pathway_order <- utils::head(pathway_order, n = min(10, length(pathway_order)))

  keep_genes <- unique(unlist(lapply(pathway_order, function(pw) {
    g <- unique(le$gene_symbol[le$pathway == pw])
    utils::head(g, n = 30)
  })))

  keep_genes <- intersect(keep_genes, rownames(expr_sym))
  if (length(keep_genes) < 2) {
    return(invisible(NULL))
  }

  mat <- expr_sym[keep_genes, , drop = FALSE]
  mat <- .zscore_rows(mat)

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  top_ann <- .build_col_annotation(meta, colnames(mat))
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

  grDevices::png(out_png, width = 1800, height = 1200, res = 200)
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    top_annotation = top_ann,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 7)
  )
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  invisible(out_png)
}

plot_score_heatmap <- function(score_mat, meta, out_png, title) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || !requireNamespace("circlize", quietly = TRUE)) {
    warning("ComplexHeatmap/circlize missing; skipping score heatmap.")
    return(invisible(NULL))
  }

  mat <- as.matrix(score_mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(invisible(NULL))
  }

  mat <- mat[rowSums(is.finite(mat)) > 0, , drop = FALSE]
  if (nrow(mat) == 0) {
    return(invisible(NULL))
  }

  mat <- .zscore_rows(mat)

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  top_ann <- .build_col_annotation(meta, colnames(mat))
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))

  grDevices::png(out_png, width = 1800, height = 1200, res = 200)
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "z-score",
    col = col_fun,
    top_annotation = top_ann,
    column_title = title,
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 7)
  )
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
  invisible(out_png)
}

sanitize_filename <- function(x) {
  out <- gsub("[^A-Za-z0-9_]+", "_", x)
  out <- gsub("_+", "_", out)
  gsub("^_|_$", "", out)
}

plot_trajectory <- function(score_mat, meta, pathway, out_png, time_col = "time", group_col = "group") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed; skipping trajectory plot.")
    return(invisible(NULL))
  }

  mat <- as.matrix(score_mat)
  if (!(pathway %in% rownames(mat))) {
    return(invisible(NULL))
  }
  if (is.null(meta) || nrow(meta) == 0 || is.null(rownames(meta)) || !(time_col %in% names(meta))) {
    return(invisible(NULL))
  }

  samples <- intersect(colnames(mat), rownames(meta))
  if (length(samples) < 3) {
    return(invisible(NULL))
  }

  group_vals <- if (group_col %in% names(meta)) as.character(meta[samples, group_col]) else rep("all", length(samples))
  df <- data.frame(
    sample = samples,
    time = meta[samples, time_col],
    group = group_vals,
    score = as.numeric(mat[pathway, samples]),
    stringsAsFactors = FALSE
  )
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 3) {
    return(invisible(NULL))
  }

  time_num <- suppressWarnings(as.numeric(as.character(df$time)))
  has_numeric_time <- !all(is.na(time_num))

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  if (has_numeric_time) {
    df$time_numeric <- time_num
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time_numeric, y = score, color = group)) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::geom_smooth(method = "loess", se = TRUE) +
      ggplot2::labs(
        title = pathway,
        x = time_col,
        y = "Pathway score",
        color = group_col
      ) +
      ggplot2::theme_bw(base_size = 11)
  } else {
    df$time_factor <- factor(df$time, levels = unique(df$time))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time_factor, y = score, color = group, group = group)) +
      ggplot2::stat_summary(fun = mean, geom = "line", linewidth = 0.9) +
      ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.15) +
      ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.05), alpha = 0.6, size = 1.8) +
      ggplot2::labs(
        title = pathway,
        x = time_col,
        y = "Pathway score",
        color = group_col
      ) +
      ggplot2::theme_bw(base_size = 11)
  }

  ggplot2::ggsave(out_png, p, width = 8, height = 5, dpi = 300)
  invisible(out_png)
}

safe_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  df <- as.data.frame(x)
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_tsv(df, file = path, na = "")
  } else if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, file = path, sep = "\t", na = "")
  } else {
    utils::write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  }
  invisible(path)
}
