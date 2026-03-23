#!/usr/bin/env Rscript

# cGAS-STING pathway map from DESeq2 CSV (Gallus gallus, Ensembl IDs)
#
# How to run in RStudio:
# source("cgast_sting_map.R")
# run_cgast_sting_map(
#   input = "/path/to/DESeq2_results.csv",
#   outdir = "cgast_sting_out",
#   title = "Sacral vs Free (cGAS-STING map)"
# )
#
# Optional RStudio file-picker mode:
# source("cgast_sting_map.R")
# run_cgast_sting_map_interactive()
#
# How to run from command line:
# Rscript cgast_sting_map.R \
#   --input /path/to/DESeq2_results.csv \
#   --outdir cgast_sting_out \
#   --title "Sacral vs Free (cGAS-STING map)"

suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Gg.eg.db)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(svglite)
  library(grid)
})

print_help <- function() {
  cat(
    "Usage:\n",
    "  Rscript cgast_sting_map.R --input <csv> [--outdir <dir>] [--title <title>]\n",
    "\n",
    "RStudio usage:\n",
    "  source(\"cgast_sting_map.R\")\n",
    "  run_cgast_sting_map(input = \"...\", outdir = \"cgast_sting_out\", title = \"...\")\n",
    "  # or\n",
    "  run_cgast_sting_map_interactive()\n",
    sep = ""
  )
}

parse_args <- function(args) {
  out <- list(
    input = NULL,
    outdir = "cgast_sting_out",
    title = "cGAS-STING Pathway Map",
    help = FALSE
  )
  
  if (length(args) == 0L) {
    return(out)
  }
  
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      out$help <- TRUE
      return(out)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    
    value <- args[[i + 1L]]
    if (key == "--input") {
      out$input <- value
    } else if (key == "--outdir") {
      out$outdir <- value
    } else if (key == "--title") {
      out$title <- value
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    i <- i + 2L
  }
  
  out
}

get_pathway_spec <- function() {
  pathway_genes <- c(
    "MB21D1", "STING1", "DDX41", "IFI16", "ZBP1",
    "TRAF3", "TRAF6",
    "TBK1", "IKBKB", "CHUK", "IKBKG",
    "IRF3", "IRF7", "NFKB1", "RELA",
    "IFNB1", "IFNAR1", "IFNAR2", "JAK1", "TYK2", "STAT1", "STAT2", "IRF9",
    "ISG15", "OASL", "MX1", "IFIT5", "RSAD2", "CXCL10", "CCL5", "IL6", "TNF"
  )
  
  edge_df <- tribble(
    ~from, ~to,
    "DDX41", "STING1",
    "IFI16", "STING1",
    "MB21D1", "STING1",
    "ZBP1", "TBK1",
    "STING1", "TRAF3",
    "STING1", "TRAF6",
    "TRAF3", "TBK1",
    "TRAF6", "IKBKB",
    "IKBKG", "IKBKB",
    "STING1", "TBK1",
    "STING1", "IKBKB",
    "CHUK", "NFKB1",
    "IKBKB", "NFKB1",
    "IKBKB", "RELA",
    "TBK1", "IRF3",
    "TBK1", "IRF7",
    "IRF3", "IFNB1",
    "IRF7", "IFNB1",
    "IFNB1", "IFNAR1",
    "IFNB1", "IFNAR2",
    "IFNAR1", "JAK1",
    "IFNAR2", "TYK2",
    "JAK1", "STAT1",
    "TYK2", "STAT2",
    "STAT1", "IRF9",
    "STAT2", "IRF9",
    "IRF9", "ISG15",
    "IRF9", "MX1",
    "IRF9", "OASL",
    "IRF9", "IFIT5",
    "IRF9", "RSAD2",
    "RELA", "IL6",
    "RELA", "TNF",
    "IRF3", "CXCL10",
    "IRF7", "CCL5"
  )
  
  layout_df <- tribble(
    ~symbol,  ~x, ~y,  ~module,
    "DDX41",   0,  2.8, "Sensors",
    "IFI16",   0,  2.2, "Sensors",
    "MB21D1",  0,  1.6, "Sensors",
    "ZBP1",    0,  0.8, "Sensors",
    "STING1",  1,  1.8, "STING",
    "TRAF3",   2,  2.7, "Trafficking",
    "TRAF6",   2,  1.1, "Trafficking",
    "TBK1",    2,  2.1, "Kinases",
    "IKBKB",   2,  0.4, "Kinases",
    "IKBKG",   2, -0.2, "Kinases",
    "CHUK",    2, -0.8, "Kinases",
    "IRF7",    3,  2.9, "TFs",
    "IRF3",    3,  2.2, "TFs",
    "NFKB1",   3,  0.5, "TFs",
    "RELA",    3, -0.2, "TFs",
    "IFNB1",   4,  2.3, "IFN",
    "IFNAR1",  5,  2.8, "IFNAR/JAK-STAT",
    "IFNAR2",  5,  2.0, "IFNAR/JAK-STAT",
    "JAK1",    6,  2.8, "IFNAR/JAK-STAT",
    "TYK2",    6,  2.0, "IFNAR/JAK-STAT",
    "STAT1",   6,  3.3, "IFNAR/JAK-STAT",
    "STAT2",   6,  1.5, "IFNAR/JAK-STAT",
    "IRF9",    6,  2.3, "IFNAR/JAK-STAT",
    "CXCL10",  7,  3.4, "ISGs/Outputs",
    "CCL5",    7,  2.8, "ISGs/Outputs",
    "ISG15",   7,  2.4, "ISGs/Outputs",
    "MX1",     7,  2.0, "ISGs/Outputs",
    "OASL",    7,  1.6, "ISGs/Outputs",
    "IFIT5",   7,  1.2, "ISGs/Outputs",
    "RSAD2",   7,  0.8, "ISGs/Outputs",
    "IL6",     7,  0.1, "ISGs/Outputs",
    "TNF",     7, -0.5, "ISGs/Outputs"
  )
  
  list(pathway_genes = pathway_genes, edge_df = edge_df, layout_df = layout_df)
}

run_cgast_sting_map <- function(input,
                                outdir = "cgast_sting_out",
                                title = "cGAS-STING Pathway Map") {
  if (missing(input) || is.null(input) || !nzchar(input)) {
    stop("`input` is required.", call. = FALSE)
  }
  if (!file.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  padj_threshold <- 0.05
  
  message("Reading DESeq2 CSV: ", input)
  de_raw <- readr::read_csv(
    file = input,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  
  if (any(is.na(names(de_raw)))) {
    idx_na <- which(is.na(names(de_raw)))
    names(de_raw)[idx_na] <- paste0("unnamed_col_", seq_along(idx_na))
  }
  
  if ("Unnamed: 0" %in% names(de_raw)) {
    de_raw <- de_raw %>% dplyr::rename(ensembl_id = `Unnamed: 0`)
  } else if ("...1" %in% names(de_raw)) {
    de_raw <- de_raw %>% dplyr::rename(ensembl_id = `...1`)
  } else if ("" %in% names(de_raw)) {
    names(de_raw)[which(names(de_raw) == "")[1]] <- "ensembl_id"
  }
  
  if (!"ensembl_id" %in% names(de_raw) && ncol(de_raw) >= 1L) {
    first_values <- de_raw[[1]]
    first_non_na <- first_values[which(!is.na(first_values))[1]]
    looks_like_ensembl <- !is.na(first_non_na) && grepl("^ENS", as.character(first_non_na))
    if (looks_like_ensembl) {
      names(de_raw)[1] <- "ensembl_id"
    }
  }
  
  if (!"ensembl_id" %in% names(de_raw)) {
    stop("Could not find Ensembl ID column. Expected 'Unnamed: 0' or 'ensembl_id'.", call. = FALSE)
  }
  
  if ("log2FC_shrunken" %in% names(de_raw)) {
    if ("log2FC" %in% names(de_raw)) {
      fc_column_used <- "log2FC_shrunken (fallback to log2FC when NA)"
      de_raw <- de_raw %>%
        mutate(log2FC_used = ifelse(!is.na(log2FC_shrunken), log2FC_shrunken, log2FC))
    } else {
      fc_column_used <- "log2FC_shrunken"
      de_raw <- de_raw %>%
        mutate(log2FC_used = log2FC_shrunken)
    }
  } else if ("log2FC" %in% names(de_raw)) {
    fc_column_used <- "log2FC"
    de_raw <- de_raw %>%
      mutate(log2FC_used = log2FC)
  } else {
    fc_column_used <- "none (all NA)"
    warning("Neither log2FC_shrunken nor log2FC found. log2FC_used will be NA for all genes.")
    de_raw <- de_raw %>%
      mutate(log2FC_used = NA_real_)
  }
  
  if (!"padj" %in% names(de_raw)) {
    warning("padj column not found. All genes will be treated as not significant.")
    de_raw <- de_raw %>% mutate(padj = NA_real_)
  }
  
  de_raw <- de_raw %>%
    mutate(
      ensembl_id = as.character(ensembl_id),
      log2FC_used = suppressWarnings(as.numeric(log2FC_used)),
      padj = suppressWarnings(as.numeric(padj))
    )
  
  ensembl_keys <- unique(stats::na.omit(de_raw$ensembl_id))
  annotation <- AnnotationDbi::select(
    x = org.Gg.eg.db,
    keys = ensembl_keys,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "ENSEMBL"
  )
  
  annotation_first <- annotation %>%
    as_tibble() %>%
    dplyr::rename(ensembl_id = ENSEMBL, symbol = SYMBOL, entrez_id = ENTREZID) %>%
    group_by(ensembl_id) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  de_annot <- de_raw %>%
    left_join(annotation_first, by = "ensembl_id") %>%
    mutate(symbol = as.character(symbol))
  
  de_by_symbol <- de_annot %>%
    filter(!is.na(symbol), symbol != "") %>%
    arrange(is.na(padj), padj, desc(abs(replace_na(log2FC_used, 0)))) %>%
    group_by(symbol) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(data_symbol = symbol, ensembl_id, log2FC_used, padj)
  
  pathway <- get_pathway_spec()
  pathway_genes <- pathway$pathway_genes
  edge_df <- pathway$edge_df
  layout_df <- pathway$layout_df
  
  if (!all(pathway_genes %in% layout_df$symbol)) {
    missing_layout <- setdiff(pathway_genes, layout_df$symbol)
    stop("Manual layout missing symbols: ", paste(missing_layout, collapse = ", "), call. = FALSE)
  }
  
  pathway_symbol_map <- tibble(symbol = pathway_genes) %>%
    mutate(map_symbol = ifelse(symbol == "STING1", "TMEM173", symbol))
  
  node_status_levels <- c(
    "Missing gene",
    "Present; padj NA",
    "Present; not significant",
    "Present; significant"
  )
  
  node_tbl <- pathway_symbol_map %>%
    left_join(layout_df, by = "symbol") %>%
    left_join(de_by_symbol, by = c("map_symbol" = "data_symbol")) %>%
    mutate(
      present_in_data = !is.na(ensembl_id),
      signif_flag = !is.na(padj) & padj < padj_threshold,
      node_status = case_when(
        !present_in_data ~ "Missing gene",
        is.na(padj) ~ "Present; padj NA",
        signif_flag ~ "Present; significant",
        TRUE ~ "Present; not significant"
      ),
      node_status = factor(node_status, levels = node_status_levels),
      log2FC_plot = ifelse(node_status == "Missing gene", NA_real_, log2FC_used),
      gene_id = ensembl_id
    ) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(
      gene_id, symbol, module, x, y, log2FC_used, log2FC_plot, padj,
      signif_flag, node_status, present_in_data
    )
  
  missing_genes <- node_tbl %>%
    filter(node_status == "Missing gene") %>%
    dplyr::select(symbol, module)
  
  output_nodes <- c("ISG15", "MX1", "OASL", "IFIT5", "RSAD2", "CXCL10", "CCL5", "IL6", "TNF")
  edge_df <- edge_df %>%
    mutate(
      edge_class = factor(ifelse(to %in% output_nodes, "output", "core"), levels = c("core", "output"))
    )
  
  g <- graph_from_data_frame(
    d = edge_df,
    directed = TRUE,
    vertices = node_tbl %>%
      dplyr::rename(name = symbol) %>%
      dplyr::select(name, dplyr::everything())
  )
  
  lim <- max(abs(node_tbl$log2FC_used), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  
  match_n <- sum(node_tbl$present_in_data, na.rm = TRUE)
  total_n <- nrow(node_tbl)
  subtitle_txt <- paste0(basename(input), " | matched: ", match_n, "/", total_n, " pathway genes")
  
  status_colours <- c(
    "Missing gene" = "grey45",
    "Present; padj NA" = "black",
    "Present; not significant" = "black",
    "Present; significant" = "black"
  )
  status_linetypes <- c(
    "Missing gene" = "22",
    "Present; padj NA" = "13",
    "Present; not significant" = "solid",
    "Present; significant" = "solid"
  )
  status_linewidths <- c(
    "Missing gene" = 0.5,
    "Present; padj NA" = 0.75,
    "Present; not significant" = 0.45,
    "Present; significant" = 1.2
  )
  
  p <- ggraph(g, layout = "manual", x = x, y = y) +
    geom_edge_link(
      aes(
        start_cap = label_rect(
          node1.name,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        end_cap = label_rect(
          node2.name,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        edge_width = edge_class
      ),
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      edge_colour = "black",
      lineend = "round",
      show.legend = FALSE
    ) +
    scale_edge_width_manual(values = c("core" = 0.72, "output" = 0.56), guide = "none") +
    geom_node_label(
      aes(
        label = name,
        fill = log2FC_plot,
        colour = node_status,
        linetype = node_status,
        linewidth = node_status
      ),
      size = 3.0,
      family = "sans",
      fontface = "bold",
      text.colour = "black",
      label.padding = unit(0.22, "lines"),
      label.r = unit(0.08, "lines"),
      lineheight = 0.95
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = squish,
      na.value = "grey90",
      name = "log2FC"
    ) +
    scale_colour_manual(
      values = status_colours,
      breaks = node_status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    scale_linetype_manual(
      values = status_linetypes,
      breaks = node_status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    scale_linewidth_manual(
      values = status_linewidths,
      breaks = node_status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      colour = guide_legend(
        order = 2,
        override.aes = list(
          fill = c("grey90", "grey85", "grey85", "grey85"),
          linetype = unname(status_linetypes[node_status_levels]),
          linewidth = unname(status_linewidths[node_status_levels]),
          colour = unname(status_colours[node_status_levels])
        )
      ),
      linetype = "none",
      linewidth = "none"
    ) +
    labs(
      title = title,
      subtitle = subtitle_txt,
      caption = "Missing gene = gray dashed outline; padj NA = dotted outline; significant = thick outline."
    ) +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.caption = element_text(size = 9, colour = "grey30")
    )
  
  svg_file <- file.path(outdir, "cgast_sting_map_SVG.svg")
  png_file <- file.path(outdir, "cgast_sting_map_PNG.png")
  nodes_file <- file.path(outdir, "cgast_sting_nodes.tsv")
  missing_file <- file.path(outdir, "cgast_sting_missing_genes.tsv")
  all_results_file <- file.path(outdir, "cgast_sting_all_results.csv")
  metadata_file <- file.path(outdir, "cgast_sting_run_metadata.yml")
  
  svglite::svglite(file = svg_file, width = 14, height = 7.5, bg = "white")
  print(p)
  invisible(dev.off())
  
  ggsave(
    filename = png_file,
    plot = p,
    width = 14,
    height = 7.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  node_tbl %>%
    dplyr::select(gene_id, symbol, log2FC_used, padj, present_in_data, node_status, signif_flag) %>%
    readr::write_tsv(nodes_file, na = "NA")
  
  missing_genes %>%
    readr::write_tsv(missing_file, na = "NA")
  
  de_annot %>%
    dplyr::select(ensembl_id, symbol, entrez_id, log2FC_used, padj, dplyr::everything()) %>%
    readr::write_csv(all_results_file, na = "NA")
  
  metadata_lines <- c(
    paste0("input_file: \"", gsub("\"", "\\\\\"", normalizePath(input, winslash = "/", mustWork = FALSE)), "\""),
    paste0("input_basename: \"", gsub("\"", "\\\\\"", basename(input)), "\""),
    paste0("timestamp: \"", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "\""),
    paste0("fc_column_used: \"", gsub("\"", "\\\\\"", fc_column_used), "\""),
    paste0("padj_threshold: ", format(padj_threshold, scientific = FALSE)),
    paste0("matched_count: ", match_n),
    paste0("total_pathway_genes: ", total_n)
  )
  writeLines(metadata_lines, metadata_file, useBytes = TRUE)
  
  message("Done.")
  message("Outputs:")
  message(" - ", svg_file)
  message(" - ", png_file)
  message(" - ", nodes_file)
  message(" - ", missing_file)
  message(" - ", all_results_file)
  message(" - ", metadata_file)
  
  invisible(
    list(
      svg = svg_file,
      png = png_file,
      nodes = nodes_file,
      missing = missing_file,
      all_results = all_results_file,
      metadata = metadata_file
    )
  )
}

prepare_de_table_for_pathways <- function(input,
                                          lfc_column_preference = c("log2FC_shrunken", "log2FC"),
                                          id_column_candidates = c("ensembl_id", "Unnamed: 0", "gene", "id")) {
  if (!file.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }
  
  de_raw <- readr::read_csv(
    file = input,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  
  if (any(is.na(names(de_raw)))) {
    idx_na <- which(is.na(names(de_raw)))
    names(de_raw)[idx_na] <- paste0("unnamed_col_", seq_along(idx_na))
  }
  
  id_candidates <- unique(c(id_column_candidates, "ensembl_id", "Unnamed: 0", "...1"))
  id_col <- id_candidates[id_candidates %in% names(de_raw)][1]
  
  if (!is.na(id_col) && nzchar(id_col) && id_col != "ensembl_id") {
    names(de_raw)[names(de_raw) == id_col] <- "ensembl_id"
  }
  
  if (!"ensembl_id" %in% names(de_raw) && "" %in% names(de_raw)) {
    names(de_raw)[which(names(de_raw) == "")[1]] <- "ensembl_id"
  }
  
  if (!"ensembl_id" %in% names(de_raw) && ncol(de_raw) >= 1L) {
    first_values <- de_raw[[1]]
    first_non_na <- first_values[which(!is.na(first_values))[1]]
    looks_like_ensembl <- !is.na(first_non_na) && grepl("^ENS", as.character(first_non_na))
    if (looks_like_ensembl) {
      names(de_raw)[1] <- "ensembl_id"
    }
  }
  
  if (!"ensembl_id" %in% names(de_raw)) {
    stop(
      "Could not find Ensembl ID column. Checked: ",
      paste(id_column_candidates, collapse = ", "),
      call. = FALSE
    )
  }
  
  available_fc <- lfc_column_preference[lfc_column_preference %in% names(de_raw)]
  if (length(available_fc) >= 1L) {
    primary_fc <- available_fc[1]
    secondary_fc <- if (length(available_fc) >= 2L) available_fc[2] else NA_character_
    if (!is.na(secondary_fc)) {
      fc_column_used <- paste0(primary_fc, " (fallback to ", secondary_fc, " when NA)")
      de_raw <- de_raw %>%
        mutate(log2FC_used = ifelse(!is.na(.data[[primary_fc]]), .data[[primary_fc]], .data[[secondary_fc]]))
    } else {
      fc_column_used <- primary_fc
      de_raw <- de_raw %>% mutate(log2FC_used = .data[[primary_fc]])
    }
  } else {
    fc_column_used <- "none (all NA)"
    warning("None of the preferred log2FC columns were found: ", paste(lfc_column_preference, collapse = ", "))
    de_raw <- de_raw %>% mutate(log2FC_used = NA_real_)
  }
  
  if (!"padj" %in% names(de_raw)) {
    warning("padj column not found. All genes will be treated as not significant.")
    de_raw <- de_raw %>% mutate(padj = NA_real_)
  }
  
  de_raw <- de_raw %>%
    mutate(
      ensembl_id = as.character(ensembl_id),
      log2FC_used = suppressWarnings(as.numeric(log2FC_used)),
      padj = suppressWarnings(as.numeric(padj))
    )
  
  ensembl_keys <- unique(stats::na.omit(de_raw$ensembl_id))
  annotation <- AnnotationDbi::select(
    x = org.Gg.eg.db,
    keys = ensembl_keys,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "ENSEMBL"
  )
  
  annotation_first <- annotation %>%
    as_tibble() %>%
    dplyr::rename(ensembl_id = ENSEMBL, symbol = SYMBOL, entrez_id = ENTREZID) %>%
    group_by(ensembl_id) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  de_annot <- de_raw %>%
    left_join(annotation_first, by = "ensembl_id") %>%
    mutate(symbol = as.character(symbol))
  
  de_by_symbol <- de_annot %>%
    filter(!is.na(symbol), symbol != "") %>%
    arrange(is.na(padj), padj, desc(abs(replace_na(log2FC_used, 0)))) %>%
    group_by(symbol) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(data_symbol = symbol, ensembl_id, log2FC_used, padj)
  
  list(
    de_raw = de_raw,
    de_annot = de_annot,
    de_by_symbol = de_by_symbol,
    fc_column_used = fc_column_used,
    input_basename = basename(input)
  )
}

get_pathway_config <- function(pathway) {
  pathway <- tolower(pathway)
  
  if (pathway == "apoptosis") {
    nodes <- tribble(
      ~node_id, ~symbol_key, ~label, ~x, ~y, ~module, ~compartment,
      "TNF", "TNF", "TNF", 0.0, 0.2, "Ligand/Receptor", NA_character_,
      "TNFRSF1A", "TNFRSF1A", "TNFRSF1A", 1.0, 0.2, "Ligand/Receptor", "Plasma membrane",
      "TRADD", "TRADD", "TRADD", 2.0, 0.2, "Adaptor", "Cytosol",
      "FADD", "FADD", "FADD", 3.0, 0.8, "Caspase axis", "Cytosol",
      "CASP8", "CASP8", "CASP8", 4.0, 0.8, "Caspase axis", "Cytosol",
      "CASP7", "CASP7", "CASP7", 5.0, 1.2, "Executioners", "Cytosol",
      "CASP3", "CASP3", "CASP3", 6.2, 0.8, "Executioners", "Cytosol",
      "BID", "BID", "BID", 4.8, 0.0, "Mitochondrial axis", "Mitochondria",
      "BAX", "BAX", "BAX", 5.8, -0.3, "Mitochondrial axis", "Mitochondria",
      "CYCS", "CYCS", "CYCS", 6.8, -0.3, "Mitochondrial axis", "Mitochondria",
      "APAF1", "APAF1", "APAF1", 7.8, -0.3, "Mitochondrial axis", "Cytosol",
      "CASP9", "CASP9", "CASP9", 8.8, -0.3, "Mitochondrial axis", "Cytosol"
    )
    
    edges <- tribble(
      ~from, ~to, ~edge_class, ~edge_pathway,
      "TNF", "TNFRSF1A", "core", "shared_death_input",
      "TNFRSF1A", "TRADD", "core", "shared_death_input",
      "TRADD", "FADD", "core", "apoptosis_branch",
      "FADD", "CASP8", "core", "apoptosis_branch",
      "CASP8", "CASP3", "output", "apoptosis_branch",
      "CASP8", "CASP7", "output", "apoptosis_branch",
      "CASP8", "BID", "core", "apoptosis_branch",
      "BID", "BAX", "core", "apoptosis_branch",
      "BAX", "CYCS", "core", "apoptosis_branch",
      "CYCS", "APAF1", "core", "apoptosis_branch",
      "APAF1", "CASP9", "core", "apoptosis_branch",
      "CASP9", "CASP3", "output", "apoptosis_branch"
    )
    
    synonyms <- list(
      TNFRSF1A = c("TNFR1", "TNFRSF1"),
      CASP3 = c("CASP3A", "CASP3-like"),
      BAX = c("BAX1")
    )
    
    return(list(name = "apoptosis", display_name = "Apoptosis", nodes = nodes, edges = edges, synonyms = synonyms))
  }
  
  if (pathway == "necroptosis") {
    nodes <- tribble(
      ~node_id, ~symbol_key, ~label, ~x, ~y, ~module, ~compartment,
      "TNF", "TNF", "TNF", 0.0, 0.2, "Ligand/Receptor", NA_character_,
      "TNFRSF1A", "TNFRSF1A", "TNFRSF1A", 1.0, 0.2, "Ligand/Receptor", "Plasma membrane",
      "TRADD", "TRADD", "TRADD", 2.0, 0.2, "Adaptor", "Cytosol",
      "RIPK1", "RIPK1", "RIPK1", 3.2, 0.2, "Core necroptosis", "Cytosol",
      "RIPK3", "RIPK3", "RIPK3", 4.2, 0.2, "Core necroptosis", "Cytosol",
      "MLKL", "MLKL", "MLKL", 5.2, 0.2, "Executioner", "Plasma membrane",
      "ZBP1", "ZBP1", "ZBP1", 3.2, 1.0, "Inputs", "Cytosol",
      "CASP8", "CASP8", "CASP8", 3.2, -0.6, "Cross-talk", "Cytosol"
    )
    
    edges <- tribble(
      ~from, ~to, ~edge_class, ~edge_pathway,
      "TNF", "TNFRSF1A", "core", "shared_death_input",
      "TNFRSF1A", "TRADD", "core", "shared_death_input",
      "TRADD", "RIPK1", "core", "necroptosis_branch",
      "RIPK1", "RIPK3", "core", "necroptosis_branch",
      "RIPK3", "MLKL", "output", "necroptosis_branch",
      "ZBP1", "RIPK3", "core", "necroptosis_branch",
      "CASP8", "RIPK1", "core", "necroptosis_branch"
    )
    
    synonyms <- list(
      TNFRSF1A = c("TNFR1", "TNFRSF1"),
      RIPK1 = c("RIP1"),
      RIPK3 = c("RIP3"),
      MLKL = c("MLKL1"),
      ZBP1 = c("DAI")
    )
    
    return(list(name = "necroptosis", display_name = "Necroptosis", nodes = nodes, edges = edges, synonyms = synonyms))
  }
  
  if (pathway == "combined_death") {
    core <- get_pathway_spec()
    core_nodes <- core$layout_df %>%
      transmute(
        node_id = symbol,
        symbol_key = symbol,
        label = symbol,
        x = x,
        y = y,
        module = module,
        compartment = case_when(
          symbol == "STING1" ~ "ER/Golgi",
          symbol %in% c("IRF3", "IRF7", "NFKB1", "RELA", "STAT1", "STAT2", "IRF9") ~ "Nucleus",
          symbol %in% c("IFNAR1", "IFNAR2") ~ "Plasma membrane",
          symbol %in% c("ISG15", "MX1", "OASL", "IFIT5", "RSAD2") ~ "Cytosol",
          TRUE ~ "Cytosol"
        )
      )
    
    core_output_nodes <- c("ISG15", "MX1", "OASL", "IFIT5", "RSAD2", "CXCL10", "CCL5", "IL6", "TNF")
    core_edges <- core$edge_df %>%
      mutate(
        edge_class = ifelse(to %in% core_output_nodes, "output", "core"),
        edge_pathway = "cgas_sting"
      )
    
    death_nodes <- tribble(
      ~node_id, ~symbol_key, ~label, ~x, ~y, ~module, ~compartment,
      "TNFRSF1A", "TNFRSF1A", "TNFRSF1A", 8.2, -0.2, "Death modules", "Plasma membrane",
      "TRADD", "TRADD", "TRADD", 10.2, -0.2, "Death modules", "Cytosol",
      "FADD", "FADD", "FADD", 11.2, 0.4, "Apoptosis", "Cytosol",
      "CASP8", "CASP8", "CASP8", 12.2, 0.4, "Apoptosis", "Cytosol",
      "CASP7", "CASP7", "CASP7", 13.2, 0.9, "Apoptosis", "Cytosol",
      "CASP3", "CASP3", "CASP3", 14.2, 0.4, "Apoptosis", "Cytosol",
      "BID", "BID", "BID", 12.2, -0.4, "Apoptosis", "Mitochondria",
      "BAX", "BAX", "BAX", 13.2, -0.8, "Apoptosis", "Mitochondria",
      "CYCS", "CYCS", "CYCS", 14.2, -0.8, "Apoptosis", "Mitochondria",
      "APAF1", "APAF1", "APAF1", 15.2, -0.8, "Apoptosis", "Cytosol",
      "CASP9", "CASP9", "CASP9", 16.2, -0.8, "Apoptosis", "Cytosol",
      "RIPK1", "RIPK1", "RIPK1", 11.2, -2.2, "Necroptosis", "Cytosol",
      "RIPK3", "RIPK3", "RIPK3", 12.2, -2.2, "Necroptosis", "Cytosol",
      "MLKL", "MLKL", "MLKL", 14.2, -2.2, "Necroptosis", "Plasma membrane"
    )
    
    death_edges <- tribble(
      ~from, ~to, ~edge_class, ~edge_pathway,
      "TNF", "TNFRSF1A", "core", "shared_death_input",
      "TNFRSF1A", "TRADD", "core", "shared_death_input",
      "TRADD", "FADD", "core", "apoptosis_branch",
      "FADD", "CASP8", "core", "apoptosis_branch",
      "CASP8", "CASP7", "output", "apoptosis_branch",
      "CASP8", "CASP3", "output", "apoptosis_branch",
      "CASP8", "BID", "core", "apoptosis_branch",
      "BID", "BAX", "core", "apoptosis_branch",
      "BAX", "CYCS", "core", "apoptosis_branch",
      "CYCS", "APAF1", "core", "apoptosis_branch",
      "APAF1", "CASP9", "core", "apoptosis_branch",
      "CASP9", "CASP3", "output", "apoptosis_branch",
      "TRADD", "RIPK1", "core", "necroptosis_branch",
      "TNFRSF1A", "RIPK1", "core", "necroptosis_branch",
      "RIPK1", "RIPK3", "core", "necroptosis_branch",
      "RIPK3", "MLKL", "output", "necroptosis_branch",
      "CASP8", "RIPK1", "core", "necroptosis_branch"
    )
    
    nodes <- bind_rows(core_nodes, death_nodes) %>%
      distinct(node_id, .keep_all = TRUE)
    edges <- bind_rows(core_edges, death_edges)
    
    synonyms <- list(
      STING1 = c("TMEM173"),
      TNFRSF1A = c("TNFR1", "TNFRSF1"),
      CASP3 = c("CASP3A", "CASP3-like"),
      BAX = c("BAX1"),
      RIPK1 = c("RIP1"),
      RIPK3 = c("RIP3"),
      MLKL = c("MLKL1"),
      ZBP1 = c("DAI")
    )
    
    compartment_boxes <- tibble::tribble(
      ~compartment,             ~xmin,  ~xmax,  ~ymin,  ~ymax,  ~node_ids,                                                                    ~label_nudge_x, ~label_nudge_y, ~label_pad,
      "Plasma membrane",         7.45,   9.00,  -0.75,   0.35,  list(c("TNFRSF1A")),                                                           0,              0,              0.30,
      "Cytosol",                 10.05,  13.45,  -2.65,   1.35,  list(c("TRADD", "FADD", "CASP8", "CASP7", "RIPK1", "RIPK3")),                 0,              0,              0.30,
      "Nucleus",                  2.15,   6.95,   1.05,   3.85,  list(c("IRF7", "IRF3", "STAT1", "STAT2", "IRF9", "NFKB1", "RELA")),           0,              0,              0.30,
      "Mitochondria",            13.05,  14.85,  -1.30,  -0.30,  list(c("BAX", "CYCS")),                                                        0,              0,              0.30,
      "ER/Golgi",                 0.15,   1.85,   1.05,   2.55,  list(c("STING1")),                                                             0,              0.30,           0.30,
      "Plasma membrane (MLKL)", 13.55,  14.85,  -2.65,  -1.80,  list(c("MLKL")),                                                               0,              0,              0.25
    )
    
    return(list(
      name = "combined_death",
      display_name = "Combined Death Modules",
      nodes = nodes,
      edges = edges,
      synonyms = synonyms,
      compartment_boxes = compartment_boxes
    ))
  }
  
  stop("Unknown pathway config: ", pathway, call. = FALSE)
}

resolve_symbol_match <- function(symbol_key, synonyms, available_symbols) {
  candidates <- unique(c(symbol_key, synonyms[[symbol_key]]))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0L || length(available_symbols) == 0L) {
    return(NA_character_)
  }
  
  available_symbols <- as.character(available_symbols)
  match_idx <- match(toupper(candidates), toupper(available_symbols), nomatch = 0L)
  match_idx <- match_idx[match_idx > 0L]
  if (length(match_idx) == 0L) {
    return(NA_character_)
  }
  available_symbols[match_idx[1]]
}

build_pathway_node_table <- function(config, de_by_symbol, p_cutoff = 0.05) {
  required_cols <- c("node_id", "symbol_key", "label", "x", "y")
  if (!all(required_cols %in% names(config$nodes))) {
    stop("Config nodes must include: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  
  nodes <- config$nodes
  if (!"module" %in% names(nodes)) {
    nodes$module <- "Pathway"
  }
  if (!"compartment" %in% names(nodes)) {
    nodes$compartment <- NA_character_
  }
  
  available_symbols <- de_by_symbol$data_symbol
  synonyms <- config$synonyms
  if (is.null(synonyms)) {
    synonyms <- list()
  }
  
  node_status_levels <- c(
    "Missing gene",
    "Present; padj NA",
    "Present; not significant",
    "Present; significant"
  )
  
  nodes <- nodes %>%
    rowwise() %>%
    mutate(map_symbol = resolve_symbol_match(symbol_key, synonyms, available_symbols)) %>%
    ungroup()
  
  nodes %>%
    left_join(de_by_symbol, by = c("map_symbol" = "data_symbol")) %>%
    mutate(
      present_in_data = !is.na(ensembl_id),
      signif_flag = !is.na(padj) & padj < p_cutoff,
      node_status = case_when(
        !present_in_data ~ "Missing gene",
        is.na(padj) ~ "Present; padj NA",
        signif_flag ~ "Present; significant",
        TRUE ~ "Present; not significant"
      ),
      node_status = factor(node_status, levels = node_status_levels),
      log2FC_plot = ifelse(node_status == "Missing gene", NA_real_, log2FC_used),
      compartment = as.character(compartment),
      gene_id = ensembl_id
    ) %>%
    distinct(node_id, .keep_all = TRUE)
}

choose_compartment_label <- function(box, nodes_in_box, pad = 0.45) {
  # Place label centered along the top edge of the box.
  # vjust=1 means the top of the text aligns with ymax, so it hangs into the box.
  tibble::tibble(
    label_x = (box$xmin + box$xmax) / 2,
    label_y = box$ymax,
    hjust   = 0.5,
    vjust   = 1,
    pos_id  = "top_center"
  )
}

infer_compartment_boxes <- function(nodes_df,
                                    pad_x = 0.6,
                                    pad_y = 0.6,
                                    min_nodes_per_box = 2,
                                    force_singleton_compartments = NULL) {
  if (!"compartment" %in% names(nodes_df)) {
    return(tibble())
  }
  
  nodes_with_comp <- nodes_df %>%
    filter(!is.na(compartment), compartment != "")
  if (nrow(nodes_with_comp) == 0) {
    return(tibble())
  }
  
  if (is.null(force_singleton_compartments)) {
    force_singleton_compartments <- character(0)
  }
  
  keep_compartments <- nodes_with_comp %>%
    count(compartment, name = "n_nodes") %>%
    filter(n_nodes >= min_nodes_per_box | compartment %in% force_singleton_compartments) %>%
    pull(compartment)
  
  if (length(keep_compartments) == 0) {
    return(tibble())
  }
  
  boxes <- nodes_with_comp %>%
    filter(compartment %in% keep_compartments) %>%
    group_by(compartment) %>%
    summarise(
      xmin = min(x, na.rm = TRUE) - pad_x,
      xmax = max(x, na.rm = TRUE) + pad_x,
      ymin = min(y, na.rm = TRUE) - pad_y,
      ymax = max(y, na.rm = TRUE) + pad_y,
      .groups = "drop"
    )
  
  label_tbl <- purrr::map_dfr(seq_len(nrow(boxes)), function(i) {
    box <- boxes[i, , drop = FALSE]
    comp <- box$compartment[[1]]
    nodes_in_box <- nodes_with_comp %>%
      filter(compartment == comp) %>%
      dplyr::select(x, y)
    choose_compartment_label(box = box, nodes_in_box = nodes_in_box, pad = 0.45)
  }) %>%
    dplyr::select(label_x, label_y, hjust, vjust, pos_id)
  
  bind_cols(boxes, label_tbl)
}

compute_manual_compartment_boxes <- function(nodes_df, compartment_spec) {
  if (is.null(compartment_spec) || nrow(compartment_spec) == 0) {
    return(tibble())
  }
  
  if (!"compartment" %in% names(compartment_spec)) {
    stop("Manual compartment spec must include a `compartment` column.", call. = FALSE)
  }
  
  spec <- as_tibble(compartment_spec)
  has_explicit_bounds <- all(c("xmin", "xmax", "ymin", "ymax") %in% names(spec))
  has_node_pad_bounds <- all(c("node_ids", "pad_x", "pad_y") %in% names(spec))
  if (!has_explicit_bounds && !has_node_pad_bounds) {
    stop(
      "Manual compartment spec must include either explicit bounds (xmin/xmax/ymin/ymax) ",
      "or node_ids + pad_x/pad_y.",
      call. = FALSE
    )
  }
  
  if (!"label_nudge_x" %in% names(spec)) {
    spec$label_nudge_x <- 0
  }
  if (!"label_nudge_y" %in% names(spec)) {
    spec$label_nudge_y <- 0
  }
  if (!"label_pad" %in% names(spec)) {
    spec$label_pad <- 0.45
  }
  if (!"node_ids" %in% names(spec)) {
    spec$node_ids <- rep(list(NULL), nrow(spec))
  }
  
  boxes <- if (has_explicit_bounds) {
    spec %>%
      dplyr::transmute(
        compartment = as.character(compartment),
        xmin = as.numeric(xmin),
        xmax = as.numeric(xmax),
        ymin = as.numeric(ymin),
        ymax = as.numeric(ymax),
        node_ids = node_ids,
        label_nudge_x = as.numeric(label_nudge_x),
        label_nudge_y = as.numeric(label_nudge_y),
        label_pad = as.numeric(label_pad)
      )
  } else {
    purrr::map_dfr(seq_len(nrow(spec)), function(i) {
      comp <- spec$compartment[[i]]
      node_ids <- spec$node_ids[[i]]
      if (is.null(node_ids)) {
        return(tibble())
      }
      node_ids <- as.character(unlist(node_ids, recursive = TRUE, use.names = FALSE))
      comp_nodes <- nodes_df %>% filter(node_id %in% node_ids)
      if (nrow(comp_nodes) == 0) {
        warning("Manual compartment '", comp, "' has no matching nodes in this map.")
        return(tibble())
      }
      tibble(
        compartment = comp,
        xmin = min(comp_nodes$x, na.rm = TRUE) - spec$pad_x[[i]],
        xmax = max(comp_nodes$x, na.rm = TRUE) + spec$pad_x[[i]],
        ymin = min(comp_nodes$y, na.rm = TRUE) - spec$pad_y[[i]],
        ymax = max(comp_nodes$y, na.rm = TRUE) + spec$pad_y[[i]],
        node_ids = list(node_ids),
        label_nudge_x = spec$label_nudge_x[[i]],
        label_nudge_y = spec$label_nudge_y[[i]],
        label_pad = spec$label_pad[[i]]
      )
    })
  }
  
  label_tbl <- purrr::map_dfr(seq_len(nrow(boxes)), function(i) {
    box <- boxes[i, , drop = FALSE]
    node_ids <- box$node_ids[[1]]
    if (is.null(node_ids) || length(node_ids) == 0) {
      nodes_in_box <- nodes_df %>%
        filter(
          !is.na(compartment),
          compartment == box$compartment[[1]],
          x >= box$xmin[[1]], x <= box$xmax[[1]],
          y >= box$ymin[[1]], y <= box$ymax[[1]]
        ) %>%
        dplyr::select(x, y)
    } else {
      node_ids <- as.character(unlist(node_ids, recursive = TRUE, use.names = FALSE))
      nodes_in_box <- nodes_df %>%
        filter(node_id %in% node_ids) %>%
        dplyr::select(x, y)
    }
    
    lab <- choose_compartment_label(
      box = box,
      nodes_in_box = nodes_in_box,
      pad = box$label_pad[[1]]
    )
    lab$label_x <- lab$label_x + box$label_nudge_x[[1]]
    lab$label_y <- lab$label_y + box$label_nudge_y[[1]]
    lab
  }) %>%
    dplyr::select(label_x, label_y, hjust, vjust, pos_id)
  
  bind_cols(
    boxes %>% dplyr::select(compartment, xmin, xmax, ymin, ymax),
    label_tbl
  )
}

build_pathway_plot <- function(node_tbl,
                               edge_df,
                               title,
                               subtitle,
                               show_compartments = FALSE,
                               compartment_pad_x = 0.6,
                               compartment_pad_y = 0.6,
                               min_nodes_per_compartment = 2,
                               force_singleton_compartments = NULL,
                               compartment_boxes_spec = NULL) {
  if (!"edge_class" %in% names(edge_df)) {
    edge_df <- edge_df %>% mutate(edge_class = "core")
  }
  if (!"edge_pathway" %in% names(edge_df)) {
    edge_df <- edge_df %>% mutate(edge_pathway = "Pathway flow")
  }
  edge_df <- edge_df %>%
    mutate(
      edge_class = as.character(edge_class),
      edge_pathway = as.character(edge_pathway),
      curvature = dplyr::case_when(
        TRUE ~ 0
      ),
      is_curved = curvature != 0
    )
  
  edges_joined <- edge_df %>%
    left_join(node_tbl %>% dplyr::select(node_id, x, y), by = c("from" = "node_id")) %>%
    dplyr::rename(x_from = x, y_from = y) %>%
    left_join(node_tbl %>% dplyr::select(node_id, x, y), by = c("to" = "node_id")) %>%
    dplyr::rename(x_to = x, y_to = y)
  
  missing_xy <- edges_joined %>%
    dplyr::filter(is.na(x_from) | is.na(x_to) | is.na(y_from) | is.na(y_to))
  if (nrow(missing_xy) > 0) {
    warning(
      "Edges dropped due to missing node coords: ",
      paste0(missing_xy$from, "->", missing_xy$to, collapse = ", ")
    )
  }
  if (all(c("TNFRSF1A", "TRADD") %in% node_tbl$node_id)) {
    stopifnot(any(edges_joined$from == "TNFRSF1A" & edges_joined$to == "TRADD"))
  }
  
  edges_joined_valid <- edges_joined %>%
    dplyr::filter(!is.na(x_from), !is.na(y_from), !is.na(x_to), !is.na(y_to))
  
  edge_overlay <- edges_joined_valid %>%
    dplyr::filter(is_curved)
  edges_main_joined <- edges_joined_valid %>%
    dplyr::filter(!is_curved)
  
  edge_df_main <- edges_main_joined %>%
    dplyr::select(from, to, edge_class, edge_pathway, curvature, is_curved)
  if (nrow(edge_df_main) == 0) {
    stop("No drawable edges after coordinate join.", call. = FALSE)
  }
  
  g <- graph_from_data_frame(
    d = edge_df_main,
    directed = TRUE,
    vertices = node_tbl %>%
      dplyr::rename(name = node_id) %>%
      dplyr::select(name, dplyr::everything())
  )
  
  lim <- max(abs(node_tbl$log2FC_used), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) {
    lim <- 1
  }
  
  status_levels <- c(
    "Missing gene",
    "Present; padj NA",
    "Present; not significant",
    "Present; significant"
  )
  status_colours <- c(
    "Missing gene" = "grey45",
    "Present; padj NA" = "black",
    "Present; not significant" = "black",
    "Present; significant" = "black"
  )
  status_linetypes <- c(
    "Missing gene" = "22",
    "Present; padj NA" = "13",
    "Present; not significant" = "solid",
    "Present; significant" = "solid"
  )
  status_linewidths <- c(
    "Missing gene" = 0.5,
    "Present; padj NA" = 0.75,
    "Present; not significant" = 0.45,
    "Present; significant" = 1.2
  )
  
  edge_width_values <- c("core" = 0.72, "output" = 0.56)
  extra_edge_classes <- setdiff(unique(edge_df_main$edge_class), names(edge_width_values))
  if (length(extra_edge_classes) > 0) {
    edge_width_values[extra_edge_classes] <- 0.72
  }
  
  edge_pathway_levels <- unique(edge_df_main$edge_pathway)
  use_edge_pathway_legend <- length(edge_pathway_levels) > 1L
  edge_pathway_palette <- c(
    "cgas_sting" = "black",
    "shared_death_input" = "#1F4E79",
    "apoptosis_branch" = "#A33A3A",
    "necroptosis_branch" = "#2F7F56",
    "pathway_flow" = "black",
    "Pathway flow" = "black"
  )
  missing_edge_levels <- setdiff(edge_pathway_levels, names(edge_pathway_palette))
  if (length(missing_edge_levels) > 0) {
    hue_vals <- scales::hue_pal()(length(missing_edge_levels))
    names(hue_vals) <- missing_edge_levels
    edge_pathway_palette <- c(edge_pathway_palette, hue_vals)
  }
  edge_pathway_palette <- edge_pathway_palette[unique(c(edge_pathway_levels, names(edge_pathway_palette)))]
  edge_pathway_labels <- c(
    cgas_sting = "cGAS-STING",
    shared_death_input = "Shared death input",
    apoptosis_branch = "Apoptosis branch",
    necroptosis_branch = "Necroptosis branch",
    pathway_flow = "Pathway flow",
    `Pathway flow` = "Pathway flow"
  )
  
  if (use_edge_pathway_legend) {
    edge_df_main <- edge_df_main %>%
      mutate(edge_pathway = factor(edge_pathway, levels = edge_pathway_levels))
  }
  
  edge_layer_straight <- if (use_edge_pathway_legend) {
    geom_edge_link(
      aes(
        start_cap = label_rect(
          node1.label,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        end_cap = label_rect(
          node2.label,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        edge_width = edge_class,
        edge_colour = edge_pathway,
        filter = !is_curved
      ),
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      lineend = "round",
      show.legend = TRUE
    )
  } else {
    geom_edge_link(
      aes(
        start_cap = label_rect(
          node1.label,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        end_cap = label_rect(
          node2.label,
          padding = margin(2.5, 3.5, 2.5, 3.5, "mm")
        ),
        edge_width = edge_class,
        filter = !is_curved
      ),
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      edge_colour = "black",
      lineend = "round",
      show.legend = FALSE
    )
  }
  
  edge_layer_curved <- if (use_edge_pathway_legend) {
    geom_edge_arc(
      aes(
        edge_width = edge_class,
        edge_colour = edge_pathway,
        filter = is_curved
      ),
      strength = 0.35,
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      lineend = "round",
      show.legend = FALSE
    )
  } else {
    geom_edge_arc(
      aes(edge_width = edge_class, filter = is_curved),
      strength = 0.35,
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      edge_colour = "black",
      lineend = "round",
      show.legend = FALSE
    )
  }
  
  compartment_boxes <- if (isTRUE(show_compartments)) {
    if (!is.null(compartment_boxes_spec) && nrow(compartment_boxes_spec) > 0) {
      compute_manual_compartment_boxes(
        nodes_df = node_tbl,
        compartment_spec = compartment_boxes_spec
      )
    } else {
      infer_compartment_boxes(
        nodes_df = node_tbl,
        pad_x = compartment_pad_x,
        pad_y = compartment_pad_y,
        min_nodes_per_box = min_nodes_per_compartment,
        force_singleton_compartments = force_singleton_compartments
      )
    }
  } else {
    tibble()
  }
  
  if (nrow(edge_overlay) > 0) {
    edge_overlay <- edge_overlay %>%
      mutate(
        dx = x_to - x_from,
        dy = y_to - y_from,
        dist = sqrt(dx^2 + dy^2),
        shrink_amt = ifelse(dist > 0, pmin(0.25, pmax((dist / 2) - 1e-6, 0)), 0),
        ux = ifelse(dist > 0, dx / dist, 0),
        uy = ifelse(dist > 0, dy / dist, 0),
        x_from_plot = x_from + ux * shrink_amt,
        y_from_plot = y_from + uy * shrink_amt,
        x_to_plot = x_to - ux * shrink_amt,
        y_to_plot = y_to - uy * shrink_amt,
        edge_colour_plot = unname(edge_pathway_palette[as.character(edge_pathway)])
      )
    if (any(is.na(edge_overlay$edge_colour_plot))) {
      edge_overlay$edge_colour_plot[is.na(edge_overlay$edge_colour_plot)] <- "black"
    }
  }
  
  x_candidates <- c(node_tbl$x)
  y_candidates <- c(node_tbl$y)
  if (nrow(compartment_boxes) > 0) {
    x_candidates <- c(x_candidates, compartment_boxes$xmin, compartment_boxes$xmax)
    y_candidates <- c(y_candidates, compartment_boxes$ymin, compartment_boxes$ymax)
  }
  x_lim <- c(min(x_candidates, na.rm = TRUE) - 1.20, max(x_candidates, na.rm = TRUE) + 0.80)
  y_lim <- c(min(y_candidates, na.rm = TRUE) - 0.80, max(y_candidates, na.rm = TRUE) + 1.10)
  
  ggraph(g, layout = "manual", x = x, y = y) +
    {
      if (nrow(compartment_boxes) > 0) {
        geom_rect(
          data = compartment_boxes,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = "grey95",
          alpha = 0.10,
          color = "grey85",
          linewidth = 0.3
        )
      } else {
        NULL
      }
    } +
    {
      if (nrow(compartment_boxes) > 0) {
        geom_label(
          data = compartment_boxes,
          aes(
            x = label_x,
            y = label_y,
            label = compartment,
            hjust = hjust,
            vjust = vjust
          ),
          inherit.aes = FALSE,
          color = "grey30",
          fill = "white",
          size = 3.0,
          fontface = "italic",
          alpha = 0.92,
          label.padding = unit(0.13, "lines"),
          label.r = unit(0.06, "lines"),
          label.size = 0.25
        )
      } else {
        NULL
      }
    } +
    edge_layer_straight +
    edge_layer_curved +
    scale_edge_width_manual(values = edge_width_values, guide = "none") +
    geom_node_label(
      aes(
        label = label,
        fill = log2FC_plot,
        colour = node_status,
        linetype = node_status,
        linewidth = node_status
      ),
      size = 3.0,
      family = "sans",
      fontface = "bold",
      text.colour = "black",
      label.padding = unit(0.22, "lines"),
      label.r = unit(0.08, "lines"),
      lineheight = 0.95
    ) +
    {
      if (nrow(edge_overlay) > 0) {
        # Render each curved edge individually so each gets its own curvature sign/magnitude.
        purrr::map(seq_len(nrow(edge_overlay)), function(i) {
          row <- edge_overlay[i, , drop = FALSE]
          geom_curve(
            data = row,
            aes(
              x = x_from_plot,
              y = y_from_plot,
              xend = x_to_plot,
              yend = y_to_plot
            ),
            inherit.aes = FALSE,
            curvature = row$curvature,
            linewidth = 0.9,
            lineend = "round",
            colour = row$edge_colour_plot,
            arrow = grid::arrow(type = "closed", length = grid::unit(0.20, "cm")),
            show.legend = FALSE
          )
        })
      } else {
        NULL
      }
    } +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = squish,
      na.value = "grey90",
      name = "log2FC"
    ) +
    scale_colour_manual(
      values = status_colours,
      breaks = status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    scale_linetype_manual(
      values = status_linetypes,
      breaks = status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    scale_linewidth_manual(
      values = status_linewidths,
      breaks = status_levels,
      drop = FALSE,
      name = "Node status"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      colour = guide_legend(
        order = 2,
        override.aes = list(
          fill = c("grey90", "grey85", "grey85", "grey85"),
          linetype = unname(status_linetypes[status_levels]),
          linewidth = unname(status_linewidths[status_levels]),
          colour = unname(status_colours[status_levels])
        )
      ),
      linetype = "none",
      linewidth = "none"
    ) +
    {
      if (use_edge_pathway_legend) {
        scale_edge_colour_manual(
          values = edge_pathway_palette,
          breaks = edge_pathway_levels,
          drop = FALSE,
          labels = unname(edge_pathway_labels[edge_pathway_levels]),
          name = "Edge pathway"
        )
      } else {
        scale_edge_colour_manual(values = c("Pathway flow" = "black"), guide = "none")
      }
    } +
    {
      if (use_edge_pathway_legend) {
        guides(edge_colour = guide_legend(order = 3))
      } else {
        guides(edge_colour = "none")
      }
    } +
    labs(
      title = title,
      subtitle = subtitle,
      caption = "Missing gene = gray dashed outline; padj NA = dotted outline; significant = thick outline."
    ) +
    coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "off") +
    theme_void(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.caption = element_text(size = 9, colour = "grey30")
    )
}

write_pathway_outputs <- function(prefix, outdir, plot_obj, node_tbl, de_annot) {
  png_file <- file.path(outdir, paste0(prefix, "_map.png"))
  nodes_file <- file.path(outdir, paste0(prefix, "_nodes.tsv"))
  missing_file <- file.path(outdir, paste0(prefix, "_missing.tsv"))
  all_results_file <- file.path(outdir, paste0(prefix, "_all_results.csv"))
  
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = 14,
    height = 7.5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  node_tbl %>%
    dplyr::transmute(
      gene_id = gene_id,
      symbol = symbol_key,
      label = label,
      log2FC_used = log2FC_used,
      padj = padj,
      present_in_data = present_in_data,
      node_status = node_status,
      signif_flag = signif_flag
    ) %>%
    readr::write_tsv(nodes_file, na = "NA")
  
  node_tbl %>%
    filter(node_status == "Missing gene") %>%
    dplyr::transmute(symbol = symbol_key, label = label, module = module) %>%
    readr::write_tsv(missing_file, na = "NA")
  
  de_annot %>%
    dplyr::select(ensembl_id, symbol, entrez_id, log2FC_used, padj, dplyr::everything()) %>%
    readr::write_csv(all_results_file, na = "NA")
  
  list(
    png = png_file,
    nodes = nodes_file,
    missing = missing_file,
    all_results = all_results_file
  )
}

make_interpretation_panel <- function(nodes_annot,
                                      edges_cfg = NULL,
                                      title = "Interpretation notes",
                                      top_n = 3,
                                      p_cutoff = 0.05) {
  boundary_text <- paste(
    "- DESeq2 overlay shows transcriptional change, not protein activation.",
    "- Edges are canonical pathway wiring, not inferred causality in this dataset.",
    "- Apoptosis/necroptosis assignment needs orthogonal validation in tissue.",
    sep = "\n"
  )
  
  validation_text <- paste(
    "Suggested validation markers:",
    "- Apoptosis: cleaved CASP3, cleaved CASP8 (plus TUNEL if available).",
    "- Necroptosis: pMLKL, RIPK3/MLKL localization, MLKL membrane translocation.",
    "- STING axis: STING1 ER-to-Golgi trafficking, IFNB1/ISGs by RNAscope.",
    "- TNF axis: TNF/TNFRSF1A/TRADD complex-localization assays.",
    sep = "\n"
  )
  
  top_text <- "Top changing genes by module: no mapped genes with finite log2FC."
  req_cols <- c("module", "label", "log2FC_used", "present_in_data", "padj")
  if (all(req_cols %in% names(nodes_annot))) {
    top_tbl <- nodes_annot %>%
      filter(present_in_data, !is.na(log2FC_used), !is.na(module), module != "") %>%
      group_by(module) %>%
      slice_max(order_by = abs(log2FC_used), n = top_n, with_ties = FALSE) %>%
      ungroup()
    
    if (nrow(top_tbl) > 0) {
      module_lines <- top_tbl %>%
        mutate(
          sig_tag = ifelse(!is.na(padj) & padj < p_cutoff, "*", ""),
          gene_txt = paste0(label, "(", format(round(log2FC_used, 2), nsmall = 2), sig_tag, ")")
        ) %>%
        group_by(module) %>%
        summarise(line = paste0(dplyr::first(module), ": ", paste(gene_txt, collapse = ", ")), .groups = "drop") %>%
        pull(line)
      
      top_text <- paste(
        "Top changing genes by module (abs log2FC; * = padj < 0.05):",
        paste0("- ", module_lines, collapse = "\n"),
        sep = "\n"
      )
    }
  }
  
  ggplot() +
    annotate("text", x = 0, y = 1.00, label = title, hjust = 0, vjust = 1, size = 5.0, fontface = "bold") +
    annotate("text", x = 0, y = 0.88, label = boundary_text, hjust = 0, vjust = 1, size = 3.25, lineheight = 1.15) +
    annotate("text", x = 0, y = 0.58, label = validation_text, hjust = 0, vjust = 1, size = 3.20, lineheight = 1.15) +
    annotate("text", x = 0, y = 0.22, label = top_text, hjust = 0, vjust = 1, size = 3.00, lineheight = 1.10) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 8, 12, 12)
    )
}

save_multipanel_figure <- function(p_left,
                                   p_right,
                                   outfile,
                                   width = 16,
                                   height = 8,
                                   dpi = 300,
                                   rel_widths = c(3.2, 1.4)) {
  grDevices::png(
    filename = outfile,
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(rel_widths, "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  print(p_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
  invisible(outfile)
}

render_pathway_map <- function(config,
                               prepared_de,
                               outdir,
                               title,
                               p_cutoff = 0.05,
                               show_compartments = FALSE,
                               compartment_pad_x = 0.6,
                               compartment_pad_y = 0.6,
                               min_nodes_per_compartment = 2,
                               force_singleton_compartments = NULL,
                               return_plot = FALSE) {
  node_tbl <- build_pathway_node_table(config, prepared_de$de_by_symbol, p_cutoff = p_cutoff)
  matched_n <- sum(node_tbl$present_in_data, na.rm = TRUE)
  total_n <- nrow(node_tbl)
  subtitle_txt <- paste0(prepared_de$input_basename, " | matched: ", matched_n, "/", total_n, " pathway genes")
  plot_title <- paste0(title, " - ", config$display_name)
  
  plot_obj <- build_pathway_plot(
    node_tbl = node_tbl,
    edge_df = config$edges,
    title = plot_title,
    subtitle = subtitle_txt,
    show_compartments = show_compartments,
    compartment_pad_x = compartment_pad_x,
    compartment_pad_y = compartment_pad_y,
    min_nodes_per_compartment = min_nodes_per_compartment,
    force_singleton_compartments = force_singleton_compartments,
    compartment_boxes_spec = if ("compartment_boxes" %in% names(config)) config$compartment_boxes else NULL
  )
  
  if (isTRUE(return_plot)) {
    return(
      list(
        plot = plot_obj,
        node_tbl = node_tbl,
        subtitle = subtitle_txt,
        title = plot_title
      )
    )
  }
  
  write_pathway_outputs(
    prefix = config$name,
    outdir = outdir,
    plot_obj = plot_obj,
    node_tbl = node_tbl,
    de_annot = prepared_de$de_annot
  )
}

run_pathway_maps <- function(input,
                             outdir,
                             title,
                             pathways = c("apoptosis", "necroptosis"),
                             p_cutoff = 0.05,
                             lfc_column_preference = c("log2FC_shrunken", "log2FC"),
                             id_column_candidates = c("ensembl_id", "Unnamed: 0", "gene", "id"),
                             include_combined = TRUE) {
  if (missing(input) || !nzchar(input)) {
    stop("`input` is required.", call. = FALSE)
  }
  if (missing(outdir) || !nzchar(outdir)) {
    stop("`outdir` is required.", call. = FALSE)
  }
  if (missing(title) || !nzchar(title)) {
    title <- "Pathway maps"
  }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  prepared_de <- prepare_de_table_for_pathways(
    input = input,
    lfc_column_preference = lfc_column_preference,
    id_column_candidates = id_column_candidates
  )
  
  show_compartments <- TRUE
  min_nodes_per_compartment <- 2
  force_singleton_compartments <- c("Plasma membrane", "ER/Golgi")
  
  requested <- unique(tolower(pathways))
  valid <- c("apoptosis", "necroptosis", "combined_death")
  invalid <- setdiff(requested, valid)
  if (length(invalid) > 0) {
    stop("Unknown pathway names: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  
  out <- list()
  for (pathway_name in requested) {
    if (pathway_name == "combined_death" && !isTRUE(include_combined)) {
      next
    }
    cfg <- get_pathway_config(pathway_name)
    out[[pathway_name]] <- render_pathway_map(
      config = cfg,
      prepared_de = prepared_de,
      outdir = outdir,
      title = title,
      p_cutoff = p_cutoff,
      show_compartments = show_compartments,
      min_nodes_per_compartment = min_nodes_per_compartment,
      force_singleton_compartments = force_singleton_compartments
    )
  }
  
  if (isTRUE(include_combined) && !("combined_death" %in% names(out))) {
    cfg <- get_pathway_config("combined_death")
    out[["combined_death"]] <- render_pathway_map(
      config = cfg,
      prepared_de = prepared_de,
      outdir = outdir,
      title = title,
      p_cutoff = p_cutoff,
      show_compartments = show_compartments,
      min_nodes_per_compartment = min_nodes_per_compartment,
      force_singleton_compartments = force_singleton_compartments
    )
  }
  
  if ("combined_death" %in% names(out)) {
    cfg <- get_pathway_config("combined_death")
    combined_plot <- render_pathway_map(
      config = cfg,
      prepared_de = prepared_de,
      outdir = outdir,
      title = title,
      p_cutoff = p_cutoff,
      show_compartments = show_compartments,
      min_nodes_per_compartment = min_nodes_per_compartment,
      force_singleton_compartments = force_singleton_compartments,
      return_plot = TRUE
    )
    
    p_interpret <- make_interpretation_panel(
      nodes_annot = combined_plot$node_tbl,
      edges_cfg = cfg$edges,
      title = "Interpretation and validation"
    )
    multipanel_png <- file.path(outdir, "combined_death_multipanel.png")
    save_multipanel_figure(
      p_left = combined_plot$plot,
      p_right = p_interpret,
      outfile = multipanel_png,
      width = 16,
      height = 8,
      dpi = 300,
      rel_widths = c(3.2, 1.4)
    )
    out[["combined_death"]]$multipanel_png <- multipanel_png
  }
  
  out
}

run_cgast_sting_map_interactive <- function() {
  message("Interactive mode: select your DESeq2 CSV.")
  input <- tryCatch(file.choose(), error = function(e) "")
  if (!nzchar(input)) {
    stop("No input file selected.", call. = FALSE)
  }
  
  outdir <- readline("Output directory [cgast_sting_out]: ")
  if (!nzchar(outdir)) {
    outdir <- "cgast_sting_out"
  }
  
  title <- readline("Plot title [cGAS-STING Pathway Map]: ")
  if (!nzchar(title)) {
    title <- "cGAS-STING Pathway Map"
  }
  
  run_cgast_sting_map(input = input, outdir = outdir, title = title)
}

main <- function() {
  parsed <- parse_args(commandArgs(trailingOnly = TRUE))
  
  if (isTRUE(parsed$help)) {
    print_help()
    return(invisible(NULL))
  }
  
  if (is.null(parsed$input) || !nzchar(parsed$input)) {
    stop(
      "Missing --input. In RStudio, run:\n",
      "  source(\"cgast_sting_map.R\")\n",
      "  run_cgast_sting_map_interactive()\n",
      call. = FALSE
    )
  }
  
  run_cgast_sting_map(
    input = parsed$input,
    outdir = parsed$outdir,
    title = parsed$title
  )
}

# Execute only when run as a script (not when sourced in RStudio).
if (sys.nframe() == 0L) {
  main()
}