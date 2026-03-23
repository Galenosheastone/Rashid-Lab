run_directional_enrichment_module <- function(
  tbl,
  up_sig,
  down_sig,
  gene_universe,
  orgDb,
  kegg_term2gene,
  kegg_term2name,
  minGSSize,
  maxGSSize,
  dirs,
  save_plot,
  register_if_exists,
  safe_slug
) {
  as_df_safe <- function(x) {
    if (is.null(x)) return(data.frame())
    out <- tryCatch(as.data.frame(x), error = function(e) data.frame())
    if (is.null(out) || nrow(out) == 0) return(data.frame())
    out
  }

  as_sheet <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(Note = "No results", stringsAsFactors = FALSE))
    }
    as.data.frame(df)
  }

  run_enrich_go <- function(genes, ont) {
    if (length(genes) == 0) return(NULL)
    tryCatch(
      clusterProfiler::enrichGO(
        gene = genes,
        OrgDb = orgDb,
        keyType = "ENTREZID",
        ont = ont,
        universe = gene_universe,
        pAdjustMethod = "BH",
        readable = TRUE,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      ),
      error = function(e) NULL
    )
  }

  run_enrich_kegg <- function(genes) {
    if (length(genes) == 0) return(NULL)
    out <- tryCatch(
      clusterProfiler::enricher(
        gene = genes,
        universe = gene_universe,
        TERM2GENE = kegg_term2gene,
        TERM2NAME = kegg_term2name,
        pAdjustMethod = "BH",
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      ),
      error = function(e) NULL
    )
    if (!is.null(out) && nrow(as_df_safe(out)) > 0) {
      out <- tryCatch(clusterProfiler::setReadable(out, orgDb, keyType = "ENTREZID"), error = function(e) out)
    }
    out
  }

  run_enrich_reactome <- function(genes) {
    if (length(genes) == 0) return(NULL)
    out <- tryCatch(
      ReactomePA::enrichPathway(
        gene = genes,
        universe = gene_universe,
        organism = "chicken",
        pAdjustMethod = "BH",
        readable = FALSE,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize
      ),
      error = function(e) NULL
    )
    if (!is.null(out) && nrow(as_df_safe(out)) > 0) {
      out <- tryCatch(clusterProfiler::setReadable(out, orgDb, keyType = "ENTREZID"), error = function(e) out)
    }
    out
  }

  make_dot2 <- function(res, stem, title) {
    df <- as_df_safe(res)
    if (nrow(df) == 0) return(invisible(NULL))
    p <- tryCatch(
      enrichplot::dotplot(res, showCategory = min(15, nrow(df))) + ggplot2::ggtitle(title),
      error = function(e) NULL
    )
    if (is.null(p)) return(invisible(NULL))
    save_plot(p, file.path(dirs$direction, paste0(safe_slug(stem), ".pdf")), 8, 6)
  }

  run_compare <- function(fun, ...) {
    tryCatch(
      clusterProfiler::compareCluster(
        geneCluster = cluster_list,
        fun = fun,
        ...
      ),
      error = function(e) NULL
    )
  }

  stack_block <- function(label, df) {
    if (is.null(df) || nrow(df) == 0) {
      return(data.frame(section = label, Note = "No results", stringsAsFactors = FALSE))
    }
    dplyr::bind_cols(data.frame(section = label, stringsAsFactors = FALSE), as.data.frame(df))
  }

  sig_up <- unique(as.character(up_sig$ENTREZID))
  sig_up <- sig_up[!is.na(sig_up) & nzchar(sig_up)]
  sig_down <- unique(as.character(down_sig$ENTREZID))
  sig_down <- sig_down[!is.na(sig_down) & nzchar(sig_down)]
  cluster_list <- list(Up = sig_up, Down = sig_down)
  cluster_list <- cluster_list[lengths(cluster_list) > 0]

  up_go_bp <- run_enrich_go(sig_up, "BP")
  down_go_bp <- run_enrich_go(sig_down, "BP")
  up_go_mf <- run_enrich_go(sig_up, "MF")
  down_go_mf <- run_enrich_go(sig_down, "MF")
  up_kegg <- run_enrich_kegg(sig_up)
  down_kegg <- run_enrich_kegg(sig_down)
  up_react <- run_enrich_reactome(sig_up)
  down_react <- run_enrich_reactome(sig_down)

  make_dot2(up_go_bp, "directional_go_bp_up", "GO BP ORA (Up)")
  make_dot2(down_go_bp, "directional_go_bp_down", "GO BP ORA (Down)")
  make_dot2(up_go_mf, "directional_go_mf_up", "GO MF ORA (Up)")
  make_dot2(down_go_mf, "directional_go_mf_down", "GO MF ORA (Down)")
  make_dot2(up_kegg, "directional_kegg_up", "KEGG ORA (Up)")
  make_dot2(down_kegg, "directional_kegg_down", "KEGG ORA (Down)")
  make_dot2(up_react, "directional_reactome_up", "Reactome ORA (Up)")
  make_dot2(down_react, "directional_reactome_down", "Reactome ORA (Down)")

  cmp_go_bp <- NULL
  cmp_go_mf <- NULL
  cmp_kegg <- NULL
  cmp_react <- NULL
  if (length(cluster_list) >= 2) {
    cmp_go_bp <- run_compare(
      "enrichGO",
      OrgDb = orgDb,
      keyType = "ENTREZID",
      ont = "BP",
      universe = gene_universe,
      pAdjustMethod = "BH",
      readable = TRUE,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    cmp_go_mf <- run_compare(
      "enrichGO",
      OrgDb = orgDb,
      keyType = "ENTREZID",
      ont = "MF",
      universe = gene_universe,
      pAdjustMethod = "BH",
      readable = TRUE,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    cmp_kegg <- run_compare(
      "enricher",
      universe = gene_universe,
      TERM2GENE = kegg_term2gene,
      TERM2NAME = kegg_term2name,
      pAdjustMethod = "BH",
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
    cmp_react <- run_compare(
      "enrichPathway",
      universe = gene_universe,
      organism = "chicken",
      pAdjustMethod = "BH",
      readable = FALSE,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    )
  }

  save_compare <- function(res, stem, title) {
    df <- as_df_safe(res)
    if (nrow(df) == 0) return(invisible(NULL))
    p <- tryCatch(enrichplot::dotplot(res, showCategory = min(20, nrow(df))) + ggplot2::ggtitle(title), error = function(e) NULL)
    if (is.null(p)) return(invisible(NULL))
    save_plot(p, file.path(dirs$direction, paste0(safe_slug(stem), ".pdf")), 9, 6)
  }

  save_compare(cmp_go_bp, "comparecluster_go_bp", "compareCluster GO BP (Up vs Down)")
  save_compare(cmp_go_mf, "comparecluster_go_mf", "compareCluster GO MF (Up vs Down)")
  save_compare(cmp_kegg, "comparecluster_kegg", "compareCluster KEGG (Up vs Down)")
  save_compare(cmp_react, "comparecluster_reactome", "compareCluster Reactome (Up vs Down)")

  go_bp_up_df <- as_df_safe(up_go_bp)
  go_bp_down_df <- as_df_safe(down_go_bp)
  go_bp_cmp_df <- as_df_safe(cmp_go_bp)
  go_mf_up_df <- as_df_safe(up_go_mf)
  go_mf_down_df <- as_df_safe(down_go_mf)
  go_mf_cmp_df <- as_df_safe(cmp_go_mf)
  kegg_up_df <- as_df_safe(up_kegg)
  kegg_down_df <- as_df_safe(down_kegg)
  kegg_cmp_df <- as_df_safe(cmp_kegg)
  react_up_df <- as_df_safe(up_react)
  react_down_df <- as_df_safe(down_react)
  react_cmp_df <- as_df_safe(cmp_react)

  go_bp_xlsx <- file.path(dirs$tables_sig, "directional_ORA_GO_BP.xlsx")
  go_mf_xlsx <- file.path(dirs$tables_sig, "directional_ORA_GO_MF.xlsx")
  kegg_xlsx <- file.path(dirs$tables_sig, "directional_ORA_KEGG.xlsx")
  react_xlsx <- file.path(dirs$tables_sig, "directional_ORA_Reactome.xlsx")

  openxlsx::write.xlsx(
    list(Up = as_sheet(go_bp_up_df), Down = as_sheet(go_bp_down_df), CompareCluster = as_sheet(go_bp_cmp_df)),
    file = go_bp_xlsx,
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    list(Up = as_sheet(go_mf_up_df), Down = as_sheet(go_mf_down_df), CompareCluster = as_sheet(go_mf_cmp_df)),
    file = go_mf_xlsx,
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    list(Up = as_sheet(kegg_up_df), Down = as_sheet(kegg_down_df), CompareCluster = as_sheet(kegg_cmp_df)),
    file = kegg_xlsx,
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    list(Up = as_sheet(react_up_df), Down = as_sheet(react_down_df), CompareCluster = as_sheet(react_cmp_df)),
    file = react_xlsx,
    overwrite = TRUE
  )

  register_if_exists("table", "directional_ORA_GO_BP", go_bp_xlsx, "Directional ORA GO BP.")
  register_if_exists("table", "directional_ORA_GO_MF", go_mf_xlsx, "Directional ORA GO MF.")
  register_if_exists("table", "directional_ORA_KEGG", kegg_xlsx, "Directional ORA KEGG.")
  register_if_exists("table", "directional_ORA_Reactome", react_xlsx, "Directional ORA Reactome.")

  stacked <- dplyr::bind_rows(
    stack_block("GO_BP_Up", go_bp_up_df),
    stack_block("GO_BP_Down", go_bp_down_df),
    stack_block("GO_BP_CompareCluster", go_bp_cmp_df),
    stack_block("GO_MF_Up", go_mf_up_df),
    stack_block("GO_MF_Down", go_mf_down_df),
    stack_block("GO_MF_CompareCluster", go_mf_cmp_df),
    stack_block("KEGG_Up", kegg_up_df),
    stack_block("KEGG_Down", kegg_down_df),
    stack_block("KEGG_CompareCluster", kegg_cmp_df),
    stack_block("Reactome_Up", react_up_df),
    stack_block("Reactome_Down", react_down_df),
    stack_block("Reactome_CompareCluster", react_cmp_df)
  )

  list(
    go_bp = list(up = go_bp_up_df, down = go_bp_down_df, compare = go_bp_cmp_df),
    go_mf = list(up = go_mf_up_df, down = go_mf_down_df, compare = go_mf_cmp_df),
    kegg = list(up = kegg_up_df, down = kegg_down_df, compare = kegg_cmp_df),
    reactome = list(up = react_up_df, down = react_down_df, compare = react_cmp_df),
    stacked = stacked
  )
}
