pick_ids_for_go <- function(gene_annotation) {
  if ("Ensembl_GeneID" %in% colnames(gene_annotation)) {
    ensembl_ids <- strip_ens_version(trimws(as.character(gene_annotation$Ensembl_GeneID)))
    ensembl_ids[is.na(ensembl_ids)] <- ""
    if (any(nzchar(ensembl_ids))) {
      return(list(ids = ensembl_ids, keyType = "ENSEMBL", source_column = "Ensembl_GeneID"))
    }
  }

  if ("NCBI_GeneID" %in% colnames(gene_annotation)) {
    entrez_ids <- trimws(as.character(gene_annotation$NCBI_GeneID))
    entrez_ids[is.na(entrez_ids)] <- ""
    if (any(nzchar(entrez_ids))) {
      return(list(ids = entrez_ids, keyType = "ENTREZID", source_column = "NCBI_GeneID"))
    }
  }

  if ("SYMBOL" %in% colnames(gene_annotation)) {
    symbol_ids <- trimws(as.character(gene_annotation$SYMBOL))
    symbol_ids[is.na(symbol_ids)] <- ""
    if (any(nzchar(symbol_ids))) {
      return(list(ids = symbol_ids, keyType = "SYMBOL", source_column = "SYMBOL"))
    }
  }

  NULL
}

collapse_ranked_statistics <- function(ranked_stats, ids) {
  ranked_stats <- as.numeric(ranked_stats)
  ids <- as.character(ids)

  keep <- is.finite(ranked_stats) & ranked_stats != 0 & !is.na(ids) & nzchar(ids)
  if (!any(keep)) {
    return(numeric())
  }

  collapsed <- tapply(ranked_stats[keep], ids[keep], function(x) x[which.max(abs(x))])
  collapsed_vec <- as.numeric(collapsed)
  names(collapsed_vec) <- names(collapsed)
  sort(collapsed_vec, decreasing = TRUE)
}

map_ranked_ids_to_entrez <- function(ranked_stats, key_type, OrgDb) {
  if (length(ranked_stats) == 0L) {
    return(numeric())
  }

  if (identical(key_type, "ENTREZID")) {
    names(ranked_stats) <- as.character(names(ranked_stats))
    return(ranked_stats)
  }

  id_map <- tryCatch(
    clusterProfiler::bitr(
      names(ranked_stats),
      fromType = key_type,
      toType = "ENTREZID",
      OrgDb = OrgDb
    ),
    error = function(e) NULL
  )

  if (is.null(id_map) || nrow(id_map) == 0L) {
    return(numeric())
  }

  id_map <- id_map[!duplicated(id_map$ENTREZID), , drop = FALSE]
  mapped_stats <- ranked_stats[id_map[[key_type]]]
  collapse_ranked_statistics(mapped_stats, id_map$ENTREZID)
}

run_contrast_gsea <- function(raw_results, gene_annotation, config, label) {
  if (!isTRUE(config$run_gsea)) {
    return(list(status = "skipped", reason = "run_gsea is FALSE", messages = character()))
  }

  if (is.null(config$OrgDb)) {
    return(list(
      status = "skipped",
      reason = paste0("OrgDb package '", config$OrgDb_pkg, "' is not available"),
      messages = character()
    ))
  }

  ranked_stats <- raw_results$stat
  names(ranked_stats) <- rownames(raw_results)
  id_info <- pick_ids_for_go(gene_annotation)

  if (is.null(id_info)) {
    message("No usable gene IDs found for GSEA in ", label, "; skipping GSEA.")
    return(list(
      status = "skipped",
      reason = "No usable Ensembl, ENTREZ, or SYMBOL IDs were available",
      messages = character()
    ))
  }

  go_ranked_stats <- collapse_ranked_statistics(ranked_stats, id_info$ids)
  if (length(go_ranked_stats) == 0L) {
    return(list(
      status = "skipped",
      reason = "No non-zero ranked statistics remained after ID collapse",
      id_info = id_info,
      messages = character()
    ))
  }

  messages <- character()

  go_result <- tryCatch(
    clusterProfiler::gseGO(
      geneList = go_ranked_stats,
      keyType = id_info$keyType,
      ont = "ALL",
      OrgDb = config$OrgDb,
      minGSSize = 3,
      maxGSSize = 800,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    ),
    error = function(e) {
      messages <<- c(messages, paste("gseGO failed:", conditionMessage(e)))
      NULL
    }
  )

  kegg_ranked_stats <- map_ranked_ids_to_entrez(go_ranked_stats, id_info$keyType, config$OrgDb)
  if (length(kegg_ranked_stats) == 0L) {
    messages <- c(messages, "No ENTREZ IDs were available for KEGG GSEA.")
  }

  kegg_result <- NULL
  if (length(kegg_ranked_stats) > 0L) {
    kegg_result <- tryCatch(
      clusterProfiler::gseKEGG(
        geneList = kegg_ranked_stats,
        organism = config$kegg_org,
        keyType = "ncbi-geneid",
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH"
      ),
      error = function(e) {
        messages <<- c(messages, paste("gseKEGG failed:", conditionMessage(e)))
        NULL
      }
    )
  }

  list(
    status = "completed",
    reason = NA_character_,
    id_info = id_info,
    go = list(result = go_result, table = if (!is.null(go_result)) as.data.frame(go_result) else NULL),
    kegg = list(
      result = kegg_result,
      table = if (!is.null(kegg_result)) as.data.frame(kegg_result) else NULL
    ),
    messages = messages
  )
}

write_gsea_outputs <- function(gsea_result, output_dir, label) {
  if (is.null(gsea_result)) {
    return(invisible(NULL))
  }

  if (!is.null(gsea_result$go$table) && nrow(gsea_result$go$table) > 0L) {
    utils::write.csv(
      gsea_result$go$table,
      file.path(output_dir, paste0("GO_GSEA_", label, ".csv")),
      row.names = FALSE
    )

    grDevices::pdf(file.path(output_dir, paste0("GO_GSEA_dotplot_", label, ".pdf")))
    print(
      enrichplot::dotplot(gsea_result$go$result, showCategory = 10, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign)
    )
    grDevices::dev.off()
  }

  if (!is.null(gsea_result$kegg$table) && nrow(gsea_result$kegg$table) > 0L) {
    utils::write.csv(
      gsea_result$kegg$table,
      file.path(output_dir, paste0("KEGG_GSEA_", label, ".csv")),
      row.names = FALSE
    )

    grDevices::pdf(file.path(output_dir, paste0("KEGG_GSEA_dotplot_", label, ".pdf")))
    print(
      enrichplot::dotplot(gsea_result$kegg$result, showCategory = 10, split = ".sign") +
        ggplot2::facet_grid(. ~ .sign)
    )
    grDevices::dev.off()
  }

  status_lines <- c(
    paste("status:", gsea_result$status %||% "unknown"),
    if (!is.null(gsea_result$reason) && !is.na(gsea_result$reason)) paste("reason:", gsea_result$reason) else NULL,
    if (length(gsea_result$messages) > 0L) c("messages:", paste("-", gsea_result$messages)) else NULL
  )

  writeLines(
    status_lines,
    file.path(output_dir, paste0("GSEA_status_", label, ".txt"))
  )

  invisible(TRUE)
}
