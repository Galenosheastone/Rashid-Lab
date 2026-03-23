run_activity_inference_module <- function(
  tbl,
  tbl_nodup,
  rank_col,
  activity_use_human_orthologs,
  dirs,
  save_plot,
  register_if_exists
) {
  out <- list(
    ortholog_map = data.frame(),
    tf_scores = data.frame(),
    progeny_scores = data.frame(),
    notes = data.frame(),
    tf_xlsx = NULL,
    progeny_xlsx = NULL,
    ortholog_xlsx = NULL
  )

  note_lines <- character(0)
  add_note <- function(msg) {
    note_lines <<- c(note_lines, as.character(msg))
    message(msg)
  }

  write_note_file <- function() {
    if (length(note_lines) == 0) return(invisible(NULL))
    note_file <- file.path(dirs$tables_act, "activity_inference_note.txt")
    writeLines(note_lines, note_file)
    register_if_exists("table", "activity_inference_note", note_file, "Activity module notes.")
    out$notes <<- data.frame(message = note_lines, stringsAsFactors = FALSE)
    invisible(note_file)
  }

  if (!"SYMBOL" %in% names(tbl_nodup)) tbl_nodup$SYMBOL <- NA_character_
  if (!"gene_label" %in% names(tbl_nodup)) tbl_nodup$gene_label <- tbl_nodup$gene
  if (!"gene" %in% names(tbl_nodup)) tbl_nodup$gene <- tbl_nodup$gene_label

  if (!(rank_col %in% names(tbl_nodup))) {
    rank_col <- if ("stat" %in% names(tbl_nodup)) "stat" else "logFC"
  }
  if (!(rank_col %in% names(tbl_nodup))) {
    add_note("Activity inference skipped: no usable rank/stat column in tbl_nodup.")
    write_note_file()
    return(out)
  }

  tbl_work <- tbl_nodup %>%
    dplyr::mutate(
      input_symbol = dplyr::coalesce(.data$SYMBOL, .data$gene_label, .data$gene),
      stat_for_activity = suppressWarnings(as.numeric(.data[[rank_col]]))
    ) %>%
    dplyr::filter(!is.na(input_symbol), nzchar(input_symbol), !is.na(stat_for_activity))

  if (nrow(tbl_work) == 0) {
    add_note("Activity inference skipped: no non-missing symbols/stat values after filtering.")
    write_note_file()
    return(out)
  }

  ortholog_map <- data.frame()
  target_symbol_col <- "input_symbol"
  if (isTRUE(activity_use_human_orthologs)) {
    if (requireNamespace("babelgene", quietly = TRUE)) {
      species_opts <- c("Gallus gallus", "chicken", "9031")
      ortho_try <- NULL
      for (sp in species_opts) {
        ortho_try <- tryCatch(
          babelgene::orthologs(
            genes = unique(tbl_work$input_symbol),
            species = sp,
            human = TRUE,
            min_support = 1,
            top = TRUE
          ),
          error = function(e) NULL
        )
        if (!is.null(ortho_try) && nrow(as.data.frame(ortho_try)) > 0) break
      }

      if (!is.null(ortho_try) && nrow(as.data.frame(ortho_try)) > 0) {
        ortholog_map <- as.data.frame(ortho_try)
        if (!"symbol" %in% names(ortholog_map)) ortholog_map$symbol <- NA_character_
        if (!"human_symbol" %in% names(ortholog_map)) ortholog_map$human_symbol <- NA_character_
        if (!"entrez" %in% names(ortholog_map)) ortholog_map$entrez <- NA_character_
        if (!"human_entrez" %in% names(ortholog_map)) ortholog_map$human_entrez <- NA_character_
        if (!"support" %in% names(ortholog_map)) ortholog_map$support <- NA_character_
        if (!"support_n" %in% names(ortholog_map)) ortholog_map$support_n <- NA_character_

        ortholog_map <- ortholog_map %>%
          dplyr::transmute(
            chicken_symbol = .data$symbol,
            chicken_entrez = as.character(.data$entrez),
            human_symbol = .data$human_symbol,
            human_entrez = as.character(.data$human_entrez),
            support = as.character(.data$support),
            support_n = as.character(.data$support_n)
          ) %>%
          dplyr::distinct()

        tbl_work <- tbl_work %>%
          dplyr::left_join(ortholog_map %>% dplyr::select(chicken_symbol, human_symbol), by = c("input_symbol" = "chicken_symbol")) %>%
          dplyr::mutate(target_symbol = dplyr::coalesce(.data$human_symbol, .data$input_symbol))
        target_symbol_col <- "target_symbol"
      } else {
        add_note("Ortholog mapping failed or returned no rows; using input symbols directly for activity scoring.")
        tbl_work$target_symbol <- tbl_work$input_symbol
        target_symbol_col <- "target_symbol"
      }
    } else {
      add_note("babelgene not installed; using input symbols directly for activity scoring.")
      tbl_work$target_symbol <- tbl_work$input_symbol
      target_symbol_col <- "target_symbol"
    }
  } else {
    tbl_work$target_symbol <- tbl_work$input_symbol
    target_symbol_col <- "target_symbol"
  }

  ortholog_xlsx <- file.path(dirs$tables_act, "chicken_to_human_orthologs.xlsx")
  openxlsx::write.xlsx(
    if (nrow(ortholog_map) > 0) ortholog_map else data.frame(Note = "No ortholog mapping generated", stringsAsFactors = FALSE),
    ortholog_xlsx,
    overwrite = TRUE
  )
  register_if_exists("table", "chicken_to_human_orthologs", ortholog_xlsx, "Chicken to human ortholog mapping.")
  out$ortholog_xlsx <- ortholog_xlsx

  stat_tbl <- tbl_work %>%
    dplyr::group_by(target_symbol = .data[[target_symbol_col]]) %>%
    dplyr::summarise(stat = mean(stat_for_activity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(!is.na(target_symbol), nzchar(target_symbol), !is.na(stat))

  if (nrow(stat_tbl) < 10) {
    add_note("Activity inference skipped: insufficient mapped symbols for ULM.")
    write_note_file()
    out$ortholog_map <- ortholog_map
    return(out)
  }

  mat_stat <- matrix(stat_tbl$stat, ncol = 1)
  rownames(mat_stat) <- stat_tbl$target_symbol
  colnames(mat_stat) <- "contrast_stat"

  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    add_note("decoupleR not installed; skipping TF and PROGENy activity inference.")
    write_note_file()
    out$ortholog_map <- ortholog_map
    return(out)
  }

  standardize_net <- function(net, default_source = NULL) {
    if (is.null(net)) return(data.frame())
    net <- as.data.frame(net)
    if (nrow(net) == 0) return(data.frame())

    if (!"source" %in% names(net)) {
      if ("tf" %in% names(net)) net$source <- net$tf
      if ("pathway" %in% names(net)) net$source <- net$pathway
      if ("signature" %in% names(net)) net$source <- net$signature
    }
    if (!"target" %in% names(net)) {
      if ("gene" %in% names(net)) net$target <- net$gene
      if ("target_genes" %in% names(net)) net$target <- net$target_genes
    }
    if (!"mor" %in% names(net)) {
      if ("weight" %in% names(net)) net$mor <- as.numeric(net$weight)
      else if ("likelihood" %in% names(net)) net$mor <- as.numeric(net$likelihood)
      else net$mor <- 1
    }

    if (!all(c("source", "target", "mor") %in% names(net))) return(data.frame())

    if (!is.null(default_source) && "source" %in% names(net)) {
      net$source <- ifelse(is.na(net$source) | net$source == "", default_source, net$source)
    }

    net <- net %>%
      dplyr::transmute(
        source = as.character(.data$source),
        target = as.character(.data$target),
        mor = suppressWarnings(as.numeric(.data$mor)),
        confidence = if ("confidence" %in% names(net)) as.character(.data$confidence) else NA_character_
      ) %>%
      dplyr::filter(!is.na(source), nzchar(source), !is.na(target), nzchar(target), !is.na(mor))

    net
  }

  net_tf <- tryCatch(decoupleR::get_collectri(organism = "human"), error = function(e) NULL)
  net_tf <- standardize_net(net_tf)

    if (nrow(net_tf) == 0 && requireNamespace("dorothea", quietly = TRUE)) {
    dorothea_net <- tryCatch(get("dorothea_hs", envir = asNamespace("dorothea")), error = function(e) NULL)
    if (!is.null(dorothea_net)) {
      dorothea_net <- as.data.frame(dorothea_net)
      if (!"source" %in% names(dorothea_net) && "tf" %in% names(dorothea_net)) dorothea_net$source <- dorothea_net$tf
      if ("confidence" %in% names(dorothea_net)) {
        dorothea_net <- dorothea_net[dorothea_net$confidence %in% c("A", "B", "C"), , drop = FALSE]
      }
      net_tf <- standardize_net(dorothea_net)
    }
  }

  net_progeny <- tryCatch(decoupleR::get_progeny(organism = "human", top = 500), error = function(e) NULL)
  net_progeny <- standardize_net(net_progeny)

  run_ulm_safe <- function(mat, net) {
    if (is.null(net) || nrow(net) == 0) return(data.frame())
    net <- net[net$target %in% rownames(mat), , drop = FALSE]
    if (nrow(net) == 0) return(data.frame())

    ulm_res <- tryCatch(
      decoupleR::run_ulm(
        mat = mat,
        net = net,
        .source = "source",
        .target = "target",
        .mor = "mor"
      ),
      error = function(e1) {
        tryCatch(
          decoupleR::run_ulm(
            mat = mat,
            network = net,
            .source = "source",
            .target = "target",
            .mor = "mor"
          ),
          error = function(e2) NULL
        )
      }
    )

    if (is.null(ulm_res)) return(data.frame())
    ulm_df <- as.data.frame(ulm_res)
    if (nrow(ulm_df) == 0) return(data.frame())

    if (!"regulator" %in% names(ulm_df) && "source" %in% names(ulm_df)) ulm_df$regulator <- ulm_df$source
    if (!"activity" %in% names(ulm_df) && "score" %in% names(ulm_df)) ulm_df$activity <- suppressWarnings(as.numeric(ulm_df$score))
    if (!"activity" %in% names(ulm_df) && "estimate" %in% names(ulm_df)) ulm_df$activity <- suppressWarnings(as.numeric(ulm_df$estimate))
    if (!"activity" %in% names(ulm_df) && "statistic" %in% names(ulm_df)) ulm_df$activity <- suppressWarnings(as.numeric(ulm_df$statistic))

    if (!"pvalue" %in% names(ulm_df) && "p_value" %in% names(ulm_df)) ulm_df$pvalue <- suppressWarnings(as.numeric(ulm_df$p_value))
    if (!"pvalue" %in% names(ulm_df) && "pval" %in% names(ulm_df)) ulm_df$pvalue <- suppressWarnings(as.numeric(ulm_df$pval))
    if (!"padj" %in% names(ulm_df) && "pvalue" %in% names(ulm_df)) ulm_df$padj <- stats::p.adjust(ulm_df$pvalue, method = "BH")

    ulm_df
  }

  tf_scores <- run_ulm_safe(mat_stat, net_tf)
  progeny_scores <- run_ulm_safe(mat_stat, net_progeny)

  tf_xlsx <- file.path(dirs$tables_act, "tf_activity_ulm.xlsx")
  progeny_xlsx <- file.path(dirs$tables_act, "progeny_activity_ulm.xlsx")

  openxlsx::write.xlsx(
    if (nrow(tf_scores) > 0) tf_scores else data.frame(Note = "No TF activity results", stringsAsFactors = FALSE),
    tf_xlsx,
    overwrite = TRUE
  )
  openxlsx::write.xlsx(
    if (nrow(progeny_scores) > 0) progeny_scores else data.frame(Note = "No PROGENy activity results", stringsAsFactors = FALSE),
    progeny_xlsx,
    overwrite = TRUE
  )
  register_if_exists("table", "tf_activity_ulm", tf_xlsx, "ULM TF activity scores.")
  register_if_exists("table", "progeny_activity_ulm", progeny_xlsx, "ULM PROGENy pathway activity scores.")

  plot_activity <- function(df, out_pdf, title, top_n = 20) {
    if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

    reg_col <- if ("regulator" %in% names(df)) "regulator" else if ("source" %in% names(df)) "source" else names(df)[1]
    act_col <- if ("activity" %in% names(df)) "activity" else if ("score" %in% names(df)) "score" else if ("estimate" %in% names(df)) "estimate" else if ("statistic" %in% names(df)) "statistic" else NULL
    if (is.null(act_col)) return(invisible(NULL))

    d <- df %>%
      dplyr::mutate(activity_value = suppressWarnings(as.numeric(.data[[act_col]]))) %>%
      dplyr::filter(!is.na(activity_value), !is.na(.data[[reg_col]]), nzchar(as.character(.data[[reg_col]]))) %>%
      dplyr::arrange(dplyr::desc(abs(activity_value))) %>%
      dplyr::slice_head(n = top_n)

    if (nrow(d) == 0) return(invisible(NULL))

    d <- d %>%
      dplyr::mutate(
        sig_label = if ("padj" %in% names(d)) dplyr::if_else(!is.na(.data$padj) & .data$padj < 0.05, "*", "") else ""
      )

    p <- ggplot2::ggplot(d, ggplot2::aes(x = stats::reorder(.data[[reg_col]], activity_value), y = activity_value, fill = activity_value > 0)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "royalblue"), guide = "none") +
      ggplot2::geom_text(ggplot2::aes(label = sig_label), hjust = -0.2, size = 3) +
      ggplot2::labs(title = title, x = NULL, y = "Activity score") +
      ggplot2::theme_bw(base_size = 11)

    save_plot(p, out_pdf, 9, 6)
  }

  plot_activity(tf_scores, file.path(dirs$activity, "tf_activity_top20.pdf"), "Top TF activities (ULM)", top_n = 20)
  plot_activity(progeny_scores, file.path(dirs$activity, "progeny_activity_top14.pdf"), "PROGENy pathway activities (ULM)", top_n = 14)

  if (nrow(tf_scores) == 0) add_note("No TF activity scores were produced (network overlap may be too low).")
  if (nrow(progeny_scores) == 0) add_note("No PROGENy activity scores were produced (network overlap may be too low).")
  write_note_file()

  out$ortholog_map <- ortholog_map
  out$tf_scores <- tf_scores
  out$progeny_scores <- progeny_scores
  out$tf_xlsx <- tf_xlsx
  out$progeny_xlsx <- progeny_xlsx
  out
}
