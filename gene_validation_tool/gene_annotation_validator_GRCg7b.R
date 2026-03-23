##############################################################################
# gene_annotation_validator_GRCg7b.R
#
# Purpose : Validate pathway gene annotations in the GRCg7b chicken genome
#           assembly (Ensembl release 109 / latest available via biomaRt).
#           Classifies each gene as PRESENT, ALIASED, or ABSENT using a
#           four-step resolution cascade.
#
# Output  : gene_annotation_validator_GRCg7b.csv
#
# Author  : Generated for Rashid Lab chicken RNA-seq project
# Date    : 2026-03-19
##############################################################################

suppressPackageStartupMessages({
  library(biomaRt)
  library(AnnotationDbi)
  library(org.Gg.eg.db)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Null-coalescing operator (base R doesn't have %||% before R 4.4)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ─── 0. GENE LIST DEFINITION ─────────────────────────────────────────────────

pathway_genes <- list(
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "FADD", "CASP8", "TNFRSF1A", "TNFRSF1B",
    "ZBP1", "TLR3", "TLR4", "IFNAR1", "IFNAR2", "DAI"
  ),
  cGAS_STING_ISG = c(
    "MB21D1", "STING1", "TBK1", "IRF3", "IRF7", "IRF9", "STAT1", "STAT2",
    "MX1", "MX2", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "RSAD2",
    "ISG15", "ISG20", "IFITM1", "IFITM3", "OAS1", "OAS2", "OAS3",
    "CCL5", "CXCL10"
  ),
  Apoptosis = c(
    "CASP3", "CASP7", "CASP9", "BAX", "BID", "CYCS", "APAF1", "CASP8", "FADD"
  ),
  NFkB_Cytokine = c(
    "NFKB1", "RELA", "IKBKB", "IKBKG", "CHUK", "TRAF3", "TRAF6",
    "TBK1", "IL6", "TNF"
  ),
  Osteoblast_differentiation = c(
    "BMP2", "BMP4", "BMP7", "SMAD1", "SMAD5", "SMAD9", "RUNX2", "SP7",
    "COL1A1", "COL1A2", "ALPL", "BGLAP", "SOST", "DKK1", "WNT3A", "WNT5A",
    "FZD1", "FZD2", "LRP5", "LRP6", "AXIN2", "CTNNB1", "TWIST1",
    "ID1", "ID2", "ID3"
  ),
  Osteoclast_differentiation = c(
    "TNFSF11", "TNFRSF11A", "TNFRSF11B", "CSF1", "CSF1R", "NFATC1",
    "CTSK", "MMP9", "ACP5", "DCSTAMP", "ATP6V0D2", "SLC4A2",
    "CALCR", "ITGB3", "ITGAV"
  ),
  Inflammasome = c(
    "NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD",
    "P2RX7", "HMGB1", "IL33", "AIM2"
  )
)

# Flatten to data.frame and collapse duplicate symbols into comma-sep modules
gene_module_df <- lapply(names(pathway_genes), function(mod) {
  data.frame(human_symbol = pathway_genes[[mod]], module = mod,
             stringsAsFactors = FALSE)
}) |> bind_rows() |>
  group_by(human_symbol) |>
  summarise(pathway_module = paste(sort(unique(module)), collapse = ", "),
            .groups = "drop")

all_genes <- gene_module_df$human_symbol
cat(sprintf("\n=== Gene Annotation Validator — GRCg7b (Ensembl) ===\n"))
cat(sprintf("Total unique gene symbols to check : %d\n\n", length(all_genes)))

# ─── 1. CONNECT TO BIOMART ───────────────────────────────────────────────────

biomart_available <- FALSE
gg_mart  <- NULL
hs_mart  <- NULL
ensembl_version_str <- "not available"

cat("Connecting to Ensembl biomaRt ...\n")
tryCatch({
  gg_mart <- useEnsembl("genes", dataset = "ggallus_gene_ensembl")
  hs_mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
  # ensemblVersion() was removed in newer biomaRt; pull version from mart object
  ensembl_version_str <- tryCatch(
    as.character(gg_mart@host),
    error = function(e2) "unknown"
  )
  # Try the slot that stores dataset version info
  ensembl_version_str <- tryCatch({
    attrs <- listAttributes(gg_mart)
    ver   <- grep("version", attrs$description, ignore.case = TRUE, value = TRUE)
    if (length(ver) == 0) "unknown" else paste(gg_mart@biomart, gg_mart@host, sep = " @ ")
  }, error = function(e2) paste(gg_mart@biomart, gg_mart@host, sep = " @ "))
  cat(sprintf("  biomaRt connected. Ensembl : %s\n\n", ensembl_version_str))
  biomart_available <- TRUE
}, error = function(e) {
  cat(sprintf("  WARNING: biomaRt unavailable (%s).\n", conditionMessage(e)))
  cat("  Falling back to org.Gg.eg.db only.\n\n")
})

# ─── HELPER: clean up a biomaRt result frame ─────────────────────────────────

clean_bm <- function(df) {
  df <- df[!is.na(df[[1]]) & df[[1]] != "", , drop = FALSE]
  df[!duplicated(df), , drop = FALSE]
}

# ─── 2. STEP 1 — DIRECT SYMBOL MATCH IN CHICKEN MART ────────────────────────

results <- list()   # accumulate per-gene result rows

direct_hits <- data.frame()
unresolved  <- all_genes

if (biomart_available) {
  cat("Step 1 — Direct symbol match in ggallus_gene_ensembl ...\n")

  # Query by external_gene_name (chicken symbol)
  bm_direct <- tryCatch(
    getBM(
      attributes = c("external_gene_name", "ensembl_gene_id",
                     "hgnc_symbol", "gene_biotype"),
      filters    = "external_gene_name",
      values     = all_genes,
      mart       = gg_mart
    ) |> clean_bm(),
    error = function(e) {
      cat("  WARNING: direct query failed:", conditionMessage(e), "\n")
      data.frame()
    }
  )

  if (nrow(bm_direct) > 0) {
    # Genes whose query symbol exactly matches the chicken external_gene_name
    matched_symbols <- intersect(toupper(all_genes),
                                 toupper(bm_direct$external_gene_name))

    direct_hits <- bm_direct |>
      filter(toupper(external_gene_name) %in% matched_symbols) |>
      rename(chicken_symbol   = external_gene_name,
             chicken_ensembl_id = ensembl_gene_id)

    # Collapse multiple rows per symbol (e.g. multiple Ensembl IDs)
    direct_summary <- direct_hits |>
      group_by(toupper_sym = toupper(chicken_symbol)) |>
      summarise(
        chicken_ensembl_id = paste(unique(chicken_ensembl_id), collapse = "|"),
        chicken_symbol     = paste(unique(chicken_symbol),     collapse = "|"),
        biotype_flag       = ifelse(
          all(grepl("pseudogene", gene_biotype, ignore.case = TRUE)),
          "pseudogene only", ""
        ),
        multi_flag = ifelse(n_distinct(chicken_ensembl_id) > 1,
                            "multiple chicken orthologs found", ""),
        .groups = "drop"
      )

    for (sym in all_genes) {
      row <- direct_summary[direct_summary$toupper_sym == toupper(sym), ]
      if (nrow(row) > 0) {
        results[[sym]] <- data.frame(
          human_symbol       = sym,
          status             = "PRESENT",
          chicken_ensembl_id = row$chicken_ensembl_id,
          chicken_symbol     = row$chicken_symbol,
          resolution_method  = "direct_symbol",
          notes              = paste(
            Filter(nchar, c(row$biotype_flag, row$multi_flag)),
            collapse = "; "
          ),
          stringsAsFactors   = FALSE
        )
      }
    }

    resolved_step1 <- names(results)
    unresolved     <- setdiff(all_genes, resolved_step1)
    cat(sprintf("  Resolved : %d  |  Unresolved : %d\n\n",
                length(resolved_step1), length(unresolved)))
  } else {
    cat("  No direct hits returned.\n\n")
  }
}

# ─── 3. STEP 2 — HUMAN ORTHOLOG MAPPING ──────────────────────────────────────

if (biomart_available && length(unresolved) > 0) {
  cat(sprintf("Step 2 — Human ortholog mapping for %d genes ...\n",
              length(unresolved)))

  # Map human symbols → human Ensembl IDs
  hs_ids <- tryCatch(
    getBM(
      attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters    = "hgnc_symbol",
      values     = unresolved,
      mart       = hs_mart
    ) |> clean_bm(),
    error = function(e) {
      cat("  WARNING: human ID lookup failed:", conditionMessage(e), "\n")
      data.frame()
    }
  )

  if (nrow(hs_ids) > 0) {
    # Map human Ensembl IDs → chicken orthologs
    gg_ortho <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id",
                       "ggallus_homolog_ensembl_gene",
                       "ggallus_homolog_associated_gene_name",
                       "ggallus_homolog_orthology_type"),
        filters    = "ensembl_gene_id",
        values     = hs_ids$ensembl_gene_id,
        mart       = hs_mart
      ) |> clean_bm() |>
        filter(nchar(ggallus_homolog_ensembl_gene) > 0),
      error = function(e) {
        cat("  WARNING: ortholog query failed:", conditionMessage(e), "\n")
        data.frame()
      }
    )

    if (nrow(gg_ortho) > 0) {
      # Join: human symbol → human ensembl → chicken ensembl
      ortho_joined <- hs_ids |>
        inner_join(gg_ortho, by = "ensembl_gene_id") |>
        rename(
          human_ensembl      = ensembl_gene_id,
          chicken_ensembl_id = ggallus_homolog_ensembl_gene,
          chicken_symbol     = ggallus_homolog_associated_gene_name,
          orthology_type     = ggallus_homolog_orthology_type
        )

      for (sym in unresolved) {
        rows <- ortho_joined[toupper(ortho_joined$hgnc_symbol) == toupper(sym), ]
        if (nrow(rows) > 0) {
          multi <- nrow(rows) > 1
          notes_parts <- character(0)
          if (multi) notes_parts <- c(notes_parts, "multiple chicken orthologs found")
          ortho_types <- unique(rows$orthology_type)
          if (any(grepl("many", ortho_types, ignore.case = TRUE)))
            notes_parts <- c(notes_parts,
                             paste("orthology_type:", paste(ortho_types, collapse = "/")))

          results[[sym]] <- data.frame(
            human_symbol       = sym,
            status             = "ALIASED",
            chicken_ensembl_id = paste(unique(rows$chicken_ensembl_id), collapse = "|"),
            chicken_symbol     = paste(unique(rows$chicken_symbol),     collapse = "|"),
            resolution_method  = "human_ortholog",
            notes              = paste(notes_parts, collapse = "; "),
            stringsAsFactors   = FALSE
          )
        }
      }
    }
  }

  resolved_step2 <- setdiff(names(results), setdiff(all_genes, unresolved))
  unresolved     <- setdiff(unresolved, names(results))
  cat(sprintf("  Newly resolved : %d  |  Unresolved : %d\n\n",
              length(resolved_step2), length(unresolved)))
}

# ─── 4. STEP 3 — org.Gg.eg.db ALIAS LOOKUP ───────────────────────────────────

if (length(unresolved) > 0) {
  cat(sprintf("Step 3 — org.Gg.eg.db alias lookup for %d genes ...\n",
              length(unresolved)))

  # Try ALIAS keytype (catches synonyms)
  orgdb_alias <- tryCatch(
    AnnotationDbi::select(
      org.Gg.eg.db,
      keys    = unresolved,
      keytype = "ALIAS",
      columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME")
    ),
    error = function(e) NULL
  )

  # Also try direct SYMBOL lookup for any that failed ALIAS
  orgdb_symbol <- tryCatch(
    AnnotationDbi::select(
      org.Gg.eg.db,
      keys    = unresolved,
      keytype = "SYMBOL",
      columns = c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME")
    ),
    error = function(e) NULL
  )

  orgdb_combined <- bind_rows(
    if (!is.null(orgdb_alias))  orgdb_alias  |> mutate(query = ALIAS  %||% SYMBOL),
    if (!is.null(orgdb_symbol)) orgdb_symbol |> mutate(query = SYMBOL)
  ) |>
    filter(!is.na(ENSEMBL) & nchar(ENSEMBL) > 0) |>
    distinct()

  # Match query back to original human symbol (case-insensitive)
  for (sym in unresolved) {
    rows <- orgdb_combined[toupper(orgdb_combined$query) == toupper(sym) |
                             toupper(orgdb_combined$SYMBOL) == toupper(sym), ]
    rows <- rows[!is.na(rows$ENSEMBL), ]
    if (nrow(rows) > 0) {
      multi <- n_distinct(rows$ENSEMBL) > 1
      # If the chicken SYMBOL exactly matches the query (case-insensitive),
      # classify as PRESENT; otherwise it is an alias
      exact_sym_match <- any(toupper(rows$SYMBOL) == toupper(sym))
      status_step3    <- ifelse(exact_sym_match, "PRESENT", "ALIASED")
      method_step3    <- ifelse(exact_sym_match, "direct_symbol", "org_db_alias")
      results[[sym]] <- data.frame(
        human_symbol       = sym,
        status             = status_step3,
        chicken_ensembl_id = paste(unique(rows$ENSEMBL), collapse = "|"),
        chicken_symbol     = paste(unique(rows$SYMBOL),  collapse = "|"),
        resolution_method  = method_step3,
        notes              = ifelse(multi, "multiple chicken orthologs found", ""),
        stringsAsFactors   = FALSE
      )
    }
  }

  newly3     <- setdiff(names(results),
                        c(setdiff(all_genes, unresolved),
                          setdiff(unresolved, names(results))))
  unresolved <- setdiff(unresolved, names(results))
  cat(sprintf("  Newly resolved : %d  |  Unresolved : %d\n\n",
              length(setdiff(unresolved, setdiff(all_genes, names(results)))),
              length(unresolved)))
}

# ─── 5. STEP 4 — MARK REMAINING AS ABSENT ────────────────────────────────────

if (length(unresolved) > 0) {
  cat(sprintf("Step 4 — Marking %d genes as ABSENT.\n\n", length(unresolved)))
  for (sym in unresolved) {
    results[[sym]] <- data.frame(
      human_symbol       = sym,
      status             = "ABSENT",
      chicken_ensembl_id = NA_character_,
      chicken_symbol     = NA_character_,
      resolution_method  = "absent",
      notes              = "",
      stringsAsFactors   = FALSE
    )
  }
}

# ─── 6. ASSEMBLE FINAL TABLE ─────────────────────────────────────────────────

final_df <- bind_rows(results) |>
  left_join(gene_module_df, by = "human_symbol") |>
  select(human_symbol, status, chicken_ensembl_id, chicken_symbol,
         resolution_method, pathway_module, notes) |>
  arrange(pathway_module, status, human_symbol)

# ─── 7. WRITE CSV ────────────────────────────────────────────────────────────

out_file <- "gene_annotation_validator_GRCg7b.csv"

# Add header comment with version info
header_lines <- c(
  paste0("# gene_annotation_validator_GRCg7b.csv"),
  paste0("# Generated  : ", Sys.time()),
  paste0("# Ensembl ver: ", ensembl_version_str),
  paste0("# biomaRt    : ", ifelse(biomart_available, "available", "UNAVAILABLE — org.Gg.eg.db fallback only")),
  paste0("# Genome     : GRCg7b (Gallus gallus)"),
  ""
)
writeLines(header_lines, out_file)
suppressWarnings(
  write.table(final_df, file = out_file, sep = ",", quote = TRUE,
              row.names = FALSE, col.names = TRUE, append = TRUE)
)
cat(sprintf("Output written to: %s\n\n", out_file))

# ─── 8. CONSOLE SUMMARY ──────────────────────────────────────────────────────

sep  <- paste(rep("─", 60), collapse = "")
sep2 <- paste(rep("═", 60), collapse = "")

cat(sep2, "\n")
cat("SUMMARY REPORT — Gene Annotation Validator (GRCg7b)\n")
cat(sprintf("Ensembl version : %s\n", ensembl_version_str))
cat(sprintf("biomaRt status  : %s\n",
            ifelse(biomart_available, "connected", "unavailable — org.Gg.eg.db only")))
cat(sep2, "\n\n")

total <- nrow(final_df)
status_counts <- table(final_df$status)

cat(sprintf("Total genes checked : %d\n\n", total))
cat("Overall status breakdown:\n")
for (s in c("PRESENT", "ALIASED", "ABSENT")) {
  n   <- as.integer(status_counts[s])
  n   <- ifelse(is.na(n), 0L, n)
  pct <- round(100 * n / total, 1)
  cat(sprintf("  %-8s : %3d  (%5.1f%%)\n", s, n, pct))
}

cat("\n", sep, "\n", sep = "")
cat("Per-module breakdown:\n")
cat(sep, "\n")

modules_all <- sort(unique(unlist(strsplit(final_df$pathway_module, ", "))))
for (mod in modules_all) {
  mod_rows <- final_df[grepl(mod, final_df$pathway_module, fixed = TRUE), ]
  n_total   <- nrow(mod_rows)
  n_present <- sum(mod_rows$status == "PRESENT")
  n_aliased <- sum(mod_rows$status == "ALIASED")
  n_absent  <- sum(mod_rows$status == "ABSENT")
  cat(sprintf("\n  %-35s (n=%d)\n", mod, n_total))
  cat(sprintf("    PRESENT  : %d\n", n_present))
  cat(sprintf("    ALIASED  : %d\n", n_aliased))
  cat(sprintf("    ABSENT   : %d\n", n_absent))
  absent_genes <- mod_rows$human_symbol[mod_rows$status == "ABSENT"]
  if (length(absent_genes) > 0)
    cat(sprintf("    >> ABSENT : %s\n", paste(absent_genes, collapse = ", ")))
}

cat("\n", sep2, "\n", sep = "")
cat("ABSENT GENES (all modules):\n")
absent_all <- final_df[final_df$status == "ABSENT", "human_symbol"]
if (length(absent_all) == 0) {
  cat("  None — all genes resolved.\n")
} else {
  cat(sprintf("  %s\n", paste(sort(absent_all), collapse = "\n  ")))
}

cat("\n", sep, "\n", sep = "")
cat("ALIASED GENES (chicken symbol differs from human symbol):\n")
aliased_df <- final_df[final_df$status == "ALIASED" &
                          !is.na(final_df$chicken_symbol) &
                          final_df$chicken_symbol != final_df$human_symbol, ]
if (nrow(aliased_df) == 0) {
  cat("  None.\n")
} else {
  for (i in seq_len(nrow(aliased_df))) {
    cat(sprintf("  %-12s → %-20s [%s] (%s)\n",
                aliased_df$human_symbol[i],
                aliased_df$chicken_symbol[i],
                aliased_df$chicken_ensembl_id[i],
                aliased_df$resolution_method[i]))
  }
}
cat("\n", sep2, "\n", sep = "")
