# ============================================================
# oxysterol_ensembl_lookup.R
# Purpose: Map oxysterol pathway human gene symbols to chicken 
#          Ensembl IDs using GRCg7b annotation.
# 
# IMPORTANT: This uses the same Ensembl ID-based logic established
#            for other pathway panels to avoid the symbol-based
#            misclassification issue (cf. RIPK3/MLKL experience).
#
# Inputs:  - oxysterol_pathway_genes.csv (master gene list)
#          - Your existing DESeq2 results or count matrix with 
#            Ensembl IDs as row identifiers
#
# Outputs: - oxysterol_pathway_chicken_verified.csv
#          - oxysterol_pathway_genesets_ensembl.gmt (ready for scoring)
#          - Console report of found/missing/ambiguous genes
# ============================================================

library(org.Gg.eg.db)
library(AnnotationDbi)
library(dplyr)
library(readr)

# ── 1. Load the master gene list ────────────────────────────
gene_table <- read_csv("oxysterol_pathway_genes.csv", show_col_types = FALSE)

# Human symbols to look up
human_symbols <- unique(gene_table$Human_Symbol)

cat("=== Oxysterol Pathway: Chicken GRCg7b Annotation Lookup ===\n")
cat(sprintf("Genes to look up: %d\n\n", length(human_symbols)))

# ── 2. Strategy: Multi-pass lookup ──────────────────────────
# Pass 1: Direct SYMBOL match in org.Gg.eg.db
# Pass 2: ALIAS match for genes that fail Pass 1
# Pass 3: Manual curation for known problem genes
#
# The final classification must use Ensembl IDs present in 
# your count matrix, NOT just the annotation database.

# Pass 1: Direct symbol lookup
symbol_to_ensembl <- tryCatch({
  res <- AnnotationDbi::select(
    org.Gg.eg.db,
    keys = human_symbols,
    keytype = "SYMBOL",
    columns = c("ENSEMBL", "GENENAME", "ENTREZID")
  )
  res
}, error = function(e) {
  cat("Pass 1 (SYMBOL) returned errors for some keys — this is expected.\n")
  # Try one by one
  results <- lapply(human_symbols, function(sym) {
    tryCatch({
      AnnotationDbi::select(
        org.Gg.eg.db,
        keys = sym,
        keytype = "SYMBOL",
        columns = c("ENSEMBL", "GENENAME", "ENTREZID")
      )
    }, error = function(e2) {
      data.frame(SYMBOL = sym, ENSEMBL = NA, GENENAME = NA, ENTREZID = NA)
    })
  })
  bind_rows(results)
})

# Classify results
found_pass1 <- symbol_to_ensembl %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

missing_pass1 <- setdiff(human_symbols, found_pass1$SYMBOL)

cat(sprintf("Pass 1 (SYMBOL match): %d found, %d missing\n",
            nrow(found_pass1), length(missing_pass1)))

if (length(missing_pass1) > 0) {
  cat("Missing after Pass 1:", paste(missing_pass1, collapse = ", "), "\n\n")
}

# Pass 2: Try ALIAS lookup for missing genes
if (length(missing_pass1) > 0) {
  alias_results <- lapply(missing_pass1, function(sym) {
    tryCatch({
      AnnotationDbi::select(
        org.Gg.eg.db,
        keys = sym,
        keytype = "ALIAS",
        columns = c("SYMBOL", "ENSEMBL", "GENENAME", "ENTREZID")
      )
    }, error = function(e) {
      data.frame(ALIAS = sym, SYMBOL = NA, ENSEMBL = NA, GENENAME = NA, ENTREZID = NA)
    })
  })
  alias_df <- bind_rows(alias_results)
  found_pass2 <- alias_df %>% filter(!is.na(ENSEMBL))
  
  cat(sprintf("Pass 2 (ALIAS match): %d additional found\n", 
              n_distinct(found_pass2$ENSEMBL)))
  
  if (nrow(found_pass2) > 0) {
    # Merge back with original human symbol
    found_pass2$Human_Symbol <- found_pass2$ALIAS
    cat("ALIAS matches:\n")
    print(found_pass2 %>% select(Human_Symbol = ALIAS, Chicken_Symbol = SYMBOL, ENSEMBL))
  }
}

# ── 3. Cross-reference against your count matrix ────────────
# CRITICAL: A gene is only "Present" if its Ensembl ID exists
# in the row names of your count matrix / DESeq2 object.
#
# UPDATE THIS PATH to point to your counts file or DESeq2 results:

# counts_file <- "path/to/your/counts_matrix.csv"  # or .tsv
# If using DESeq2 results CSV:
# deseq_file <- "path/to/your/deseq2_results.csv"

# Example (uncomment and modify for your setup):
# count_ensembl_ids <- rownames(read.csv(counts_file, row.names = 1))
# -- OR --
# deseq_res <- read.csv(deseq_file)
# count_ensembl_ids <- deseq_res$Ensembl_ID  # adjust column name

# For now, create the lookup table without count matrix verification:
all_found <- found_pass1 %>%
  select(Human_Symbol = SYMBOL, Ensembl_ID = ENSEMBL, Gene_Name = GENENAME)

# ── 4. Build verified gene table ────────────────────────────
verified <- gene_table %>%
  left_join(all_found, by = "Human_Symbol") %>%
  mutate(
    Annotation_Status = case_when(
      !is.na(Ensembl_ID) ~ "Found_in_org.Gg.eg.db",
      TRUE ~ "NOT_FOUND — requires manual curation"
    )
  )

# ── 5. Report ───────────────────────────────────────────────
cat("\n=== ANNOTATION SUMMARY ===\n")
cat(sprintf("Total genes: %d\n", nrow(verified)))
cat(sprintf("Found (Ensembl ID): %d\n", sum(!is.na(verified$Ensembl_ID))))
cat(sprintf("Missing: %d\n", sum(is.na(verified$Ensembl_ID))))

cat("\nBy module:\n")
verified %>%
  group_by(Module) %>%
  summarise(
    Total = n(),
    Found = sum(!is.na(Ensembl_ID)),
    Missing = sum(is.na(Ensembl_ID)),
    .groups = "drop"
  ) %>%
  print(n = 20)

if (any(is.na(verified$Ensembl_ID))) {
  cat("\n*** GENES REQUIRING MANUAL CURATION ***\n")
  cat("These genes were not found by symbol or alias in org.Gg.eg.db.\n")
  cat("Check GRCg7b GTF/GFF directly, or search Ensembl BioMart for orthologs.\n")
  cat("(Recall: MB21D1/cGAS and RIPK3 had this same issue.)\n\n")
  verified %>%
    filter(is.na(Ensembl_ID)) %>%
    select(Module, Human_Symbol, Chicken_Symbol, Full_Name, Chicken_Annotation_Risk) %>%
    print(n = 50)
}

# ── 6. Save outputs ─────────────────────────────────────────
write_csv(verified, "oxysterol_pathway_chicken_verified.csv")
cat("\nSaved: oxysterol_pathway_chicken_verified.csv\n")

# ── 7. Build Ensembl-ID-based GMT for scoring ───────────────
# Only include genes with verified Ensembl IDs

# Define module-level gene sets
gmt_sets <- list(
  OXYSTEROL_FULL = verified %>% filter(!is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_BIOSYNTHESIS = verified %>% filter(Module == "1_Biosynthesis_Enzymatic", !is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_HEDGEHOG = verified %>% filter(Module == "2A_Hedgehog_Osteogenic", !is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_CHEMOTAXIS = verified %>% filter(Module == "2B_Immune_Chemotaxis" | Human_Symbol %in% c("CH25H", "CYP7B1"), !is.na(Ensembl_ID)) %>% pull(Ensembl_ID) %>% unique(),
  OXYSTEROL_LXR_EFFLUX = verified %>% filter(Human_Symbol %in% c("NR1H3", "NR1H2", "ABCA1", "ABCG1"), !is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_CHOLESTEROL_HOMEOSTASIS = verified %>% filter(Human_Symbol %in% c("SREBF2", "SCAP", "HMGCR", "INSIG1", "INSIG2", "LDLR", "HMGCS1"), !is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_TRANSPORT = verified %>% filter(Module == "4_Transport", !is.na(Ensembl_ID)) %>% pull(Ensembl_ID),
  OXYSTEROL_ISG_BRIDGE = verified %>% filter(Human_Symbol %in% c("CH25H", "CYP7B1", "GPR183"), !is.na(Ensembl_ID)) %>% pull(Ensembl_ID)
)

# Write GMT
gmt_lines <- sapply(names(gmt_sets), function(name) {
  genes <- gmt_sets[[name]]
  if (length(genes) == 0) return(NULL)
  paste(c(name, paste0("Oxysterol_", name), genes), collapse = "\t")
})
gmt_lines <- gmt_lines[!sapply(gmt_lines, is.null)]

writeLines(gmt_lines, "oxysterol_pathway_genesets_ensembl.gmt")
cat("Saved: oxysterol_pathway_genesets_ensembl.gmt\n")

cat("\nGene set sizes (Ensembl-verified only):\n")
for (name in names(gmt_sets)) {
  cat(sprintf("  %s: %d genes\n", name, length(gmt_sets[[name]])))
}

# ── 8. Integration with existing count matrix ───────────────
# Once you have your count matrix loaded, run this block to 
# classify genes as Present vs Missing:
#
# count_ids <- rownames(dds)  # or from your counts file
# verified$In_Count_Matrix <- verified$Ensembl_ID %in% count_ids
# 
# cat("\n=== COUNT MATRIX CROSS-REFERENCE ===\n")
# cat(sprintf("In count matrix: %d / %d\n", 
#     sum(verified$In_Count_Matrix, na.rm = TRUE),
#     sum(!is.na(verified$Ensembl_ID))))
# 
# # Genes with Ensembl IDs but NOT in count matrix
# # (may indicate annotation version mismatch)
# missing_from_counts <- verified %>%
#   filter(!is.na(Ensembl_ID), !In_Count_Matrix)
# if (nrow(missing_from_counts) > 0) {
#   cat("\n*** ANNOTATION MISMATCH: Found in org.Gg.eg.db but NOT in counts ***\n")
#   print(missing_from_counts %>% select(Human_Symbol, Ensembl_ID))
# }

cat("\n=== DONE ===\n")
cat("Next steps:\n")
cat("1. Run this script in your R environment with org.Gg.eg.db loaded\n")
cat("2. Manually curate any 'NOT_FOUND' genes via Ensembl BioMart\n")
cat("3. Cross-reference against your count matrix (uncomment Section 8)\n")
cat("4. Use the Ensembl GMT file for GSVA/ssGSEA/singscore scoring\n")
