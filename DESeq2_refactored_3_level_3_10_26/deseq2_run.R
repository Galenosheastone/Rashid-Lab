## -----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)
set.seed(666)

# ==============================================================================
# Configuration (EDIT THESE FOR EACH PROJECT, or override via params)
# ==============================================================================
output_dir            <- "deseq2_outputs"
counts_file           <- "Sacral24_realigned_GRCg7b_counts_matchedIDs_MAIN_forDESeq2.csv"  # CSV/TSV with 'Geneid' + sample columns
metadata_file         <- "GRCg7b_metadata_CLEAN.csv"                                      # Optional CSV with rownames = sample IDs; must include 'condition'
ref_level_opt         <- NA_character_                                                    # set to a level name or NA to auto-pick Free if present
ssrnaseq_output_file  <- "DESeq2_full_combined_ssRNAseq.csv"
results_alpha         <- 0.05
min_count             <- 10
min_reps_within_group <- 3
contrast_plan  <- data.frame(
  group1 = c("Sacral", "Pygo", "Sacral"),
  group0 = c("Free",   "Free", "Pygo"),
  stringsAsFactors = FALSE
)
filter_mode    <- "strict"   # "legacy" (= ≥1 read anywhere) or "strict" (= ≥ min_count reads in ≥ min_reps_within_group samples within at least one condition)
use_rlog       <- FALSE      # if TRUE use rlog; else VST
run_lrt        <- TRUE       # Likelihood Ratio Test vs ~1
run_gsea       <- TRUE       # GO + KEGG GSEA using clusterProfiler
OrgDb_pkg      <- "org.Gg.eg.db" # Annotation package (species-specific)
kegg_org       <- "gga"          # KEGG organism code (e.g., "hsa","mmu","gga")
shrink_method  <- "apeglm"       # "apeglm" or "ashr"
n_workers      <- 4              # parallel workers

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Package management
# ==============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

core_pkgs <- c("DESeq2","apeglm","ashr","AnnotationDbi",
               "clusterProfiler","enrichplot","pathview")
qc_pkgs   <- c("ggplot2","pheatmap","vsn")
para_pkgs <- c("BiocParallel")

for (pkg in c(core_pkgs, qc_pkgs, para_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = TRUE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Try to load OrgDb if present (used only if run_gsea == TRUE)
OrgDb <- NULL
if (run_gsea) {
  if (!requireNamespace(OrgDb_pkg, quietly = TRUE)) {
    BiocManager::install(OrgDb_pkg, ask = FALSE, update = TRUE)
  }
  ok <- suppressWarnings(suppressMessages(
    require(OrgDb_pkg, character.only = TRUE)
  ))
  if (ok) {
    # More robust: org.*.eg.db packages usually export an OrgDb object with the same name as the package
    OrgDb <- tryCatch(get(OrgDb_pkg, envir = asNamespace(OrgDb_pkg)),
                      error = function(e) {
                        tryCatch(get(OrgDb_pkg, envir = .GlobalEnv), error = function(e2) NULL)
                      })
    if (is.null(OrgDb)) message("Could not retrieve OrgDb object from ", OrgDb_pkg, "; GSEA will be skipped.")
  } else {
    message("Annotation package ", OrgDb_pkg, " not available; GSEA will be skipped.")
  }
}

# ==============================================================================
# Parallel registration (portable)
# ==============================================================================
if (.Platform$OS.type == "windows") {
  BiocParallel::register(BiocParallel::SnowParam(workers = n_workers))
} else {
  BiocParallel::register(BiocParallel::MulticoreParam(workers = n_workers))
}

# ==============================================================================
# Helpers for sample/condition harmonization
# ==============================================================================
clean_names <- function(x) {
  x <- basename(x)
  # remove common file suffixes/extensions
  x <- sub("\\.(bam|sam|counts(?:\\.txt)?|txt|fq|fastq|gz)$", "", x, ignore.case = TRUE)
  x <- sub("\\.sorted$", "", x, ignore.case = TRUE)
  # remove trailing lane/read tokens: _S#, _R1/_R2, _001 (and optional combos)
  x <- sub("(_S\\d+)?(_R[12])?(_00[12])?$", "", x, ignore.case = TRUE)
  x <- sub("(_S\\d+)?(_R[12])?(_001)?$",  "", x, ignore.case = TRUE)
  # optional aligner tokens at the very end
  x <- sub("(\\.aligned)?$", "", x, ignore.case = TRUE)
  x
}

canonicalize_condition <- function(x) {
  x_chr <- trimws(as.character(x))
  x_chr[is.na(x_chr)] <- ""
  x_key <- tolower(gsub("[^a-z0-9]+", "", x_chr))

  out <- x_chr
  out[grepl("^sacral", x_key) | grepl("^exp", x_key)] <- "Sacral"
  out[grepl("^free",   x_key)] <- "Free"
  out[grepl("^pygo",   x_key)] <- "Pygo"
  out
}

infer_condition_from_sample <- function(x) {
  prefixes <- sub("_.*$", "", x)
  canonicalize_condition(prefixes)
}

keep_by_group_abundance <- function(dds_obj, min_count = 10, min_reps = 3) {
  stopifnot(min_count >= 1, min_reps >= 1)
  count_mat <- counts(dds_obj)
  keep_list <- lapply(levels(dds_obj$condition), function(cond) {
    idx <- dds_obj$condition == cond
    rowSums(count_mat[, idx, drop = FALSE] >= min_count) >= min_reps
  })
  Reduce("|", keep_list)
}

# ==============================================================================
# I/O — counts and metadata
# ==============================================================================
# Read counts flexibly (csv/tsv)
ext <- tools::file_ext(counts_file)
Expression <- if (tolower(ext) %in% c("csv")) {
  read.csv(counts_file, check.names = FALSE)
} else {
  read.delim(counts_file, check.names = FALSE)
}

stopifnot("Geneid" %in% colnames(Expression))

# ---- FIX: sanitize Geneid to avoid blank/duplicate rownames ----
# Behavior:
#   1) Trim whitespace in Geneid
#   2) If blank/NA, fallback to Ensembl_GeneID (first token before ';') if present
#   3) If still blank, assign UNLABELED_000001 ...
#   4) Make everything unique with make.unique()
first_token <- function(x) sub(";.*$", "", x)

geneid_original <- as.character(Expression$Geneid)
geneid_clean <- trimws(geneid_original)
geneid_clean[is.na(geneid_clean)] <- ""

# fallback to Ensembl if Geneid missing/blank
if ("Ensembl_GeneID" %in% colnames(Expression)) {
  ens <- as.character(Expression$Ensembl_GeneID)
  ens <- trimws(ens)
  ens[is.na(ens)] <- ""
  ens_first <- first_token(ens)
  missing_geneid <- (geneid_clean == "")
  geneid_clean[missing_geneid & ens_first != ""] <- ens_first[missing_geneid & ens_first != ""]
}

# final fallback: UNLABELED ids
still_missing <- (geneid_clean == "")
if (any(still_missing)) {
  geneid_clean[still_missing] <- sprintf("UNLABELED_%06d", which(still_missing))
}

# ensure unique
geneid_clean <- make.unique(geneid_clean)

# Keep a copy of annotation columns if present
anno_cols <- intersect(colnames(Expression),
  c("SYMBOL","NCBI_GeneID","Ensembl_GeneID","gene_biotype","Chr","Start","End","Strand","Length")
)

# Gene annotation table (safe rownames now)
gene_annot <- Expression[, c("Geneid", anno_cols), drop = FALSE]
gene_annot$Geneid_original <- geneid_original
gene_annot$Geneid_clean <- geneid_clean
rownames(gene_annot) <- geneid_clean

# ---- CRITICAL FIX: ensure SYMBOL exists so downstream subsetting never fails ----
if (!"SYMBOL" %in% colnames(gene_annot)) gene_annot$SYMBOL <- NA_character_
if (!"Ensembl_GeneID" %in% colnames(gene_annot)) gene_annot$Ensembl_GeneID <- NA_character_

# Optional: if SYMBOL missing but we have OrgDb + Ensembl IDs, map ENSEMBL -> SYMBOL
# (This is optional, but improves *_by_SYMBOL exports and sometimes GSEA ID handling.)
strip_ens_version <- function(x) sub("\\.[0-9]+$", "", x)
if (all(is.na(gene_annot$SYMBOL)) && !is.null(OrgDb) &&
    "Ensembl_GeneID" %in% colnames(gene_annot) && any(nzchar(gene_annot$Ensembl_GeneID))) {

  ens_ids <- strip_ens_version(as.character(gene_annot$Ensembl_GeneID))
  ens_ids[is.na(ens_ids)] <- ""
  ok_ens <- nzchar(ens_ids)

  if (any(ok_ens)) {
    sym_map <- tryCatch(
      AnnotationDbi::mapIds(OrgDb,
                            keys = unique(ens_ids[ok_ens]),
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first"),
      error = function(e) NULL
    )
    if (!is.null(sym_map)) {
      gene_annot$SYMBOL[ok_ens] <- unname(sym_map[ens_ids[ok_ens]])
    }
  }
}

# ---- Build counts matrix (exclude annotation; keep only numeric sample columns)
rownames(Expression) <- geneid_clean
sample_cols <- setdiff(colnames(Expression), c("Geneid", anno_cols))

is_num <- vapply(Expression[, sample_cols, drop = FALSE],
                 function(x) is.numeric(x) || is.integer(x),
                 logical(1))
if (!all(is_num)) {
  msg <- sprintf("Non-numeric columns found among putative sample columns: %s",
                 paste(sample_cols[!is_num], collapse = ", "))
  warning(msg)
}
sample_cols <- sample_cols[is_num]

cts <- as.matrix(Expression[, sample_cols, drop = FALSE])
storage.mode(cts) <- "integer"

# ==============================================================================
# Harmonize sample names (paths/extensions) between counts and metadata
# ==============================================================================
colnames(cts) <- clean_names(colnames(cts))

if (file.exists(metadata_file)) {
  coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  stopifnot("condition" %in% colnames(coldata))
  rownames(coldata) <- clean_names(rownames(coldata))
} else {
  message("Metadata file '", metadata_file, "' not found; inferring condition labels from sample names.")
  coldata <- data.frame(
    condition = infer_condition_from_sample(colnames(cts)),
    row.names = colnames(cts),
    stringsAsFactors = FALSE
  )
}

# Quick sanity preview (helps debugging if something still mismatches)
cat("Counts IDs (first 12):\n"); print(sort(colnames(cts))[1:min(12, ncol(cts))])
cat("Metadata IDs:\n"); print(sort(rownames(coldata)))

# Check set equality, then enforce identical order
if (!setequal(colnames(cts), rownames(coldata))) {
  stop(
    "Counts columns and metadata sample IDs don't match even after cleaning.\n",
    "Not in metadata: ", paste(setdiff(colnames(cts), rownames(coldata)), collapse = ", "), "\n",
    "Not in counts: ",   paste(setdiff(rownames(coldata), colnames(cts)), collapse = ", ")
  )
}
coldata <- coldata[colnames(cts), , drop = FALSE]

# ==============================================================================
# Factor setup — sanitize labels, normalize aliases, and choose reference
# ==============================================================================
coldata$condition_original <- as.character(coldata$condition)
coldata$condition <- canonicalize_condition(coldata$condition_original)
coldata$condition <- factor(coldata$condition)
levels(coldata$condition) <- make.names(levels(coldata$condition))

contrast_plan$group1 <- make.names(canonicalize_condition(contrast_plan$group1))
contrast_plan$group0 <- make.names(canonicalize_condition(contrast_plan$group0))
contrast_plan <- unique(contrast_plan)

required_groups <- unique(c(contrast_plan$group1, contrast_plan$group0))
missing_groups <- setdiff(required_groups, levels(coldata$condition))
if (length(missing_groups) > 0) {
  stop("Required contrast groups are missing from the condition levels: ",
       paste(missing_groups, collapse = ", "),
       ". Available levels: ", paste(levels(coldata$condition), collapse = ", "))
}

if (is.na(ref_level_opt)) {
  ref_level <- if ("Free" %in% levels(coldata$condition)) "Free" else levels(coldata$condition)[1]
} else {
  ref_level_opt <- make.names(canonicalize_condition(ref_level_opt))
  if (!ref_level_opt %in% levels(coldata$condition)) {
    stop("Requested ref_level '", ref_level_opt, "' not found in condition levels: ",
         paste(levels(coldata$condition), collapse = ", "))
  }
  ref_level <- ref_level_opt
}
coldata$condition <- relevel(coldata$condition, ref = ref_level)

# ==============================================================================
# DESeq2 dataset, filtering, and fit
# ==============================================================================
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

if (tolower(filter_mode) == "legacy") {
  keep <- rowSums(counts(dds) > 0) > 0
} else {
  keep <- keep_by_group_abundance(dds,
                                  min_count = min_count,
                                  min_reps = min_reps_within_group)
}
dds <- dds[keep, ]

dds <- DESeq(dds, parallel = TRUE)
coef_names <- resultsNames(dds)

# Optional LRT
if (run_lrt) {
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1, parallel = TRUE)
  res_lrt <- results(dds_lrt, alpha = results_alpha)
  write.csv(as.data.frame(res_lrt), file.path(output_dir, "LRT_condition_vs_null.csv"))
}

# ==============================================================================
# Transform for QC
# ==============================================================================
vsd_fun <- if (use_rlog) rlog else vst
vsd <- vsd_fun(dds, blind = FALSE)

# ==============================================================================
# QC plots
# ==============================================================================
# PCA
p_pca <- plotPCA(vsd, intgroup = "condition") + ggplot2::ggtitle("PCA — variance-stabilised counts")
ggplot2::ggsave(file.path(output_dir, "QC_PCA_vst_or_rlog.pdf"), p_pca, width = 6.5, height = 5)

# Size factors
pdf(file.path(output_dir, "QC_size_factors.pdf"))
barplot(sizeFactors(dds), las=2, ylab="sizeFactor", main="Size factors")
dev.off()

# Dispersion trend
pdf(file.path(output_dir, "QC_dispersion_trend.pdf")); plotDispEsts(dds); dev.off()

# Sample distance heatmap
mat <- as.matrix(dist(t(assay(vsd))))
pdf(file.path(output_dir, "QC_sample_distance_heatmap.pdf"))
pheatmap(mat, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean",
         main = "Sample–sample distance (VST/rlog)")
dev.off()

# Cook's distance
pdf(file.path(output_dir, "QC_cooks_distance_boxplot.pdf"))
boxplot(log10(assays(dds)[["cooks"]] + 1e-8), range=0, outline=FALSE, las=2,
        ylab="log10 Cook's distance", main="Influential observations")
dev.off()

# mean–SD plot
pdf(file.path(output_dir, "QC_meanSD_vst_or_rlog.pdf")); vsn::meanSdPlot(assay(vsd)); dev.off()

# p-value histogram on first available coefficient
first_coef <- grep("^condition", resultsNames(dds), value = TRUE)[1]
if (!is.na(first_coef)) {
  res_tmp <- results(dds, name = first_coef, alpha = results_alpha, parallel = TRUE)
  pdf(file.path(output_dir, "QC_pvalue_histogram.pdf"))
  hist(res_tmp$pvalue, breaks = 50, col = "grey",
       main = paste("P-value distribution —", first_coef), xlab = "p-value")
  dev.off()
}

# ==============================================================================
# Helpers for robust shrinkage, exporting contrasts, GSEA ID selection, and combined export
# ==============================================================================
find_coef_name <- function(dds_obj, grp1, grp0) {
  rn <- resultsNames(dds_obj)
  g1 <- gsub("[^A-Za-z0-9.]+", ".", grp1)
  g0 <- gsub("[^A-Za-z0-9.]+", ".", grp0)
  pat  <- paste0("^condition[_:]*", g1, ".*_vs_.*", g0, "$")
  hits <- grep(pat, rn, value = TRUE)
  if (length(hits) >= 1) return(hits[1])
  pat2 <- paste0("^condition_", make.names(grp1), "_vs_", make.names(grp0), "$")
  hits2 <- grep(pat2, rn, value = TRUE)
  if (length(hits2) >= 1) return(hits2[1])
  stop("Could not find coefficient for ", grp1, " vs ", grp0, " in resultsNames().")
}

shrink_contrast <- function(dds_obj, contrast_vec, shrink_method = "apeglm") {
  grp1 <- as.character(contrast_vec[2]); grp0 <- as.character(contrast_vec[3])
  if (tolower(shrink_method) == "ashr") {
    return(lfcShrink(dds_obj, contrast = contrast_vec, type = "ashr", parallel = TRUE))
  } else {
    try_coef <- try(find_coef_name(dds_obj, grp1, grp0), silent = TRUE)
    if (!inherits(try_coef, "try-error")) {
      return(lfcShrink(dds_obj, coef = try_coef, type = "apeglm", parallel = TRUE))
    } else {
      dds2 <- dds_obj
      dds2$condition <- relevel(dds2$condition, ref = grp0)
      dds2 <- DESeq(dds2, parallel = TRUE)
      coef2 <- find_coef_name(dds2, grp1, grp0)
      return(lfcShrink(dds2, coef = coef2, type = "apeglm", parallel = TRUE))
    }
  }
}

labelize <- function(a, b) gsub("[^A-Za-z0-9._-]+", "_", paste0(a, "_vs_", b))

pick_ids_for_go <- function(gann) {
  if ("Ensembl_GeneID" %in% colnames(gann) && any(nzchar(gann$Ensembl_GeneID))) {
    return(list(ids = gann$Ensembl_GeneID, keyType = "ENSEMBL"))
  }
  if ("NCBI_GeneID" %in% colnames(gann) && any(!is.na(gann$NCBI_GeneID))) {
    return(list(ids = as.character(gann$NCBI_GeneID), keyType = "ENTREZID"))
  }
  if ("SYMBOL" %in% colnames(gann) && any(nzchar(gann$SYMBOL))) {
    return(list(ids = gann$SYMBOL, keyType = "SYMBOL"))
  }
  return(NULL)
}

to_entrez <- function(vec_ids, keyType_in, OrgDb) {
  if (keyType_in == "ENTREZID") {
    ents <- unique(stats::na.omit(vec_ids))
    return(as.character(ents))
  }
  suppressWarnings({
    df <- tryCatch(
      clusterProfiler::bitr(vec_ids, fromType = keyType_in, toType = "ENTREZID", OrgDb = OrgDb),
      error = function(e) NULL
    )
  })
  if (is.null(df) || nrow(df) == 0) return(character(0))
  df <- df[!duplicated(df$ENTREZID), ]
  return(as.character(df$ENTREZID))
}

make_de_call <- function(log2fc, padj, group1, group0, alpha = 0.05) {
  out <- rep("not_sig", length(log2fc))
  sig <- !is.na(padj) & padj < alpha
  out[sig & log2fc > 0] <- paste0(group1, "_up")
  out[sig & log2fc < 0] <- paste0(group0, "_up")
  out
}

build_ssrnaseq_combined_output <- function(dds_obj, gene_annot_df, contrast_tables, contrast_plan_df, results_alpha = 0.05) {
  gene_ids <- rownames(dds_obj)

  annot_df <- gene_annot_df[gene_ids, , drop = FALSE]
  annot_df <- annot_df[, setdiff(colnames(annot_df), "Geneid_clean"), drop = FALSE]

  raw_counts <- as.data.frame(counts(dds_obj, normalized = FALSE))
  raw_counts <- raw_counts[gene_ids, , drop = FALSE]
  colnames(raw_counts) <- paste0("raw_", colnames(raw_counts))

  norm_counts <- as.data.frame(counts(dds_obj, normalized = TRUE))
  norm_counts <- norm_counts[gene_ids, , drop = FALSE]
  colnames(norm_counts) <- paste0("norm_", colnames(norm_counts))

  norm_mat <- counts(dds_obj, normalized = TRUE)
  mean_norm_list <- lapply(levels(dds_obj$condition), function(cond) {
    idx <- dds_obj$condition == cond
    vals <- rowMeans(norm_mat[, idx, drop = FALSE])
    out <- data.frame(vals, row.names = gene_ids, check.names = FALSE)
    colnames(out) <- paste0("mean_norm_", cond)
    out
  })
  mean_norm_df <- do.call(cbind, mean_norm_list)

  first_lbl <- labelize(contrast_plan_df$group1[1], contrast_plan_df$group0[1])
  if (is.null(contrast_tables[[first_lbl]])) {
    stop("Could not find the first contrast table needed to build the combined ssRNAseq export.")
  }
  base_mean <- contrast_tables[[first_lbl]][gene_ids, "baseMean", drop = TRUE]

  combined <- data.frame(
    Geneid_clean = gene_ids,
    annot_df,
    baseMean = base_mean,
    raw_counts,
    norm_counts,
    mean_norm_df,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  for (k in seq_len(nrow(contrast_plan_df))) {
    group1 <- contrast_plan_df$group1[k]
    group0 <- contrast_plan_df$group0[k]
    lbl <- labelize(group1, group0)
    out <- contrast_tables[[lbl]]
    if (is.null(out)) next

    out <- out[gene_ids, , drop = FALSE]
    stat_cols <- intersect(c("log2FC", "log2FC_shrunken", "lfcSE_shrunken", "stat", "pvalue", "padj"), colnames(out))
    block <- out[, stat_cols, drop = FALSE]
    colnames(block) <- paste0(lbl, "__", colnames(block))
    block[[paste0(lbl, "__DE_call")]] <- make_de_call(out$log2FC, out$padj, group1, group0,
                                                      alpha = results_alpha)

    combined <- cbind(combined, block)
  }

  rownames(combined) <- gene_ids
  combined
}

run_contrast <- function(dds_obj, group1, group0, do_gsea = TRUE, shrink_method = "apeglm",
                         results_alpha = 0.05) {
  lbl <- labelize(group1, group0)
  contrast_vec <- c("condition", group1, group0)

  # Raw Wald results
  res_raw <- results(dds_obj, contrast = contrast_vec, alpha = results_alpha, parallel = TRUE)
  res_df  <- as.data.frame(res_raw)

  # Shrunken LFCs
  res_shr <- tryCatch(
    shrink_contrast(dds_obj, contrast_vec, shrink_method = shrink_method),
    error = function(e) { message("Shrinkage failed (", lbl, "): ", e$message); NULL }
  )

  # Compose export table
  if (!is.null(res_shr)) {
    out <- data.frame(
      baseMean        = res_raw$baseMean,
      log2FC          = res_raw$log2FoldChange,
      log2FC_shrunken = res_shr$log2FoldChange,
      lfcSE_shrunken  = res_shr$lfcSE,
      stat            = res_raw$stat,
      pvalue          = res_raw$pvalue,
      padj            = res_raw$padj,
      row.names       = rownames(res_raw)
    )
  } else {
    out <- data.frame(
      baseMean = res_raw$baseMean,
      log2FC   = res_raw$log2FoldChange,
      stat     = res_raw$stat,
      pvalue   = res_raw$pvalue,
      padj     = res_raw$padj,
      row.names = rownames(res_raw)
    )
  }

  # Main export
  write.csv(out, file.path(output_dir, paste0("DESeq2_", lbl, ".csv")))

  # ---- Extra exports with alternative row labels (robust if SYMBOL absent) ----
  gann_sub <- gene_annot[rownames(out), , drop = FALSE]
  if (!"SYMBOL" %in% colnames(gann_sub)) gann_sub$SYMBOL <- NA_character_
  if (!"Ensembl_GeneID" %in% colnames(gann_sub)) gann_sub$Ensembl_GeneID <- NA_character_

  # Ensembl-labeled
  ens_lab <- gann_sub$Ensembl_GeneID
  ens_lab[is.na(ens_lab) | !nzchar(ens_lab)] <- rownames(out)   # fallback to Geneid_clean
  ens_lab <- make.unique(as.character(ens_lab))
  out_by_ens <- out
  rownames(out_by_ens) <- ens_lab
  write.csv(out_by_ens, file.path(output_dir, paste0("DESeq2_", lbl, "_by_ENS.csv")))

  # SYMBOL-labeled
  sym_lab <- gann_sub$SYMBOL
  sym_lab[is.na(sym_lab) | !nzchar(sym_lab)] <- rownames(out)   # fallback to Geneid_clean
  sym_lab <- make.unique(as.character(sym_lab))
  out_by_sym <- out
  rownames(out_by_sym) <- sym_lab
  write.csv(out_by_sym, file.path(output_dir, paste0("DESeq2_", lbl, "_by_SYMBOL.csv")))

  # Sig counts at the configured FDR cutoff
  sig <- !is.na(res_raw$padj) & res_raw$padj < results_alpha
  n_up   <- sum(sig & res_raw$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(sig & res_raw$log2FoldChange < 0, na.rm = TRUE)
  write.table(data.frame(up = n_up, down = n_down),
              file = file.path(output_dir, paste0("sig_gene_counts_", lbl, ".txt")),
              sep = "\t", row.names = FALSE)

  # MA plot
  pdf(file.path(output_dir, paste0("MA_", lbl, ".pdf"))); plotMA(res_raw, ylim = c(-5,5)); dev.off()

  # Volcano
  vdf <- transform(res_df,
                   neglog10padj = -log10(pmin(padj, 1)),
                   log2FC       = log2FoldChange)
  vdf <- vdf[is.finite(vdf$neglog10padj) & is.finite(vdf$log2FC), ]
  pdf(file.path(output_dir, paste0("volcano_", lbl, ".pdf")))
  plot(vdf$log2FC, vdf$neglog10padj, pch = 20,
       xlab = "log2 fold-change", ylab = "-log10(FDR)",
       main = paste0("Volcano — ", group1, " vs ", group0))
  abline(h = -log10(results_alpha), lty = 2)
  dev.off()

  # -----------------------------  GSEA  ---------------------------------------
  if (do_gsea && !is.null(OrgDb)) {
    ranked <- res_raw$stat
    names(ranked) <- rownames(res_raw)

    gann <- gene_annot[rownames(res_raw), , drop = FALSE]
    sel <- pick_ids_for_go(gann)

    if (!is.null(sel)) {
      ids_vec  <- sel$ids
      keyType  <- sel$keyType

      keep <- !is.na(ranked) & ranked != 0 & !is.na(ids_vec) & nzchar(as.character(ids_vec))
      ranked2 <- ranked[keep]
      ids2    <- as.character(ids_vec[keep])

      # collapse duplicates by max |stat|
      rk_tbl <- tapply(ranked2, ids2, function(x) x[which.max(abs(x))])
      ranked_named <- sort(unlist(rk_tbl), decreasing = TRUE)

      # GO GSEA
      gse_go <- tryCatch(
        clusterProfiler::gseGO(geneList = ranked_named, keyType = keyType, ont = "ALL",
                               OrgDb = OrgDb, nPerm = 10000,
                               minGSSize = 3, maxGSSize = 800,
                               pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE),
        error = function(e) { message("gseGO failed (", lbl, "): ", e$message); NULL }
      )
      if (!is.null(gse_go) && nrow(as.data.frame(gse_go)) > 0) {
        pdf(file.path(output_dir, paste0("GO_GSEA_dotplot_", lbl, ".pdf")))
        print(enrichplot::dotplot(gse_go, showCategory = 10, split = ".sign") + ggplot2::facet_grid(. ~ .sign))
        dev.off()
      }

      # KEGG GSEA (requires ENTREZ)
      map_df <- tryCatch(
        clusterProfiler::bitr(names(ranked_named), fromType = keyType, toType = "ENTREZID", OrgDb = OrgDb),
        error = function(e) NULL
      )
      if (!is.null(map_df) && nrow(map_df) > 0) {
        map_df <- map_df[!duplicated(map_df$ENTREZID), ]
        rk2 <- ranked_named[map_df[[keyType]]]
        names(rk2) <- as.character(map_df$ENTREZID)
        rk2 <- sort(rk2[!is.na(rk2)], decreasing = TRUE)

        kk2 <- tryCatch(
          clusterProfiler::gseKEGG(geneList = rk2, organism = kegg_org, keyType = "ncbi-geneid",
                                   minGSSize = 3, maxGSSize = 800,
                                   pvalueCutoff = 0.05, pAdjustMethod = "BH"),
          error = function(e) { message("gseKEGG failed (", lbl, "): ", e$message); NULL }
        )
        if (!is.null(kk2) && nrow(as.data.frame(kk2)) > 0) {
          pdf(file.path(output_dir, paste0("KEGG_GSEA_dotplot_", lbl, ".pdf")))
          print(enrichplot::dotplot(kk2, showCategory = 10, split = ".sign") + ggplot2::facet_grid(. ~ .sign))
          dev.off()
        }
      }
    } else {
      message("No usable gene IDs found for GSEA in gene_annot; skipping GSEA for ", lbl)
    }
  }

  out
}

# ==============================================================================
# Run only the requested contrasts
# ==============================================================================
contrast_results <- vector("list", nrow(contrast_plan))
names(contrast_results) <- vapply(
  seq_len(nrow(contrast_plan)),
  function(k) labelize(contrast_plan$group1[k], contrast_plan$group0[k]),
  character(1)
)

for (k in seq_len(nrow(contrast_plan))) {
  lbl <- labelize(contrast_plan$group1[k], contrast_plan$group0[k])
  contrast_results[[lbl]] <- run_contrast(dds,
                                          group1 = contrast_plan$group1[k],
                                          group0 = contrast_plan$group0[k],
                                          do_gsea = run_gsea,
                                          shrink_method = shrink_method,
                                          results_alpha = results_alpha)
}

combined_ssrnaseq <- build_ssrnaseq_combined_output(dds, gene_annot, contrast_results, contrast_plan,
                                                    results_alpha = results_alpha)
write.csv(combined_ssrnaseq, file.path(output_dir, ssrnaseq_output_file), row.names = FALSE)

# ==============================================================================
# Save normalized counts, contrast plan, and session info
# ==============================================================================
write.csv(counts(dds, normalized = TRUE), file.path(output_dir, "normalized_counts.csv"))
write.csv(as.data.frame(colData(dds)),          file.path(output_dir, "coldata_used.csv"))
write.csv(as.data.frame(gene_annot[rownames(dds), , drop = FALSE]),
          file.path(output_dir, "gene_annotation_used.csv"))
write.csv(contrast_plan, file.path(output_dir, "contrast_plan_used.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

# Console summary
cat("\n=== Run summary ===\n")
cat("Samples:", ncol(dds), " | Genes (after filter):", nrow(dds), "\n")
cat("Condition levels:", paste(levels(dds$condition), collapse = ", "), "\n")
cat("Reference level:", ref_level, "\n")
cat("Abundance filter:", if (tolower(filter_mode) == "legacy") {
      "legacy (>0 reads in at least 1 sample)"
    } else {
      paste0(">=", min_count, " reads in >=", min_reps_within_group, " samples within at least 1 condition")
    }, "\n")
cat("Results alpha:", results_alpha, "\n")
cat("Contrasts:", paste(paste(contrast_plan$group1, "vs", contrast_plan$group0), collapse = "; "), "\n")
cat("Combined ssRNAseq export:", file.path(output_dir, ssrnaseq_output_file), "\n")
cat("Outputs in:", normalizePath(output_dir), "\n")


