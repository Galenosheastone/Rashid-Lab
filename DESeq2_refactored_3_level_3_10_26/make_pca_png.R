library(DESeq2)
library(ggplot2)

wd <- "/Users/galen2/Documents/Documents_Folder/Rashid_Lab/Code_MAIN/R_Scripts/DESeq2_refactored_3_level_3_10_26"
setwd(wd)

cd <- read.csv("deseq2_outputs/coldata_used.csv", row.names = 1, check.names = FALSE)
cd$condition <- factor(cd$condition, levels = c("Free", "Pygo", "Sacral"))

Expression <- read.csv("Sacral24_realigned_GRCg7b_counts_matchedIDs_MAIN_forDESeq2.csv",
                       check.names = FALSE)
anno_cols <- intersect(colnames(Expression),
  c("SYMBOL","NCBI_GeneID","Ensembl_GeneID","gene_biotype","Chr","Start","End","Strand","Length"))
sample_cols <- setdiff(colnames(Expression), c("Geneid", anno_cols))

cts <- as.matrix(Expression[, sample_cols])
rownames(cts) <- make.unique(trimws(as.character(Expression$Geneid)))
storage.mode(cts) <- "integer"
cts <- cts[, rownames(cd)]

dds <- DESeqDataSetFromMatrix(countData = cts, colData = cd, design = ~ condition)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)

p <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA - variance-stabilised counts") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("deseq2_outputs/QC_PCA_vst_or_rlog.png", p, width = 7, height = 5.5, dpi = 150)
cat("Saved: deseq2_outputs/QC_PCA_vst_or_rlog.png\n")
