# Rashid Lab — RNA-seq Analysis Pipeline

This repository contains R and Python scripts used in the Rashid Lab for bulk and single-cell RNA-seq analysis. Our work focuses on comparing gene expression across tissues and developmental states, identifying differentially expressed genes, and translating those changes into pathway- and mechanism-level interpretation.

---

## Repository Structure

| Folder | Description |
|--------|-------------|
| `DESeq2_updated` | Current DESeq2 differential expression pipeline with the latest parameter tuning and output structure |
| `DESeq2_refactored_3_level_3_10_26` | Refactored DESeq2 pipeline supporting multi-level (3-group) contrasts |
| `DESeq2_Clemson_output` | DESeq2 analysis outputs for the Clemson dataset |
| `Gene_enrichment` | Gene ontology (GO) and pathway enrichment analysis scripts |
| `RNAseq_pathway_pipeline` | End-to-end pipeline from count matrices to pathway-level results |
| `ssRNAseq_analysis` | Single-cell RNA-seq analysis scripts and outputs |
| `HOX_focus` | Scripts focused on HOX gene expression analysis |
| `cGAST:STING_codex` | cGAS-STING pathway-focused analysis and visualization |
| `3D_network map` | Scripts for generating 3D gene network visualizations |
| `cnet_plot_tool` | Concept network (cnet) plot generation tool |
| `sPLSDA` | Sparse Partial Least Squares Discriminant Analysis scripts |
| `gene_validation_tool` | Tools for validating gene lists and cross-referencing annotations |
| `ggallus_conversion_tool` | Gene ID conversion utilities for *Gallus gallus* (chicken) genome |
| `Project_specific_analysis` | One-off or project-specific analysis scripts |
| `Claude_RNAseq` | Experimental scripts developed with AI assistance |

---

## Dependencies

Most scripts require **R (>= 4.0)** with the following commonly used packages:

- DESeq2
- clusterProfiler
- ggplot2
- enrichplot
- mixOmics (for sPLSDA)
- dplyr, tidyr, readr (tidyverse)

Python scripts require Python 3 with pandas, numpy, matplotlib.

To install Bioconductor packages in R:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "enrichplot"))

---

## General Usage

1. Prepare your count matrix - Scripts expect a gene x sample count matrix (CSV or TSV) and a corresponding sample metadata table.
2. Run DESeq2 - Start with the DESeq2_updated folder for the most current differential expression workflow.
3. Pathway enrichment - Feed DESeq2 output into Gene_enrichment or RNAseq_pathway_pipeline for GO/KEGG enrichment analysis.
4. Visualization - Use cnet_plot_tool and 3D_network map for network-level visualization of results.

Each folder contains its own scripts and, where applicable, example input/output files.

---

## Contact

For questions about this repository, please contact Galen O'Shea-Stone @ galenoshea@gmail.com
