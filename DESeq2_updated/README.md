# DESeq2 Bulk RNA-seq Pipeline

This project refactors the original monolithic R Markdown analysis into a modular DESeq2-based bulk RNA-seq workflow. It preserves the same core scientific intent: count import, metadata handling, DESeq2 fitting, optional LRT, QC plots, targeted contrasts, LFC shrinkage, optional GO/KEGG GSEA, and downstream exports.

## Project layout

```text
.
├── run_deseq2_analysis.R
├── config/
│   └── analysis_config.yaml
├── R/
│   ├── io.R
│   ├── validation.R
│   ├── cleaning.R
│   ├── deseq_model.R
│   ├── qc_plots.R
│   ├── contrasts.R
│   ├── gsea.R
│   └── exports.R
├── README.md
└── outputs/
```

## What the pipeline does

The pipeline:

- loads a bulk RNA-seq count matrix and metadata
- cleans sample IDs and gene IDs with explicit logging
- validates inputs before model fitting
- fits a DESeq2 Wald model with design `~ condition`
- optionally fits a separate LRT model against `~ 1`
- generates standard QC plots from `vst` or `rlog`
- runs requested pairwise contrasts
- performs LFC shrinkage with `apeglm` or `ashr`
- optionally runs GO and KEGG GSEA
- exports normalized counts, metadata, annotation, per-contrast result tables, and a combined gene-level export

## Expected counts file format

The counts file must be CSV or TSV and must contain a `Geneid` column.

Expected structure:

- one row per gene or feature
- one `Geneid` column
- zero or more annotation columns such as `SYMBOL`, `Ensembl_GeneID`, `NCBI_GeneID`, `gene_biotype`, `Chr`, `Start`, `End`, `Strand`, `Length`
- one numeric sample column per library

Notes:

- sample columns must be numeric
- non-numeric sample candidates are dropped and listed in the repair log
- blank `Geneid` values are repaired explicitly and logged
- if `Geneid` is blank, the pipeline first tries `Ensembl_GeneID`
- if both are blank, the pipeline generates `UNLABELED_000001` style IDs

## Expected metadata file format

The preferred metadata format is a CSV or TSV with sample IDs as row names and at least one column named `condition`.

Example:

```text
sample_id,condition,batch
Sacral_1,Sacral,1
Sacral_2,Sacral,1
Free_1,Free,1
Pygo_1,Pygo,2
```

Two metadata layouts are supported:

1. First column contains sample IDs and is read as row names.
2. A regular table with an explicit `sample_id` column and a `condition` column.

Required metadata columns:

- `condition`

Additional metadata columns are preserved and written back out in `coldata_used.csv`.

## How condition labels are handled

Condition normalization is intentionally conservative.

- Exact case-insensitive matches for `Sacral`, `Free`, and `Pygo` are accepted.
- If `allow_condition_aliases: true`, limited numeric suffix aliases such as `Sacral_1` or `Pygo-2` can be canonicalized after normalization.
- Risky guesses from the original script, such as mapping `exp` to `Sacral`, have been removed.
- If a condition label cannot be mapped confidently, the pipeline stops instead of guessing.

Metadata inference from sample names is optional and discouraged.

- Set `infer_metadata_from_sample_names: true` only when you do not have a metadata file.
- If inference is enabled, the pipeline uses the leading sample-name token before `.`, `_`, or `-`.
- Any inferred or canonicalized labels are recorded in `input_repair_log.csv`.

## Configuration

Edit [`config/analysis_config.yaml`](/Users/galen2/Documents/Documents_Folder/Rashid_Lab/Code_MAIN/R_Scripts/DESeq2_updated/config/analysis_config.yaml) before running.

Key fields:

- `output_dir`: output directory
- `counts_file`: counts matrix path
- `metadata_file`: metadata table path
- `ref_level`: reference level for DESeq2; `null` auto-picks `Free` when present
- `filter_mode`: `legacy` or `strict`
- `min_count`, `min_reps_within_group`: strict filter parameters
- `use_rlog`: use `rlog` instead of `vst` for QC
- `run_lrt`: fit `test = "LRT"` against `~ 1`
- `run_gsea`: enable GO and KEGG GSEA
- `OrgDb_pkg`: species-specific OrgDb package name
- `kegg_org`: KEGG organism code
- `shrink_method`: `apeglm` or `ashr`
- `n_workers`: BiocParallel worker count
- `infer_metadata_from_sample_names`: optional metadata fallback
- `allow_condition_aliases`: optional conservative alias support
- `contrast_plan`: list of requested pairwise contrasts

## How to configure contrasts

Use a YAML list of `group1` / `group0` entries:

```yaml
contrast_plan:
  - group1: Sacral
    group0: Free
  - group1: Pygo
    group0: Free
  - group1: Sacral
    group0: Pygo
```

These labels are canonicalized with the same conservative logic used for metadata. If a requested group is missing from the cleaned metadata `condition` levels, the run stops before fitting.

## How to run

From the project root:

```bash
Rscript run_deseq2_analysis.R
```

Or pass an explicit config path:

```bash
Rscript run_deseq2_analysis.R config/analysis_config.yaml
```

## Required packages

Do not install packages inside the pipeline. Install them beforehand in your R environment.

Core packages used by the pipeline:

- `yaml`
- `DESeq2`
- `BiocParallel`
- `ggplot2`
- `pheatmap`
- `vsn`
- `apeglm` or `ashr`

Additional packages required when `run_gsea: true`:

- `AnnotationDbi`
- `clusterProfiler`
- `enrichplot`
- `pathview`
- the configured OrgDb package such as `org.Gg.eg.db`

## Output files

Main run-level outputs written to `output_dir`:

- `normalized_counts.csv`
- `coldata_used.csv`
- `gene_annotation_used.csv`
- `contrast_plan_used.csv`
- `sessionInfo.txt`
- `input_repair_log.csv`
- `DESeq2_combined_contrast_export.csv`
- `LRT_condition_vs_null.csv` if `run_lrt: true`

QC outputs:

- `QC_PCA_vst_or_rlog.pdf`
- `QC_size_factors.pdf`
- `QC_dispersion_trend.pdf`
- `QC_sample_distance_heatmap.pdf`
- `QC_cooks_distance_boxplot.pdf`
- `QC_meanSD_vst_or_rlog.pdf`
- `QC_pvalue_histogram.pdf`

Per-contrast outputs for each requested comparison:

- `DESeq2_<group1>_vs_<group0>.csv`
- `DESeq2_<group1>_vs_<group0>_by_ENS.csv`
- `DESeq2_<group1>_vs_<group0>_by_SYMBOL.csv`
- `sig_gene_counts_<group1>_vs_<group0>.txt`
- `MA_<group1>_vs_<group0>.pdf`
- `volcano_<group1>_vs_<group0>.pdf`

Optional GSEA outputs:

- `GO_GSEA_<contrast>.csv`
- `GO_GSEA_dotplot_<contrast>.pdf`
- `KEGG_GSEA_<contrast>.csv`
- `KEGG_GSEA_dotplot_<contrast>.pdf`
- `GSEA_status_<contrast>.txt`

## Combined export contents

`DESeq2_combined_contrast_export.csv` contains:

- `Geneid_clean`
- annotation columns carried through from the counts file
- original raw counts for all retained genes
- normalized counts for all retained genes
- mean normalized counts per condition
- per-contrast statistics with names like `<contrast>__log2FC`
- per-contrast `DE_call` values such as `Sacral_up`, `Free_up`, or `not_sig`

## Behavior changes from original script

Compared with the monolithic R Markdown version:

- runtime package installation has been removed
- risky condition guessing such as `exp -> Sacral` has been removed
- metadata inference is now opt-in instead of an automatic fallback
- all input repairs are written to `input_repair_log.csv`
- validation is stricter and stops earlier on ambiguous inputs
- the combined export is renamed from `DESeq2_full_combined_ssRNAseq.csv` to `DESeq2_combined_contrast_export.csv`
- `apeglm` refits are still supported, but the refit is now explicit and messaged
- GSEA is isolated so annotation issues do not crash the full pipeline

## Important caveats

- This is a bulk RNA-seq DESeq2 pipeline, not single-cell RNA-seq.
- Logged auto-cleaning is for traceability, not a substitute for proper upstream preprocessing.
- Metadata inference from sample names should be treated as a fallback, not a default operating mode.
- If your metadata uses condition names outside `Sacral`, `Free`, and `Pygo`, update the labels upstream or extend the canonicalization rules deliberately.
