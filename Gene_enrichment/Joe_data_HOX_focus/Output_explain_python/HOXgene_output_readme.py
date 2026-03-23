
# Create a comprehensive README describing each output file.
import textwrap, os, json, datetime

readme_path = "/mnt/data/README_HOX_GDF11_outputs.md"

today = datetime.date.today().isoformat()

content = f"""# HOX/GDF11 Post-hoc Analysis — Output Guide
_Date generated: {today}_

This README explains **every file** output by `posthoc_HOX_GDF11_analysis.R`, what it contains, and how to interpret it.  
Files are grouped by **tables (CSV)** and **figures (PDF/PNG/SVG)**. PDF/PNG/SVG are identical visual content in different formats.

> Tip: PDFs (vector) are best for print; SVGs (vector) are ideal for Illustrator/Inkscape; PNGs are 600 dpi for slides.

---

## Tables (CSV)

### 1) Gene/Expression Mappings & Long Data
- **`HOX_GDF11_gene_map.csv`**  
  Mapping for each row used in HOX/GDF11 analyses. Columns:
  - `Geneid`: row ID used in your DESeq2 matrices.
  - `SYMBOL_display`: preferred human-readable label (SYMBOL if available, otherwise fallback to Geneid).
  - `SYMBOL_raw`: gene symbol as provided by annotation (may be missing).
  - `HOX_cluster`: \"A\"/\"B\"/\"C\"/\"D\"/\"Other\" (based on `SYMBOL_raw`).
  - `HOX_num`: HOX paralog group (1–13) if detected.
  - `Ensembl_GeneID`, `NCBI_GeneID`: IDs if present.

- **`HOX_GDF11_expression_long.csv`**  
  Long/tidy expression table for HOX (A–D) + GDF11 across all samples. Columns:
  - `Geneid`, `SYMBOL_display`, `HOX_cluster`, `HOX_num`
  - `sample`: sample ID (matches `coldata_used.csv`).
  - `norm_count`: DESeq2 **normalized** count (not log-transformed).
  - `log2_norm`: `log2(norm_count + 1)` (used in most plots).
  - `condition`: sample condition from your metadata.

### 2) DE Results (per contrast) filtered to HOX/GDF11
- **`HOX_GDF11_all_contrasts_tidy.csv`**  
  Row-bind of all DE result files from `deseq2_outputs/` filtered to HOX A–D and GDF11. Includes (where available):
  - `baseMean`, `log2FC`, `log2FC_shrunken`, `lfcSE_shrunken`, `stat`, `pvalue`, `padj`, `contrast`, `Geneid`, `SYMBOL_display`.

- **`HOX_GDF11_*.csv`** (e.g.,  
  `HOX_GDF11_anterior.sacral_vs_posterior.sacral.csv`,  
  `HOX_GDF11_posterior.sacral_vs_sacralized.caudal.csv`, and their `_by_SYMBOL` / `_by_ENS` variants)  
  These are **contrast-specific** subsets (HOX/GDF11 only) mirroring the parent DESeq2 tables. Variants:
  - `_by_SYMBOL`: rows relabeled to symbols (fallback to Geneid if missing).
  - `_by_ENS`: rows relabeled to Ensembl IDs.
  Each table contains fold changes and multiple-testing columns where available.

### 3) Sacral-oriented Summary Scores
- **`HOX_AP_scores_per_sample.csv`**  
  Per-sample summary used to relate to **sacral somites**:
  - `posterior_score`: mean `log2_norm` across HOX PG10–13 (A–D), a proxy for sacral signal.
  - `anterior_score`: mean `log2_norm` across HOX PG1–3.
  - `AP_index`: `posterior_score − anterior_score` (positive ⇒ posterior-enriched).
  - `condition`

- **GDF11 summary tables**  
  - `GDF11_log2norm_per_sample.csv`: per-sample GDF11 `log2_norm` with condition.  
  - `GDF11_prevalence_by_condition.csv`: fraction of samples with GDF11 `log2_norm > 1` per condition (threshold adjustable in script).  
  - `GDF11_summary_overall.csv`: overall mean, SD, and sample count for GDF11.

---

## Figures

> Unless specified, heatmaps are **z-scored per gene** (`log2(norm+1)` then row-scaled, capped at ±3 SD) and include column annotation by `condition`.

### 1) Global and Cluster-Specific Heatmaps
- **`HOX_GDF11_heatmap_ALL_clusters.(pdf|png|svg)`**  
  Heatmap of HOX A–D genes plus GDF11 across all samples.  
  - Rows ordered by cluster (A→D) and paralog group (1→13).  
  - **Row gaps** show cluster boundaries.  
  - Use this to verify A→P gradients and GDF11 co-expression at the posterior.

- **`HOX[ABCD]_heatmap.(pdf|png|svg)`** (e.g., `HOXA_heatmap.pdf`)  
  One heatmap per HOX cluster (A/B/C/D), rows ordered by paralog number (1→13).  
  - Highlights cluster-specific A→P progression independent of other clusters.

- **`HOX_GDF11_heatmap_samples_sorted_by_posteriorScore.(pdf|png|svg)`**  
  Same gene set as the global heatmap, but **samples are ordered** by the posterior HOX score (PG10–13).  
  - Useful to see whether posterior-rich samples cluster together; aligns with sacral expectations.

### 2) Schematic Grids (A–D × PG1–13)
- **`HOX_schematic_grid.(pdf|png|svg)`**  
  4×13 tile plot showing **mean** `log2(norm+1)` per HOX cluster (rows: HOXA at top) and paralog group (columns: PG1→PG13).  
  - Each tile labeled (e.g., A10, B5).  
  - Readout resembles in situ atlases: posterior (PG9–13) brightness suggests sacral enrichment.

- **Category splits**  
  - `HOX_schematic_grid_anterior.(pdf|png|svg)` → PG1–3 only.  
  - `HOX_schematic_grid_thoracic.(pdf|png|svg)` → PG4–8.  
  - `HOX_schematic_grid_posterior_sacral.(pdf|png|svg)` → PG9–13.  
  These zoom into A/P segments to compare categories with prior in situ patterns.

### 3) Condition-Level and Per-Cluster A→P Profiles
- **`HOX_AP_profiles_by_condition.(pdf|png|svg)`**  
  Line plots of **mean** `log2(norm+1)` vs paralog group (PG1→PG13) **per condition**, faceted by cluster (A/B/C/D).  
  - Expect rising curves toward PG9–13 for sacral-like conditions.

- **`HOX[A-D]_AP_profile_by_condition.(pdf|png|svg)`** (e.g., `HOXA_AP_profile_by_condition.pdf`)  
  Same as above but **one cluster per file**.  
  - Helps inspect cluster-specific behavior (e.g., whether HOXC/D are more posterior-weighted).

### 4) Sacral Score & A→P Index Summaries
- **`posterior_HOX_score_by_condition.(pdf|png|svg)`**  
  Box/jitter plot of the **posterior HOX score** (PG10–13) by condition.  
  - Higher values imply posteriorized/sacral identity.

- **`AP_index_by_condition.(pdf|png|svg)`**  
  Box/jitter plot of **A→P index** = (posterior − anterior).  
  - Values > 0 indicate posterior bias; compare against in situ expectations for sacral somites.

### 5) Correlations & Multivariate Structure
- **`Corr_genes_PG9_13_plus_GDF11.(pdf|png|svg)`**  
  Gene–gene Pearson correlation heatmap among **HOX PG9–13 (A–D)** and **GDF11**.  
  - Strong positive blocks across clusters suggest coordinated posterior programs; GDF11 correlation indicates coupling.

- **`HOX_GDF11_sample_correlation.(pdf|png|svg)`**  
  Sample–sample Pearson correlation based on HOX+GDF11 log2 data (variable genes only).  
  - Good for checking whether sacral/caudal samples cluster together.

- **`HOX_GDF11_PCA.(pdf|png|svg)`**  
  PCA on HOX+GDF11 (row-scaled).  
  - PC1/PC2 often track A→P state; color by condition reveals separation along posteriorization.

### 6) GDF11-specific Checks
- **`GDF11_by_condition.(pdf|png|svg)`**  
  Box/jitter plot of **GDF11** `log2(norm+1)` by condition; dashed line shows prevalence threshold used in the summary table.

- **`GDF11_vs_posteriorHOX_scatter.(pdf|png|svg)`**  
  Scatter of per-sample **GDF11** vs mean **posterior HOX (PG9–13)** with linear fit and Pearson *r*.  
  - Tests whether GDF11 tracks posteriorization expected from in situ.

### 7) LFC/Lollipop (per contrast)
- **`HOX_GDF11_LFC_<contrast>.(pdf|png|svg)`** (and `_by_SYMBOL`, `_by_ENS` variants)  
  Lollipop plots of HOX/GDF11 **log2 fold changes** for each DE contrast (shrunken if available).  
  - Rightward lollipops = up in first group of the contrast label; leftward = down.

---

## How to read these in the context of sacral somites

- Sacral identity is typically reflected by **elevated PG10–13** across HOX clusters; look for:  
  1) Bright PG9–13 tiles in schematic grids (especially **C/D** clusters).  
  2) Higher **posterior HOX score** and **AP index** in sacral-like conditions.  
  3) Positive correlations among HOX PG9–13 and a positive association with **GDF11**.
- **GDF11 presence**: Check `GDF11_by_condition` and `GDF11_*` tables; presence and correlation with posterior HOX supports classical in situ patterns.

---

## File naming conventions

- Most plots are emitted in **three formats**: `.pdf`, `.png`, `.svg` with identical content.  
- DE contrast files mirror the base names from your `deseq2_outputs/` directory.  
- “`_by_SYMBOL`” and “`_by_ENS`” tables relabel rows for easier human reading.

---

## Provenance

- Inputs pulled from the **previous DESeq2 run**:  
  `normalized_counts.csv`, `coldata_used.csv`, `gene_annotation_used.csv`.  
- Normalization displayed in plots is **log2(normalized + 1)**; heatmaps use **per-gene z-scoring** (±3 cap).

---

If you’d like this README exported as **PDF** for sharing, say the word and I’ll produce a print-styled PDF version.
"""

os.makedirs(os.path.dirname(readme_path), exist_ok=True)
with open(readme_path, "w") as f:
    f.write(content)

readme_path
