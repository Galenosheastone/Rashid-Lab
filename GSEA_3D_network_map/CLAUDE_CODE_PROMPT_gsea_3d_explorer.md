# GSEA 3D Network Explorer — Build Prompt for Claude Code

## Project Summary

Build a new Python script called `build_gsea_3d_network.py` that takes a DESeq2 results CSV as input, runs a GSEA-style preranked enrichment analysis against gene set databases, and produces a self-contained interactive 3D force-directed network HTML file — similar in spirit and UI quality to the existing `build_3d_network.py` in this directory, but designed as an **exploratory** tool rather than a curated-pathway tool.

The existing `build_3d_network.py` and `pathway_network_3d.html` are your reference implementations for the 3D visualization, HTML template structure, control panel UI, and general architecture. Study them carefully — the new tool should feel like a sibling of the existing one, reusing the same visual language (dark theme, Three.js + 3d-force-graph, control panel on left, detail panel on right, same tooltip style, same export features).

---

## Organism & Data Context

- **Organism**: Gallus gallus (chicken). Gene symbols in the DESeq2 CSV are mostly HGNC-style symbols (e.g., TNF, CASP8, RUNX2) with some chicken-specific identifiers (LOC######). There are also some Excel date-mangled gene names (e.g., "1-Mar", "2-Sep") which should be handled gracefully (skip or warn).
- **Input format**: DESeq2 results CSV with columns: `""` (unnamed first col = gene symbol), `baseMean`, `log2FC`, `log2FC_shrunken`, `lfcSE_shrunken`, `stat`, `pvalue`, `padj`. See `EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv` for the exact format (~14,000 rows).
- **Ranking metric for GSEA prerank**: Use `stat` column (the Wald statistic) as the ranking metric. This is standard for DESeq2 → GSEA prerank. Fall back to `log2FC` if `stat` is missing/NA.

---

## Architecture Overview

The script should follow this pipeline:

```
DESeq2 CSV
    ↓
[1] Read & rank genes (by stat column)
    ↓
[2] Load gene set database (GMT file — bundled or user-supplied)
    ↓
[3] Run GSEA prerank (using gseapy)
    ↓
[4] Filter to significant gene sets (FDR < cutoff)
    ↓
[5] Cluster related gene sets (by gene overlap, Jaccard similarity)
    ↓
[6] Build 3D network: gene-set nodes + gene nodes + edges
    ↓
[7] Overlay DE data (log2FC, padj) onto gene nodes
    ↓
[8] Generate self-contained HTML (same approach as build_3d_network.py)
```

---

## Detailed Requirements

### Step 1: Read DESeq2 CSV

- Reuse the CSV reading logic from `build_3d_network.py` (`read_deseq2` function) but also extract the `stat` column.
- Build a ranked gene list: `{gene_symbol: stat_value}`, sorted descending by stat.
- Filter out genes with NA stat values.
- Warn (don't crash) on suspected Excel date-mangled gene names (regex: `^\d+-[A-Z][a-z]{2}$`).

### Step 2: Gene Set Database

- Use `gseapy` to handle GMT files. The script should support a `--gene-sets` CLI argument that accepts either:
  - A path to a local `.gmt` file
  - A shorthand name for a built-in MSigDB collection (e.g., `"H"` for Hallmark, `"C2_CP_KEGG"` for KEGG, `"C5_GO_BP"` for GO Biological Process)
- **Default**: Use MSigDB Hallmark gene sets (`H`) as the default since these are well-curated and produce a manageable number of results.
- **IMPORTANT chicken gene compatibility note**: MSigDB gene sets use human gene symbols. Most chicken genes in this dataset use the same symbols (chicken orthologs share HGNC symbols). This is an acceptable approximation — the user understands this mapping isn't perfect. Do NOT try to do ortholog conversion. Just do case-insensitive matching of gene symbols against gene set members.
- `gseapy.prerank()` can download MSigDB GMT files automatically if you pass the collection name to `gene_sets` parameter.

### Step 3: Run GSEA Prerank

- Use `gseapy.prerank()` with:
  - `rnk`: the ranked gene list (pandas Series or dict, gene symbol → stat)
  - `gene_sets`: GMT file path or MSigDB collection name
  - `outdir`: a temp directory (use `tempfile.mkdtemp()`)
  - `permutation_num`: 1000 (default, allow CLI override with `--permutations`)
  - `seed`: 42 for reproducibility
  - `min_size`: 15 (gene sets smaller than this are skipped)
  - `max_size`: 500 (gene sets larger than this are skipped)
  - `no_plot`: True (we don't need the individual GSEA plots)
- Capture the results dataframe from `gseapy.prerank()` — it contains: `Term`, `ES`, `NES`, `NOM p-val`, `FDR q-val`, `FWER p-val`, `Lead_genes`, etc.

### Step 4: Filter Significant Gene Sets

- Filter to gene sets with `FDR q-val < fdr_cutoff` (default 0.25, standard GSEA threshold; allow CLI override with `--fdr-cutoff`).
- If fewer than 3 gene sets pass, relax to the top 20 by |NES| and warn the user.
- If more than 40 gene sets pass, take the top 40 by |NES| and note this in the output.
- Print a summary: number of gene sets tested, number significant, top 10 by NES.

### Step 5: Cluster Gene Sets into Super-Groups

This is the key step that replaces the hand-curated pathway Z-layers. Use Jaccard similarity of gene membership to cluster related gene sets:

- For each pair of significant gene sets, compute Jaccard index of their lead-edge gene lists.
- Build a similarity matrix.
- Use agglomerative clustering (scipy `fcluster` with a distance threshold, or a simple approach: cut the dendrogram to get 4–10 clusters). Target ~5–8 clusters for a clean 3D layout.
- Assign each cluster a Z-layer in 3D space (spaced 150 units apart, same as the existing tool).
- Assign each cluster a color from a curated palette (see below).
- Name each cluster after its most significant (lowest FDR) member gene set — this becomes the "pathway layer" label.

**Cluster color palette** (expand as needed):
```python
CLUSTER_COLORS = [
    '#6A5ACD',  # slate blue
    '#A33A3A',  # dark red
    '#2F7F56',  # forest green
    '#FF8C00',  # dark orange
    '#1F4E79',  # navy
    '#1E90FF',  # dodger blue
    '#8B6914',  # goldenrod
    '#CC5599',  # mauve
    '#2EAAAA',  # teal
    '#8855CC',  # purple
]
```

### Step 6: Build the 3D Network

The network has two types of nodes and two types of edges:

**Node types:**
1. **Gene-set nodes** (the enriched pathways/terms): Larger, diamond-shaped or distinctly styled. Positioned at the center of their cluster's Z-layer. Labeled with the gene set name (truncated to ~40 chars for display). Store: `term_name`, `NES`, `FDR`, `cluster_id`, `cluster_name`, `es_direction` (up/down based on NES sign).
2. **Gene nodes** (individual genes that are in the lead-edge of significant gene sets): Sized/colored by DE data (same scheme as existing tool — color by log2FC blue→white→red, size by significance). Store: `symbol`, `log2FC`, `padj`, `node_status`, `member_of` (list of gene set names).

**Edge types:**
1. **Gene-set ↔ Gene edges**: Connect each gene-set node to its lead-edge genes. These are the primary structural edges.
2. **Gene-set ↔ Gene-set edges** (optional but valuable): Connect gene sets that share significant gene overlap (Jaccard > 0.1). Style these as dashed/faded to distinguish from membership edges.

**3D Layout strategy:**
- Gene-set nodes: placed at their cluster's Z-layer, spread in X/Y within the layer.
- Gene nodes: positioned near the centroid of all the gene-set nodes they belong to (averaged X/Y/Z). Genes belonging to multiple gene sets will naturally bridge between clusters — these are the equivalent of "hub nodes" in the existing tool.
- Use 3d-force-graph's force simulation to refine positions, but pin gene-set nodes to their Z-layers (same technique as the existing tool uses `z_target`).

### Step 7: Overlay DE Data

- For each gene node, look up log2FC and padj from the DESeq2 data.
- Use the same node_status logic as the existing tool: `"Present; significant"`, `"Present; not significant"`, `"Missing gene"`.
- Compute `max_lfc` for color scale (same as existing tool).

### Step 8: Generate HTML

- Use the same approach as `build_3d_network.py`: embed a giant `HTML_TEMPLATE` string with `{placeholder}` substitution for the DATA JSON.
- **Reuse as much of the existing HTML/CSS/JS as possible** from `build_3d_network.py`. The key differences in the UI:

**Control panel changes:**
- "Pathways" section → "Gene Set Clusters" section: checkboxes for each cluster (named after representative gene set), colored dots.
- Add a "Gene Sets" subsection: expandable list of all significant gene sets, grouped by cluster. Clicking a gene set highlights it and its member genes.
- Filters section: keep "Significant only" and "Hub nodes only" (hubs = genes in 2+ gene sets). Add "Upregulated sets only" / "Downregulated sets only" toggle (filter by NES sign).
- Add an NES color bar to the legend (for gene-set nodes).

**Tooltip changes:**
- Gene nodes: same as existing (gene symbol, log2FC, padj, pathways → gene sets).
- Gene-set nodes: show term name, NES, FDR q-val, number of lead-edge genes, cluster name.

**Detail panel changes:**
- When clicking a gene-set node: show full term name, NES, ES, FDR, FWER, lead-edge gene list (with their log2FC), link-outs (if KEGG/GO, construct URL).
- When clicking a gene node: show symbol, log2FC, padj, list of gene sets it belongs to (with their NES).

**Keep these features from existing tool (copy/adapt):**
- Search box (search genes AND gene set names)
- Dark/light mode toggle
- Colorblind mode
- Layer labels in 3D space
- Screenshot export (PNG)
- CSV export (nodes and edges)
- Camera reset, reheat simulation
- Keyboard shortcuts
- Responsive hamburger menu
- Loading spinner

**Legend updates:**
- log2FC color bar (same as existing)
- NES color bar (for gene-set nodes): negative NES = blue, positive NES = red
- Node type key: ◆ Gene Set, ● Significant Gene, ○ Non-significant Gene
- Cluster color key (dynamic, based on discovered clusters)

---

## CLI Interface

```
python3 build_gsea_3d_network.py \
  --deseq2-csv DESeq2_results.csv \
  --gene-sets H \
  --output gsea_network_3d.html \
  --title "Free vs Pygo — GSEA 3D Explorer" \
  --fdr-cutoff 0.25 \
  --p-cutoff 0.05 \
  --permutations 1000 \
  --max-sets 40 \
  --min-size 15 \
  --max-size 500
```

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--deseq2-csv` | Yes | — | Path to DESeq2 results CSV |
| `--gene-sets` | No | `H` | MSigDB collection name or path to .gmt file |
| `--output` | No | `gsea_network_3d.html` | Output HTML path |
| `--title` | No | `GSEA 3D Network Explorer` | Title in browser |
| `--fdr-cutoff` | No | `0.25` | FDR threshold for significant gene sets |
| `--p-cutoff` | No | `0.05` | padj threshold for gene significance coloring |
| `--permutations` | No | `1000` | Number of GSEA permutations |
| `--max-sets` | No | `40` | Maximum gene sets to display |
| `--min-size` | No | `15` | Minimum gene set size |
| `--max-size` | No | `500` | Maximum gene set size |

---

## Dependencies

- `gseapy` (pip install gseapy) — GSEA prerank engine + MSigDB GMT downloads
- `pandas` (pip install pandas) — data handling
- `scipy` (pip install scipy) — hierarchical clustering
- `numpy` (pip install numpy) — distance matrices
- Standard library: `json`, `csv`, `argparse`, `pathlib`, `collections`, `tempfile`, `re`

---

## Files in the Working Directory for Reference

These files are already present and should be studied but NOT modified:

| File | Purpose |
|------|---------|
| `build_3d_network.py` | **REFERENCE** — existing 3D network builder. Study the HTML template, CSS, JS, node rendering, control panel, force-graph config. Reuse patterns heavily. |
| `pathway_network_3d.html` | **REFERENCE** — example output of existing tool. Open in browser to see the target UX quality. |
| `pathway_mapper_v4.R` | Not needed for the new tool (this is the R pathway registry). |
| `EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv` | **TEST DATA** — use this to test the new script. ~14,000 genes, chicken DESeq2 output. |
| `README.md` | Documents the existing tool. |

---

## Implementation Notes & Gotchas

1. **gseapy prerank output**: The `.res2d` attribute of the prerank result object is a pandas DataFrame with columns: `Name` (alias `Term`), `ES`, `NES`, `NOM p-val`, `FDR q-val`, `FWER p-val`, `Tag %`, `Gene %`, `Lead_genes`. `Lead_genes` is a semicolon-separated string of gene symbols.

2. **Gene symbol case matching**: MSigDB gene sets use UPPERCASE human symbols. The chicken DESeq2 data uses mixed case. Always compare UPPERCASE to UPPERCASE.

3. **The HTML template in build_3d_network.py uses Python f-string-style `{{` and `}}` for literal JS braces inside a `.format()` call.** Follow the same pattern — all JS `{` must be doubled to `{{` in the template string except for the `{placeholder}` substitutions.

4. **3d-force-graph node custom objects**: The existing tool creates custom Three.js objects for nodes (sprite labels, colored spheres, gold hub rings). Study the `buildNodeObject` function in the HTML template carefully. Adapt it to handle the two node types (gene-set diamonds vs gene-node circles).

5. **Force simulation Z-pinning**: The existing tool pins nodes to their pathway Z-layer on each simulation tick. Do the same for gene-set nodes (pin to cluster Z-layer). Let gene nodes float freely in Z, pulled by force toward their connected gene-set nodes.

6. **File size**: The HTML output with embedded DATA JSON should be under 5 MB. If the gene node count gets large (>500), consider only including lead-edge genes rather than all DE genes.

7. **The existing tool loads Three.js and 3d-force-graph from CDN (unpkg)**. Use the same CDN URLs and versions.

8. **Test with the example data**: After building, run `python3 build_gsea_3d_network.py --deseq2-csv EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv` and verify the output HTML opens and renders correctly in a browser.

---

## Deliverables

1. `build_gsea_3d_network.py` — the main build script (single file, like the existing tool)
2. Updated `README.md` — add a section documenting the new GSEA explorer tool
3. A `requirements.txt` with the pip dependencies (gseapy, pandas, scipy, numpy)

---

## Quality Bar

- The 3D visualization should feel as polished as the existing `pathway_network_3d.html`
- Print informative progress messages during the build (same style as existing tool: `[1/N] step description`)
- Handle edge cases gracefully: no significant gene sets, too many gene sets, missing columns in CSV
- The output HTML should work offline after first load (CDN scripts get cached by browser)
