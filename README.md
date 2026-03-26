# 3D Pathway Network Explorer

> **To run the GSEA 3D Network Explorer, use this command:**
> ```bash
> conda activate rnaseq && python3 build_gsea_3d_network.py
> ```

An interactive 3D force-directed network visualization of seven signaling pathways, colored by DESeq2 differential expression results. Runs entirely in the browser — no server required.

---

## Files

| File | Description |
|------|-------------|
| `build_3d_network.py` | Python script that generates the HTML visualization |
| `pathway_mapper_v4.R` | Input: R pathway definitions (nodes, edges, synonyms) |
| `DESeq2_Free_vs_Pygostyle_by_ENS_2023.csv` | Input: DESeq2 results (gene symbols + log2FC + padj) |
| `pathway_network_3d.html` | Output: self-contained interactive HTML file |

---

## Usage

Run from the terminal:

```bash
cd "/path/to/3D_network map"

python3 build_3d_network.py
```

By default, the script uses `pathway_mapper_v4.R` and `DESeq2_FILE_PATH_HERE.csv` from the same folder as `build_3d_network.py`.

If you get a missing-file error, the message now tells you exactly what to edit in `build_3d_network.py`: `default_r_script` for the R file or `default_deseq2_csv` for the DESeq2 CSV.

To override those inputs:

```bash
python3 build_3d_network.py \
  --r-script pathway_mapper_v4.R \
  --deseq2-csv DESeq2_FILE_PATH_HERE.csv \
  --output pathway_network_3d.html \
  --title "Sacral vs Free — 3D Pathway Network" \
  --p-cutoff 0.05
```

Then open `pathway_network_3d.html` in any modern browser (Chrome, Firefox, Safari). An internet connection is required on first open to load the Three.js and 3d-force-graph libraries from CDN.

### Arguments

| Argument | Required | Default | Description |
|----------|----------|---------|-------------|
| `--r-script` | No | local `pathway_mapper_v4.R` | Path to `pathway_mapper_v4.R` |
| `--deseq2-csv` | No | local `DESeq2_FILE_PATH_HERE.csv` | Path to DESeq2 results CSV |
| `--output` | No | `pathway_network_3d.html` | Output HTML filename |
| `--title` | No | `Sacral vs Free — 3D Pathway Network` | Title shown in the browser |
| `--p-cutoff` | No | `0.05` | Significance threshold for padj |

### Dependencies

Python standard library only — no `pip install` needed:
- `re`, `json`, `csv`, `argparse`, `pathlib`, `collections`

---

## Pathways

Seven signaling pathways are included, each assigned a z-layer in 3D space:

| Pathway | Display Name | Color | Z-layer |
|---------|-------------|-------|---------|
| `cgas_sting` | cGAS-STING / Type I IFN Signaling | Slate blue | 0 |
| `apoptosis` | Apoptosis | Dark red | 150 |
| `necroptosis` | Necroptosis | Forest green | 300 |
| `inflammasome` | Inflammasome (NLRP3 Canonical) | Dark orange | 450 |
| `osteoclast` | Osteoclast Differentiation (RANKL-RANK) | Navy | 600 |
| `osteoblast` | Osteoblast Differentiation (BMP/Wnt) | Dodger blue | 750 |
| `oxysterol` | Oxysterol Signaling in Bone Fusion | Goldenrod | 900 |

---

## Network Summary (current dataset)

- **122 nodes** (unique gene symbols across all pathways)
- **140 edges**
- **13 hub nodes** appearing in 2+ pathways (marked with gold rings)
- **66 significant** (padj < 0.05), **37 not significant**, **19 missing genes**

### Key hub nodes

| Gene | Pathways |
|------|----------|
| TNF | cgas_sting, apoptosis, necroptosis |
| TNFRSF1A | apoptosis, necroptosis |
| TRADD | apoptosis, necroptosis |
| CASP8 | apoptosis, necroptosis |
| NFKB1 | cgas_sting, inflammasome, osteoclast |
| RELA | cgas_sting, inflammasome |
| TRAF6 | cgas_sting, inflammasome, osteoclast |
| RUNX2 | osteoblast, oxysterol |
| SP7 | osteoblast, oxysterol |
| TNFSF11 | osteoclast, oxysterol |
| ZBP1 | cgas_sting, necroptosis |
| IFNB1 | cgas_sting, oxysterol |
| IFNAR1 | cgas_sting, oxysterol |

### Synonym resolution

Some genes are listed under alternative symbols in the DESeq2 results. The script checks synonyms defined in `pathway_mapper_v4.R` automatically. The most important case:

- **STING1** is matched via its synonym **TMEM173** in the CSV

### Expected missing genes (~19)

These genes are absent from the DESeq2 dataset and render as small grey nodes. This is biologically informative, not an error — they likely represent chicken genome gaps or genes below detection threshold:

`ZBP1, IFI16, IKBKG, IFNB1, CXCL10, ISG15, TNF, BAX, IRF3, NLRP3, PYCARD, AIM2, NLRC4, IL1B, GSDMD, WNT3A, DVL2, HSD3B7, SULT2B1`

---

## Visualization Guide

### Node appearance

| Style | Meaning |
|-------|---------|
| Large sphere | Present; significant (padj < 0.05) |
| Medium sphere | Present; not significant or padj NA |
| Small sphere | Missing gene (not in DESeq2 dataset) |
| Red color | Upregulated (log2FC > 0) |
| Blue color | Downregulated (log2FC < 0) |
| White/grey | Near-zero FC or missing |
| Gold ring | Hub node (appears in 2+ pathways) |
| White ring | Current search match |

Color scale is capped at ±3 log2FC so that biologically meaningful fold-changes (1–2) show clear color contrast.

### Edge appearance

- **Thick edges** = "core" signaling edges; **thin** = output/downstream edges
- Edge width also scales with the |log2FC| of the source node
- Edge color indicates which pathway the edge belongs to
- Moving particles show signal flow direction

### Controls panel (left sidebar)

| Control | Function |
|---------|----------|
| Search box | Find genes by symbol or label; matching nodes pulse and get a white ring |
| Pathway checkboxes | Show/hide individual pathways |
| **only** button | Isolate a single pathway (hides all others, flies camera to that layer) |
| Show All | Restore all pathways |
| Significant only | Filter to padj < 0.05 nodes |
| Hub nodes only | Show only multi-pathway genes |
| Show labels | Toggle gene name labels |
| Show arrows | Toggle directional arrowheads |
| Light background | Switch to light theme |
| Reset Camera | Return to default view |
| Reheat Simulation | Restart the force layout |
| Screenshot PNG | Save current view as PNG |
| Export Node CSV | Download visible nodes as CSV |

### Keyboard shortcuts

| Key | Action |
|-----|--------|
| `R` | Reset camera |
| `1` – `7` | Fly to pathway layer (1=cGAS-STING … 7=Oxysterol) |
| `L` | Toggle labels |
| `H` | Toggle hub-only filter |
| `Space` | Pause / resume force simulation |

### Interactions

- **Hover** a node → tooltip with gene, log2FC, padj, −log10(padj), status
- **Click** a node → highlight neighborhood, fly camera to node, open detail panel
- **Click background** → reset highlight
- **Legend** header → click to collapse/expand

---

## Organism

*Gallus gallus* (chicken). The DESeq2 CSV uses gene symbols in the first column despite the filename suggesting Ensembl IDs — no ID conversion is needed.

---

---

## GSEA 3D Network Explorer (`build_gsea_3d_network.py`)

A companion tool that performs GSEA preranked enrichment analysis against gene set databases and visualizes the results as an exploratory 3D network. Unlike the curated-pathway tool above, this one is fully data-driven: the network structure is determined by the enrichment results.

### How it works

```
DESeq2 CSV → rank genes by Wald stat → GSEA prerank (gseapy)
  → filter significant gene sets (FDR < 0.25)
  → cluster gene sets by Jaccard similarity of lead-edge genes
  → build 3D network: gene-set nodes (diamonds) + gene nodes (spheres)
  → overlay log2FC / padj → generate self-contained HTML
```

### Quick start

```bash
# Activate an environment with gseapy, pandas, scipy, numpy
conda activate rnaseq

python3 build_gsea_3d_network.py \
  --deseq2-csv EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv \
  --gene-sets H \
  --output gsea_network_3d.html \
  --title "Free vs Pygo — GSEA 3D Explorer"
```

Open `gsea_network_3d.html` in any modern browser. Internet connection needed on first load for CDN libraries.

### Gene set library options

| `--gene-sets` value | Library | Source |
|---------------------|---------|--------|
| `H` | MSigDB Hallmark (default — 50 well-curated sets) | Enrichr |
| `KEGG` | KEGG 2019 Human | Enrichr |
| `KEGG_GGA` | KEGG Gallus gallus (chicken-specific) | KEGG REST API |
| `REACTOME` | Reactome 2022 | Enrichr |
| `REACTOME_GGA` | Reactome Gallus gallus (chicken-specific) | Reactome bulk mapping + UniProt |
| `C5_GO_BP` | GO Biological Process | Enrichr |
| `C5_GO_MF` | GO Molecular Function | Enrichr |
| `C5_GO_CC` | GO Cellular Component | Enrichr |
| `WIKIPATHWAYS` | WikiPathways Human | Enrichr |
| `WP_GGA` | WikiPathways Gallus gallus (chicken-specific) | WikiPathways GMT + NCBI gene_info |
| `/path/to/file.gmt` | Any local GMT file | Local |

Gene sets are downloaded automatically and cached in `~/.gsea3d_cache/` — subsequent runs load instantly from cache. The three `_GGA` libraries use chicken-native gene IDs converted to HGNC-style symbols, which match directly against the gene column in a chicken DESeq2 CSV.

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--deseq2-csv` | required | Path to DESeq2 results CSV |
| `--gene-sets` | `H` | Gene set library or GMT file path |
| `--output` | `gsea_network_3d.html` | Output HTML path |
| `--title` | `GSEA 3D Network Explorer` | Title in browser tab |
| `--fdr-cutoff` | `0.25` | FDR threshold for significant gene sets |
| `--p-cutoff` | `0.05` | padj threshold for gene significance coloring |
| `--permutations` | `1000` | GSEA permutations (use 100 for quick testing) |
| `--max-sets` | `40` | Maximum gene sets to show |
| `--min-size` | `15` | Minimum gene set size |
| `--max-size` | `500` | Maximum gene set size |

### Node appearance

| Style | Meaning |
|-------|---------|
| Diamond (◆) | Gene set node; size ∝ \|NES\|; color = NES (blue=down, red=up) |
| Colored ring on diamond | Cluster membership color |
| Large sphere | Significant gene (padj < threshold) |
| Small sphere | Non-significant or missing gene |
| Red sphere | Upregulated (log2FC > 0) |
| Blue sphere | Downregulated (log2FC < 0) |
| Gold ring | Hub gene (in ≥2 gene sets) |
| White ring | Search match |

### Controls unique to GSEA tool

- **Gene Set Clusters** section: checkboxes to show/hide each cluster; "only" button to isolate one cluster
- **Gene Sets** section: expandable list of all significant gene sets grouped by cluster; click any to fly to it
- **Upregulated/Downregulated filter**: show only gene sets with NES > 0 or NES < 0
- **Hide overlap edges**: toggle the faded dashed edges between gene sets with shared lead-edge genes
- **Legend**: includes NES color bar (blue → red) in addition to the log2FC color bar

### Dependencies

```
pip install gseapy pandas scipy numpy
```

Or via conda: `conda install -c bioconda gseapy scipy numpy pandas`

---

## Notes for re-running

- The script is **read-only** with respect to `pathway_mapper_v4.R` — it only reads the file, never modifies it.
- To use a different DESeq2 result file, pass it with `--deseq2-csv`. The first (unnamed) column must contain gene symbols.
- To change the significance cutoff, use `--p-cutoff 0.01` (or any value).
- The output HTML is fully self-contained — share it by sending the single `.html` file. Recipients need only a browser and internet access for the CDN libraries.
