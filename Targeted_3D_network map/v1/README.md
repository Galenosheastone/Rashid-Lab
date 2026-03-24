# 3D Pathway Network Explorer

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

## Notes for re-running

- The script is **read-only** with respect to `pathway_mapper_v4.R` — it only reads the file, never modifies it.
- To use a different DESeq2 result file, pass it with `--deseq2-csv`. The first (unnamed) column must contain gene symbols.
- To change the significance cutoff, use `--p-cutoff 0.01` (or any value).
- The output HTML is fully self-contained — share it by sending the single `.html` file. Recipients need only a browser and internet access for the CDN libraries.
