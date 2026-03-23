"""
Per-pathway gene expression figures — one PDF per pathway.
Shows per-gene normalized expression boxplots (Sacral / Free / Pygo),
with DESeq2 significance annotations and log2FC arrows.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA  = Path("/mnt/user-data/uploads")
OUT   = Path("/mnt/user-data/outputs/pathway_gene_plots")
OUT.mkdir(exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
expr  = pd.read_csv(DATA / "expr_pathway_genes.csv")
gsets = pd.read_csv(DATA / "curated_gene_sets_populated.csv")
cov   = pd.read_csv(DATA / "03_gene_set_coverage.csv")

# ── Constants ──────────────────────────────────────────────────────────────────
GROUPS    = ["Sacral", "Free", "Pygo"]
GRP_COL   = {"Sacral": "#1B6CA8", "Free": "#D95F02", "Pygo": "#7570B3"}
GRP_LIGHT = {"Sacral": "#AECDE8", "Free": "#FDBE85", "Pygo": "#BCBDDC"}

NORM_COLS = {
    g: [c for c in expr.columns if c.startswith("norm_") and
        c.replace("norm_", "").rsplit("_", 1)[0] == g]
    for g in GROUPS
}

PATHWAY_LABELS = {
    "necroptosis": "Necroptosis",
    "TNF_NFkB": "TNF–NF-κB",
    "TLR_cGAS_STING": "TLR / cGAS-STING",
    "complement": "Complement",
    "heterophil_myeloid": "Heterophil / Myeloid",
    "angiogenesis": "Angiogenesis",
    "osteogenesis": "Osteogenesis",
    "ECM_remodeling": "ECM Remodeling",
    "sterol_cholesterol_fatty_acid_metabolism": "Sterol / FA Metabolism",
    "fracture_healing": "Fracture Healing",
}

DE_COMPARISONS = [
    ("Sacral_vs_Free",  "Sacral", "Free"),
    ("Pygo_vs_Free",    "Pygo",   "Free"),
    ("Sacral_vs_Pygo",  "Sacral", "Pygo"),
]

def sig_stars(p):
    if pd.isna(p):     return ""
    if p < 0.001:      return "***"
    if p < 0.01:       return "**"
    if p < 0.05:       return "*"
    if p < 0.1:        return "·"
    return ""

def pval_color(p):
    if pd.isna(p):  return "#CBD5E1"
    if p < 0.05:    return "#EF4444"
    if p < 0.1:     return "#F59E0B"
    return "#CBD5E1"

# Set gene symbol as index for fast lookup
expr_idx = expr.set_index("SYMBOL")

pathways = list(gsets.set_name.unique())

for pw_name in pathways:
    pw_label = PATHWAY_LABELS.get(pw_name, pw_name)

    # Genes in this pathway (from gene set definition)
    pw_genes_df = gsets[gsets.set_name == pw_name].copy()
    pw_genes_all = pw_genes_df.gene_symbol.tolist()

    # Keep only genes present in expression data
    pw_genes = [g for g in pw_genes_all if g in expr_idx.index]
    pw_missing = [g for g in pw_genes_all if g not in expr_idx.index]

    if not pw_genes:
        print(f"  Skipping {pw_name}: no matched genes in expression data")
        continue

    # Sort genes by mean expression across all samples (descending)
    norm_cols_all = [c for c in expr.columns if c.startswith("norm_")]
    mean_expr = expr_idx.loc[pw_genes, norm_cols_all].mean(axis=1)
    pw_genes_sorted = mean_expr.sort_values(ascending=False).index.tolist()

    n_genes = len(pw_genes_sorted)

    # ── Figure layout ────────────────────────────────────────────────────────
    # Each gene gets one subplot; arrange in rows of up to 6
    ncols = min(n_genes, 6)
    nrows = int(np.ceil(n_genes / ncols))
    fig_w = max(14, ncols * 2.4)
    fig_h = nrows * 3.6 + 2.2   # extra for title + footer

    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(fig_w, fig_h),
                              squeeze=False)
    fig.patch.set_facecolor("#F8FAFC")

    for gi, gene in enumerate(pw_genes_sorted):
        row, col = divmod(gi, ncols)
        ax = axes[row][col]
        ax.set_facecolor("#F8FAFC")

        row_data = expr_idx.loc[gene]

        # ── Boxplot + jitter ─────────────────────────────────────────────────
        bp_data = [row_data[NORM_COLS[g]].values.astype(float) for g in GROUPS]
        np.random.seed(gi)

        bp = ax.boxplot(bp_data, positions=range(3), patch_artist=True,
                        widths=0.42, showfliers=False,
                        medianprops=dict(color="white", linewidth=2.0),
                        whiskerprops=dict(color="#94A3B8", linewidth=1.1),
                        capprops=dict(color="#94A3B8", linewidth=1.1),
                        boxprops=dict(linewidth=0))
        for patch, g in zip(bp["boxes"], GROUPS):
            patch.set_facecolor(GRP_LIGHT[g])
            patch.set_alpha(0.75)

        for gi2, g in enumerate(GROUPS):
            vals = row_data[NORM_COLS[g]].values.astype(float)
            jit = np.random.uniform(-0.13, 0.13, len(vals))
            ax.scatter(gi2 + jit, vals, color=GRP_COL[g], s=28, zorder=4,
                       edgecolors="white", linewidths=0.5, alpha=0.92)

        # ── DE significance brackets ─────────────────────────────────────────
        all_vals = np.concatenate(bp_data)
        y_max = np.nanmax(all_vals)
        y_min = np.nanmin(all_vals)
        y_range = y_max - y_min if y_max != y_min else y_max * 0.2
        step = y_range * 0.20

        bracket_level = 0
        for comp, g1, g2 in DE_COMPARISONS:
            padj_col = f"{comp}__padj"
            lfc_col  = f"{comp}__log2FC_shrunken"
            if padj_col not in row_data.index:
                continue
            padj = row_data[padj_col]
            lfc  = row_data.get(lfc_col, np.nan)
            stars = sig_stars(padj)
            if not stars:
                continue

            xi = GROUPS.index(g1)
            xj = GROUPS.index(g2)
            yb = y_max + step * (bracket_level + 0.7)

            color = pval_color(padj)
            ax.plot([xi, xi, xj, xj],
                    [yb - step*0.18, yb, yb, yb - step*0.18],
                    color=color, lw=0.9)
            label = stars
            if not pd.isna(lfc):
                direction = "↑" if lfc > 0 else "↓"
                label = f"{direction}{abs(lfc):.1f} {stars}"
            ax.text((xi + xj) / 2, yb + step * 0.06, label,
                    ha="center", va="bottom", fontsize=7.5,
                    color=color, fontweight="bold")
            bracket_level += 1

        # ── Axis formatting ──────────────────────────────────────────────────
        ax.set_xticks(range(3))
        ax.set_xticklabels(GROUPS, fontsize=7.5, rotation=30, ha="right")
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.axhline(0, color="#E2E8F0", lw=0.6, ls="--", zorder=1)

        # Gene title — bold + DE call indicator
        de_call = ""
        for comp, g1, g2 in DE_COMPARISONS:
            call_col = f"{comp}__DE_call"
            if call_col in row_data.index:
                call = str(row_data[call_col])
                if call != "not_sig" and call != "nan":
                    de_call = "  ✦"
                    break

        ax.set_title(f"{gene}{de_call}", fontsize=9, fontweight="bold",
                     color="#1E293B", pad=3)

        if col == 0:
            ax.set_ylabel("Norm. expression", fontsize=7.5)

    # Hide unused subplots
    for gi in range(n_genes, nrows * ncols):
        row, col = divmod(gi, ncols)
        axes[row][col].set_visible(False)

    # ── Figure-level title & annotations ────────────────────────────────────
    # Coverage info from coverage file
    cov_row = cov[cov.set_name == pw_name]
    if not cov_row.empty:
        cr = cov_row.iloc[0]
        cov_str = (f"{cr.matched_gene_count}/{cr.input_gene_set_size} genes matched "
                   f"({cr.percent_coverage:.0f}%)")
        if pw_missing:
            cov_str += f"  |  Not in expression data: {', '.join(pw_missing)}"
    else:
        cov_str = ""

    fig.suptitle(pw_label, fontsize=16, fontweight="bold",
                 color="#1E293B", y=0.99)
    fig.text(0.5, 0.965, cov_str, ha="center", fontsize=8.5,
             color="#64748B", style="italic")

    # Legend row
    grp_handles = [mpatches.Patch(color=GRP_COL[g], label=g) for g in GROUPS]
    de_handle   = mlines.Line2D([], [], color="white", marker="*",
                                 markerfacecolor="#EF4444", markersize=9,
                                 label="✦ = any DESeq2 sig. (padj<0.05)")
    star_note   = mpatches.Patch(color="none",
                                  label="Brackets: ↑/↓log2FC  *p<.05  **p<.01  ***p<.001  ·p<.1")
    fig.legend(handles=grp_handles + [de_handle, star_note],
               loc="lower center", ncol=6, fontsize=8.5,
               framealpha=0.95, edgecolor="#CBD5E1",
               bbox_to_anchor=(0.5, 0.0))

    fig.tight_layout(rect=[0, 0.045, 1, 0.96])

    fname = pw_name.replace("/", "_")
    fig.savefig(OUT / f"{fname}.pdf", bbox_inches="tight",
                facecolor=fig.get_facecolor())
    fig.savefig(OUT / f"{fname}.png", dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Done: {pw_label}  ({n_genes} genes)")

print(f"\nAll figures saved to {OUT}")
