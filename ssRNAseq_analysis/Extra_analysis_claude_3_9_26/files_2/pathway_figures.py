"""
Pathway scoring figure suite
Produces 5 publication-quality figures as individual PDFs + one combined summary.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as mticker
import seaborn as sns
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA = Path("/mnt/user-data/uploads")
OUT  = Path("/mnt/user-data/outputs")
OUT.mkdir(exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
cov   = pd.read_csv(DATA / "03_gene_set_coverage.csv")
gsva  = pd.read_csv(DATA / "04_scores_gsva.csv").set_index("pathway")
ssg   = pd.read_csv(DATA / "05_scores_ssgsea.csv").set_index("pathway")
sing  = pd.read_csv(DATA / "06_scores_singscore.csv").set_index("pathway")
st_g  = pd.read_csv(DATA / "10_stats_gsva.csv")
st_s  = pd.read_csv(DATA / "11_stats_ssgsea.csv")
st_k  = pd.read_csv(DATA / "12_stats_singscore.csv")

# ── Shared palette & helpers ──────────────────────────────────────────────────
GROUPS  = ["Sacral", "Free", "Pygo"]
GRP_COL = {"Sacral": "#1B6CA8", "Free": "#D95F02", "Pygo": "#7570B3"}
GRP_LIGHT = {"Sacral": "#AECDE8", "Free": "#FDBE85", "Pygo": "#BCBDDC"}
METHOD_COL = {"GSVA": "#0369A1", "ssGSEA": "#059669", "singscore": "#7C3AED"}

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

def pwy_label(s):
    return PATHWAY_LABELS.get(s, s)

def parse_group(col):
    # "norm_Sacral_1" → "Sacral"
    return col.replace("norm_", "").rsplit("_", 1)[0]

def score_long(df, method):
    rows = []
    for pathway, row in df.iterrows():
        for col, val in row.items():
            rows.append({"pathway": pathway, "sample": col,
                         "group": parse_group(col), "score": val, "method": method})
    return pd.DataFrame(rows)

long_gsva = score_long(gsva, "GSVA")
long_ssg  = score_long(ssg,  "ssGSEA")
long_sing = score_long(sing, "singscore")
long_all  = pd.concat([long_gsva, long_ssg, long_sing], ignore_index=True)

def sig_stars(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    if p < 0.1:   return "·"
    return ""

def omnibus_pvals(st):
    return st[st.test_scope == "omnibus"].set_index("pathway")[["p_value","p_adj"]]

def pairwise_pvals(st):
    return st[st.test_scope == "pairwise"].copy()

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1a — Gene set coverage bar chart
# ══════════════════════════════════════════════════════════════════════════════
fig1a, ax = plt.subplots(figsize=(9, 5.5))
fig1a.patch.set_facecolor("#F8FAFC")
ax.set_facecolor("#F8FAFC")

cov_sorted = cov.sort_values("percent_coverage", ascending=True)
pathways_sorted = [pwy_label(p) for p in cov_sorted.set_name]
pct   = cov_sorted.percent_coverage.values
n_in  = cov_sorted.matched_gene_count.values
n_out = cov_sorted.unmatched_gene_count.values
total = n_in + n_out

bar_colors = ["#0EA5E9" if p >= 80 else "#F59E0B" if p >= 60 else "#EF4444"
              for p in pct]

bars = ax.barh(pathways_sorted, pct, color=bar_colors, height=0.55,
               zorder=3, alpha=0.9)
ax.barh(pathways_sorted, 100 - pct, left=pct, color="#E2E8F0",
        height=0.55, zorder=2)

for i, (b, ni, tot) in enumerate(zip(bars, n_in, total)):
    w = b.get_width()
    ax.text(w + 1.5, i, f"{ni}/{tot} genes", va="center", ha="left",
            fontsize=9, color="#334155")
    ax.text(min(w - 2, 96), i, f"{w:.0f}%", va="center", ha="right",
            fontsize=9, fontweight="bold",
            color="white" if w > 30 else "#334155")

ax.set_xlim(0, 120)
ax.set_xlabel("Gene set coverage (%)", fontsize=12)
ax.set_title("Gene Set Coverage", fontsize=14, fontweight="bold", pad=12,
             color="#1E293B")
ax.axvline(80, color="#94A3B8", lw=0.8, ls="--", zorder=1)
ax.axvline(60, color="#CBD5E1", lw=0.8, ls="--", zorder=1)
ax.tick_params(axis="y", labelsize=11)
ax.tick_params(axis="x", labelsize=10)
ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.spines[["top","right"]].set_visible(False)

legend_patches = [
    mpatches.Patch(color="#0EA5E9", label="≥80%"),
    mpatches.Patch(color="#F59E0B", label="60–79%"),
    mpatches.Patch(color="#EF4444", label="<60%"),
    mpatches.Patch(color="#E2E8F0", label="Not matched"),
]
ax.legend(handles=legend_patches, loc="lower right", fontsize=9,
          framealpha=0.9, edgecolor="#CBD5E1")

fig1a.tight_layout()
fig1a.savefig(OUT / "fig1a_coverage.pdf", bbox_inches="tight",
              facecolor=fig1a.get_facecolor())
fig1a.savefig(OUT / "fig1a_coverage.png", dpi=180, bbox_inches="tight",
              facecolor=fig1a.get_facecolor())
plt.close(fig1a)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1b — Gene membership matrix (wide/landscape, PPT-friendly)
#   Genes on X-axis (columns), Pathways on Y-axis (rows)
# ══════════════════════════════════════════════════════════════════════════════

# Build gene × pathway membership
all_genes = set()
pw_genes  = {}
for _, row in cov.iterrows():
    genes = [g.strip() for g in row.matched_genes.split(";") if g.strip()]
    pw_genes[row.set_name] = set(genes)
    all_genes.update(genes)

# Sort genes: most shared first, then alphabetical within tier
gene_counts = {g: sum(1 for pw in pw_genes.values() if g in pw)
               for g in all_genes}
genes_ordered = sorted(all_genes, key=lambda g: (-gene_counts[g], g))

pathways_order = list(cov.set_name)

n_genes    = len(genes_ordered)
n_pathways = len(pathways_order)

# Build matrix — now: pathways (rows) × genes (cols)
mat = np.zeros((n_pathways, n_genes))
for i, pw in enumerate(pathways_order):
    for j, g in enumerate(genes_ordered):
        mat[i, j] = 1 if g in pw_genes[pw] else 0

membership = np.sum(mat, axis=0)   # per-gene column sum

# ── Figure sizing: wide landscape ──────────────────────────────────────────
col_width  = 0.30   # inches per gene column
fig_width  = max(20, n_genes * col_width + 4.5)   # +4.5 for y-labels + margins
fig_height = 7.5    # fixed — fits a PPT slide comfortably

fig1b, ax2 = plt.subplots(figsize=(fig_width, fig_height))
fig1b.patch.set_facecolor("#F8FAFC")
ax2.set_facecolor("#F8FAFC")

# Alternating column bands (one per gene) for easy eye tracking
for j in range(n_genes):
    color = "#F1F5F9" if j % 2 == 0 else "#FFFFFF"
    ax2.axvspan(j - 0.5, j + 0.5, color=color, zorder=1)

# Horizontal pathway dividers
for i in range(n_pathways):
    ax2.axhline(i, color="#E2E8F0", lw=0.6, zorder=2)

# Tier boundaries — vertical dashed lines between gene tiers
tier_boundaries_x = []
prev_count = gene_counts[genes_ordered[0]]
for j, g in enumerate(genes_ordered):
    c = gene_counts[g]
    if c != prev_count:
        tier_boundaries_x.append(j - 0.5)
        prev_count = c

for xb in tier_boundaries_x:
    ax2.axvline(xb, color="#475569", lw=1.0, ls="--", zorder=4, alpha=0.55)

# Tier labels — placed just below the bottom pathway row in data coords
tier_starts = [0] + [int(b + 0.5) for b in tier_boundaries_x]
tier_ends   = [int(b - 0.5) for b in tier_boundaries_x] + [n_genes - 1]
for ts, te in zip(tier_starts, tier_ends):
    g_example = genes_ordered[ts]
    count = gene_counts[g_example]
    mid_x = (ts + te) / 2
    label = f"In {count} pathway{'s' if count > 1 else ''}"
    ax2.text(mid_x, -0.75, label, fontsize=8.5, color="#64748B",
             va="top", ha="center", style="italic")

# Dots
dot_size = max(40, int(col_width * 240))
cmap_dots = {1: "#94A3B8", 2: "#F59E0B", 3: "#EF4444"}

for i, pw in enumerate(pathways_order):
    for j, g in enumerate(genes_ordered):
        if mat[i, j] == 1:
            m = int(membership[j])
            c = cmap_dots.get(m, "#7C3AED")
            ax2.scatter(j, i, s=dot_size, color=c, zorder=3,
                        linewidths=0, alpha=0.88)

# X-axis — gene labels (rotated)
gene_fontsize = min(9.5, max(6.5, col_width * 28))
ax2.set_xticks(range(n_genes))
ax2.set_xticklabels(genes_ordered, rotation=60, ha="right",
                     fontsize=gene_fontsize, fontfamily="monospace")
ax2.xaxis.set_tick_params(length=0)

# Secondary x-axis on top
ax2_top = ax2.twiny()
ax2_top.set_xlim(ax2.get_xlim())
ax2_top.set_xticks(range(n_genes))
ax2_top.set_xticklabels(genes_ordered, rotation=60, ha="left",
                          fontsize=gene_fontsize, fontfamily="monospace")
ax2_top.xaxis.set_tick_params(length=0)

# Y-axis — pathway labels
ax2.set_yticks(range(n_pathways))
ax2.set_yticklabels([pwy_label(p) for p in pathways_order],
                     fontsize=11, fontweight="bold")
ax2.tick_params(axis="y", length=0)

ax2.set_xlim(-0.5, n_genes - 0.5)
ax2.set_ylim(-1.4, n_pathways - 0.5)   # extra room at bottom for tier labels
ax2.spines[:].set_visible(False)

# Legend — top-right, outside plot area
legend2 = [
    mpatches.Patch(color="#94A3B8", label="In 1 pathway  (unique)"),
    mpatches.Patch(color="#F59E0B", label="In 2 pathways (shared)"),
    mpatches.Patch(color="#EF4444", label="In 3 pathways (shared)"),
    mpatches.Patch(color="#7C3AED", label="In 4+ pathways (core)"),
]
ax2.legend(handles=legend2, fontsize=9,
           framealpha=0.95, edgecolor="#CBD5E1",
           title="Gene pathway membership", title_fontsize=9.5,
           loc="upper left", bbox_to_anchor=(1.01, 1.0),
           borderaxespad=0)

ax2.set_title("Gene Membership Across Pathways", fontsize=14,
               fontweight="bold", pad=14, color="#1E293B")

fig1b.tight_layout()
fig1b.savefig(OUT / "fig1b_gene_membership.pdf", bbox_inches="tight",
              facecolor=fig1b.get_facecolor())
fig1b.savefig(OUT / "fig1b_gene_membership.png", dpi=150, bbox_inches="tight",
              facecolor=fig1b.get_facecolor())
plt.close(fig1b)
print("Figure 1 done")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Pathway score heatmap (all 3 methods side by side)
# ══════════════════════════════════════════════════════════════════════════════
def make_heatmap_mat(df):
    """Sort columns by group order"""
    cols = []
    for g in GROUPS:
        cols += [c for c in df.columns if parse_group(c) == g]
    return df[cols]

gsva_m  = make_heatmap_mat(gsva)
ssg_m   = make_heatmap_mat(ssg)
sing_m  = make_heatmap_mat(sing)

row_labels = [pwy_label(p) for p in gsva.index]
col_labels = [parse_group(c) + "\n" + c.split("_")[-1] for c in gsva_m.columns]

fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6.5),
                            gridspec_kw={"wspace": 0.06})
fig2.patch.set_facecolor("#F8FAFC")

cmap_heat = LinearSegmentedColormap.from_list(
    "rwb", ["#053061", "#FFFFFF", "#67001F"], N=256)

for ax, df_h, title, st in zip(axes2,
                                [gsva_m, ssg_m, sing_m],
                                ["GSVA", "ssGSEA", "singscore"],
                                [st_g, st_s, st_k]):
    mat = df_h.values
    vabs = np.nanmax(np.abs(mat))
    im = ax.imshow(mat, cmap=cmap_heat, aspect="auto",
                   vmin=-vabs, vmax=vabs)

    # Group dividers and labels
    n_per = 4
    for gi, grp in enumerate(GROUPS):
        start = gi * n_per
        mid   = start + n_per / 2 - 0.5
        ax.text(mid, -0.8, grp, ha="center", va="bottom", fontsize=10,
                fontweight="bold", color=list(GRP_COL.values())[gi])
        if gi > 0:
            ax.axvline(start - 0.5, color="white", lw=2)
        # Color band at top
        ax.add_patch(mpatches.Rectangle(
            (start - 0.5, -1.5), n_per, 0.8,
            color=list(GRP_COL.values())[gi], clip_on=False))

    # p-value stars on right
    omni = omnibus_pvals(st)
    for ri, pw in enumerate(df_h.index):
        stars = sig_stars(omni.loc[pw, "p_value"]) if pw in omni.index else ""
        if stars:
            ax.text(df_h.shape[1] + 0.1, ri, stars, va="center",
                    fontsize=11, color="#EF4444", fontweight="bold")

    ax.set_xticks(range(df_h.shape[1]))
    ax.set_xticklabels(
        [c.split("_")[-1] for c in df_h.columns],
        fontsize=8, rotation=0)
    ax.set_yticks(range(len(row_labels)))
    ax.tick_params(length=0)

    if ax is axes2[0]:
        ax.set_yticklabels(row_labels, fontsize=9.5)
    else:
        ax.set_yticklabels([])

    ax.set_title(title, fontsize=13, fontweight="bold",
                 color=METHOD_COL[title], pad=14)
    ax.spines[:].set_visible(False)

    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.ax.tick_params(labelsize=7.5)
    cb.set_label("Score", fontsize=8)

fig2.suptitle("Pathway Scores by Group — All Three Methods\n"
              "Stars indicate nominal p < 0.05 (Kruskal-Wallis)",
              fontsize=14, fontweight="bold", y=1.04, color="#1E293B")
fig2.tight_layout(rect=[0, 0, 1, 0.97])
fig2.savefig(OUT / "fig2_heatmap_all_methods.pdf", bbox_inches="tight",
             facecolor=fig2.get_facecolor())
fig2.savefig(OUT / "fig2_heatmap_all_methods.png", dpi=180,
             bbox_inches="tight", facecolor=fig2.get_facecolor())
plt.close(fig2)
print("Figure 2 done")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Pathway significance overview (dotplot/lollipop across methods)
# ══════════════════════════════════════════════════════════════════════════════
def median_by_group(df):
    """Return dict of {pathway: {group: median_score}}"""
    result = {}
    for pw in df.index:
        result[pw] = {}
        for g in GROUPS:
            cols = [c for c in df.columns if parse_group(c) == g]
            result[pw][g] = np.median(df.loc[pw, cols])
    return result

medians_gsva = median_by_group(gsva)
medians_ssg  = median_by_group(ssg)
medians_sing = median_by_group(sing)

fig3, axes3 = plt.subplots(1, 3, figsize=(20, 7), sharey=True,
                            gridspec_kw={"wspace": 0.04})
fig3.patch.set_facecolor("#F8FAFC")

pw_order = list(cov.set_name)
pw_labels = [pwy_label(p) for p in pw_order]
y_pos = {pw: i for i, pw in enumerate(pw_order)}
n_pw  = len(pw_order)

for ax, medians, st, title in zip(
        axes3,
        [medians_gsva, medians_ssg, medians_sing],
        [st_g, st_s, st_k],
        ["GSVA", "ssGSEA", "singscore"]):

    ax.set_facecolor("#F8FAFC")
    pw_pairs = pairwise_pvals(st)

    for pw in pw_order:
        yi = y_pos[pw]
        vals = [medians[pw][g] for g in GROUPS]
        vmin, vmax = min(vals), max(vals)

        # Range line
        ax.plot([vmin, vmax], [yi, yi], color="#CBD5E1", lw=2.5, zorder=2)

        for gi, g in enumerate(GROUPS):
            ax.scatter(medians[pw][g], yi,
                       color=GRP_COL[g], s=100, zorder=4,
                       edgecolors="white", linewidths=0.8)

        # Significance annotation — find best pairwise
        pw_sub = pw_pairs[pw_pairs.pathway == pw]
        if not pw_sub.empty:
            best_p = pw_sub.p_value.min()
            best_row = pw_sub.loc[pw_sub.p_value.idxmin()]
            stars = sig_stars(best_p)
            if stars:
                xann = vmax + (vmax - vmin) * 0.08 + 0.02
                ax.text(xann, yi, stars, va="center", ha="left",
                        fontsize=11, color="#EF4444", fontweight="bold")
                # Arrow indicating direction
                g1, g2 = best_row.group1, best_row.group2
                if g1 in medians[pw] and g2 in medians[pw]:
                    diff = medians[pw][g1] - medians[pw][g2]

    # Alternating row shading
    for i, pw in enumerate(pw_order):
        if i % 2 == 0:
            ax.axhspan(i - 0.45, i + 0.45, color="#F1F5F9", zorder=1)

    ax.axvline(0, color="#94A3B8", lw=0.8, ls="--", zorder=2)
    ax.set_yticks(range(n_pw))
    ax.set_yticklabels([], fontsize=10)
    ax.set_ylim(-0.7, n_pw - 0.3)
    if ax is axes3[0]:
        for yi, label in enumerate(pw_labels):
            ax.text(-0.01, yi, label, transform=ax.get_yaxis_transform(),
                    fontsize=9.5, va="center", ha="right", color="#1E293B")
    ax.set_xlabel("Median pathway score", fontsize=10)
    ax.set_title(title, fontsize=13, fontweight="bold",
                 color=METHOD_COL[title], pad=10)
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

# Legend
handles = [mpatches.Patch(color=GRP_COL[g], label=g) for g in GROUPS]
axes3[1].legend(handles=handles, loc="lower right", fontsize=9,
                framealpha=0.95, edgecolor="#CBD5E1")

# Stars legend
axes3[2].text(1.05, 0.98, "Stars = nominal p\n* < .05  ** < .01  *** < .001\n· < .10",
              transform=axes3[2].transAxes, fontsize=8, va="top",
              bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                        edgecolor="#CBD5E1", alpha=0.95))

fig3.suptitle("Pathway Score Range by Group — Median per Group\n"
              "Stars = most significant pairwise comparison",
              fontsize=14, fontweight="bold", y=1.02, color="#1E293B")
fig3.subplots_adjust(left=0.22, right=0.96, top=0.88, bottom=0.1)
fig3.savefig(OUT / "fig3_score_lollipop.pdf", bbox_inches="tight",
             facecolor=fig3.get_facecolor())
fig3.savefig(OUT / "fig3_score_lollipop.png", dpi=180, bbox_inches="tight",
             facecolor=fig3.get_facecolor())
plt.close(fig3)
print("Figure 3 done")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Per-pathway score boxplots with individual points (GSVA focused)
# ══════════════════════════════════════════════════════════════════════════════
n_pw = len(pw_order)
ncols = 5
nrows = 2
fig4, axes4 = plt.subplots(nrows, ncols, figsize=(18, 7.5))
fig4.patch.set_facecolor("#F8FAFC")

pw_pairs_g = pairwise_pvals(st_g)

for idx, pw in enumerate(pw_order):
    row, col = divmod(idx, ncols)
    ax = axes4[row][col]
    ax.set_facecolor("#F8FAFC")

    scores_by_group = {}
    for g in GROUPS:
        cols = [c for c in gsva.columns if parse_group(c) == g]
        scores_by_group[g] = gsva.loc[pw, cols].values

    # Box plots
    bp_data = [scores_by_group[g] for g in GROUPS]
    bp = ax.boxplot(bp_data, positions=range(len(GROUPS)), patch_artist=True,
                    widths=0.45, showfliers=False,
                    medianprops=dict(color="white", linewidth=2),
                    whiskerprops=dict(color="#94A3B8", linewidth=1.2),
                    capprops=dict(color="#94A3B8", linewidth=1.2),
                    boxprops=dict(linewidth=0))
    for patch, g in zip(bp["boxes"], GROUPS):
        patch.set_facecolor(GRP_LIGHT[g])
        patch.set_alpha(0.7)

    # Jitter points
    np.random.seed(42)
    for gi, g in enumerate(GROUPS):
        jitter = np.random.uniform(-0.12, 0.12, len(scores_by_group[g]))
        ax.scatter(gi + jitter, scores_by_group[g],
                   color=GRP_COL[g], s=32, zorder=4,
                   edgecolors="white", linewidths=0.5, alpha=0.9)

    # Significance brackets
    pw_sub = pw_pairs_g[pw_pairs_g.pathway == pw].copy()
    pw_sub = pw_sub[pw_sub.p_value < 0.1]
    pw_sub = pw_sub.sort_values("p_value")

    y_max = max(v for vals in bp_data for v in vals)
    y_range = y_max - min(v for vals in bp_data for v in vals)
    bracket_step = y_range * 0.18

    for bi, (_, brow) in enumerate(pw_sub.iterrows()):
        g1i = GROUPS.index(brow.group1)
        g2i = GROUPS.index(brow.group2)
        stars = sig_stars(brow.p_value)
        if not stars:
            continue
        yb = y_max + bracket_step * (bi + 0.6)
        ax.plot([g1i, g1i, g2i, g2i],
                [yb - bracket_step*0.2, yb, yb, yb - bracket_step*0.2],
                color="#475569", lw=1)
        ax.text((g1i + g2i) / 2, yb + bracket_step * 0.05, stars,
                ha="center", va="bottom", fontsize=9.5,
                color="#EF4444" if brow.p_value < 0.05 else "#94A3B8",
                fontweight="bold")

    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GROUPS, fontsize=8, rotation=30, ha="right")
    ax.set_title(pwy_label(pw), fontsize=9.5, fontweight="bold",
                 color="#1E293B", pad=4)
    ax.tick_params(axis="y", labelsize=7.5)
    ax.spines[["top","right"]].set_visible(False)
    ax.axhline(0, color="#CBD5E1", lw=0.7, ls="--", zorder=1)

    if col == 0:
        ax.set_ylabel("GSVA score", fontsize=8.5)

fig4.suptitle("GSVA Pathway Scores by Group\n"
              "Boxes = IQR; points = individual samples; brackets = nominal significance",
              fontsize=13, fontweight="bold", y=1.02, color="#1E293B")
fig4.tight_layout(rect=[0, 0, 1, 0.97])
fig4.savefig(OUT / "fig4_boxplots_gsva.pdf", bbox_inches="tight",
             facecolor=fig4.get_facecolor())
fig4.savefig(OUT / "fig4_boxplots_gsva.png", dpi=180, bbox_inches="tight",
             facecolor=fig4.get_facecolor())
plt.close(fig4)
print("Figure 4 done")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 5 — Multi-method consensus: 2×2 layout, spacious
# ══════════════════════════════════════════════════════════════════════════════
from scipy import stats as scipy_stats

def method_medians(df):
    rows = []
    for pw in df.index:
        for g in GROUPS:
            cols = [c for c in df.columns if parse_group(c) == g]
            rows.append({"pathway": pw, "group": g,
                         "median": np.median(df.loc[pw, cols])})
    return pd.DataFrame(rows)

mm_gsva = method_medians(gsva).rename(columns={"median": "GSVA"})
mm_ssg  = method_medians(ssg).rename(columns={"median": "ssGSEA"})
mm_sing = method_medians(sing).rename(columns={"median": "singscore"})
mm = mm_gsva.merge(mm_ssg, on=["pathway","group"]).merge(mm_sing, on=["pathway","group"])

# Per-pathway average across groups (for single label per pathway dot)
mm_avg = mm.groupby("pathway")[["GSVA","ssGSEA","singscore"]].mean().reset_index()

fig5, axes5 = plt.subplots(2, 2, figsize=(16, 14))
fig5.patch.set_facecolor("#F8FAFC")

pairs = [("GSVA", "ssGSEA"), ("GSVA", "singscore"), ("ssGSEA", "singscore")]
scatter_axes = [axes5[0,0], axes5[0,1], axes5[1,0]]

for ax, (m1, m2) in zip(scatter_axes, pairs):
    ax.set_facecolor("#F8FAFC")

    # Draw points per group
    for grp in GROUPS:
        sub = mm[mm.group == grp]
        ax.scatter(sub[m1], sub[m2], color=GRP_COL[grp], s=110, zorder=4,
                   edgecolors="white", linewidths=1.0, alpha=0.88, label=grp)

    # Single pathway label per pathway — placed at the mean position across groups
    # with a small horizontal offset to avoid sitting on dots
    for _, row in mm_avg.iterrows():
        xv, yv = row[m1], row[m2]
        # Nudge label right if on right half, left if on left half
        xrange = mm[m1].max() - mm[m1].min()
        yrange = mm[m2].max() - mm[m2].min()
        dx = xrange * 0.04
        dy = yrange * 0.015
        ax.annotate(pwy_label(row.pathway),
                    xy=(xv, yv), xytext=(xv + dx, yv + dy),
                    fontsize=8.5, color="#334155",
                    va="center", ha="left",
                    arrowprops=dict(arrowstyle="-", color="#CBD5E1",
                                   lw=0.6, shrinkA=5, shrinkB=0))

    # Regression line over full data range
    r, _ = scipy_stats.pearsonr(mm[m1], mm[m2])
    slope, intercept, *_ = scipy_stats.linregress(mm[m1], mm[m2])
    x_lo, x_hi = mm[m1].min() - 0.05, mm[m1].max() + 0.05
    xs = np.linspace(x_lo, x_hi, 200)
    ax.plot(xs, slope * xs + intercept, color="#94A3B8", lw=1.5,
            ls="--", zorder=2)

    ax.axhline(0, color="#E2E8F0", lw=0.8, zorder=1)
    ax.axvline(0, color="#E2E8F0", lw=0.8, zorder=1)

    ax.set_xlabel(m1, fontsize=13, fontweight="bold", color=METHOD_COL[m1],
                  labelpad=8)
    ax.set_ylabel(m2, fontsize=13, fontweight="bold", color=METHOD_COL[m2],
                  labelpad=8)
    ax.set_title(f"{m1} vs {m2}", fontsize=13, fontweight="bold",
                 color="#1E293B", pad=10)
    ax.tick_params(labelsize=10)
    ax.spines[["top","right"]].set_visible(False)

    # r annotation — top-left corner, away from labels
    ax.text(0.04, 0.96, f"r = {r:.2f}",
            transform=ax.transAxes, fontsize=11, va="top", color="#1E293B",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#CBD5E1", alpha=0.95))

    # Expand axes limits to give label room
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    ax.set_xlim(cur_xlim[0] - 0.06, cur_xlim[1] + 0.28)
    ax.set_ylim(cur_ylim[0] - 0.04, cur_ylim[1] + 0.06)

# Shared group legend on the first scatter
handles5 = [mpatches.Patch(color=GRP_COL[g], label=g) for g in GROUPS]
scatter_axes[0].legend(handles=handles5, loc="lower right", fontsize=10,
                       framealpha=0.95, edgecolor="#CBD5E1", title="Group",
                       title_fontsize=10)

# ── Bottom-right: significance panel ─────────────────────────────────────────
ax_sig = axes5[1, 1]
ax_sig.set_facecolor("#F8FAFC")

for st, method, color in [(st_g, "GSVA", METHOD_COL["GSVA"]),
                            (st_s, "ssGSEA", METHOD_COL["ssGSEA"]),
                            (st_k, "singscore", METHOD_COL["singscore"])]:
    omni = omnibus_pvals(st)
    xs = [-np.log10(omni.loc[pw, "p_value"]) if pw in omni.index else 0
          for pw in pw_order]
    ys = list(range(len(pw_order)))
    ax_sig.plot(xs, ys, "o-", color=color, lw=2, ms=9,
                markeredgecolor="white", markeredgewidth=0.8,
                label=method, alpha=0.88, zorder=3)

# Alternating row shading
for i in range(len(pw_order)):
    if i % 2 == 0:
        ax_sig.axhspan(i - 0.45, i + 0.45, color="#F1F5F9", zorder=1)

ax_sig.axvline(-np.log10(0.05), color="#EF4444", lw=1.2, ls="--", zorder=2)
ax_sig.axvline(-np.log10(0.10), color="#F59E0B", lw=1.0, ls=":",  zorder=2)

# Threshold labels at top
ax_sig.text(-np.log10(0.05) + 0.03, len(pw_order) - 0.3,
            "p = .05", fontsize=9, color="#EF4444", va="top")
ax_sig.text(-np.log10(0.10) + 0.03, len(pw_order) - 0.3,
            "p = .10", fontsize=9, color="#F59E0B", va="top")

ax_sig.set_yticks(range(len(pw_order)))
ax_sig.set_yticklabels(pw_labels, fontsize=10)
ax_sig.set_xlabel("–log₁₀(Kruskal-Wallis p)", fontsize=12, labelpad=8)
ax_sig.set_title("Omnibus Significance per Method", fontsize=13,
                  fontweight="bold", color="#1E293B", pad=10)
ax_sig.spines[["top","right"]].set_visible(False)
ax_sig.legend(fontsize=10, framealpha=0.95, edgecolor="#CBD5E1",
              loc="lower right")
ax_sig.tick_params(axis="y", length=0, labelsize=10)
ax_sig.set_ylim(-0.7, len(pw_order) - 0.3)

fig5.suptitle("Method Comparison — Median Scores & Omnibus Significance",
              fontsize=16, fontweight="bold", y=1.01, color="#1E293B")
fig5.tight_layout(rect=[0, 0, 1, 0.99])
fig5.subplots_adjust(hspace=0.38, wspace=0.32)
fig5.savefig(OUT / "fig5_method_comparison.pdf", bbox_inches="tight",
             facecolor=fig5.get_facecolor())
fig5.savefig(OUT / "fig5_method_comparison.png", dpi=180, bbox_inches="tight",
             facecolor=fig5.get_facecolor())
plt.close(fig5)
print("Figure 5 done")

print("\nAll figures written to:", OUT)
