#!/usr/bin/env python3
"""
build_gsea_3d_network.py  —  GSEA 3D Network Explorer
DESeq2 CSV  →  GSEA prerank  →  self-contained interactive 3D HTML

Usage:
  python3 build_gsea_3d_network.py --deseq2-csv results.csv --gene-sets H --output out.html

Dependencies: gseapy, pandas, scipy, numpy  (pip install gseapy pandas scipy numpy)
"""

import re
import json
import csv
import math
import argparse
import sys
import tempfile
import shutil
from pathlib import Path
from collections import defaultdict

# ============================================================
# Constants
# ============================================================

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

# Okabe-Ito colorblind-safe palette
CLUSTER_COLORS_CB = [
    '#E69F00', '#56B4E9', '#009E73', '#F0E442',
    '#0072B2', '#D55E00', '#CC79A7', '#999999',
    '#44AA99', '#882255',
]

EXCEL_DATE_RE = re.compile(r'^\d+-[A-Z][a-z]{2}$')

# ============================================================
# Step 1: Read DESeq2 CSV
# ============================================================

def read_deseq2(csv_path):
    """Read DESeq2 CSV. Returns (de_table, rnk_dict)."""
    de_table = {}
    rnk_dict = {}
    n_mangled = 0

    def safe_float(v):
        if v is None:
            return None
        s = str(v).strip()
        if s.lower() in ('na', 'nan', ''):
            return None
        try:
            f = float(s)
            return None if math.isnan(f) else f
        except (ValueError, TypeError):
            return None

    with open(csv_path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            symbol = ''
            for cand in ('', 'Unnamed: 0', '...1', 'gene', 'Gene', 'symbol'):
                v = row.get(cand, '').strip()
                if v:
                    symbol = v
                    break
            if not symbol:
                keys = list(row.keys())
                if keys:
                    symbol = row[keys[0]].strip()
            if not symbol:
                continue

            if EXCEL_DATE_RE.match(symbol):
                n_mangled += 1
                continue

            lfc = safe_float(row.get('log2FC_shrunken')) or safe_float(row.get('log2FC'))
            padj = safe_float(row.get('padj'))
            stat = safe_float(row.get('stat'))

            key = symbol.upper()
            de_table[key] = {'symbol': symbol, 'log2FC': lfc, 'padj': padj, 'stat': stat}
            rank_val = stat if stat is not None else lfc
            if rank_val is not None:
                rnk_dict[key] = rank_val

    if n_mangled:
        print(f"      WARNING: {n_mangled} Excel date-mangled gene names skipped (e.g. '1-Mar')")

    rnk_dict = dict(sorted(rnk_dict.items(), key=lambda x: x[1], reverse=True))
    return de_table, rnk_dict


# ============================================================
# Step 2: Run GSEA prerank
# ============================================================

# Map user-friendly shorthand names to Enrichr library names (gseapy 1.x).
# Special sentinel 'KEGG_GGA' is handled separately via the KEGG REST API.
_GENE_SET_ALIASES = {
    'H':              'MSigDB_Hallmark_2020',
    'HALLMARK':       'MSigDB_Hallmark_2020',
    'C2_CP_KEGG':     'KEGG_2019_Human',
    'KEGG':           'KEGG_2019_Human',
    'KEGG_GGA':       'KEGG_GGA',           # Gallus gallus — fetched from KEGG REST API
    'C5_GO_BP':       'GO_Biological_Process_2021',
    'C5_GO_MF':       'GO_Molecular_Function_2021',
    'C5_GO_CC':       'GO_Cellular_Component_2021',
    'REACTOME':       'Reactome_2022',
    'C2_CP_REACTOME': 'Reactome_2022',
    'BIOCARTA':       'BioCarta_2016',
    'WIKIPATHWAYS':   'WikiPathways_2019_Human',
    'WP':             'WikiPathways_2019_Human',
}


def _fetch_kegg_gga():
    """
    Fetch KEGG pathway gene sets for Gallus gallus (gga) via the KEGG REST API.

    Strategy: use KO (KEGG Orthology) numbers as the primary bridge to HGNC symbols.
    The KEGG /list/gga gene names are unreliable — many unannotated entries appear
    as 'CDS', 'LOC123456', etc. which don't match the HGNC-style symbols in a
    chicken DESeq2 CSV.  By going gga_gene → KO → human_gene → HGNC_symbol we
    get gene names that match the DESeq2 ranked list.

    Four bulk API requests, no per-pathway loops:
      1. link/ko/gga        gga gene ID → KO number
      2. link/hsa/ko        KO number   → human gene ID
      3. list/hsa           human gene ID → HGNC symbol
      4. link/gga/pathway   all pathway–gene associations
      5. list/pathway/gga   pathway ID  → human-readable name
      (+6. list/gga         direct gga symbol as last-resort fallback)
    """
    import urllib.request

    def kegg_get(endpoint):
        url = f'https://rest.kegg.jp/{endpoint}'
        with urllib.request.urlopen(url, timeout=60) as resp:
            return resp.read().decode('utf-8')

    def is_valid_symbol(sym):
        """True if sym looks like a real HGNC gene symbol."""
        _BAD = {'CDS', 'ORF', 'NA', 'PREDICTED', 'UNKNOWN', 'HYPOTHETICAL'}
        return bool(sym and
                    not sym.startswith('[') and
                    sym.upper() not in _BAD and
                    re.match(r'^[A-Za-z][A-Za-z0-9\-_.]{1,}$', sym))

    # 1. gga gene ID → KO  (primary path to HGNC symbols)
    print('      Fetching gga → KO orthology mappings…')
    gid2ko = {}
    for line in kegg_get('link/ko/gga').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        gid2ko[parts[0].strip()] = parts[1].strip()

    # 2. KO → human gene IDs
    print('      Fetching KO → human gene mappings…')
    ko2hsa = defaultdict(list)
    for line in kegg_get('link/hsa/ko').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        ko2hsa[parts[0].strip()].append(parts[1].strip())

    # 3. Human gene ID → HGNC symbol
    # Format: hsa:ID \t CDS \t CHROMOSOME \t SYMBOL, aliases; description
    # Gene symbol is in column index 3, not 1.
    print('      Fetching human HGNC gene symbols…')
    hsa2sym = {}
    for line in kegg_get('list/hsa').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 4:
            continue
        sym = parts[3].split(';')[0].split(',')[0].strip()
        if is_valid_symbol(sym):
            hsa2sym[parts[0].strip()] = sym.upper()

    # 4. Direct gga symbol as last-resort fallback
    # Format: gga:ID \t CDS \t CHROMOSOME \t SYMBOL, aliases; description
    # Gene symbol is in column index 3, not 1.
    print('      Fetching direct Gallus gallus gene symbols (fallback)…')
    id2sym_direct = {}
    for line in kegg_get('list/gga').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 4:
            continue
        sym = parts[3].split(';')[0].split(',')[0].strip()
        if is_valid_symbol(sym):
            id2sym_direct[parts[0].strip()] = sym.upper()

    # Build final gid → HGNC symbol: KO-based first, direct fallback second
    id2sym = {}
    n_ko = n_direct = 0
    for gid in set(gid2ko) | set(id2sym_direct):
        ko = gid2ko.get(gid)
        sym = None
        if ko:
            for hid in ko2hsa.get(ko, []):
                sym = hsa2sym.get(hid)
                if sym:
                    break
        if sym:
            id2sym[gid] = sym
            n_ko += 1
        elif gid in id2sym_direct:
            id2sym[gid] = id2sym_direct[gid]
            n_direct += 1

    print(f'      Mapped {len(id2sym)} genes: '
          f'{n_ko} via KO→human, {n_direct} direct gga symbol only')

    # 5. Pathway names
    print('      Fetching KEGG pathway list…')
    pw_names = {}
    for line in kegg_get('list/pathway/gga').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        pw_id = parts[0].strip().replace('path:', '')
        pw_names[pw_id] = parts[1].strip()

    # 6. All pathway–gene links in one request
    print('      Fetching pathway–gene associations…')
    pw_genes = defaultdict(set)
    for line in kegg_get('link/gga/pathway').strip().split('\n'):
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        pw_id = parts[0].strip().replace('path:', '')
        sym   = id2sym.get(parts[1].strip())
        if sym:
            pw_genes[pw_id].add(sym)

    gs_dict = {pw_names.get(pid, pid): list(genes)
               for pid, genes in pw_genes.items() if genes}
    print(f'      Built {len(gs_dict)} KEGG Gallus gallus gene sets')
    return gs_dict


def _resolve_gene_sets(gene_sets_arg):
    """
    Resolve gene_sets argument to a value gseapy.prerank() accepts.
    - If it's a local .gmt file path: return the path string.
    - If it's a known shorthand (H, KEGG, …): download from Enrichr, cache locally, and return dict.
    - Otherwise: return as-is (let gseapy handle it).
    """
    import gseapy as gp
    import json as _json

    # Local file takes priority
    p = Path(gene_sets_arg)
    if p.exists() and p.suffix.lower() == '.gmt':
        print(f"      Using local GMT file: {p}")
        return str(p)

    # Resolve alias
    resolved = _GENE_SET_ALIASES.get(gene_sets_arg.upper(), gene_sets_arg)
    if resolved != gene_sets_arg:
        print(f"      Resolved '{gene_sets_arg}' → '{resolved}'")

    # Cache path (Feature 4)
    cache_dir = Path.home() / '.gsea3d_cache'
    cache_dir.mkdir(exist_ok=True)
    safe_name = re.sub(r'[^A-Za-z0-9_\-]', '_', resolved)
    cache_file = cache_dir / f'{safe_name}.json'

    if cache_file.exists():
        try:
            gs_dict = _json.loads(cache_file.read_text(encoding='utf-8'))
            print(f"      Loaded {len(gs_dict)} gene sets from cache: {cache_file}")
            return gs_dict
        except Exception:
            pass  # corrupted cache, re-download

    # KEGG Gallus gallus — fetched from KEGG REST API, not Enrichr
    if resolved == 'KEGG_GGA':
        try:
            gs_dict = _fetch_kegg_gga()
            if gs_dict:
                try:
                    cache_file.write_text(_json.dumps(gs_dict, ensure_ascii=False), encoding='utf-8')
                    print(f"      Cached to {cache_file}")
                except Exception as e:
                    print(f"      Warning: could not write cache: {e}")
                return gs_dict
        except Exception as e:
            sys.exit(f"ERROR: Could not fetch KEGG gga from KEGG REST API: {e}")

    try:
        print(f"      Downloading gene sets from Enrichr: {resolved} …")
        gs_dict = gp.get_library(name=resolved, organism='Human')
        if gs_dict:
            print(f"      Downloaded {len(gs_dict)} gene sets")
            try:
                serializable = {k: list(v) for k, v in gs_dict.items()}
                cache_file.write_text(_json.dumps(serializable, ensure_ascii=False), encoding='utf-8')
                print(f"      Cached to {cache_file}")
            except Exception as e:
                print(f"      Warning: could not write cache: {e}")
            return gs_dict
    except Exception as e:
        print(f"      Could not download via get_library ('{e}'); passing name directly to prerank")

    # Fall back: pass name directly, let gseapy try
    return gene_sets_arg


def run_gsea(rnk_dict, gene_sets_arg, permutations, min_size, max_size):
    """Run gseapy.prerank. Returns normalized results DataFrame."""
    try:
        import gseapy as gp
        import pandas as pd
    except ImportError as e:
        print(f"\nERROR: Could not import a required package: {e}")
        print(f"  Python executable: {sys.executable}")
        print(f"  This Python does NOT have gseapy/pandas installed.")
        print(f"  Fix: open a terminal and run:")
        print(f"    conda activate rnaseq   # (or whichever env has gseapy)")
        print(f"    python3 {Path(__file__).name}")
        print(f"  OR install into the Python Spyder is using:")
        print(f"    {sys.executable} -m pip install gseapy pandas scipy numpy")
        sys.exit(1)

    rnk = pd.Series(rnk_dict)
    gene_sets_resolved = _resolve_gene_sets(gene_sets_arg)

    tmpdir = tempfile.mkdtemp(prefix='gsea3d_')
    print(f"      Ranked list: {len(rnk)} genes  |  Running {permutations} permutations…")

    # Build kwargs — gseapy 1.x deprecated 'processes', use 'threads'
    import inspect
    prerank_sig = inspect.signature(gp.prerank)
    thread_kwarg = 'threads' if 'threads' in prerank_sig.parameters else 'processes'

    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets_resolved,
            outdir=tmpdir,
            permutation_num=permutations,
            seed=42,
            min_size=min_size,
            max_size=max_size,
            no_plot=True,
            **{thread_kwarg: 1},
            verbose=False,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # Extract DataFrame
    if hasattr(pre_res, 'res2d') and pre_res.res2d is not None:
        df = pre_res.res2d.copy()
    elif hasattr(pre_res, 'results') and pre_res.results is not None:
        df = pre_res.results.copy()
    else:
        df = pre_res.copy() if hasattr(pre_res, 'copy') else pre_res

    if hasattr(df, 'index') and str(df.index.name) in ('Term', 'Name'):
        df = df.reset_index()

    # Normalize column names across gseapy versions
    # In gseapy 1.x: columns are ['Name', 'Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', ...]
    # where 'Name' == 'prerank' (input label, not the gene set name) and 'Term' is the gene set name.
    # We must NOT rename 'Name' to 'Term' when 'Term' already exists.
    has_term_col = 'Term' in df.columns
    rename = {}
    for col in df.columns:
        cl = col.lower().strip().replace('-', '_').replace(' ', '_')
        if cl == 'term':
            rename[col] = 'Term'
        elif cl == 'name' and not has_term_col:
            rename[col] = 'Term'  # Only if no 'Term' column present
        elif cl == 'es':
            rename[col] = 'ES'
        elif cl == 'nes':
            rename[col] = 'NES'
        elif 'nom' in cl and 'p' in cl:
            rename[col] = 'NOM_pval'
        elif 'fdr' in cl:
            rename[col] = 'FDR_qval'
        elif 'fwer' in cl:
            rename[col] = 'FWER_pval'
        elif cl in ('lead_genes', 'leading_edge', 'leadingedge', 'genes'):
            rename[col] = 'Lead_genes'
    df = df.rename(columns=rename)

    for col in ('ES', 'NES', 'NOM_pval', 'FDR_qval', 'FWER_pval'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    return df, pre_res


# ============================================================
# Step 3: Filter significant gene sets
# ============================================================

def filter_genesets(df, fdr_cutoff, max_sets):
    """Filter and rank gene sets. Returns filtered DataFrame."""
    import pandas as pd

    if 'Term' not in df.columns:
        print(f"      WARNING: 'Term' column not found. Columns: {list(df.columns)}")

    df = df.copy()
    if 'NES' in df.columns:
        df['_abs_NES'] = df['NES'].abs()
        df = df.sort_values('_abs_NES', ascending=False)

    n_total = len(df)

    if 'FDR_qval' in df.columns:
        sig = df[df['FDR_qval'] < fdr_cutoff].copy()
        n_sig = len(sig)
        print(f"      {n_total} gene sets tested  |  {n_sig} pass FDR < {fdr_cutoff}")
        if n_sig < 3:
            print(f"      Relaxing: taking top 20 by |NES|")
            sig = df.head(20).copy()
    else:
        print(f"      WARNING: FDR_qval column not found; taking top 20 by |NES|")
        sig = df.head(20).copy()
        n_sig = len(sig)

    if len(sig) > max_sets:
        print(f"      Capping at {max_sets} gene sets (top by |NES|)")
        sig = sig.head(max_sets).copy()

    print(f"\n      Top gene sets by NES:")
    for i, (_, row) in enumerate(sig.head(10).iterrows()):
        term = str(row.get('Term', '?'))
        nes = row.get('NES', float('nan'))
        fdr = row.get('FDR_qval', float('nan'))
        try:
            print(f"        {i+1:2d}. {term:<45s}  NES={nes:+.3f}  FDR={fdr:.4f}")
        except (ValueError, TypeError):
            print(f"        {i+1:2d}. {term}")

    return sig


# ============================================================
# Step 4: Cluster gene sets by Jaccard similarity
# ============================================================

def cluster_genesets(sig_df, target_min=4, target_max=8):
    """Cluster gene sets by Jaccard similarity. Returns (assignments, gene_sets_genes)."""
    try:
        import numpy as np
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage, fcluster
    except ImportError:
        sys.exit("ERROR: pip install scipy numpy")

    terms = [str(t) for t in sig_df['Term']]
    n = len(terms)

    gene_sets_genes = {}
    for _, row in sig_df.iterrows():
        term = str(row['Term'])
        raw = row.get('Lead_genes', '') or ''
        if isinstance(raw, str):
            genes = frozenset(g.strip().upper() for g in raw.split(';') if g.strip())
        else:
            genes = frozenset()
        gene_sets_genes[term] = genes

    if n <= target_max:
        return {t: i for i, t in enumerate(terms)}, gene_sets_genes

    dist_matrix = np.ones((n, n))
    np.fill_diagonal(dist_matrix, 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            g1, g2 = gene_sets_genes[terms[i]], gene_sets_genes[terms[j]]
            inter = len(g1 & g2)
            union = len(g1 | g2)
            d = 1.0 - inter / union if union > 0 else 1.0
            dist_matrix[i, j] = dist_matrix[j, i] = d

    condensed = squareform(dist_matrix)
    Z = linkage(condensed, method='average')

    best_labels = None
    for t in [0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]:
        labels = fcluster(Z, t=t, criterion='distance')
        nc = len(set(labels))
        if target_min <= nc <= target_max:
            best_labels = labels
            break

    if best_labels is None:
        best_labels = fcluster(Z, t=min(target_max, n), criterion='maxclust')

    assignments = {terms[i]: int(best_labels[i]) - 1 for i in range(n)}
    print(f"      Clustered {n} gene sets into {len(set(assignments.values()))} clusters")
    return assignments, gene_sets_genes


# ============================================================
# Feature 7: Extract Running Enrichment Score curves
# ============================================================

def extract_res_curves(pre_res, sig_df):
    """Extract RES curves for significant gene sets (downsampled to max 150 points).
    Returns dict: {term: [float, ...]}
    """
    curves = {}
    terms = set(str(t) for t in sig_df['Term'])
    results_attr = None
    if hasattr(pre_res, 'results') and isinstance(pre_res.results, dict):
        results_attr = pre_res.results
    if results_attr is None:
        return curves
    for term_key, term_data in results_attr.items():
        if str(term_key) not in terms:
            continue
        res = None
        if isinstance(term_data, dict):
            res = term_data.get('RES') or term_data.get('res')
        elif hasattr(term_data, 'RES'):
            res = term_data.RES
        if res is None:
            continue
        try:
            res_list = [float(x) for x in res]
        except (TypeError, ValueError):
            continue
        n = len(res_list)
        if n > 150:
            step = n / 150
            res_list = [res_list[int(i * step)] for i in range(150)]
        curves[str(term_key)] = [round(x, 4) for x in res_list]
    return curves


# ============================================================
# Feature 5: Compute Jaccard similarity matrix
# ============================================================

def compute_jaccard_matrix(sig_df, gene_sets_genes):
    """Compute N×N Jaccard similarity matrix for all significant gene sets.
    Returns (terms_list, flat_matrix) where flat_matrix is NxN as a flat list (distance = 1 - jaccard).
    """
    terms = [str(t) for t in sig_df['Term']]
    n = len(terms)
    flat = [0.0] * (n * n)
    for i in range(n):
        for j in range(n):
            if i == j:
                flat[i * n + j] = 0.0
            elif j > i:
                g1 = gene_sets_genes.get(terms[i], frozenset())
                g2 = gene_sets_genes.get(terms[j], frozenset())
                inter = len(g1 & g2)
                union = len(g1 | g2)
                dist = 1.0 - inter / union if union > 0 else 1.0
                flat[i * n + j] = round(dist, 4)
                flat[j * n + i] = round(dist, 4)
    return terms, flat


# ============================================================
# Step 5: Build 3D network
# ============================================================

def _pretty_term(term):
    """Convert HALLMARK_APOPTOSIS → Apoptosis, or leave human-readable names as-is."""
    # If the term already has spaces and no underscores it's already human-readable
    if ' ' in term and '_' not in term:
        return term.title() if term == term.upper() else term
    s = term
    for pfx in ('HALLMARK_', 'KEGG_', 'REACTOME_', 'GO_BP_', 'GO_MF_', 'GO_CC_',
                'BIOCARTA_', 'WP_', 'PID_', 'NABA_', 'C2_CP_', 'HP_', 'GOCC_',
                'GOBP_', 'GOMF_', 'GTRD_'):
        s = s.replace(pfx, '')
    return s.replace('_', ' ').title()


def _safe_id(term):
    """Create a safe string ID from a gene set term name (for use in JS/HTML)."""
    return re.sub(r'[^A-Za-z0-9_\-]', '_', term)


def build_network(sig_df, cluster_assignments, gene_sets_genes, de_table, p_cutoff,
                  res_curves=None, nes_b_lk=None, fdr_b_lk=None, de_table_b=None):
    """Build nodes + links. Returns (nodes, links, clusters_meta).

    Optional args for Features 1 and 7:
      res_curves  – dict {term: [float]} RES curve data
      nes_b_lk    – dict {term: float} NES values from comparison 2
      fdr_b_lk    – dict {term: float} FDR values from comparison 2
      de_table_b  – dict {gene_upper: {...}} DE data for comparison 2
    """
    if res_curves is None:
        res_curves = {}
    if nes_b_lk is None:
        nes_b_lk = {}
    if fdr_b_lk is None:
        fdr_b_lk = {}
    if de_table_b is None:
        de_table_b = {}
    # Lookup tables from sig_df
    fdr_lk = {}; nes_lk = {}; es_lk = {}; nom_lk = {}; fwer_lk = {}; lead_lk = {}
    for _, row in sig_df.iterrows():
        t = str(row['Term'])
        fdr_lk[t] = float(row.get('FDR_qval', 1) or 1)
        nes_lk[t] = float(row.get('NES', 0) or 0)
        es_lk[t] = float(row.get('ES', 0) or 0)
        nom_lk[t] = float(row.get('NOM_pval', 1) or 1)
        fwer_lk[t] = float(row.get('FWER_pval', 1) or 1)
        raw = row.get('Lead_genes', '') or ''
        if isinstance(raw, str):
            lead_lk[t] = [g.strip().upper() for g in raw.split(';') if g.strip()]
        else:
            lead_lk[t] = []

    # Cluster → Z position
    cluster_terms = defaultdict(list)
    for term, cid in cluster_assignments.items():
        cluster_terms[cid].append(term)
    unique_cids = sorted(cluster_terms.keys())
    cluster_z_map = {cid: i * 150 for i, cid in enumerate(unique_cids)}

    # Cluster representative = lowest FDR member
    cluster_repr = {cid: min(terms, key=lambda t: fdr_lk.get(t, 1))
                    for cid, terms in cluster_terms.items()}

    # Build clusters metadata
    clusters_meta = {}
    for i, cid in enumerate(unique_cids):
        repr_term = cluster_repr[cid]
        name = _pretty_term(repr_term)
        if len(name) > 32:
            name = name[:30] + '…'
        clusters_meta[f'cluster_{cid}'] = {
            'display_name': name,
            'color': CLUSTER_COLORS[i % len(CLUSTER_COLORS)],
            'z': cluster_z_map[cid],
            'repr_term': repr_term,
            'n_sets': len(cluster_terms[cid]),
        }

    # Layout gene-set nodes: circular within each cluster Z-layer
    gs_positions = {}
    for cid, terms in cluster_terms.items():
        z = cluster_z_map[cid]
        n = len(terms)
        radius = 160 if n > 1 else 0
        terms_sorted = sorted(terms, key=lambda t: abs(nes_lk.get(t, 0)), reverse=True)
        for i, term in enumerate(terms_sorted):
            angle = 2 * math.pi * i / max(n, 1)
            x = round(radius * math.cos(angle), 1)
            y = round(radius * math.sin(angle), 1)
            gs_positions[term] = (x, y, z)

    # Accumulate gene node positions (centroid of connected gene-set nodes)
    gene_pos_acc = defaultdict(lambda: [0.0, 0.0, 0.0, 0])
    gene_membership = defaultdict(set)
    for term, genes in lead_lk.items():
        if term not in gs_positions:
            continue
        gx, gy, gz = gs_positions[term]
        for gene in genes:
            gene_membership[gene].add(term)
            a = gene_pos_acc[gene]
            a[0] += gx; a[1] += gy; a[2] += gz; a[3] += 1

    # Build nodes
    nodes = []
    all_ids = set()

    # Gene-set nodes
    for term, (gx, gy, gz) in gs_positions.items():
        cid = cluster_assignments[term]
        gs_id = 'gs_' + _safe_id(term)   # sanitize for JS/JSON safety
        lbl_short = _pretty_term(term)
        if len(lbl_short) > 36:
            lbl_short = lbl_short[:34] + '…'
        lead_genes_list = lead_lk.get(term, [])

        # Feature 1: comparison NES/FDR
        nes_b = nes_b_lk.get(term)
        fdr_b = fdr_b_lk.get(term)
        # Feature 7: RES curve
        res_curve = res_curves.get(term, [])
        nodes.append({
            'id': gs_id,
            'type': 'geneset',
            'label': term.replace('_', ' '),
            'label_short': lbl_short,
            'term': term,
            'NES': round(nes_lk.get(term, 0), 4),
            'FDR': round(fdr_lk.get(term, 1), 6),
            'ES': round(es_lk.get(term, 0), 4),
            'NOM_pval': round(nom_lk.get(term, 1), 6),
            'FWER_pval': round(fwer_lk.get(term, 1), 6),
            'NES_B': round(nes_b, 4) if nes_b is not None else None,
            'FDR_B': round(fdr_b, 6) if fdr_b is not None else None,
            'res_curve': res_curve,
            'cluster_id': f'cluster_{cid}',
            'cluster_name': clusters_meta[f'cluster_{cid}']['display_name'],
            'n_lead': len(lead_genes_list),
            'lead_genes': lead_genes_list[:60],
            'x': gx, 'y': gy, 'z': gz,
            'x_initial': gx, 'y_initial': gy, 'z_target': gz,
            'is_hub': False,
        })
        all_ids.add(gs_id)

    # Gene nodes
    for gene, acc in gene_pos_acc.items():
        count = acc[3]
        if count == 0:
            continue
        base_x = acc[0] / count
        base_y = acc[1] / count
        base_z = acc[2] / count
        # Small deterministic jitter to prevent exact overlap
        jx = round(base_x + (hash(gene) % 60 - 30), 1)
        jy = round(base_y + (hash(gene[::-1]) % 60 - 30), 1)
        jz = round(base_z, 1)

        de = de_table.get(gene.upper())
        if de:
            lfc = de['log2FC']
            padj = de['padj']
            status = ('Present; significant' if padj is not None and padj < p_cutoff
                      else 'Present; padj NA' if padj is None
                      else 'Present; not significant')
        else:
            lfc = padj = None
            status = 'Missing gene'

        # Feature 1: comparison 2 DE data
        de_b = de_table_b.get(gene.upper()) if de_table_b else None
        lfc_b = de_b['log2FC'] if de_b else None
        padj_b = de_b['padj'] if de_b else None

        terms_list = sorted(gene_membership[gene])
        # Use safe IDs for gene_sets (matches gs_* node IDs in the network)
        gene_sets_ids = [_safe_id(t) for t in terms_list]
        gene_sets_labels = terms_list  # human-readable display names
        clusters_list = list(set(f'cluster_{cluster_assignments[t]}'
                                 for t in terms_list if t in cluster_assignments))
        nodes.append({
            'id': gene,
            'type': 'gene',
            'label': gene,
            'gene_sets': gene_sets_ids,      # safe IDs (used for nodeMap lookup)
            'gene_sets_labels': gene_sets_labels,  # human-readable (for display)
            'clusters': clusters_list,
            'log2FC': round(lfc, 4) if lfc is not None else None,
            'padj': padj,
            'log2FC_B': round(lfc_b, 4) if lfc_b is not None else None,
            'padj_B': padj_b,
            'node_status': status,
            'is_hub': len(terms_list) > 1,
            'x': jx, 'y': jy, 'z': jz,
            'x_initial': round(base_x, 1),
            'y_initial': round(base_y, 1),
            'z_target': jz,
        })
        all_ids.add(gene)

    # Build links
    links = []
    seen = set()

    for term, genes in lead_lk.items():
        if term not in cluster_assignments:
            continue
        gs_id = 'gs_' + _safe_id(term)
        cid = cluster_assignments[term]
        for gene in genes:
            if gene not in all_ids:
                continue
            key = (gs_id, gene)
            if key not in seen:
                seen.add(key)
                links.append({
                    'source': gs_id,
                    'target': gene,
                    'edge_type': 'membership',
                    'cluster_id': f'cluster_{cid}',
                })

    # Gene-set ↔ Gene-set overlap (Jaccard > 0.1)
    terms_list = list(gs_positions.keys())
    for i in range(len(terms_list)):
        for j in range(i + 1, len(terms_list)):
            t1, t2 = terms_list[i], terms_list[j]
            g1 = gene_sets_genes.get(t1, frozenset())
            g2 = gene_sets_genes.get(t2, frozenset())
            inter = len(g1 & g2)
            union = len(g1 | g2)
            if union > 0 and inter / union > 0.1:
                links.append({
                    'source': 'gs_' + _safe_id(t1),
                    'target': 'gs_' + _safe_id(t2),
                    'edge_type': 'overlap',
                    'jaccard': round(inter / union, 3),
                    'cluster_id': None,
                })

    return nodes, links, clusters_meta


def compute_max_lfc(nodes):
    vals = [abs(n['log2FC']) for n in nodes
            if n.get('type') == 'gene' and n.get('log2FC') is not None]
    return round(max(vals), 2) if vals else 3.0


# ============================================================
# Step 7: Generate HTML
# ============================================================

GSEA_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<script src="https://unpkg.com/three@0.160.0/build/three.min.js"></script>
<script src="https://unpkg.com/3d-force-graph@1"
  onerror="document.getElementById('loading-msg').textContent='ERROR: Failed to load 3d-force-graph from CDN. Check internet connection.';"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    --bg: #0a0a1a; --text: #e0e0f0; --panel-bg: rgba(14,14,35,0.92);
    --border: rgba(100,100,200,0.2); --heading: #7080c0; --btn-bg: rgba(40,50,120,0.7);
    --btn-text: #c0ccff; --input-bg: rgba(30,30,60,0.8); --tooltip-bg: rgba(10,10,30,0.95);
    --legend-bg: rgba(12,12,32,0.90); --muted: #888; --title-color: #c0c8ff;
    --title-bar-bg: rgba(10,10,40,0.95); --detail-bg: rgba(10,10,30,0.96);
    background: var(--bg); color: var(--text);
    font-family: 'Segoe UI', system-ui, sans-serif; font-size: 13px; overflow: hidden;
  }}
  body.light-mode {{
    --bg: #f4f4f8; --text: #1a1a2e; --panel-bg: rgba(240,240,248,0.95);
    --border: rgba(100,100,160,0.25); --heading: #4050a0; --btn-bg: rgba(200,210,240,0.8);
    --btn-text: #2a3a7a; --input-bg: rgba(220,225,240,0.9); --tooltip-bg: rgba(245,245,255,0.97);
    --legend-bg: rgba(240,240,248,0.95); --muted: #666; --title-color: #2a3a7a;
    --title-bar-bg: rgba(230,230,248,0.97); --detail-bg: rgba(240,240,248,0.97);
  }}
  #title-bar {{
    position: fixed; top: 0; left: 0; right: 0; height: 38px;
    background: var(--title-bar-bg); border-bottom: 1px solid var(--border);
    display: flex; align-items: center; padding: 0 16px; z-index: 100; gap: 16px;
  }}
  #title-bar h1 {{ font-size: 14px; font-weight: 600; color: var(--title-color); white-space: nowrap; }}
  #title-bar .stats {{ font-size: 11px; color: var(--muted); white-space: nowrap; }}
  #layout {{ display: flex; position: fixed; top: 38px; left: 0; right: 0; bottom: 0; }}
  #control-panel {{
    width: 248px; min-width: 248px; background: var(--panel-bg);
    border-right: 1px solid var(--border); overflow-y: auto; padding: 10px; z-index: 50;
    scrollbar-width: thin; scrollbar-color: #333 transparent;
  }}
  #graph-container {{ flex: 1; position: relative; overflow: hidden; }}
  .panel-section {{ margin-bottom: 14px; }}
  .panel-section h3 {{
    font-size: 10px; text-transform: uppercase; letter-spacing: 1px; color: var(--heading);
    margin-bottom: 7px; border-bottom: 1px solid var(--border); padding-bottom: 4px;
  }}
  .check-row {{
    display: flex; align-items: center; gap: 7px; margin-bottom: 5px; cursor: pointer;
  }}
  .check-row:hover {{ background: rgba(100,120,255,0.08); border-radius: 3px; }}
  .check-row input[type=checkbox] {{ cursor: pointer; accent-color: #6a8aff; flex-shrink:0; }}
  .pathway-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
  .check-label {{ font-size: 12px; color: var(--text); line-height: 1.3; }}
  #search-box {{
    width: 100%; background: var(--input-bg); border: 1px solid var(--border);
    border-radius: 4px; color: var(--text); padding: 5px 8px; font-size: 12px;
    outline: none; margin-bottom: 6px;
  }}
  #search-box:focus {{ border-color: rgba(100,150,255,0.6); }}
  #search-results {{ font-size: 11px; color: var(--muted); min-height: 16px; }}
  .btn {{
    display: block; width: 100%; margin-bottom: 5px; padding: 5px 8px;
    background: var(--btn-bg); border: 1px solid var(--border);
    border-radius: 4px; color: var(--btn-text); font-size: 11px; cursor: pointer;
    text-align: left; transition: background 0.15s;
  }}
  .btn:hover {{ background: rgba(60,80,160,0.8); }}
  .kbd {{
    display: inline-block; background: rgba(60,60,100,0.6);
    border: 1px solid rgba(100,100,180,0.4); border-radius: 3px;
    padding: 1px 5px; font-size: 10px; font-family: monospace; color: #aabf;
  }}
  #tooltip {{
    position: fixed; background: var(--tooltip-bg); border: 1px solid rgba(100,120,255,0.4);
    border-radius: 6px; padding: 10px 13px; pointer-events: none; z-index: 200; display: none;
    min-width: 200px; max-width: 300px; font-size: 12px; line-height: 1.6;
    box-shadow: 0 4px 20px rgba(0,0,0,0.6); color: var(--text);
  }}
  #tooltip .tt-gene {{ font-size: 14px; font-weight: 700; color: var(--title-color); margin-bottom: 6px; }}
  #tooltip .tt-row {{ display: flex; gap: 8px; }}
  #tooltip .tt-key {{ color: var(--heading); min-width: 90px; flex-shrink:0; }}
  #tooltip .tt-val {{ color: var(--text); word-break: break-word; }}
  #tooltip .tt-sig {{ color: #ff9966; font-weight: 600; }}
  #tooltip .tt-up {{ color: #ff8888; }}
  #tooltip .tt-down {{ color: #88aaff; }}
  #detail-panel {{
    position: fixed; right: 0; top: 38px; bottom: 0; width: 290px;
    background: var(--detail-bg); border-left: 1px solid var(--border);
    padding: 12px; overflow-y: auto; z-index: 80; display: none;
    font-size: 12px; line-height: 1.6; color: var(--text); pointer-events: auto;
  }}
  #detail-panel h2 {{ font-size: 15px; color: var(--title-color); margin-bottom: 10px; word-break: break-word; }}
  #detail-panel .dp-key {{ color: var(--heading); }}
  #detail-panel .dp-val {{ color: var(--text); }}
  #detail-panel .dp-section {{
    margin-top: 10px; color: var(--heading); font-size: 10px;
    text-transform: uppercase; letter-spacing: 1px;
    border-bottom: 1px solid var(--border); padding-bottom: 3px; margin-bottom: 6px;
  }}
  #detail-panel .neighbor-row {{ display: flex; justify-content: space-between; margin-bottom: 3px; }}
  #close-detail {{ float: right; cursor: pointer; color: var(--muted); font-size: 16px; line-height: 1; }}
  #legend {{
    position: fixed; bottom: 16px; right: 16px; background: var(--legend-bg);
    border: 1px solid var(--border); border-radius: 8px; padding: 10px 14px;
    z-index: 60; min-width: 175px; font-size: 11px; color: var(--text);
  }}
  #legend-header {{ display: flex; align-items: center; justify-content: space-between; cursor: pointer; user-select: none; }}
  #legend h4 {{ font-size: 10px; text-transform: uppercase; letter-spacing: 1px; color: var(--heading); margin: 0; }}
  #legend-toggle-icon {{ font-size: 12px; color: var(--heading); line-height: 1; margin-left: 8px; }}
  #legend-body {{ margin-top: 8px; }}
  #legend-body.collapsed {{ display: none; }}
  .lfc-bar {{
    width: 140px; height: 12px;
    background: linear-gradient(to right, #2166ac, #f7f7f7, #b2182b);
    border-radius: 3px; margin: 4px 0 2px;
  }}
  .nes-bar {{
    width: 140px; height: 12px;
    background: linear-gradient(to right, #2166ac, #f7f7f7, #b2182b);
    border-radius: 3px; margin: 4px 0 2px;
  }}
  .lfc-labels {{ display: flex; justify-content: space-between; color: var(--muted); font-size: 10px; }}
  .legend-row {{ display: flex; align-items: center; gap: 7px; margin-bottom: 4px; }}
  #loading {{
    position: fixed; inset: 0; background: var(--bg);
    display: flex; flex-direction: column; align-items: center; justify-content: center;
    z-index: 999; gap: 16px;
  }}
  .spinner {{
    width: 40px; height: 40px; border: 3px solid rgba(100,130,255,0.2);
    border-top-color: #6a8aff; border-radius: 50%; animation: spin 0.8s linear infinite;
  }}
  @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
  #loading p {{ color: var(--muted); font-size: 13px; }}
  /* Gene-sets subsection */
  #gene-sets-list {{
    max-height: 220px; overflow-y: auto; margin-top: 4px;
    scrollbar-width: thin; scrollbar-color: #333 transparent;
  }}
  .gs-item {{
    padding: 3px 6px; margin-bottom: 2px; border-radius: 3px; cursor: pointer;
    font-size: 11px; border-left: 3px solid transparent;
    transition: background 0.1s;
  }}
  .gs-item:hover {{ background: rgba(100,120,255,0.12); }}
  .gs-item.active {{ background: rgba(100,120,255,0.2); }}
  .gs-item .gs-nes {{ float: right; font-size: 10px; }}
  .gs-cluster-header {{
    font-size: 10px; color: var(--heading); text-transform: uppercase;
    letter-spacing: 0.5px; margin: 6px 0 3px; padding-left: 4px;
  }}
  @media (max-width: 768px) {{
    #control-panel {{ display: none; }}
    #control-panel.open {{
      display: block; position: fixed; z-index: 150; top: 38px; bottom: 0; left: 0;
    }}
    #hamburger {{ display: block !important; }}
  }}
  #hamburger {{
    display: none; background: none; border: none; color: var(--title-color);
    font-size: 20px; cursor: pointer; padding: 0 6px; line-height: 1;
  }}
  /* Feature 2: FDR/NES sliders */
  input[type=range] {{ width: 100%; accent-color: #6a8aff; }}
  .slider-row {{ margin-bottom: 8px; }}
  .slider-row label {{ display: flex; justify-content: space-between; font-size: 11px; color: var(--text); margin-bottom: 2px; }}
  .slider-row label span {{ color: #6a8aff; font-weight: 600; }}
  /* Feature 3: Bar chart panel */
  #bar-chart-panel {{
    position: fixed; right: 0; top: 38px; bottom: 0; width: 260px;
    background: var(--detail-bg); border-left: 1px solid var(--border);
    z-index: 79; display: none; overflow-y: auto; padding: 10px;
  }}
  #bar-chart-panel h3 {{ font-size: 11px; color: var(--heading); text-transform: uppercase;
    letter-spacing: 1px; margin-bottom: 8px; }}
  #btn-bar-chart {{
    background: var(--btn-bg); border: 1px solid var(--border); border-radius: 4px;
    color: var(--btn-text); font-size: 11px; cursor: pointer; padding: 4px 10px;
    margin-left: auto; white-space: nowrap;
  }}
  #btn-bar-chart:hover {{ background: rgba(60,80,160,0.8); }}
  /* Feature 1: subtitle for comparison */
  #title-subtitle {{ font-size: 10px; color: var(--muted); white-space: nowrap; }}
</style>
</head>
<body>

<div id="loading">
  <div class="spinner"></div>
  <p id="loading-msg">Initializing GSEA 3D Network…</p>
</div>

<div id="title-bar">
  <button id="hamburger" onclick="document.getElementById('control-panel').classList.toggle('open')" title="Toggle controls">&#9776;</button>
  <h1>{title}</h1>
  <span id="title-subtitle" style="display:none;"></span>
  <span class="stats" id="graph-stats"></span>
  <span class="stats">Drag to rotate &nbsp;|&nbsp; Scroll to zoom &nbsp;|&nbsp; Right-drag to pan</span>
  <button id="btn-bar-chart" title="Toggle bar chart panel">&#128202; Bar Chart</button>
</div>

<div id="layout">
  <div id="control-panel">
    <div class="panel-section">
      <h3>Search</h3>
      <input id="search-box" type="text" placeholder="Gene symbol or pathway name…" autocomplete="off">
      <div id="search-results"></div>
    </div>

    <div class="panel-section">
      <h3>Gene Set Clusters</h3>
      <div id="cluster-checkboxes"></div>
      <button class="btn" id="btn-show-all" style="margin-top:4px;">Show All</button>
      <button class="btn" id="btn-isolate-active" style="color:#ffcc66;">Focus on checked ↑</button>
    </div>

    <div class="panel-section">
      <h3>Gene Sets &nbsp;<span id="gs-toggle-btn" style="cursor:pointer;font-size:10px;color:var(--muted);">▶ expand</span></h3>
      <div id="gene-sets-list" style="display:none;"></div>
    </div>

    <div class="panel-section">
      <h3>Filters</h3>
      <div class="check-row">
        <input type="checkbox" id="filter-significant">
        <label for="filter-significant" class="check-label">Significant genes only (padj &lt; {p_cutoff})</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="filter-hubs">
        <label for="filter-hubs" class="check-label">Hub genes only (≥2 gene sets) <span class="kbd">H</span></label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="filter-up-sets">
        <label for="filter-up-sets" class="check-label">Upregulated sets only (NES &gt; 0)</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="filter-down-sets">
        <label for="filter-down-sets" class="check-label">Downregulated sets only (NES &lt; 0)</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="filter-hide-overlap" checked>
        <label for="filter-hide-overlap" class="check-label">Hide overlap edges</label>
      </div>
      <div class="slider-row" style="margin-top:8px;">
        <label>FDR threshold <span id="fdr-slider-val">1.00</span></label>
        <input type="range" id="fdr-slider" min="0" max="1" step="0.01" value="1">
      </div>
      <div class="slider-row">
        <label>Min |NES| <span id="nes-slider-val">0.0</span></label>
        <input type="range" id="nes-slider" min="0" max="3" step="0.1" value="0">
      </div>
    </div>

    <div class="panel-section">
      <h3>Display</h3>
      <div class="check-row">
        <input type="checkbox" id="toggle-labels" checked>
        <label for="toggle-labels" class="check-label">Show labels <span class="kbd">L</span></label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="toggle-light">
        <label for="toggle-light" class="check-label">Light background</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="toggle-colorblind">
        <label for="toggle-colorblind" class="check-label">Colorblind-safe colors</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="toggle-layer-labels" checked>
        <label for="toggle-layer-labels" class="check-label">Show cluster labels</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="toggle-color-cluster">
        <label for="toggle-color-cluster" class="check-label">Color sets by cluster</label>
      </div>
    </div>

    <div class="panel-section">
      <h3>Camera &nbsp;<span class="kbd">R</span>=reset</h3>
      <button class="btn" id="btn-reset-camera">Reset Camera</button>
      <button class="btn" id="btn-reheat">Reheat Simulation <span class="kbd">Space</span></button>
    </div>

    <div class="panel-section">
      <h3>Export</h3>
      <button class="btn" id="btn-screenshot">Screenshot PNG</button>
      <button class="btn" id="btn-export-csv">Export Node CSV</button>
      <button class="btn" id="btn-export-edge-csv">Export Edge CSV</button>
    </div>

    <div class="panel-section" style="color:#555;font-size:10px;line-height:1.7;">
      Organism: <em style="color:#666">Gallus gallus</em><br>
      p-cutoff: {p_cutoff}<br>
      ◆ Gene Set node &nbsp; ● Gene node<br>
      Double-click node to reset view<br>
      <span class="kbd">H</span> hub-only &nbsp;<span class="kbd">L</span> labels &nbsp;<span class="kbd">Space</span> pause
    </div>
  </div>
  <div id="graph-container"></div>
</div>

<div id="tooltip"></div>

<div id="detail-panel">
  <span id="close-detail">✕</span>
  <div id="detail-content"></div>
</div>

<div id="bar-chart-panel">
  <h3>&#128202; Gene Set Rankings (NES)</h3>
  <canvas id="bar-chart-canvas"></canvas>
</div>

<div id="legend">
  <div id="legend-header" onclick="toggleLegend()">
    <h4>Legend</h4>
    <span id="legend-toggle-icon">▲</span>
  </div>
  <div id="legend-body">
    <div style="margin-top:8px;color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:4px;">Gene log2FC</div>
    <div class="lfc-bar"></div>
    <div class="lfc-labels"><span>−{max_lfc}</span><span>0</span><span>+{max_lfc}</span></div>
    <div style="margin-top:8px;color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:4px;">Gene Set NES</div>
    <div class="nes-bar"></div>
    <div class="lfc-labels"><span>Down</span><span>0</span><span>Up</span></div>
    <div style="margin-top:10px;">
      <div class="legend-row">
        <div style="width:14px;height:14px;transform:rotate(45deg);background:#aaa;border:2px solid #fff;flex-shrink:0;"></div>
        <span>Gene Set (◆ size ∝ |NES|)</span>
      </div>
      <div class="legend-row">
        <div style="width:16px;height:16px;border-radius:50%;background:#aaa;border:2px solid #fff;flex-shrink:0;"></div>
        <span>Sig. gene (padj &lt; {p_cutoff})</span>
      </div>
      <div class="legend-row">
        <div style="width:10px;height:10px;border-radius:50%;background:#666;flex-shrink:0;"></div>
        <span>Non-sig. / missing gene</span>
      </div>
      <div class="legend-row">
        <div style="width:14px;height:14px;border-radius:50%;background:#aaa;border:2px solid #ffdd88;flex-shrink:0;"></div>
        <span>Hub gene (≥2 gene sets)</span>
      </div>
    </div>
    <div style="margin-top:10px;" id="legend-clusters"></div>
  </div>
</div>

<script>
// ============================================================
// Embedded Data
// ============================================================
const DATA = {data_json};

// ============================================================
// Constants & State
// ============================================================
const CLUSTER_COLORS = {{}};
const CLUSTER_COLORS_CB = {{}};
const CB_PALETTE = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#999999','#44AA99','#882255'];
DATA.cluster_keys.forEach((k, i) => {{
  CLUSTER_COLORS[k] = DATA.clusters[k].color;
  CLUSTER_COLORS_CB[k] = CB_PALETTE[i % CB_PALETTE.length];
}});
const CLUSTER_Z = {{}};
DATA.cluster_keys.forEach(k => {{ CLUSTER_Z[k] = DATA.clusters[k].z; }});

const state = {{
  visibleClusters: new Set(DATA.cluster_keys),
  showLabels: true,
  colorblindMode: false,
  showLayerLabels: true,
  filterSignificant: false,
  filterHubs: false,
  filterUpSets: false,
  filterDownSets: false,
  hideOverlapEdges: true,
  highlightNodes: new Set(),
  highlightLinks: new Set(),
  selectedNode: null,
  paused: false,
  searchMatches: new Set(),
  fdrThreshold: 1.0,    // Feature 2
  nesThreshold: 0.0,    // Feature 2
  colorByCluster: false, // Feature 6
  barChartVisible: false, // Feature 3
}};

const _layerLabelSprites = [];

function getClusterColor(clusterId) {{
  return (state.colorblindMode ? CLUSTER_COLORS_CB : CLUSTER_COLORS)[clusterId] || '#888';
}}

// ============================================================
// Color helpers
// ============================================================
const LFC_CAP = {max_lfc};

function lfc2color(lfc) {{
  if (lfc === null || lfc === undefined || isNaN(lfc)) return '#555555';
  const t = Math.max(-1, Math.min(1, lfc / LFC_CAP));
  if (t <= 0) {{
    const s = -t;
    return `rgb(${{Math.round(33  + (247-33) *(1-s))}},${{Math.round(102 + (247-102)*(1-s))}},${{Math.round(172 + (247-172)*(1-s))}})`;
  }} else {{
    return `rgb(${{Math.round(247 + (178-247)*t)}},${{Math.round(247 + (24 -247)*t)}},${{Math.round(247 + (43 -247)*t)}})`;
  }}
}}

function nes2color(nes) {{
  if (nes === null || nes === undefined || isNaN(nes)) return '#888888';
  const t = Math.max(-1, Math.min(1, nes / 3));
  if (t <= 0) {{
    const s = -t;
    return `rgb(${{Math.round(247 + (33-247)*s)}},${{Math.round(247 + (102-247)*s)}},${{Math.round(247 + (172-247)*s)}})`;
  }} else {{
    return `rgb(${{Math.round(247 + (178-247)*t)}},${{Math.round(247 + (24-247)*t)}},${{Math.round(247 + (43-247)*t)}})`;
  }}
}}

function nodeVal(node) {{
  if (node.type === 'geneset') {{
    const nesAbs = Math.min(Math.abs(node.NES || 0), 3);
    return 10 + nesAbs * 4;
  }}
  switch (node.node_status) {{
    case 'Present; significant':     return 12;
    case 'Present; not significant': return 5;
    case 'Present; padj NA':         return 5;
    default:                         return 2;
  }}
}}

// ============================================================
// Visibility
// ============================================================
function isGeneSetVisible(node) {{
  if (!state.visibleClusters.has(node.cluster_id)) return false;
  if (state.filterUpSets   && (node.NES || 0) <= 0) return false;
  if (state.filterDownSets && (node.NES || 0) >= 0) return false;
  if (state.fdrThreshold < 1 && node.FDR > state.fdrThreshold) return false;
  if (state.nesThreshold > 0 && Math.abs(node.NES || 0) < state.nesThreshold) return false;
  return true;
}}

function isGeneNodeVisible(node) {{
  if (state.filterSignificant && node.node_status !== 'Present; significant') return false;
  if (state.filterHubs && !node.is_hub) return false;
  return (node.clusters || []).some(c => state.visibleClusters.has(c)) &&
         (node.gene_sets || []).some(gsId => {{
           const gsNode = nodeMap['gs_' + gsId];
           return gsNode && isGeneSetVisible(gsNode);
         }});
}}

function isNodeVisible(node) {{
  return node.type === 'geneset' ? isGeneSetVisible(node) : isGeneNodeVisible(node);
}}

function fmt(v, d) {{
  d = d !== undefined ? d : 3;
  if (v === null || v === undefined || (typeof v === 'number' && isNaN(v))) return 'NA';
  return (+v).toFixed(d);
}}

function negLog10(p) {{
  if (p === null || p === undefined || p <= 0) return 'NA';
  return (-Math.log10(p)).toFixed(2);
}}

// Pre-build lookup maps
const nodeMap = Object.fromEntries(DATA.nodes.map(n => [n.id, n]));
const neighborNodes = {{}};
const nodeLinks = {{}};
DATA.links.forEach(lk => {{
  const s = lk.source, t = lk.target;
  (neighborNodes[s] = neighborNodes[s] || new Set()).add(t);
  (neighborNodes[t] = neighborNodes[t] || new Set()).add(s);
  (nodeLinks[s] = nodeLinks[s] || []).push(lk);
  (nodeLinks[t] = nodeLinks[t] || []).push(lk);
}});

// ============================================================
// Graph data filter
// ============================================================
function getGraphData() {{
  const visNodes = DATA.nodes.filter(isNodeVisible);
  const visIds = new Set(visNodes.map(n => n.id));
  const visLinks = DATA.links.filter(lk => {{
    if (state.hideOverlapEdges && lk.edge_type === 'overlap') return false;
    const s = typeof lk.source === 'object' ? lk.source.id : lk.source;
    const t = typeof lk.target === 'object' ? lk.target.id : lk.target;
    return visIds.has(s) && visIds.has(t);
  }});
  return {{ nodes: visNodes, links: visLinks }};
}}

// ============================================================
// Three.js node builder
// ============================================================
function makeTextSprite(T, text, fontSize) {{
  const canvas = document.createElement('canvas');
  const scale = 3;
  canvas.width = 512 * scale;
  canvas.height = 64 * scale;
  const ctx = canvas.getContext('2d');
  ctx.scale(scale, scale);
  ctx.clearRect(0, 0, 512, 64);
  ctx.font = `bold ${{fontSize}}px 'Segoe UI', sans-serif`;
  const lightMode = document.body.classList.contains('light-mode');
  ctx.fillStyle = lightMode ? 'rgba(26,26,46,0.95)' : 'rgba(220,225,255,0.95)';
  ctx.fillText(text, 3, fontSize + 8);
  const texture = new T.CanvasTexture(canvas);
  texture.needsUpdate = true;
  const mat = new T.SpriteMaterial({{ map: texture, transparent: true, depthWrite: false, depthTest: false }});
  const sprite = new T.Sprite(mat);
  sprite.scale.set(90, 22, 1);
  return sprite;
}}

function getThree() {{ return (typeof THREE !== 'undefined') ? THREE : null; }}

function buildNodeObject(node) {{
  const T = getThree();
  if (!T) return null;

  const group = new T.Group();

  if (node.type === 'geneset') {{
    // --- Gene-set node: diamond (octahedron) ---
    const nesAbs = Math.min(Math.abs(node.NES || 0), 3);
    const size = 8 + nesAbs * 3;

    // Feature 6: color by cluster vs NES
    if (state.colorByCluster) {{
      const clusterFillColor = getClusterColor(node.cluster_id);
      const matC = new T.MeshLambertMaterial({{
        color: new T.Color(clusterFillColor), transparent: true, opacity: 0.92
      }});
      group.add(new T.Mesh(new T.OctahedronGeometry(size), matC));
    }} else if (DATA.has_comparison && node.NES_B !== null && node.NES_B !== undefined) {{
      // Feature 1: split upper/lower halves by NES_A / NES_B
      const geo = new T.OctahedronGeometry(size);
      const nonIdx = geo.toNonIndexed();
      const posArr = nonIdx.attributes.position.array;
      const nVerts = posArr.length / 3;
      const colors = new Float32Array(nVerts * 3);
      const colorA = new T.Color(nes2color(node.NES));
      const colorB = new T.Color(nes2color(node.NES_B));
      for (let i = 0; i < nVerts; i += 3) {{
        // centroid Y of triangle
        const cy = (posArr[i*3+1] + posArr[(i+1)*3+1] + posArr[(i+2)*3+1]) / 3;
        const c = cy >= 0 ? colorA : colorB;
        for (let k = 0; k < 3; k++) {{
          colors[(i+k)*3]   = c.r;
          colors[(i+k)*3+1] = c.g;
          colors[(i+k)*3+2] = c.b;
        }}
      }}
      nonIdx.setAttribute('color', new T.BufferAttribute(colors, 3));
      const matSplit = new T.MeshLambertMaterial({{
        vertexColors: true, transparent: true, opacity: 0.92
      }});
      group.add(new T.Mesh(nonIdx, matSplit));
    }} else {{
      const nesColor = nes2color(node.NES);
      const mat = new T.MeshLambertMaterial({{
        color: new T.Color(nesColor), transparent: true, opacity: 0.92
      }});
      group.add(new T.Mesh(new T.OctahedronGeometry(size), mat));
    }}

    // Cluster color outline ring
    const clusterColor = getClusterColor(node.cluster_id);
    const ringMat = new T.MeshBasicMaterial({{
      color: new T.Color(clusterColor), transparent: true, opacity: 0.70
    }});
    group.add(new T.Mesh(new T.TorusGeometry(size + 1.5, 0.8, 8, 24), ringMat));

    // Search match ring
    if (state.searchMatches.has(node.id)) {{
      const sRingMat = new T.MeshBasicMaterial({{ color: 0xffffff, transparent: true, opacity: 0.95 }});
      group.add(new T.Mesh(new T.TorusGeometry(size + 4, 1.0, 8, 24), sRingMat));
    }}

    if (state.showLabels) {{
      const sprite = makeTextSprite(T, node.label_short || node.label, 34);
      sprite.position.set(size + 2, size + 2, 0);
      group.add(sprite);
    }}
  }} else {{
    // --- Gene node: sphere (same as reference tool) ---
    const radius = nodeVal(node) * 0.7;
    const colorStr = lfc2color(node.log2FC);
    const opacity = node.node_status === 'Missing gene' ? 0.45 : 0.88;
    const mat = new T.MeshLambertMaterial({{
      color: new T.Color(colorStr), transparent: true, opacity
    }});
    group.add(new T.Mesh(new T.SphereGeometry(radius, 12, 12), mat));

    if (node.is_hub) {{
      const ringMat = new T.MeshBasicMaterial({{ color: 0xffdd44, transparent: true, opacity: 0.75 }});
      group.add(new T.Mesh(new T.TorusGeometry(radius + 2, 0.6, 8, 24), ringMat));
    }}

    if (state.searchMatches.has(node.id)) {{
      const sRingMat = new T.MeshBasicMaterial({{ color: 0xffffff, transparent: true, opacity: 0.95 }});
      group.add(new T.Mesh(new T.TorusGeometry(radius + 3.5, 0.8, 8, 24), sRingMat));
    }}

    if (state.showLabels) {{
      const sprite = makeTextSprite(T, node.label || node.id, node.is_hub ? 36 : 28);
      sprite.position.set(radius + 1, radius + 1, 0);
      group.add(sprite);

      if (node.node_status === 'Present; significant' && node.log2FC !== null && node.log2FC !== undefined) {{
        const ac = document.createElement('canvas');
        const sc2 = 3;
        ac.width = 64 * sc2; ac.height = 64 * sc2;
        const actx = ac.getContext('2d');
        actx.scale(sc2, sc2);
        actx.clearRect(0, 0, 64, 64);
        const glyph = node.log2FC > 0 ? '\u25b2' : '\u25bc';
        const glyphColor = node.log2FC > 0 ? '#ff6666' : '#6699ff';
        actx.font = 'bold 36px serif';
        actx.fillStyle = glyphColor;
        actx.textAlign = 'center';
        actx.fillText(glyph, 32, 44);
        const arrowTex = new T.CanvasTexture(ac);
        arrowTex.needsUpdate = true;
        const arrowMat = new T.SpriteMaterial({{ map: arrowTex, transparent: true, depthWrite: false, depthTest: false }});
        const arrowSprite = new T.Sprite(arrowMat);
        arrowSprite.scale.set(20, 20, 1);
        arrowSprite.position.set(radius + 1, -(radius + 4), 0);
        group.add(arrowSprite);
      }}
    }}
  }}

  return group;
}}

// ============================================================
// Graph initialization
// ============================================================
let Graph = null;

function buildGraph() {{
  const container = document.getElementById('graph-container');

  Graph = ForceGraph3D()(container)
    .backgroundColor('#0a0a1a')
    .graphData(getGraphData())
    .nodeId('id')
    .nodeVal(n => nodeVal(n))
    .nodeColor(n => {{
      if (state.highlightNodes.size > 0 && !state.highlightNodes.has(n.id))
        return 'rgba(60,60,80,0.2)';
      if (n.type === 'geneset')
        return state.colorByCluster ? getClusterColor(n.cluster_id) : nes2color(n.NES);
      return lfc2color(n.log2FC);
    }})
    .nodeOpacity(0.9)
    .nodeThreeObject(buildNodeObject)
    .nodeThreeObjectExtend(false)
    .linkSource('source')
    .linkTarget('target')
    .linkColor(lk => {{
      if (state.highlightLinks.size > 0 && !state.highlightLinks.has(lk))
        return 'rgba(40,40,60,0.10)';
      if (lk.edge_type === 'overlap') return 'rgba(130,130,180,0.30)';
      return getClusterColor(lk.cluster_id) || '#666';
    }})
    .linkWidth(lk => {{
      if (lk.edge_type === 'overlap') return 0.7;
      const tgtId = typeof lk.target === 'object' ? lk.target.id : lk.target;
      const tgt = nodeMap[tgtId];
      const lfc = (tgt && tgt.log2FC !== null) ? Math.abs(tgt.log2FC) : 0;
      const w = 1 + Math.min(lfc / LFC_CAP, 1);
      return state.highlightLinks.size > 0 && !state.highlightLinks.has(lk) ? w * 0.2 : w;
    }})
    .linkDirectionalArrowLength(0)
    .linkDirectionalParticles(lk => lk.edge_type === 'overlap' ? 0 : 2)
    .linkDirectionalParticleWidth(1.5)
    .linkDirectionalParticleColor(lk => getClusterColor(lk.cluster_id) || '#888')
    .linkDirectionalParticleSpeed(0.004)
    .onNodeClick(onNodeClick)
    .onNodeHover(onNodeHover)
    .onBackgroundClick(resetHighlight)
    .d3AlphaDecay(0.02)
    .d3VelocityDecay(0.3)
    .warmupTicks(80)
    .cooldownTicks(200)
    .onEngineStop(() => {{
      document.getElementById('loading').style.display = 'none';
      updateStats();
    }});

  try {{
    const chargeFn = Graph.d3Force('charge');
    if (chargeFn) chargeFn.strength(-150);
    const linkFn = Graph.d3Force('link');
    if (linkFn) linkFn.distance(60).strength(0.5);
    const zFn = Graph.d3Force('z');
    if (zFn) zFn.z(n => n.z_target).strength(n => n.type === 'geneset' ? 0.5 : 0.05);
    const xFn = Graph.d3Force('x');
    if (xFn) xFn.x(n => n.x_initial).strength(n => n.type === 'geneset' ? 0.3 : 0.03);
    const yFn = Graph.d3Force('y');
    if (yFn) yFn.y(n => n.y_initial).strength(n => n.type === 'geneset' ? 0.3 : 0.03);
  }} catch(e) {{ console.warn('Force setup error:', e); }}

  try {{
    const T2 = getThree();
    if (T2) {{
      const scene = Graph.scene();
      scene.add(new T2.AmbientLight(0xffffff, 0.6));
      const dirLight = new T2.DirectionalLight(0xffffff, 0.8);
      dirLight.position.set(300, 300, 300);
      scene.add(dirLight);

      _layerLabelSprites.length = 0;
      DATA.cluster_keys.forEach(ck => {{
        const cdata = DATA.clusters[ck];
        if (!cdata) return;
        const zPos = cdata.z;
        const color = getClusterColor(ck);
        const lCanvas = document.createElement('canvas');
        lCanvas.width = 1024; lCanvas.height = 96;
        const lctx = lCanvas.getContext('2d');
        lctx.clearRect(0, 0, 1024, 96);
        lctx.font = 'bold 52px "Segoe UI", sans-serif';
        lctx.fillStyle = color;
        lctx.globalAlpha = 0.40;
        lctx.fillText(cdata.display_name, 8, 66);
        lctx.globalAlpha = 1.0;
        const lTex = new T2.CanvasTexture(lCanvas);
        lTex.needsUpdate = true;
        const lMat = new T2.SpriteMaterial({{ map: lTex, transparent: true, depthWrite: false, depthTest: false }});
        const lSprite = new T2.Sprite(lMat);
        lSprite.scale.set(420, 60, 1);
        lSprite.position.set(-600, 420, zPos);
        lSprite.visible = state.showLayerLabels;
        scene.add(lSprite);
        _layerLabelSprites.push(lSprite);
      }});
    }}
  }} catch(e) {{ console.warn('Scene setup error:', e); }}

  updateStats();
}}

// ============================================================
// Tooltip
// ============================================================
const tooltipEl = document.getElementById('tooltip');

function onNodeHover(node) {{
  if (!node) {{ tooltipEl.style.display = 'none'; return; }}

  if (node.type === 'geneset') {{
    const nesSign = (node.NES || 0) > 0 ? 'tt-up' : 'tt-down';
    const cmpRows = DATA.has_comparison ? `
      <div class="tt-row"><span class="tt-key">NES (A):</span><span class="tt-val ${{nesSign}}">${{fmt(node.NES, 3)}}</span></div>
      <div class="tt-row"><span class="tt-key">NES (B):</span><span class="tt-val ${{(node.NES_B||0)>0?'tt-up':'tt-down'}}">${{node.NES_B !== null && node.NES_B !== undefined ? fmt(node.NES_B, 3) : 'NA'}}</span></div>
      <div class="tt-row"><span class="tt-key">FDR (A):</span><span class="tt-val tt-sig">${{fmt(node.FDR, 4)}}</span></div>
      <div class="tt-row"><span class="tt-key">FDR (B):</span><span class="tt-val">${{node.FDR_B !== null && node.FDR_B !== undefined ? fmt(node.FDR_B, 4) : 'NA'}}</span></div>
    ` : `
      <div class="tt-row"><span class="tt-key">NES:</span><span class="tt-val ${{nesSign}}">${{fmt(node.NES, 3)}}</span></div>
      <div class="tt-row"><span class="tt-key">ES:</span><span class="tt-val">${{fmt(node.ES, 3)}}</span></div>
      <div class="tt-row"><span class="tt-key">FDR q-val:</span><span class="tt-val tt-sig">${{fmt(node.FDR, 4)}}</span></div>
    `;
    tooltipEl.innerHTML = `
      <div class="tt-gene">${{node.label_short || node.label}}</div>
      ${{cmpRows}}
      <div class="tt-row"><span class="tt-key">Lead genes:</span><span class="tt-val">${{node.n_lead}}</span></div>
      <div class="tt-row"><span class="tt-key">Cluster:</span><span class="tt-val">${{node.cluster_name}}</span></div>
    `;
  }} else {{
    const lfc = fmt(node.log2FC);
    const padj = fmt(node.padj, 4);
    const nl10 = negLog10(node.padj);
    const lfcClass = node.log2FC > 0 ? 'tt-up' : (node.log2FC < 0 ? 'tt-down' : '');
    const sigClass = node.node_status === 'Present; significant' ? 'tt-sig' : '';
    const geneSets = (node.gene_sets_labels || node.gene_sets || []).map(g => typeof g === 'string' ? g.replace(/_/g, ' ') : g).join(', ');
    tooltipEl.innerHTML = `
      <div class="tt-gene">${{node.label || node.id}}</div>
      <div class="tt-row"><span class="tt-key">log2FC:</span><span class="tt-val ${{lfcClass}}">${{lfc}}</span></div>
      <div class="tt-row"><span class="tt-key">padj:</span><span class="tt-val ${{sigClass}}">${{padj}}</span></div>
      <div class="tt-row"><span class="tt-key">−log10(padj):</span><span class="tt-val ${{sigClass}}">${{nl10}}</span></div>
      <div class="tt-row"><span class="tt-key">Status:</span><span class="tt-val">${{node.node_status}}</span></div>
      <div class="tt-row"><span class="tt-key">Gene sets:</span><span class="tt-val">${{geneSets || 'none'}}</span></div>
    `;
  }}
  tooltipEl.style.display = 'block';
}}

document.addEventListener('mousemove', e => {{
  if (tooltipEl.style.display !== 'none') {{
    let left = e.clientX + 16;
    let top = Math.min(e.clientY + 16, window.innerHeight - 240);
    const w = tooltipEl.offsetWidth || 300;
    if (left + w > window.innerWidth - 8) left = e.clientX - w - 16;
    tooltipEl.style.left = left + 'px';
    tooltipEl.style.top = top + 'px';
  }}
}});

// ============================================================
// Node click & detail panel
// ============================================================
function onNodeClick(node) {{
  if (!node) return;
  state.selectedNode = node;
  const neighbors = neighborNodes[node.id] || new Set();
  state.highlightNodes = new Set([node.id, ...neighbors]);
  state.highlightLinks = new Set(nodeLinks[node.id] || []);
  Graph.graphData(getGraphData());
  const dist = 220;
  Graph.cameraPosition(
    {{ x: (node.x||0) + dist*0.7, y: (node.y||0) + dist*0.3, z: (node.z||0) + dist }},
    node, 800
  );
  showDetailPanel(node);
}}

function resetHighlight() {{
  state.highlightNodes.clear();
  state.highlightLinks.clear();
  state.selectedNode = null;
  if (Graph) Graph.graphData(getGraphData());
  document.getElementById('detail-panel').style.display = 'none';
  // Feature 3: restore bar chart panel position
  const bcp = document.getElementById('bar-chart-panel');
  if (bcp && state.barChartVisible) bcp.style.right = '0';
}}

function showDetailPanel(node) {{
  const panel = document.getElementById('detail-panel');
  const content = document.getElementById('detail-content');
  panel.style.display = 'block';
  // Feature 3: shift bar chart panel if open
  const bcp = document.getElementById('bar-chart-panel');
  if (bcp && state.barChartVisible) bcp.style.right = '290px';

  if (node.type === 'geneset') {{
    // Gene-set detail
    const nesClass = (node.NES||0) > 0 ? 'color:#ff8888' : 'color:#88aaff';
    const leadHtml = (node.lead_genes || []).map(g => {{
      const gn = nodeMap[g];
      if (!gn) return `<div class="neighbor-row"><span style="color:#aaa">${{g}}</span><span style="color:#555">not in graph</span></div>`;
      const lfc = fmt(gn.log2FC);
      const col = gn.log2FC > 0 ? '#ff8888' : (gn.log2FC < 0 ? '#88aaff' : '#bbb');
      return `<div class="neighbor-row">
        <span style="cursor:pointer;text-decoration:underline;color:#88aaff;"
              onclick="onNodeClick(nodeMap['${{g}}'])">${{g}}</span>
        <span style="color:${{col}}">${{lfc}}</span>
      </div>`;
    }}).join('');
    const msigdbUrl = `https://www.gsea-msigdb.org/gsea/msigdb/cards/${{encodeURIComponent(node.term)}}`;
    const cmpHtml = DATA.has_comparison ? `
      <div><span class="dp-key">NES (${{DATA.title}}):</span> <span style="${{nesClass}};font-weight:700">${{fmt(node.NES,3)}}</span></div>
      <div><span class="dp-key">NES (${{DATA.title2}}):</span> <span style="${{(node.NES_B||0)>0?'color:#ff8888':'color:#88aaff'}};font-weight:700">${{node.NES_B !== null && node.NES_B !== undefined ? fmt(node.NES_B,3) : 'NA'}}</span></div>
      <div><span class="dp-key">FDR (${{DATA.title}}):</span> <span class="dp-val">${{fmt(node.FDR,4)}}</span></div>
      <div><span class="dp-key">FDR (${{DATA.title2}}):</span> <span class="dp-val">${{node.FDR_B !== null && node.FDR_B !== undefined ? fmt(node.FDR_B,4) : 'NA'}}</span></div>
    ` : `
      <div><span class="dp-key">NES:</span> <span style="${{nesClass}};font-weight:700">${{fmt(node.NES,3)}}</span></div>
      <div><span class="dp-key">ES:</span> <span class="dp-val">${{fmt(node.ES,3)}}</span></div>
      <div><span class="dp-key">FDR q-val:</span> <span class="dp-val">${{fmt(node.FDR,4)}}</span></div>
      <div><span class="dp-key">NOM p-val:</span> <span class="dp-val">${{fmt(node.NOM_pval,4)}}</span></div>
      <div><span class="dp-key">FWER p-val:</span> <span class="dp-val">${{fmt(node.FWER_pval,4)}}</span></div>
    `;
    content.innerHTML = `
      <h2>${{node.label_short || node.label}}</h2>
      <div style="font-size:10px;color:var(--muted);margin-bottom:8px;word-break:break-all;">${{node.term}}</div>
      ${{cmpHtml}}
      <div><span class="dp-key">Lead-edge genes:</span> <span class="dp-val">${{node.n_lead}}</span></div>
      <div><span class="dp-key">Cluster:</span> <span class="dp-val">${{node.cluster_name}}</span></div>
      <div style="margin-top:6px;"><a href="${{msigdbUrl}}" target="_blank" style="color:#88aaff;font-size:11px;">MSigDB page ↗</a></div>
      ${{leadHtml ? `<div class="dp-section">Lead-edge genes (log2FC)</div>${{leadHtml}}` : ''}}
    `;
    // Feature 7: RES sparkline
    const resCurve = node.res_curve || [];
    if (resCurve.length > 0) {{
      content.innerHTML += `<div class="dp-section">Running Enrichment Score</div>
        <canvas id="res-sparkline" width="260" height="60" style="width:100%;height:60px;display:block;"></canvas>`;
      requestAnimationFrame(() => {{
        const canvas = document.getElementById('res-sparkline');
        if (!canvas) return;
        drawSparkline(canvas, resCurve, node.NES);
      }});
    }}
    // Feature 5: Jaccard similar gene sets
    const similar = getTopSimilar(node.term, 5);
    if (similar.length > 0) {{
      let simHtml = '<div class="dp-section">Most Similar Gene Sets</div>';
      similar.forEach(s => {{
        const simNode = DATA.nodes.find(n => n.type === 'geneset' && n.term === s.term);
        const barW = Math.round(s.jaccard * 100);
        const label = s.term.replace(/_/g,' ').length > 28 ? s.term.replace(/_/g,' ').slice(0,26)+'…' : s.term.replace(/_/g,' ');
        simHtml += `<div style="margin-bottom:4px;">
          <div style="display:flex;align-items:center;gap:4px;margin-bottom:2px;">
            <span style="cursor:pointer;color:#88aaff;font-size:11px;text-decoration:underline;"
                  onclick="if(nodeMap['${{simNode ? simNode.id : ''}}'])onNodeClick(nodeMap['${{simNode ? simNode.id : ''}}'])">${{label}}</span>
            <span style="color:var(--muted);font-size:10px;margin-left:auto;">${{s.jaccard.toFixed(2)}}</span>
          </div>
          <div style="height:5px;background:rgba(100,130,255,0.2);border-radius:2px;">
            <div style="height:5px;width:${{barW}}px;max-width:180px;background:#6a8aff;border-radius:2px;"></div>
          </div>
        </div>`;
      }});
      content.innerHTML += simHtml;
    }}
  }} else {{
    // Gene node detail
    const lfc = fmt(node.log2FC);
    const padj = fmt(node.padj, 4);
    const inLinks = DATA.links.filter(l => {{
      const t = typeof l.target === 'object' ? l.target.id : l.target;
      return t === node.id && l.edge_type === 'membership';
    }});
    const gsRows = (node.gene_sets || []).map(gsId => {{
      const gsn = nodeMap['gs_' + gsId];
      if (!gsn) return '';
      const nesStr = fmt(gsn.NES, 3);
      const col = (gsn.NES||0) > 0 ? '#ff8888' : '#88aaff';
      return `<div class="neighbor-row">
        <span style="cursor:pointer;text-decoration:underline;color:#88aaff;"
              onclick="onNodeClick(nodeMap['gs_${{gsId}}'])">${{gsn.label_short || gsId}}</span>
        <span style="color:${{col}}">NES ${{nesStr}}</span>
      </div>`;
    }}).join('');
    const lfcBHtml = DATA.has_comparison ? `
      <div><span class="dp-key">log2FC (B):</span> <span class="dp-val">${{fmt(node.log2FC_B)}}</span></div>
      <div><span class="dp-key">padj (B):</span> <span class="dp-val">${{fmt(node.padj_B, 4)}}</span></div>
    ` : '';
    content.innerHTML = `
      <h2>${{node.label || node.id}}</h2>
      <div><span class="dp-key">log2FC:</span> <span class="dp-val">${{lfc}}</span></div>
      <div><span class="dp-key">padj:</span> <span class="dp-val">${{padj}}</span></div>
      <div><span class="dp-key">−log10(padj):</span> <span class="dp-val">${{negLog10(node.padj)}}</span></div>
      ${{lfcBHtml}}
      <div><span class="dp-key">Status:</span> <span class="dp-val">${{node.node_status}}</span></div>
      <div><span class="dp-key">Hub gene:</span> <span class="dp-val">${{node.is_hub ? 'Yes (≥2 gene sets)' : 'No'}}</span></div>
      ${{gsRows ? `<div class="dp-section">Gene Sets (NES)</div>${{gsRows}}` : ''}}
    `;
  }}
}}

// ============================================================
// Update & stats
// ============================================================
let _recalcTimer = null;
function updateGraph() {{
  if (!Graph) return;
  const statsEl = document.getElementById('graph-stats');
  statsEl.textContent = 'Recalculating…';
  statsEl.style.color = '#ffcc66';
  clearTimeout(_recalcTimer);
  _recalcTimer = setTimeout(() => {{
    statsEl.style.color = '';
    Graph.graphData(getGraphData());
    updateStats();
  }}, 80);
}}

function updateStats() {{
  if (!Graph) return;
  const gd = Graph.graphData();
  document.getElementById('graph-stats').textContent =
    `${{(gd.nodes||[]).length}} nodes · ${{(gd.links||[]).length}} edges`;
}}

// ============================================================
// Controls
// ============================================================
function buildControls() {{
  // Cluster checkboxes
  const cbBox = document.getElementById('cluster-checkboxes');
  DATA.cluster_keys.forEach(ck => {{
    const cdata = DATA.clusters[ck];
    const color = getClusterColor(ck);
    const row = document.createElement('div');
    row.className = 'check-row';
    row.style.justifyContent = 'space-between';
    row.innerHTML = `
      <div style="display:flex;align-items:center;gap:7px;flex:1;min-width:0;">
        <input type="checkbox" id="ck-${{ck}}" checked>
        <div class="pathway-dot" data-cluster="${{ck}}" style="background:${{color}}"></div>
        <label for="ck-${{ck}}" class="check-label" style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"
               title="${{cdata.display_name}}">${{cdata.display_name}}</label>
      </div>
      <button class="isolate-btn" data-key="${{ck}}"
        style="font-size:9px;padding:2px 5px;margin-left:4px;flex-shrink:0;
               background:rgba(40,50,120,0.6);border:1px solid rgba(100,120,200,0.3);
               border-radius:3px;color:#aac;cursor:pointer;">only</button>
    `;
    row.querySelector('input').addEventListener('change', e => {{
      if (e.target.checked) state.visibleClusters.add(ck);
      else state.visibleClusters.delete(ck);
      updateGraph();
    }});
    row.querySelector('.isolate-btn').addEventListener('click', () => isolateCluster(ck));
    cbBox.appendChild(row);
  }});

  // Gene sets expandable list
  const gsListEl = document.getElementById('gene-sets-list');
  const gsToggleBtn = document.getElementById('gs-toggle-btn');
  let gsExpanded = false;
  gsToggleBtn.addEventListener('click', () => {{
    gsExpanded = !gsExpanded;
    gsListEl.style.display = gsExpanded ? 'block' : 'none';
    gsToggleBtn.textContent = gsExpanded ? '▼ collapse' : '▶ expand';
  }});

  // Populate gene sets grouped by cluster
  DATA.cluster_keys.forEach(ck => {{
    const cdata = DATA.clusters[ck];
    const color = getClusterColor(ck);
    const clusterNodes = DATA.nodes.filter(n => n.type === 'geneset' && n.cluster_id === ck);
    clusterNodes.sort((a, b) => Math.abs(b.NES||0) - Math.abs(a.NES||0));
    if (clusterNodes.length === 0) return;
    const header = document.createElement('div');
    header.className = 'gs-cluster-header';
    header.style.color = color;
    header.textContent = cdata.display_name;
    gsListEl.appendChild(header);
    clusterNodes.forEach(gsNode => {{
      const item = document.createElement('div');
      item.className = 'gs-item';
      item.style.borderLeftColor = color;
      const nesStr = (gsNode.NES||0) > 0 ? `<span style="color:#ff8888">+${{fmt(gsNode.NES,2)}}</span>`
                                          : `<span style="color:#88aaff">${{fmt(gsNode.NES,2)}}</span>`;
      item.innerHTML = `${{gsNode.label_short}} <span class="gs-nes">${{nesStr}}</span>`;
      item.addEventListener('click', () => {{
        const n = nodeMap[gsNode.id];
        if (n) onNodeClick(n);
        document.querySelectorAll('.gs-item').forEach(el => el.classList.remove('active'));
        item.classList.add('active');
      }});
      gsListEl.appendChild(item);
    }});
  }});

  document.getElementById('btn-show-all').addEventListener('click', () => {{
    DATA.cluster_keys.forEach(k => {{
      state.visibleClusters.add(k);
      const cb = document.getElementById(`ck-${{k}}`);
      if (cb) cb.checked = true;
    }});
    updateGraph();
  }});

  document.getElementById('btn-isolate-active').addEventListener('click', () => {{
    if (Graph) {{
      Graph.d3ReheatSimulation();
      const vis = Graph.graphData().nodes;
      if (vis.length > 0) {{
        const cx = vis.reduce((s,n) => s+(n.x||0), 0)/vis.length;
        const cy = vis.reduce((s,n) => s+(n.y||0), 0)/vis.length;
        const cz = vis.reduce((s,n) => s+(n.z||0), 0)/vis.length;
        Graph.cameraPosition({{x:cx+500,y:cy+200,z:cz+600}},{{x:cx,y:cy,z:cz}}, 800);
      }}
    }}
  }});

  // Legend clusters
  const legClusters = document.getElementById('legend-clusters');
  legClusters.innerHTML = '<div style="color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:6px;">Clusters</div>';
  DATA.cluster_keys.forEach(ck => {{
    const color = getClusterColor(ck);
    const cdata = DATA.clusters[ck];
    legClusters.innerHTML += `<div class="legend-row">
      <div class="cluster-legend-dot-${{ck}}" style="width:10px;height:10px;border-radius:50%;background:${{color}};flex-shrink:0;"></div>
      <span>${{cdata.display_name}}</span></div>`;
  }});

  // Filter checkboxes
  document.getElementById('filter-significant').addEventListener('change', e => {{
    state.filterSignificant = e.target.checked; updateGraph();
  }});
  document.getElementById('filter-hubs').addEventListener('change', e => {{
    state.filterHubs = e.target.checked; updateGraph();
  }});
  document.getElementById('filter-up-sets').addEventListener('change', e => {{
    state.filterUpSets = e.target.checked;
    if (e.target.checked) {{
      document.getElementById('filter-down-sets').checked = false;
      state.filterDownSets = false;
    }}
    updateGraph();
  }});
  document.getElementById('filter-down-sets').addEventListener('change', e => {{
    state.filterDownSets = e.target.checked;
    if (e.target.checked) {{
      document.getElementById('filter-up-sets').checked = false;
      state.filterUpSets = false;
    }}
    updateGraph();
  }});
  document.getElementById('filter-hide-overlap').addEventListener('change', e => {{
    state.hideOverlapEdges = e.target.checked; updateGraph();
  }});

  // Display toggles
  document.getElementById('toggle-labels').addEventListener('change', e => {{
    state.showLabels = e.target.checked;
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
  document.getElementById('toggle-light').addEventListener('change', e => {{
    document.body.classList.toggle('light-mode', e.target.checked);
    if (Graph) Graph.backgroundColor(e.target.checked ? '#f4f4f8' : '#0a0a1a');
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
  document.getElementById('toggle-colorblind').addEventListener('change', e => {{
    state.colorblindMode = e.target.checked;
    // Update cluster dots
    document.querySelectorAll('[data-cluster]').forEach(el => {{
      const ck = el.dataset.cluster;
      el.style.background = getClusterColor(ck);
    }});
    DATA.cluster_keys.forEach(ck => {{
      const dot = document.querySelector(`.cluster-legend-dot-${{ck}}`);
      if (dot) dot.style.background = getClusterColor(ck);
    }});
    if (Graph) {{
      Graph.linkColor(lk => {{
        if (lk.edge_type === 'overlap') return 'rgba(130,130,180,0.30)';
        return getClusterColor(lk.cluster_id) || '#666';
      }});
      Graph.linkDirectionalParticleColor(lk => getClusterColor(lk.cluster_id) || '#888');
      Graph.nodeThreeObject(buildNodeObject);
      Graph.refresh();
    }}
  }});
  document.getElementById('toggle-layer-labels').addEventListener('change', e => {{
    state.showLayerLabels = e.target.checked;
    _layerLabelSprites.forEach(s => {{ s.visible = state.showLayerLabels; }});
    if (Graph) Graph.refresh();
  }});

  document.getElementById('btn-reset-camera').addEventListener('click', resetCamera);
  document.getElementById('btn-reheat').addEventListener('click', () => {{
    if (Graph) Graph.d3ReheatSimulation();
  }});
  document.getElementById('close-detail').addEventListener('click', () => {{
    document.getElementById('detail-panel').style.display = 'none';
    const bcp = document.getElementById('bar-chart-panel');
    if (bcp && state.barChartVisible) bcp.style.right = '0';
  }});
  document.getElementById('btn-screenshot').addEventListener('click', doScreenshot);
  document.getElementById('btn-export-csv').addEventListener('click', exportCSV);
  document.getElementById('btn-export-edge-csv').addEventListener('click', exportEdgeCSV);
  document.getElementById('btn-bar-chart').addEventListener('click', toggleBarChart);

  // Feature 2: FDR/NES sliders
  document.getElementById('fdr-slider').addEventListener('input', function() {{
    state.fdrThreshold = parseFloat(this.value);
    document.getElementById('fdr-slider-val').textContent = state.fdrThreshold.toFixed(2);
    updateGraph();
  }});
  document.getElementById('nes-slider').addEventListener('input', function() {{
    state.nesThreshold = parseFloat(this.value);
    document.getElementById('nes-slider-val').textContent = state.nesThreshold.toFixed(1);
    updateGraph();
  }});

  // Feature 6: color by cluster toggle
  document.getElementById('toggle-color-cluster').addEventListener('change', e => {{
    state.colorByCluster = e.target.checked;
    if (Graph) {{
      Graph.nodeColor(n => {{
        if (state.highlightNodes.size > 0 && !state.highlightNodes.has(n.id))
          return 'rgba(60,60,80,0.2)';
        if (n.type === 'geneset')
          return state.colorByCluster ? getClusterColor(n.cluster_id) : nes2color(n.NES);
        return lfc2color(n.log2FC);
      }});
      Graph.nodeThreeObject(buildNodeObject);
      Graph.refresh();
    }}
  }});

  // Feature 1: show comparison subtitle
  if (DATA.has_comparison) {{
    const sub = document.getElementById('title-subtitle');
    if (sub) {{
      sub.textContent = `A: ${{DATA.title}} vs B: ${{DATA.title2}}`;
      sub.style.display = '';
    }}
  }}

  // Search (genes + gene set names)
  document.getElementById('search-box').addEventListener('input', function() {{
    const q = this.value.trim().toLowerCase();
    const results = document.getElementById('search-results');
    if (!q) {{
      results.textContent = '';
      state.searchMatches = new Set();
      if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
      return;
    }}
    const matches = DATA.nodes.filter(n =>
      n.id.toLowerCase().includes(q) ||
      (n.label && n.label.toLowerCase().includes(q)) ||
      (n.label_short && n.label_short.toLowerCase().includes(q))
    );
    state.searchMatches = new Set(matches.map(n => n.id));
    results.textContent = matches.length
      ? `${{matches.length}} match${{matches.length !== 1 ? 'es' : ''}}`
      : 'No matches';
    if (matches.length > 0 && Graph) {{
      const first = matches[0];
      const x = first.x || first.x_initial || 0;
      const y = first.y || first.y_initial || 0;
      const z = first.z || first.z_target || 0;
      Graph.cameraPosition({{x:x+200,y:y+80,z:z+350}},{{x,y,z}},800);
      pulseNodes(matches);
    }}
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
}}

// ============================================================
// Camera & Keyboard
// ============================================================
function resetCamera() {{
  if (Graph) Graph.cameraPosition({{x:0,y:0,z:1800}},{{x:0,y:0,z:0}},1000);
}}

document.addEventListener('keydown', e => {{
  if (e.target.tagName === 'INPUT') return;
  const key = e.key.toUpperCase();
  if (key === 'R') {{ resetCamera(); return; }}
  if (key === 'L') {{
    const cb = document.getElementById('toggle-labels');
    cb.checked = !cb.checked; state.showLabels = cb.checked;
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
    return;
  }}
  if (key === 'H') {{
    const cb = document.getElementById('filter-hubs');
    cb.checked = !cb.checked; state.filterHubs = cb.checked; updateGraph(); return;
  }}
  if (key === ' ') {{
    e.preventDefault();
    if (!Graph) return;
    if (state.paused) {{ Graph.resumeAnimation(); Graph.d3ReheatSimulation(); }}
    else Graph.pauseAnimation();
    state.paused = !state.paused; return;
  }}
  const digit = parseInt(key);
  if (digit >= 1 && digit <= DATA.cluster_keys.length && Graph) {{
    const ck = DATA.cluster_keys[digit - 1];
    const z = CLUSTER_Z[ck] || 0;
    Graph.cameraPosition({{x:600,y:200,z:z+700}},{{x:0,y:0,z}},800);
  }}
}});

// ============================================================
// Isolate cluster
// ============================================================
function isolateCluster(key) {{
  DATA.cluster_keys.forEach(k => {{
    const cb = document.getElementById(`ck-${{k}}`);
    if (k === key) {{ state.visibleClusters.add(k); if (cb) cb.checked = true; }}
    else {{ state.visibleClusters.delete(k); if (cb) cb.checked = false; }}
  }});
  updateGraph();
  const z = CLUSTER_Z[key] || 0;
  if (Graph) Graph.cameraPosition({{x:500,y:150,z:z+700}},{{x:0,y:0,z}},800);
}}

// ============================================================
// Search pulse animation
// ============================================================
function pulseNodes(matchNodes) {{
  if (!Graph) return;
  const T = getThree();
  if (!T) return;
  let frame = 0;
  const totalFrames = 40;
  const ids = new Set(matchNodes.map(n => n.id));
  function animate() {{
    frame++;
    const scale = 1 + 0.6 * Math.sin(frame / totalFrames * Math.PI);
    Graph.graphData().nodes.forEach(n => {{
      if (ids.has(n.id) && n.__threeObj) n.__threeObj.scale.setScalar(scale);
    }});
    if (frame < totalFrames) requestAnimationFrame(animate);
    else Graph.graphData().nodes.forEach(n => {{
      if (n.__threeObj) n.__threeObj.scale.setScalar(1);
    }});
  }}
  requestAnimationFrame(animate);
}}

// ============================================================
// Legend toggle
// ============================================================
function toggleLegend() {{
  const body = document.getElementById('legend-body');
  const icon = document.getElementById('legend-toggle-icon');
  icon.textContent = body.classList.toggle('collapsed') ? '▼' : '▲';
}}

// ============================================================
// Feature 7: RES Sparkline
// ============================================================
function drawSparkline(canvas, curve, nes) {{
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0,0,W,H);
  const lightMode = document.body.classList.contains('light-mode');
  const bgColor = lightMode ? '#f0f0f8' : '#0a0a1a';
  ctx.fillStyle = bgColor;
  ctx.fillRect(0,0,W,H);
  const mn = Math.min(...curve), mx = Math.max(...curve);
  const range = mx - mn || 1;
  const toY = v => H - 4 - ((v - mn) / range) * (H - 8);
  const zeroY = toY(0);
  ctx.strokeStyle = 'rgba(150,150,200,0.3)';
  ctx.lineWidth = 0.5;
  ctx.beginPath(); ctx.moveTo(0, zeroY); ctx.lineTo(W, zeroY); ctx.stroke();
  const fillColor = (nes||0) > 0 ? 'rgba(178,24,43,0.25)' : 'rgba(33,102,172,0.25)';
  ctx.fillStyle = fillColor;
  ctx.beginPath();
  ctx.moveTo(0, toY(curve[0]));
  curve.forEach((v, i) => ctx.lineTo(i / (curve.length-1) * W, toY(v)));
  ctx.lineTo(W, zeroY); ctx.lineTo(0, zeroY); ctx.closePath();
  ctx.fill();
  const lineColor = (nes||0) > 0 ? '#ff6666' : '#6699ff';
  ctx.strokeStyle = lineColor;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  curve.forEach((v, i) => {{
    const x = i / (curve.length-1) * W;
    i === 0 ? ctx.moveTo(x, toY(v)) : ctx.lineTo(x, toY(v));
  }});
  ctx.stroke();
  const peakIdx = curve.indexOf((nes||0) > 0 ? mx : mn);
  const px = peakIdx / (curve.length-1) * W;
  const py = toY((nes||0) > 0 ? mx : mn);
  ctx.fillStyle = lineColor;
  ctx.beginPath(); ctx.arc(px, py, 3, 0, Math.PI*2); ctx.fill();
}}

// ============================================================
// Feature 5: Jaccard similarity lookup
// ============================================================
function getTopSimilar(term, N) {{
  if (!DATA.jaccard_terms || DATA.jaccard_terms.length === 0) return [];
  const idx = DATA.jaccard_terms.indexOf(term);
  if (idx < 0) return [];
  const n = DATA.jaccard_terms.length;
  const sims = DATA.jaccard_terms.map((t, i) => ({{
    term: t,
    jaccard: 1 - DATA.jaccard_matrix[idx * n + i]
  }})).filter((x, i) => i !== idx && x.jaccard > 0);
  sims.sort((a, b) => b.jaccard - a.jaccard);
  return sims.slice(0, N || 5);
}}

// ============================================================
// Feature 3: Bar chart panel
// ============================================================
let _barChartBars = [];

function drawBarChart() {{
  const canvas = document.getElementById('bar-chart-canvas');
  if (!canvas) return;
  const gsNodes = DATA.nodes.filter(n => n.type === 'geneset' && isGeneSetVisible(n));
  gsNodes.sort((a, b) => (b.NES||0) - (a.NES||0));
  const barH = 18, gap = 2, labelW = 135, valW = 38, barAreaW = 100, pad = 8;
  const totalH = gsNodes.length * (barH + gap) + pad * 2;
  canvas.width = 230;
  canvas.height = Math.max(totalH, 40);
  canvas.style.height = canvas.height + 'px';
  const ctx = canvas.getContext('2d');
  const lightMode = document.body.classList.contains('light-mode');
  ctx.fillStyle = lightMode ? '#f0f0f8' : '#0a0a1a';
  ctx.fillRect(0,0,canvas.width,canvas.height);
  const maxNES = 3;
  _barChartBars = [];
  gsNodes.forEach((node, i) => {{
    const y = pad + i * (barH + gap);
    const nes = node.NES || 0;
    const barW = Math.round(Math.abs(nes) / maxNES * barAreaW);
    const barColor = nes > 0 ? 'rgba(178,24,43,0.75)' : 'rgba(33,102,172,0.75)';
    ctx.fillStyle = barColor;
    ctx.fillRect(labelW, y + 2, barW, barH - 4);
    const label = (node.label_short || node.label || '').slice(0,22);
    ctx.fillStyle = lightMode ? '#1a1a2e' : '#c0ccff';
    ctx.font = '10px "Segoe UI", sans-serif';
    ctx.fillText(label, 2, y + barH - 5);
    const nesStr = nes.toFixed(2);
    ctx.fillStyle = nes > 0 ? '#ff8888' : '#88aaff';
    ctx.font = 'bold 10px monospace';
    ctx.fillText(nesStr, labelW + barAreaW + 4, y + barH - 5);
    _barChartBars.push({{ node, y, h: barH }});
  }});
}}

function toggleBarChart() {{
  const panel = document.getElementById('bar-chart-panel');
  const dp = document.getElementById('detail-panel');
  state.barChartVisible = !state.barChartVisible;
  panel.style.display = state.barChartVisible ? 'block' : 'none';
  if (state.barChartVisible) {{
    // If detail panel is open, shift bar chart left
    const dpVisible = dp.style.display === 'block';
    panel.style.right = dpVisible ? '290px' : '0';
    drawBarChart();
    const canvas = document.getElementById('bar-chart-canvas');
    if (canvas) {{
      canvas.onclick = (e) => {{
        const rect = canvas.getBoundingClientRect();
        const my = e.clientY - rect.top;
        for (const bar of _barChartBars) {{
          if (my >= bar.y && my <= bar.y + bar.h) {{
            const n = nodeMap[bar.node.id];
            if (n) onNodeClick(n);
            break;
          }}
        }}
      }};
    }}
  }}
}}

// ============================================================
// Export
// ============================================================
function doScreenshot() {{
  if (!Graph) return;
  Graph.renderer().render(Graph.scene(), Graph.camera());
  const src = Graph.renderer().domElement;
  const W = src.width, H = src.height;
  const offscreen = document.createElement('canvas');
  offscreen.width = W; offscreen.height = H;
  const ctx = offscreen.getContext('2d');
  ctx.drawImage(src, 0, 0);

  const lineH = 18, pad = 10;
  const rows = DATA.cluster_keys.length + 7;
  const legH = rows * lineH + pad * 2 + 30;
  const legW = 210;
  const lx = W - legW - 12, ly = H - legH - 12;
  ctx.fillStyle = 'rgba(10,10,30,0.88)';
  ctx.beginPath();
  ctx.roundRect ? ctx.roundRect(lx,ly,legW,legH,6) : ctx.rect(lx,ly,legW,legH);
  ctx.fill();
  ctx.font = 'bold 10px sans-serif'; ctx.fillStyle = '#7080c0'; ctx.textAlign = 'left';
  let cy2 = ly + pad + 12;
  ctx.fillText('LEGEND', lx+pad, cy2); cy2 += lineH;
  ctx.font = '9px sans-serif'; ctx.fillStyle = '#7080c0';
  ctx.fillText('Gene log2FC', lx+pad, cy2); cy2 += 4;
  const barGrad = ctx.createLinearGradient(lx+pad, 0, lx+pad+120, 0);
  barGrad.addColorStop(0,'#2166ac'); barGrad.addColorStop(0.5,'#f7f7f7'); barGrad.addColorStop(1,'#b2182b');
  ctx.fillStyle = barGrad; ctx.fillRect(lx+pad, cy2, 120, 8); cy2 += 12;
  ctx.fillStyle='#888'; ctx.font='8px sans-serif';
  ctx.fillText('−',lx+pad,cy2); ctx.textAlign='center';
  ctx.fillText('0',lx+pad+60,cy2); ctx.textAlign='right';
  ctx.fillText('+',lx+pad+120,cy2); ctx.textAlign='left'; cy2+=lineH;
  ctx.font='9px sans-serif'; ctx.fillStyle='#aaa';
  ctx.fillText('◆ Gene Set (NES)   ● Gene (log2FC)', lx+pad, cy2); cy2+=lineH;
  DATA.cluster_keys.forEach(ck => {{
    const cdata = DATA.clusters[ck];
    const color = getClusterColor(ck);
    ctx.beginPath(); ctx.arc(lx+pad+5,cy2-4,5,0,Math.PI*2);
    ctx.fillStyle=color; ctx.fill();
    ctx.fillStyle='#dde'; ctx.font='9px sans-serif';
    ctx.fillText(cdata.display_name.slice(0,26), lx+pad+14, cy2); cy2+=lineH;
  }});
  const url = offscreen.toDataURL('image/png');
  const a = document.createElement('a');
  a.href=url; a.download='gsea_network_3d.png'; a.click();
}}

function exportCSV() {{
  if (!Graph) return;
  const vis = Graph.graphData().nodes || [];
  const header = 'id,type,label,gene_sets,clusters,log2FC,padj,node_status,is_hub,NES,FDR,cluster_id,cluster_name';
  const rows = vis.map(n => [
    n.id, n.type||'', n.label||'',
    (n.gene_sets||[]).join(';'),
    (n.clusters||[]).join(';'),
    n.log2FC!==null&&n.log2FC!==undefined?n.log2FC:'',
    n.padj!==null&&n.padj!==undefined?n.padj:'',
    n.node_status||'',
    n.is_hub?'TRUE':'FALSE',
    n.NES!==undefined?n.NES:'',
    n.FDR!==undefined?n.FDR:'',
    n.cluster_id||'', n.cluster_name||''
  ].map(v=>`"${{String(v).replace(/"/g,'""')}}"`).join(','));
  const csv=[header,...rows].join('\n');
  const blob=new Blob([csv],{{type:'text/csv'}});
  const url=URL.createObjectURL(blob);
  const a=document.createElement('a'); a.href=url; a.download='gsea_network_nodes.csv'; a.click();
  URL.revokeObjectURL(url);
}}

function exportEdgeCSV() {{
  if (!Graph) return;
  const vis = Graph.graphData().links || [];
  const header = 'source,target,edge_type,cluster_id,jaccard';
  const rows = vis.map(lk => {{
    const src = typeof lk.source==='object'?lk.source.id:lk.source;
    const tgt = typeof lk.target==='object'?lk.target.id:lk.target;
    return [src,tgt,lk.edge_type||'',lk.cluster_id||'',lk.jaccard||'']
      .map(v=>`"${{String(v).replace(/"/g,'""')}}"`).join(',');
  }});
  const csv=[header,...rows].join('\n');
  const blob=new Blob([csv],{{type:'text/csv'}});
  const url=URL.createObjectURL(blob);
  const a=document.createElement('a'); a.href=url; a.download='gsea_network_edges.csv'; a.click();
  URL.revokeObjectURL(url);
}}

// ============================================================
// Boot
// ============================================================
window.addEventListener('load', () => {{
  const loadMsg = document.getElementById('loading-msg');
  loadMsg.textContent = `Loading ${{DATA.nodes.length}} nodes, ${{DATA.links.length}} edges…`;

  if (typeof ForceGraph3D === 'undefined') {{
    loadMsg.textContent = 'ERROR: ForceGraph3D library not loaded. Check internet connection.';
    loadMsg.style.color = '#ff6666';
    return;
  }}

  try {{
    buildControls();
    buildGraph();
  }} catch(err) {{
    loadMsg.textContent = 'ERROR: ' + err.message;
    loadMsg.style.color = '#ff6666';
    console.error('buildGraph failed:', err);
  }}
}});
</script>
</body>
</html>
"""


def generate_html(nodes, links, clusters_meta, title, p_cutoff, max_lfc, output_path,
                  jaccard_terms=None, jaccard_matrix=None,
                  has_comparison=False, title2=None):
    """Render and write the HTML file."""
    cluster_keys = list(clusters_meta.keys())

    data_obj = {
        'title': title,
        'p_cutoff': p_cutoff,
        'max_lfc': max_lfc,
        'clusters': clusters_meta,
        'cluster_keys': cluster_keys,
        'nodes': nodes,
        'links': links,
        'has_comparison': has_comparison,
        'title2': title2 or '',
        'jaccard_terms': jaccard_terms or [],
        'jaccard_matrix': jaccard_matrix or [],
    }
    data_json = json.dumps(data_obj, ensure_ascii=False, separators=(',', ':'))

    html = GSEA_HTML_TEMPLATE.format(
        title=title,
        p_cutoff=p_cutoff,
        max_lfc=round(max_lfc, 2),
        data_json=data_json,
    )

    Path(output_path).write_text(html, encoding='utf-8')
    size_mb = Path(output_path).stat().st_size / 1_000_000
    return size_mb


# ============================================================
# Interactive gene-set library prompt
# ============================================================

_LIBRARY_MENU = [
    ('H',            'MSigDB Hallmark          — 50 curated hallmark gene sets (recommended)'),
    ('KEGG',         'KEGG 2019 Human          — ~200 metabolic & signaling pathways'),
    ('KEGG_GGA',     'KEGG Gallus gallus (gga) — chicken-specific pathways via KEGG REST API'),
    ('REACTOME',     'Reactome 2022            — ~1500 detailed pathway reactions'),
    ('C5_GO_BP',     'GO Biological Process    — ~7500 gene ontology BP terms'),
    ('C5_GO_MF',     'GO Molecular Function    — GO molecular function terms'),
    ('C5_GO_CC',     'GO Cellular Component    — GO cellular component terms'),
    ('WIKIPATHWAYS', 'WikiPathways Human       — ~300 community-curated pathways'),
]


def _prompt_gene_sets(default):
    """Interactively ask the user to choose a gene set library. Returns the chosen key."""
    print("\n┌─ Gene Set Library ──────────────────────────────────────────────────┐")
    for i, (key, desc) in enumerate(_LIBRARY_MENU, 1):
        marker = ' ◀ default' if key == default else ''
        print(f"│  {i})  {desc}{marker}")
    print("│  9)  Custom .gmt file — enter path when prompted")
    print(f"└─────────────────────────────────────────────────────────────────────┘")

    while True:
        try:
            raw = input(f"Choose [1–9, or press Enter for default '{default}']: ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return default

        if raw == '':
            print(f"      Using default: {default}")
            return default

        if raw.isdigit() and 1 <= int(raw) <= len(_LIBRARY_MENU):
            chosen = _LIBRARY_MENU[int(raw) - 1][0]
            print(f"      Selected: {chosen}")
            return chosen

        if raw == '9':
            try:
                path = input("      Path to .gmt file: ").strip()
            except (EOFError, KeyboardInterrupt):
                print()
                return default
            if Path(path).exists():
                return path
            print(f"      File not found: {path}. Try again.")
            continue

        # Accept a valid shorthand typed directly (e.g. 'KEGG')
        if raw.upper() in {k for k, _ in _LIBRARY_MENU}:
            print(f"      Selected: {raw.upper()}")
            return raw.upper()

        print("      Invalid choice — enter a number 1–9 or press Enter.")


# ============================================================
# Main
# ============================================================

def main():
    script_path = Path(__file__).resolve()
    script_dir  = script_path.parent

    # ================================================================
    # DEFAULT SETTINGS — edit these lines to change what runs when
    # you hit Run/F5 in Spyder (or run without --args in terminal).
    # ================================================================

    default_deseq2_csv = script_dir / 'EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv'
    default_output     = 'gsea_network_3d.html'
    default_title      = 'GSEA 3D Network Explorer'

    # Gene set library to use for GSEA enrichment.
    # Pick one of the shorthand names below, or supply a full path to a .gmt file.
    #
    #   'H'             — MSigDB Hallmark (50 curated sets) ← default, good starting point
    #   'KEGG'          — KEGG 2019 Human (~200 pathways)
    #   'KEGG_GGA'      — KEGG Gallus gallus (chicken-specific, via KEGG REST API)
    #   'REACTOME'      — Reactome 2022 (~1500 pathways, very detailed)
    #   'C5_GO_BP'      — GO Biological Process (~7500 terms, broad coverage)
    #   'C5_GO_MF'      — GO Molecular Function
    #   'C5_GO_CC'      — GO Cellular Component
    #   'WIKIPATHWAYS'  — WikiPathways Human (~300 community pathways)
    #   '/path/to/file.gmt'  — any local GMT file
    #
    # NOTE: libraries are cached in ~/.gsea3d_cache/ after the first download,
    # so switching between them is fast from the second run onward.
    default_gene_sets = 'H'

    # ================================================================

    parser = argparse.ArgumentParser(
        description='GSEA 3D Network Explorer — DESeq2 CSV → GSEA prerank → interactive HTML',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--deseq2-csv', default=str(default_deseq2_csv),
                        help='Path to DESeq2 results CSV')
    parser.add_argument('--deseq2-csv2', default=None,
                        help='(Optional) Second DESeq2 CSV for side-by-side comparison (Feature 1)')
    parser.add_argument('--gene-sets', default=None,
                        help='MSigDB collection name (H, KEGG, REACTOME, C5_GO_BP, C5_GO_MF, '
                             'C5_GO_CC, WIKIPATHWAYS) or path to a local .gmt file. '
                             'If omitted and running interactively, a menu is shown.')
    parser.add_argument('--output', default=default_output,
                        help='Output HTML path')
    parser.add_argument('--title', default=default_title,
                        help='Visualization title')
    parser.add_argument('--title2', default='Comparison 2',
                        help='Title label for the second comparison (Feature 1)')
    parser.add_argument('--fdr-cutoff', type=float, default=0.25,
                        help='FDR threshold for significant gene sets')
    parser.add_argument('--p-cutoff', type=float, default=0.05,
                        help='padj threshold for gene significance coloring')
    parser.add_argument('--permutations', type=int, default=1000,
                        help='Number of GSEA permutations')
    parser.add_argument('--max-sets', type=int, default=40,
                        help='Maximum gene sets to display')
    parser.add_argument('--min-size', type=int, default=15,
                        help='Minimum gene set size')
    parser.add_argument('--max-size', type=int, default=500,
                        help='Maximum gene set size')
    args = parser.parse_args()

    # If --gene-sets was not passed, either prompt (interactive terminal) or use default
    if args.gene_sets is None:
        if sys.stdin.isatty():
            args.gene_sets = _prompt_gene_sets(default_gene_sets)
        else:
            args.gene_sets = default_gene_sets

    csv_path = Path(args.deseq2_csv)
    if not csv_path.exists():
        using_default = csv_path == default_deseq2_csv
        if using_default:
            parser.error(
                f"DESeq2 CSV not found: {csv_path}\n"
                f"You are using the built-in default path.\n"
                f"Either edit `default_deseq2_csv` in main() of {script_path},\n"
                f"or run with:  --deseq2-csv \"/full/path/to/your_results.csv\"\n"
                f"Spyder:  %runfile \"{script_path}\" --args --deseq2-csv \"/full/path/to/your_results.csv\""
            )
        parser.error(f"DESeq2 CSV not found: {csv_path}")

    # Step 1
    print(f"\n[1/6] Reading DESeq2 data from: {args.deseq2_csv}")
    de_table, rnk_dict = read_deseq2(args.deseq2_csv)
    print(f"      {len(de_table)} genes loaded  |  {len(rnk_dict)} rankable")

    # Feature 1: optional second comparison
    de_table_b = {}
    rnk_dict_b = {}
    if args.deseq2_csv2:
        csv_path_b = Path(args.deseq2_csv2)
        if not csv_path_b.exists():
            print(f"WARNING: --deseq2-csv2 not found: {csv_path_b}. Skipping comparison.")
            args.deseq2_csv2 = None
        else:
            print(f"      Reading comparison CSV: {args.deseq2_csv2}")
            de_table_b, rnk_dict_b = read_deseq2(args.deseq2_csv2)
            print(f"      Comparison 2: {len(de_table_b)} genes  |  {len(rnk_dict_b)} rankable")

    # Step 2
    print(f"\n[2/6] Running GSEA prerank")
    results_df, pre_res = run_gsea(rnk_dict, args.gene_sets, args.permutations,
                                   args.min_size, args.max_size)
    print(f"      {len(results_df)} gene sets tested")

    # Feature 1: run GSEA on comparison 2 if provided
    nes_b_lk = {}
    fdr_b_lk = {}
    if args.deseq2_csv2 and rnk_dict_b:
        print(f"      Running GSEA on comparison 2…")
        try:
            results_df_b, _ = run_gsea(rnk_dict_b, args.gene_sets, args.permutations,
                                       args.min_size, args.max_size)
            import pandas as pd
            # Build lookup keyed by term
            if 'Term' in results_df_b.columns:
                for _, row in results_df_b.iterrows():
                    t = str(row['Term'])
                    nes_val = row.get('NES')
                    fdr_val = row.get('FDR_qval')
                    if nes_val is not None and not (isinstance(nes_val, float) and math.isnan(nes_val)):
                        nes_b_lk[t] = float(nes_val)
                    if fdr_val is not None and not (isinstance(fdr_val, float) and math.isnan(fdr_val)):
                        fdr_b_lk[t] = float(fdr_val)
            print(f"      Comparison 2: {len(nes_b_lk)} gene sets with NES")
        except Exception as e:
            print(f"      WARNING: Comparison 2 GSEA failed: {e}. Proceeding without comparison.")
            args.deseq2_csv2 = None

    # Step 3
    print(f"\n[3/6] Filtering significant gene sets (FDR < {args.fdr_cutoff})")
    sig_df = filter_genesets(results_df, args.fdr_cutoff, args.max_sets)
    print(f"\n      {len(sig_df)} gene sets will be visualized")

    if len(sig_df) == 0:
        print("ERROR: No gene sets to visualize. Try relaxing --fdr-cutoff or changing --gene-sets.")
        sys.exit(1)

    # Step 4
    print(f"\n[4/6] Clustering gene sets by Jaccard similarity")
    cluster_assignments, gene_sets_genes = cluster_genesets(sig_df)

    # Feature 5: compute Jaccard matrix
    print(f"      Computing Jaccard similarity matrix…")
    jaccard_terms, jaccard_matrix = compute_jaccard_matrix(sig_df, gene_sets_genes)

    # Feature 7: extract RES curves
    print(f"      Extracting RES curves…")
    res_curves = extract_res_curves(pre_res, sig_df)
    print(f"      RES curves extracted for {len(res_curves)} gene sets")

    # Step 5
    print(f"\n[5/6] Building 3D network")
    nodes, links, clusters_meta = build_network(
        sig_df, cluster_assignments, gene_sets_genes, de_table, args.p_cutoff,
        res_curves=res_curves,
        nes_b_lk=nes_b_lk,
        fdr_b_lk=fdr_b_lk,
        de_table_b=de_table_b,
    )

    gene_nodes = [n for n in nodes if n['type'] == 'gene']
    gs_nodes   = [n for n in nodes if n['type'] == 'geneset']
    mem_links  = [l for l in links if l['edge_type'] == 'membership']
    ovl_links  = [l for l in links if l['edge_type'] == 'overlap']

    print(f"      Gene-set nodes: {len(gs_nodes)}  |  Gene nodes: {len(gene_nodes)}")
    print(f"      Membership edges: {len(mem_links)}  |  Overlap edges: {len(ovl_links)}")

    sig_genes = sum(1 for n in gene_nodes if n['node_status'] == 'Present; significant')
    hub_genes = sum(1 for n in gene_nodes if n['is_hub'])
    print(f"      Significant gene nodes: {sig_genes}  |  Hub genes (≥2 sets): {hub_genes}")

    max_lfc = compute_max_lfc(nodes)
    print(f"      Max |log2FC|: {max_lfc}")

    has_comparison = bool(args.deseq2_csv2 and nes_b_lk)

    # Step 6
    print(f"\n[6/6] Generating HTML → {args.output}")
    size_mb = generate_html(
        nodes, links, clusters_meta, args.title,
        args.p_cutoff, max_lfc, args.output,
        jaccard_terms=jaccard_terms,
        jaccard_matrix=jaccard_matrix,
        has_comparison=has_comparison,
        title2=args.title2,
    )
    print(f"      Done! {args.output} ({size_mb:.1f} MB)")
    print(f"\nOpen in browser:\n  file://{Path(args.output).absolute()}")
    print(f"\nExample command to re-run:")
    print(f"  python3 {Path(__file__).name} \\")
    print(f"    --deseq2-csv {args.deseq2_csv} \\")
    print(f"    --gene-sets {args.gene_sets} \\")
    print(f"    --output {args.output} \\")
    print(f"    --title \"{args.title}\"")


if __name__ == '__main__':
    main()
