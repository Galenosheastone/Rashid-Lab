#!/usr/bin/env python3
"""
build_3d_network.py — 3D Force-Directed Pathway Network Explorer
Reads pathway_mapper_v4.R + DESeq2 results CSV → standalone HTML
"""

import re
import json
import argparse
import sys
from pathlib import Path
from collections import defaultdict

# ---------------------------------------------------------------------------
# 1. Parse the R pathway registry
# ---------------------------------------------------------------------------

UNICODE_ESCAPES = re.compile(r'\\u([0-9a-fA-F]{4})')

def decode_r_string(s):
    """Decode R string: strip quotes, resolve \\uXXXX escapes."""
    s = s.strip().strip('"').strip("'")
    s = UNICODE_ESCAPES.sub(lambda m: chr(int(m.group(1), 16)), s)
    return s

def parse_r_value(token):
    """Parse a single R token into a Python value."""
    token = token.strip()
    if token in ('NA_character_', 'NA'):
        return None
    if (token.startswith('"') and token.endswith('"')) or \
       (token.startswith("'") and token.endswith("'")):
        return decode_r_string(token)
    try:
        return float(token)
    except ValueError:
        return token  # bare symbol

def extract_tribble_content(func_body, var_name):
    """
    Extract the content inside tribble() for a given variable (nodes or edges).
    Uses balanced-paren counting, respecting string literals and # comments.
    This avoids false stops on ) inside comments (e.g., # CH25H (ISG induction)).
    Returns the inner content string, or None if not found.
    """
    pattern = re.compile(var_name + r'\s*<-\s*tribble\s*\(')
    m = pattern.search(func_body)
    if not m:
        return None
    pos = m.end()  # position right after opening (
    depth = 1
    while pos < len(func_body) and depth > 0:
        ch = func_body[pos]
        if ch in ('"', "'"):
            # Skip over string literal
            quote = ch
            pos += 1
            while pos < len(func_body):
                c2 = func_body[pos]
                if c2 == '\\':
                    pos += 2  # skip escaped char
                    continue
                if c2 == quote:
                    pos += 1
                    break
                pos += 1
            continue
        elif ch == '#':
            # Skip to end of line (comment)
            while pos < len(func_body) and func_body[pos] != '\n':
                pos += 1
            pos += 1
            continue
        elif ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        pos += 1
    return func_body[m.end():pos - 1]


def parse_tribble(block_text):
    """
    Parse a tribble() block into a list of dicts.
    block_text is the content between tribble( and the closing ).
    Handles comment lines starting with #.
    """
    lines = block_text.split('\n')
    # Remove comment lines and strip whitespace
    clean_lines = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('#') or not stripped:
            continue
        # Remove inline comments
        # Be careful: # inside a quoted string shouldn't be removed
        # Simple approach: only strip from # that's outside quotes
        result = []
        in_quote = False
        quote_char = None
        for i, ch in enumerate(stripped):
            if in_quote:
                if ch == quote_char:
                    in_quote = False
                result.append(ch)
            else:
                if ch in ('"', "'"):
                    in_quote = True
                    quote_char = ch
                    result.append(ch)
                elif ch == '#':
                    break  # rest is comment
                else:
                    result.append(ch)
        clean_lines.append(''.join(result).strip())

    # Join into one big string and tokenize
    full = ' '.join(clean_lines)

    # Extract header columns (~col_name pattern)
    headers = re.findall(r'~(\w+)', full)
    if not headers:
        return []

    # Remove header tokens from the string
    data_str = re.sub(r'~\w+', '', full)

    # Tokenize: quoted strings or bare tokens, separated by commas
    tokens = []
    token_re = re.compile(r'"[^"]*"|\'[^\']*\'|[^,\s]+')
    for m in token_re.finditer(data_str):
        t = m.group().strip()
        if t and t != ',':
            tokens.append(t)

    n_cols = len(headers)
    rows = []
    for i in range(0, len(tokens) - n_cols + 1, n_cols):
        chunk = tokens[i:i + n_cols]
        if len(chunk) != n_cols:
            break
        row = {}
        for h, v in zip(headers, chunk):
            row[h] = parse_r_value(v)
        rows.append(row)

    return rows

def parse_synonyms_block(block_text):
    """Parse synonyms list(...) into {symbol: [syn1, syn2, ...]}."""
    synonyms = {}
    # Match KEY = c("val1", "val2") patterns
    pattern = re.compile(r'(\w+)\s*=\s*c\(([^)]*)\)')
    for m in pattern.finditer(block_text):
        key = m.group(1)
        values_str = m.group(2)
        values = [decode_r_string(v.strip()) for v in values_str.split(',') if v.strip()]
        synonyms[key] = values
    return synonyms

def parse_validation_notes(block_text):
    """Parse c("note1", "note2", ...) into list of strings."""
    notes = []
    # Find the c(...) content
    m = re.search(r'c\(([\s\S]*)\)', block_text)
    if not m:
        return notes
    content = m.group(1)
    # Extract quoted strings
    for s in re.findall(r'"((?:[^"\\]|\\.)*)"|\'((?:[^\'\\]|\\.)*)\'', content):
        raw = s[0] or s[1]
        decoded = UNICODE_ESCAPES.sub(lambda mm: chr(int(mm.group(1), 16)), raw)
        notes.append(decoded)
    return notes

def parse_pathway_registry(r_script_path):
    """Parse all pathway definitions from pathway_mapper_v4.R."""
    text = Path(r_script_path).read_text(encoding='utf-8')

    pathways = {}

    # Find each pathway_registry[["name"]] <- function() { ... }
    # We need to match balanced braces
    registry_pattern = re.compile(
        r'pathway_registry\[\["(\w+)"\]\]\s*<-\s*function\s*\(\s*\)\s*\{'
    )

    for m in registry_pattern.finditer(text):
        pathway_name = m.group(1)
        start = m.end()  # position after opening {

        # Find the matching closing brace
        depth = 1
        pos = start
        while pos < len(text) and depth > 0:
            if text[pos] == '{':
                depth += 1
            elif text[pos] == '}':
                depth -= 1
            pos += 1
        func_body = text[start:pos - 1]

        # --- Parse nodes tribble ---
        nodes_content = extract_tribble_content(func_body, 'nodes')
        nodes = parse_tribble(nodes_content) if nodes_content else []

        # --- Parse edges tribble ---
        edges_content = extract_tribble_content(func_body, 'edges')
        edges = parse_tribble(edges_content) if edges_content else []

        # --- Parse synonyms ---
        synonyms_match = re.search(
            r'synonyms\s*<-\s*list\s*\(([\s\S]*?)\)\s*\n', func_body
        )
        synonyms = {}
        if synonyms_match:
            synonyms = parse_synonyms_block(synonyms_match.group(1))

        # --- Parse validation_notes ---
        val_match = re.search(
            r'validation_notes\s*<-\s*c\s*\(([\s\S]*?)\)\s*\n', func_body
        )
        validation_notes = []
        if val_match:
            validation_notes = parse_validation_notes('c(' + val_match.group(1) + ')')

        # --- Parse display_name ---
        dn_match = re.search(r'display_name\s*=\s*"([^"]+)"', func_body)
        display_name = decode_r_string(dn_match.group(1)) if dn_match else pathway_name

        pathways[pathway_name] = {
            'name': pathway_name,
            'display_name': display_name,
            'nodes': nodes,
            'edges': edges,
            'synonyms': synonyms,
            'validation_notes': validation_notes,
        }

    return pathways


# ---------------------------------------------------------------------------
# 2. Read DESeq2 CSV and build gene lookup table
# ---------------------------------------------------------------------------

def read_deseq2(csv_path):
    """
    Read DESeq2 CSV. Returns dict: {symbol_upper: {log2FC, padj, symbol}}.
    """
    import csv
    de_table = {}
    with open(csv_path, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # First column is unnamed → key is ''
            symbol = row.get('') or row.get('Unnamed: 0') or row.get('...1')
            if not symbol:
                # Try first key
                first_key = list(row.keys())[0]
                symbol = row[first_key]
            if not symbol:
                continue

            try:
                lfc = float(row.get('log2FC_shrunken') or row.get('log2FC') or 'nan')
                if lfc != lfc:  # NaN
                    lfc_raw = float(row.get('log2FC') or 'nan')
                    lfc = lfc_raw if lfc_raw == lfc_raw else None
            except (ValueError, TypeError):
                lfc = None

            try:
                padj_str = row.get('padj', '')
                padj = float(padj_str) if padj_str and padj_str.lower() != 'na' else None
            except (ValueError, TypeError):
                padj = None

            de_table[symbol.upper()] = {
                'symbol': symbol,
                'log2FC': lfc,
                'padj': padj,
            }
    return de_table


# ---------------------------------------------------------------------------
# 3. Map DE data to nodes
# ---------------------------------------------------------------------------

def map_de_to_node(symbol_key, synonyms_for_key, de_table, p_cutoff):
    """
    Look up a gene in the DE table by symbol_key then synonyms.
    Returns dict with log2FC, padj, node_status, matched_as.
    """
    # Try direct match
    hit = de_table.get(symbol_key.upper())
    matched_as = symbol_key if hit else None

    # Try synonyms
    if not hit:
        for syn in (synonyms_for_key or []):
            hit = de_table.get(syn.upper())
            if hit:
                matched_as = syn
                break

    if not hit:
        return {
            'log2FC': None,
            'padj': None,
            'node_status': 'Missing gene',
            'matched_as': None,
        }

    lfc = hit['log2FC']
    padj = hit['padj']

    if padj is None:
        status = 'Present; padj NA'
    elif padj < p_cutoff:
        status = 'Present; significant'
    else:
        status = 'Present; not significant'

    return {
        'log2FC': lfc,
        'padj': padj,
        'node_status': status,
        'matched_as': matched_as,
    }


# ---------------------------------------------------------------------------
# 4. Build unified network
# ---------------------------------------------------------------------------

PATHWAY_Z = {
    'cgas_sting':   0,
    'apoptosis':    150,
    'necroptosis':  300,
    'inflammasome': 450,
    'osteoclast':   600,
    'osteoblast':   750,
    'oxysterol':    900,
}

XY_SCALE = 90  # Scale R config x,y to 3D space

def build_network(pathways, de_table, p_cutoff):
    """
    Merge all pathways into a single node+edge graph.
    Returns (global_nodes, global_edges).
    """
    # --- Collect per-pathway data ---
    # key: symbol_key → accumulate data across pathways
    node_accumulator = defaultdict(lambda: {
        'label': None,
        'pathways': [],
        'modules': [],
        'compartments': [],
        'x_sum': 0.0,
        'y_sum': 0.0,
        'z_sum': 0.0,
        'count': 0,
        'log2FC': None,
        'padj': None,
        'node_status': None,
        'matched_as': None,
    })

    # Build a per-pathway node_id → symbol_key map for edge resolution
    pathway_node_id_to_sym = {}  # {pathway_name: {node_id: symbol_key}}

    for pathway_name, pdata in pathways.items():
        z = PATHWAY_Z.get(pathway_name, 0)
        id_to_sym = {}

        for node in pdata['nodes']:
            nid = node.get('node_id', '')
            sym = node.get('symbol_key', nid)
            id_to_sym[nid] = sym

            acc = node_accumulator[sym]
            # Prefer label without trailing markup (e.g. "NLRP3↑" → keep it)
            if acc['label'] is None:
                acc['label'] = node.get('label') or sym
            # Track pathways
            if pathway_name not in acc['pathways']:
                acc['pathways'].append(pathway_name)
            mod = node.get('module')
            if mod and mod not in acc['modules']:
                acc['modules'].append(mod)
            comp = node.get('compartment')
            if comp and comp not in acc['compartments']:
                acc['compartments'].append(comp)

            x = (node.get('x') or 0) * XY_SCALE
            y = (node.get('y') or 0) * XY_SCALE
            acc['x_sum'] += x
            acc['y_sum'] += y
            acc['z_sum'] += z
            acc['count'] += 1

        pathway_node_id_to_sym[pathway_name] = id_to_sym

    # Apply DE mapping (do this once per unique symbol_key)
    for sym, acc in node_accumulator.items():
        # Collect all synonyms across pathways
        all_synonyms = set()
        for pathway_name in acc['pathways']:
            syns = pathways[pathway_name]['synonyms'].get(sym, [])
            all_synonyms.update(syns)

        de_result = map_de_to_node(sym, list(all_synonyms), de_table, p_cutoff)
        acc.update(de_result)

    # Build final node list
    global_nodes = []
    for sym, acc in node_accumulator.items():
        n = acc['count']
        xi = round(acc['x_sum'] / n, 2)
        yi = round(acc['y_sum'] / n, 2)
        zt = round(acc['z_sum'] / n, 2)
        global_nodes.append({
            'id': sym,
            'label': acc['label'] or sym,
            'pathways': acc['pathways'],
            'modules': acc['modules'],
            'compartments': acc['compartments'],
            # x,y,z used by 3d-force-graph as initial positions
            'x': xi,
            'y': yi,
            'z': zt,
            'x_initial': xi,
            'y_initial': yi,
            'z_target': zt,
            'log2FC': acc['log2FC'],
            'padj': acc['padj'],
            'node_status': acc['node_status'],
            'matched_as': acc['matched_as'],
            'is_hub': len(acc['pathways']) > 1,
        })

    # --- Build edges ---
    seen_edges = set()
    global_edges = []

    for pathway_name, pdata in pathways.items():
        id_to_sym = pathway_node_id_to_sym[pathway_name]
        for edge in pdata['edges']:
            src_id = edge.get('from', '')
            tgt_id = edge.get('to', '')
            # Resolve node_id → symbol_key
            src = id_to_sym.get(src_id, src_id)
            tgt = id_to_sym.get(tgt_id, tgt_id)
            edge_class = edge.get('edge_class', 'core')
            edge_pathway = edge.get('edge_pathway', pathway_name)
            edge_sign = edge.get('edge_sign', 'activating')

            edge_key = (src, tgt, edge_class, pathway_name)
            if edge_key in seen_edges:
                continue
            seen_edges.add(edge_key)

            global_edges.append({
                'source': src,
                'target': tgt,
                'edge_class': edge_class,
                'edge_pathway': edge_pathway,
                'edge_sign': edge_sign,
                'pathway': pathway_name,
            })

    return global_nodes, global_edges


# ---------------------------------------------------------------------------
# 5. Compute max |log2FC| for color scale
# ---------------------------------------------------------------------------

def compute_max_lfc(nodes):
    vals = [abs(n['log2FC']) for n in nodes if n.get('log2FC') is not None]
    if not vals:
        return 3.0
    return round(max(vals), 2)


# ---------------------------------------------------------------------------
# 6. Generate HTML
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
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
    background: var(--bg);
    color: var(--text);
    font-family: 'Segoe UI', system-ui, sans-serif;
    font-size: 13px;
    overflow: hidden;
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
    background: var(--title-bar-bg);
    border-bottom: 1px solid var(--border);
    display: flex; align-items: center; padding: 0 16px; z-index: 100; gap: 16px;
  }}
  #title-bar h1 {{ font-size: 14px; font-weight: 600; color: var(--title-color); white-space: nowrap; }}
  #title-bar .stats {{ font-size: 11px; color: var(--muted); white-space: nowrap; }}
  #layout {{ display: flex; position: fixed; top: 38px; left: 0; right: 0; bottom: 0; }}
  #control-panel {{
    width: 240px; min-width: 240px;
    background: var(--panel-bg);
    border-right: 1px solid var(--border);
    overflow-y: auto; padding: 10px; z-index: 50;
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
  .check-row input[type=checkbox] {{ cursor: pointer; accent-color: #6a8aff; }}
  .pathway-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
  .check-label {{ font-size: 12px; color: var(--text); line-height: 1.3; }}
  #search-box {{
    width: 100%;
    background: var(--input-bg);
    border: 1px solid var(--border);
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
  .toggle-row {{ display: flex; align-items: center; justify-content: space-between; margin-bottom: 6px; }}
  .kbd {{
    display: inline-block; background: rgba(60,60,100,0.6);
    border: 1px solid rgba(100,100,180,0.4); border-radius: 3px;
    padding: 1px 5px; font-size: 10px; font-family: monospace; color: #aabf;
  }}
  #tooltip {{
    position: fixed;
    background: var(--tooltip-bg);
    border: 1px solid rgba(100,120,255,0.4);
    border-radius: 6px; padding: 10px 13px;
    pointer-events: none; z-index: 200; display: none;
    min-width: 200px; max-width: 280px; font-size: 12px; line-height: 1.6;
    box-shadow: 0 4px 20px rgba(0,0,0,0.6); color: var(--text);
  }}
  #tooltip .tt-gene {{ font-size: 14px; font-weight: 700; color: var(--title-color); margin-bottom: 6px; }}
  #tooltip .tt-row {{ display: flex; gap: 8px; }}
  #tooltip .tt-key {{ color: var(--heading); min-width: 90px; }}
  #tooltip .tt-val {{ color: var(--text); }}
  #tooltip .tt-sig {{ color: #ff9966; font-weight: 600; }}
  #tooltip .tt-up {{ color: #ff8888; }}
  #tooltip .tt-down {{ color: #88aaff; }}
  #detail-panel {{
    position: fixed; right: 0; top: 38px; bottom: 0; width: 280px;
    background: var(--detail-bg);
    border-left: 1px solid var(--border);
    padding: 12px; overflow-y: auto; z-index: 80; display: none;
    font-size: 12px; line-height: 1.6; color: var(--text); pointer-events: auto;
  }}
  #detail-panel h2 {{ font-size: 15px; color: var(--title-color); margin-bottom: 10px; }}
  #detail-panel .dp-key {{ color: var(--heading); }}
  #detail-panel .dp-val {{ color: var(--text); }}
  #detail-panel .dp-section {{
    margin-top: 10px; color: var(--heading); font-size: 10px;
    text-transform: uppercase; letter-spacing: 1px;
    border-bottom: 1px solid var(--border); padding-bottom: 3px; margin-bottom: 6px;
  }}
  #detail-panel .neighbor-row {{ display: flex; justify-content: space-between; margin-bottom: 3px; }}
  #detail-panel .note-item {{
    font-size: 11px; color: var(--muted); margin-bottom: 4px;
    border-left: 2px solid var(--border); padding-left: 6px;
  }}
  #close-detail {{ float: right; cursor: pointer; color: var(--muted); font-size: 16px; line-height: 1; }}
  #legend {{
    position: fixed; bottom: 16px; right: 16px;
    background: var(--legend-bg);
    border: 1px solid var(--border);
    border-radius: 8px; padding: 10px 14px; z-index: 60; min-width: 170px; font-size: 11px;
    color: var(--text);
  }}
  #legend-header {{
    display: flex; align-items: center; justify-content: space-between;
    cursor: pointer; user-select: none;
  }}
  #legend h4 {{
    font-size: 10px; text-transform: uppercase; letter-spacing: 1px;
    color: var(--heading); margin: 0;
  }}
  #legend-toggle-icon {{
    font-size: 12px; color: var(--heading); line-height: 1; margin-left: 8px;
  }}
  #legend-body {{ margin-top: 8px; }}
  #legend-body.collapsed {{ display: none; }}
  .lfc-bar {{
    width: 140px; height: 12px;
    background: linear-gradient(to right, #2166ac, #f7f7f7, #b2182b);
    border-radius: 3px; margin: 4px 0 2px;
  }}
  .lfc-labels {{ display: flex; justify-content: space-between; color: var(--muted); font-size: 10px; }}
  .legend-row {{ display: flex; align-items: center; gap: 7px; margin-bottom: 4px; }}
  #loading {{
    position: fixed; inset: 0; background: var(--bg);
    display: flex; flex-direction: column; align-items: center;
    justify-content: center; z-index: 999; gap: 16px;
  }}
  .spinner {{
    width: 40px; height: 40px;
    border: 3px solid rgba(100,130,255,0.2);
    border-top-color: #6a8aff; border-radius: 50%;
    animation: spin 0.8s linear infinite;
  }}
  @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
  #loading p {{ color: var(--muted); font-size: 13px; }}
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
</style>
</head>
<body>

<div id="loading">
  <div class="spinner"></div>
  <p id="loading-msg">Initializing 3D network...</p>
</div>

<div id="title-bar">
  <button id="hamburger" onclick="document.getElementById('control-panel').classList.toggle('open')" title="Toggle controls">&#9776;</button>
  <h1>{title}</h1>
  <span class="stats" id="graph-stats"></span>
  <span class="stats">Drag to rotate &nbsp;|&nbsp; Scroll to zoom &nbsp;|&nbsp; Right-drag to pan</span>
</div>

<div id="layout">
  <div id="control-panel">
    <div class="panel-section">
      <h3>Search</h3>
      <input id="search-box" type="text" placeholder="Gene symbol or label…" autocomplete="off">
      <div id="search-results"></div>
    </div>
    <div class="panel-section">
      <h3>Pathways</h3>
      <div id="pathway-checkboxes"></div>
      <button class="btn" id="btn-show-all" style="margin-top:4px;">Show All</button>
      <button class="btn" id="btn-isolate-active" style="color:#ffcc66;">Focus on checked ↑</button>
    </div>
    <div class="panel-section">
      <h3>Filters</h3>
      <div class="check-row">
        <input type="checkbox" id="filter-significant">
        <label for="filter-significant" class="check-label">Significant only (padj &lt; {p_cutoff})</label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="filter-hubs">
        <label for="filter-hubs" class="check-label">Hub nodes only (multi-pathway)</label>
      </div>
    </div>
    <div class="panel-section">
      <h3>Display</h3>
      <div class="check-row">
        <input type="checkbox" id="toggle-labels" checked>
        <label for="toggle-labels" class="check-label">Show labels <span class="kbd">L</span></label>
      </div>
      <div class="check-row">
        <input type="checkbox" id="toggle-arrows" checked>
        <label for="toggle-arrows" class="check-label">Show arrows</label>
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
        <label for="toggle-layer-labels" class="check-label">Show layer labels</label>
      </div>
    </div>
    <div class="panel-section">
      <h3>Camera &nbsp;<span class="kbd">R</span>=reset &nbsp;<span class="kbd">1-7</span>=layer</h3>
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

<div id="legend">
  <div id="legend-header" onclick="toggleLegend()">
    <h4>Legend</h4>
    <span id="legend-toggle-icon">▲</span>
  </div>
  <div id="legend-body">
    <div style="margin-top:8px;color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:4px;">log2 Fold-Change</div>
    <div class="lfc-bar"></div>
    <div class="lfc-labels"><span>−3</span><span>0</span><span>+3</span></div>
    <div style="margin-top:10px;">
      <div class="legend-row">
        <div style="width:18px;height:18px;border-radius:50%;background:#aaa;border:2px solid #fff;flex-shrink:0;"></div>
        <span>Significant (padj &lt; {p_cutoff})</span>
      </div>
      <div class="legend-row">
        <div style="width:11px;height:11px;border-radius:50%;background:#666;flex-shrink:0;"></div>
        <span>Not significant / padj NA</span>
      </div>
      <div class="legend-row">
        <div style="width:7px;height:7px;border-radius:50%;background:#333;flex-shrink:0;"></div>
        <span>Missing gene</span>
      </div>
      <div class="legend-row">
        <div style="width:14px;height:14px;border-radius:50%;background:#aaa;border:2px solid #ffdd88;flex-shrink:0;"></div>
        <span>Hub node (multi-pathway)</span>
      </div>
    </div>
    <div style="margin-top:10px;">
      <div style="color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:6px;">Edge types</div>
      <div class="legend-row">
        <svg width="30" height="10" style="flex-shrink:0;"><line x1="0" y1="5" x2="30" y2="5" stroke="#6a8aff" stroke-width="2" marker-end="url(#arr)"/><defs><marker id="arr" markerWidth="4" markerHeight="4" refX="3" refY="2" orient="auto"><path d="M0,0 L4,2 L0,4 Z" fill="#6a8aff"/></marker></defs></svg>
        <span>Activating</span>
      </div>
      <div class="legend-row">
        <svg width="30" height="10" style="flex-shrink:0;"><line x1="0" y1="5" x2="30" y2="5" stroke="rgba(255,68,68,0.8)" stroke-width="2" stroke-dasharray="4,3"/></svg>
        <span style="color:#ff8888;">Inhibitory</span>
      </div>
    </div>
    <div style="margin-top:10px;" id="legend-pathways"></div>
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
const PATHWAY_COLORS = {pathway_colors_json};
const PATHWAY_COLORS_CB = {{
  cgas_sting:   '#E69F00',
  apoptosis:    '#D55E00',
  necroptosis:  '#009E73',
  inflammasome: '#F0E442',
  osteoclast:   '#0072B2',
  osteoblast:   '#56B4E9',
  oxysterol:    '#CC79A7',
}};
const LAYER_Z = {layer_z_json};
const PATHWAY_KEYS = {pathway_keys_json};

const state = {{
  visiblePathways: new Set(Object.keys(DATA.pathways)),
  showLabels: true,
  showArrows: true,
  colorblindMode: false,
  showLayerLabels: true,
  filterSignificant: false,
  filterHubs: false,
  highlightNodes: new Set(),
  highlightLinks: new Set(),
  selectedNode: null,
  paused: false,
  searchMatches: new Set(),
}};

const _layerLabelSprites = [];

function getPathwayColor(key) {{
  return (state.colorblindMode ? PATHWAY_COLORS_CB : PATHWAY_COLORS)[key] || '#888';
}}

// ============================================================
// Helpers
// ============================================================
// Cap color scale at ±3 so biologically meaningful FCs (1-2) aren't washed out
const LFC_CAP = 3.0;

function lfc2color(lfc) {{
  if (lfc === null || lfc === undefined || isNaN(lfc)) return '#555555';
  const t = Math.max(-1, Math.min(1, lfc / LFC_CAP));
  if (t <= 0) {{
    const s = -t;
    return `rgb(${{Math.round(33  + (247-33) *(1-s))}},`
         + `${{Math.round(102 + (247-102)*(1-s))}},`
         + `${{Math.round(172 + (247-172)*(1-s))}})`;
  }} else {{
    return `rgb(${{Math.round(247 + (178-247)*t)}},`
         + `${{Math.round(247 + (24 -247)*t)}},`
         + `${{Math.round(247 + (43 -247)*t)}})`;
  }}
}}

function nodeVal(node) {{
  switch (node.node_status) {{
    case 'Present; significant':     return 12;
    case 'Present; not significant': return 5;
    case 'Present; padj NA':         return 5;
    default:                         return 2;
  }}
}}

function isNodeVisible(node) {{
  if (state.filterSignificant && node.node_status !== 'Present; significant') return false;
  if (state.filterHubs && !node.is_hub) return false;
  return node.pathways.some(p => state.visiblePathways.has(p));
}}

function fmt(v, d=3) {{
  if (v === null || v === undefined || (typeof v === 'number' && isNaN(v))) return 'NA';
  return (+v).toFixed(d);
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
// Graph Data Filtering
// ============================================================
function getGraphData() {{
  const visNodes = DATA.nodes.filter(isNodeVisible);
  const visIds = new Set(visNodes.map(n => n.id));
  const visLinks = DATA.links.filter(lk => {{
    const s = typeof lk.source === 'object' ? lk.source.id : lk.source;
    const t = typeof lk.target === 'object' ? lk.target.id : lk.target;
    return visIds.has(s) && visIds.has(t) && state.visiblePathways.has(lk.pathway);
  }});
  return {{ nodes: visNodes, links: visLinks }};
}}

// ============================================================
// THREE.js node object builder
// ============================================================
function makeTextSprite(T, text, fontSize) {{
  const canvas = document.createElement('canvas');
  const scale = 3;
  // Wide enough for longest labels (TNFRSF11A, β-catenin, etc.)
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

function getThree() {{
  return (typeof THREE !== 'undefined') ? THREE : null;
}}

function buildNodeObject(node) {{
  const T = getThree();
  if (!T) return null;  // fall back to default sphere rendering

  const group = new T.Group();
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

  // Search match: bright white outer ring
  if (state.searchMatches.has(node.id)) {{
    const searchRingGeo = new T.TorusGeometry(radius + 3.5, 0.8, 8, 24);
    const searchRingMat = new T.MeshBasicMaterial({{ color: 0xffffff, transparent: true, opacity: 0.95 }});
    group.add(new T.Mesh(searchRingGeo, searchRingMat));
  }}

  if (state.showLabels) {{
    const sprite = makeTextSprite(T, node.label || node.id, node.is_hub ? 36 : 32);
    sprite.position.set(radius + 1, radius + 1, 0);
    group.add(sprite);

    // log2FC direction glyph for significant nodes
    if (node.node_status === 'Present; significant' && node.log2FC !== null && node.log2FC !== undefined) {{
      const arrowCanvas = document.createElement('canvas');
      const sc2 = 3;
      arrowCanvas.width = 64 * sc2;
      arrowCanvas.height = 64 * sc2;
      const actx = arrowCanvas.getContext('2d');
      actx.scale(sc2, sc2);
      actx.clearRect(0, 0, 64, 64);
      const glyph = node.log2FC > 0 ? '\u25b2' : '\u25bc';
      const glyphColor = node.log2FC > 0 ? '#ff6666' : '#6699ff';
      actx.font = 'bold 36px serif';
      actx.fillStyle = glyphColor;
      actx.textAlign = 'center';
      actx.fillText(glyph, 32, 44);
      const arrowTex = new T.CanvasTexture(arrowCanvas);
      arrowTex.needsUpdate = true;
      const arrowMat = new T.SpriteMaterial({{ map: arrowTex, transparent: true, depthWrite: false, depthTest: false }});
      const arrowSprite = new T.Sprite(arrowMat);
      arrowSprite.scale.set(20, 20, 1);
      arrowSprite.position.set(radius + 1, -(radius + 4), 0);
      group.add(arrowSprite);
    }}
  }}

  return group;
}}

// ============================================================
// Graph Initialization
// ============================================================
let Graph = null;

function buildGraph() {{
  const container = document.getElementById('graph-container');

  // 3d-force-graph bundles Three.js internally. After init we grab it from the renderer
  // so nodeThreeObject uses the SAME Three.js instance (no duplicate conflict).
  Graph = ForceGraph3D()(container)
    .backgroundColor('#0a0a1a')
    .graphData(getGraphData())
    .nodeId('id')
    .nodeVal(n => nodeVal(n))
    .nodeColor(n => {{
      if (state.highlightNodes.size > 0 && !state.highlightNodes.has(n.id))
        return 'rgba(60,60,80,0.2)';
      return lfc2color(n.log2FC);
    }})
    .nodeOpacity(0.9)
    .nodeThreeObject(buildNodeObject)
    .nodeThreeObjectExtend(false)
    .linkSource('source')
    .linkTarget('target')
    .linkColor(lk => {{
      if (state.highlightLinks.size > 0 && !state.highlightLinks.has(lk))
        return 'rgba(40,40,60,0.12)';
      if (lk.edge_sign === 'inhibitory') return 'rgba(255,68,68,0.62)';
      return getPathwayColor(lk.pathway) || '#666';
    }})
    .linkWidth(lk => {{
      const srcId = typeof lk.source === 'object' ? lk.source.id : lk.source;
      const srcNode = nodeMap[srcId];
      const lfc = srcNode && srcNode.log2FC !== null ? Math.abs(srcNode.log2FC) : 0;
      // Inhibitory edges: slightly thinner base; activating: core=2, output=1
      const base = lk.edge_sign === 'inhibitory' ? 1.5 : (lk.edge_class === 'core' ? 2 : 1);
      const boost = 1 + Math.min(lfc / LFC_CAP, 1);
      const w = base * boost;
      return state.highlightLinks.size > 0 && !state.highlightLinks.has(lk) ? w * 0.2 : w;
    }})
    .linkDirectionalArrowLength(lk => {{
      if (!state.showArrows) return 0;
      if (lk.edge_sign === 'inhibitory') return 0;
      return 4;
    }})
    .linkDirectionalArrowRelPos(1)
    .linkDirectionalParticles(2)
    .linkDirectionalParticleWidth(lk => lk.edge_class === 'core' ? 2 : 1)
    .linkDirectionalParticleColor(lk => getPathwayColor(lk.pathway) || '#888')
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

  // Tune default forces via built-in d3Force API
  // (3d-force-graph creates charge, link, x, y, z forces automatically)
  try {{
    const chargeFn = Graph.d3Force('charge');
    if (chargeFn) chargeFn.strength(-120);

    const linkFn = Graph.d3Force('link');
    if (linkFn) linkFn.distance(50).strength(0.65);

    // Gently bias nodes toward their layer z and initial x,y
    const zFn = Graph.d3Force('z');
    if (zFn) zFn.z(n => n.z_target).strength(0.08);

    const xFn = Graph.d3Force('x');
    if (xFn) xFn.x(n => n.x_initial).strength(0.04);

    const yFn = Graph.d3Force('y');
    if (yFn) yFn.y(n => n.y_initial).strength(0.04);
  }} catch(e) {{ console.warn('Force setup error:', e); }}

  // Add ambient + directional lights
  try {{
    const T2 = getThree();
    if (T2) {{
      const scene = Graph.scene();
      scene.add(new T2.AmbientLight(0xffffff, 0.6));
      const dirLight = new T2.DirectionalLight(0xffffff, 0.8);
      dirLight.position.set(300, 300, 300);
      scene.add(dirLight);

      // Add floating pathway layer labels
      _layerLabelSprites.length = 0;
      PATHWAY_KEYS.forEach(pw => {{
        const pdata = DATA.pathways[pw];
        if (!pdata) return;
        const zPos = LAYER_Z[pw] || 0;
        const color = getPathwayColor(pw);

        const lCanvas = document.createElement('canvas');
        lCanvas.width = 1024;
        lCanvas.height = 96;
        const lctx = lCanvas.getContext('2d');
        lctx.clearRect(0, 0, 1024, 96);
        // Parse hex color to rgba at 40% opacity
        lctx.font = 'bold 52px "Segoe UI", sans-serif';
        lctx.fillStyle = color;
        lctx.globalAlpha = 0.40;
        lctx.fillText(pdata.display_name, 8, 66);
        lctx.globalAlpha = 1.0;

        const lTex = new T2.CanvasTexture(lCanvas);
        lTex.needsUpdate = true;
        const lMat = new T2.SpriteMaterial({{ map: lTex, transparent: true, depthWrite: false, depthTest: false }});
        const lSprite = new T2.Sprite(lMat);
        lSprite.scale.set(400, 60, 1);
        lSprite.position.set(-600, 400, zPos);
        lSprite.visible = state.showLayerLabels;
        scene.add(lSprite);
        _layerLabelSprites.push(lSprite);
      }});
    }}
  }} catch(e) {{ console.warn('Scene setup error:', e); }}

  updateStats();
}}

// ============================================================
// Interaction
// ============================================================
const tooltipEl = document.getElementById('tooltip');

function negLog10(p) {{
  if (p === null || p === undefined || p <= 0) return 'NA';
  return (-Math.log10(p)).toFixed(2);
}}

function onNodeHover(node) {{
  if (!node) {{ tooltipEl.style.display = 'none'; return; }}
  const lfc = fmt(node.log2FC);
  const padj = fmt(node.padj, 4);
  const nl10 = negLog10(node.padj);
  const matchNote = (node.matched_as && node.matched_as !== node.id)
    ? `<span style="color:#888"> (matched as ${{node.matched_as}})</span>` : '';
  const lfcClass = node.log2FC > 0 ? 'tt-up' : (node.log2FC < 0 ? 'tt-down' : '');
  const sigClass = node.node_status === 'Present; significant' ? 'tt-sig' : '';
  tooltipEl.innerHTML = `
    <div class="tt-gene">${{node.label || node.id}}${{matchNote}}</div>
    <div class="tt-row"><span class="tt-key">Symbol:</span><span class="tt-val">${{node.id}}</span></div>
    <div class="tt-row"><span class="tt-key">Pathways:</span><span class="tt-val">${{node.pathways.join(', ')}}</span></div>
    <div class="tt-row"><span class="tt-key">Module:</span><span class="tt-val">${{(node.modules||[]).join(', ')||'NA'}}</span></div>
    <div class="tt-row"><span class="tt-key">Compartment:</span><span class="tt-val">${{(node.compartments||[]).join(', ')||'NA'}}</span></div>
    <div class="tt-row"><span class="tt-key">log2FC:</span><span class="tt-val ${{lfcClass}}">${{lfc}}</span></div>
    <div class="tt-row"><span class="tt-key">padj:</span><span class="tt-val ${{sigClass}}">${{padj}}</span></div>
    <div class="tt-row"><span class="tt-key">−log10(padj):</span><span class="tt-val ${{sigClass}}">${{nl10}}</span></div>
    <div class="tt-row"><span class="tt-key">Status:</span><span class="tt-val">${{node.node_status}}</span></div>
  `;
  tooltipEl.style.display = 'block';
}}

document.addEventListener('mousemove', e => {{
  if (tooltipEl.style.display !== 'none') {{
    let left = e.clientX + 16;
    let top = Math.min(e.clientY + 16, window.innerHeight - 220);
    const tooltipWidth = tooltipEl.offsetWidth || 280;
    if (left + tooltipWidth > window.innerWidth - 8) {{
      left = e.clientX - tooltipWidth - 16;
    }}
    tooltipEl.style.left = left + 'px';
    tooltipEl.style.top = top + 'px';
  }}
}});

function onNodeClick(node) {{
  if (!node) return;
  state.selectedNode = node;
  const neighbors = neighborNodes[node.id] || new Set();
  state.highlightNodes = new Set([node.id, ...neighbors]);
  state.highlightLinks = new Set(nodeLinks[node.id] || []);
  Graph.graphData(getGraphData());
  // Fly camera toward node
  const dist = 180;
  Graph.cameraPosition(
    {{ x: (node.x||0) + dist * 0.7, y: (node.y||0) + dist * 0.3, z: (node.z||0) + dist }},
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
}}

function showDetailPanel(node) {{
  const panel = document.getElementById('detail-panel');
  const content = document.getElementById('detail-content');
  panel.style.display = 'block';
  const lfc = fmt(node.log2FC);
  const padj = fmt(node.padj, 4);
  const matchNote = (node.matched_as && node.matched_as !== node.id)
    ? `<div><span class="dp-key">Matched as:</span> <span class="dp-val">${{node.matched_as}}</span></div>` : '';
  const inLinks = DATA.links.filter(l => {{
    const t = typeof l.target === 'object' ? l.target.id : l.target;
    return t === node.id;
  }});
  const outLinks = DATA.links.filter(l => {{
    const s = typeof l.source === 'object' ? l.source.id : l.source;
    return s === node.id;
  }});
  function neighborRow(link, isInput) {{
    const oid = isInput
      ? (typeof link.source === 'object' ? link.source.id : link.source)
      : (typeof link.target === 'object' ? link.target.id : link.target);
    const other = nodeMap[oid];
    if (!other) return '';
    const olfc = fmt(other.log2FC);
    const col = other.log2FC > 0 ? '#ff8888' : (other.log2FC < 0 ? '#88aaff' : '#bbb');
    return `<div class="neighbor-row">
      <span style="cursor:pointer;text-decoration:underline;color:#88aaff;"
            onclick="onNodeClick(nodeMap['${{oid}}'])">${{other.label || oid}}</span>
      <span style="color:${{col}}">${{olfc}} <span style="color:#555;font-size:10px">[${{link.edge_pathway}}]</span></span>
    </div>`;
  }}
  const notes = [];
  node.pathways.forEach(p => {{
    (DATA.pathways[p]?.validation_notes || []).forEach(n => notes.push(n));
  }});
  const notesHtml = [...new Set(notes)].slice(0,6).map(n =>
    `<div class="note-item">${{n}}</div>`).join('');
  content.innerHTML = `
    <h2>${{node.label || node.id}}</h2>
    <div><span class="dp-key">Symbol:</span> <span class="dp-val">${{node.id}}</span></div>
    ${{matchNote}}
    <div><span class="dp-key">Pathways:</span> <span class="dp-val">${{node.pathways.map(p=>DATA.pathways[p]?.display_name||p).join(', ')}}</span></div>
    <div><span class="dp-key">Modules:</span> <span class="dp-val">${{(node.modules||[]).join(', ')||'NA'}}</span></div>
    <div><span class="dp-key">Compartments:</span> <span class="dp-val">${{(node.compartments||[]).join(', ')||'NA'}}</span></div>
    <div><span class="dp-key">log2FC:</span> <span class="dp-val">${{lfc}}</span></div>
    <div><span class="dp-key">padj:</span> <span class="dp-val">${{padj}}</span></div>
    <div><span class="dp-key">−log10(padj):</span> <span class="dp-val">${{negLog10(node.padj)}}</span></div>
    <div><span class="dp-key">Status:</span> <span class="dp-val">${{node.node_status}}</span></div>
    ${{inLinks.length ? `<div class="dp-section">Upstream inputs (${{inLinks.length}})</div>${{inLinks.map(l=>neighborRow(l,true)).join('')}}` : ''}}
    ${{outLinks.length ? `<div class="dp-section">Downstream targets (${{outLinks.length}})</div>${{outLinks.map(l=>neighborRow(l,false)).join('')}}` : ''}}
    ${{notesHtml ? `<div class="dp-section">Validation notes</div>${{notesHtml}}` : ''}}
  `;
}}

// ============================================================
// Update & Refresh
// ============================================================
let _recalcTimer = null;
function updateGraph() {{
  if (!Graph) return;
  // Flash "Recalculating..." in the stats bar
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
  // Pathway toggles
  const pwBox = document.getElementById('pathway-checkboxes');
  Object.entries(DATA.pathways).forEach(([key, pdata]) => {{
    const color = getPathwayColor(key);
    const row = document.createElement('div');
    row.className = 'check-row';
    row.style.justifyContent = 'space-between';
    row.innerHTML = `
      <div style="display:flex;align-items:center;gap:7px;flex:1;min-width:0;">
        <input type="checkbox" id="pw-${{key}}" checked>
        <div class="pathway-dot" style="background:${{color}}"></div>
        <label for="pw-${{key}}" class="check-label" style="overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">${{pdata.display_name}}</label>
      </div>
      <button class="isolate-btn" data-key="${{key}}"
        style="font-size:9px;padding:2px 5px;margin-left:4px;flex-shrink:0;
               background:rgba(40,50,120,0.6);border:1px solid rgba(100,120,200,0.3);
               border-radius:3px;color:#aac;cursor:pointer;">only</button>
    `;
    row.querySelector('input').addEventListener('change', e => {{
      if (e.target.checked) state.visiblePathways.add(key);
      else state.visiblePathways.delete(key);
      updateGraph();
    }});
    row.querySelector('.isolate-btn').addEventListener('click', () => isolatePathway(key));
    pwBox.appendChild(row);
  }});

  document.getElementById('btn-show-all').addEventListener('click', () => {{
    Object.keys(DATA.pathways).forEach(k => {{
      state.visiblePathways.add(k);
      const cb = document.getElementById(`pw-${{k}}`);
      if (cb) cb.checked = true;
    }});
    updateGraph();
  }});

  document.getElementById('btn-isolate-active').addEventListener('click', () => {{
    if (Graph) {{
      Graph.d3ReheatSimulation();
      const visNodes = Graph.graphData().nodes;
      if (visNodes.length > 0) {{
        const cx = visNodes.reduce((s, n) => s + (n.x || 0), 0) / visNodes.length;
        const cy = visNodes.reduce((s, n) => s + (n.y || 0), 0) / visNodes.length;
        const cz = visNodes.reduce((s, n) => s + (n.z || 0), 0) / visNodes.length;
        Graph.cameraPosition(
          {{ x: cx + 500, y: cy + 200, z: cz + 600 }},
          {{ x: cx, y: cy, z: cz }}, 800
        );
      }}
    }}
  }});

  // Legend pathways
  const legPW = document.getElementById('legend-pathways');
  legPW.innerHTML = '<div style="color:#7080c0;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:6px;">Pathways</div>';
  Object.entries(DATA.pathways).forEach(([key, pdata]) => {{
    const color = getPathwayColor(key);
    legPW.innerHTML += `<div class="legend-row pw-legend-dot-${{key}}">
      <div class="pw-dot-${{key}}" style="width:10px;height:10px;border-radius:50%;background:${{color}};flex-shrink:0;"></div>
      <span>${{pdata.display_name}}</span></div>`;
  }});

  document.getElementById('filter-significant').addEventListener('change', e => {{
    state.filterSignificant = e.target.checked; updateGraph();
  }});
  document.getElementById('filter-hubs').addEventListener('change', e => {{
    state.filterHubs = e.target.checked; updateGraph();
  }});
  document.getElementById('toggle-labels').addEventListener('change', e => {{
    state.showLabels = e.target.checked;
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
  document.getElementById('toggle-arrows').addEventListener('change', e => {{
    state.showArrows = e.target.checked;
    if (Graph) Graph.linkDirectionalArrowLength(lk => {{
      if (!state.showArrows) return 0;
      if (lk.edge_sign === 'inhibitory') return 0;
      return 4;
    }});
  }});
  document.getElementById('toggle-light').addEventListener('change', e => {{
    document.body.classList.toggle('light-mode', e.target.checked);
    if (Graph) Graph.backgroundColor(e.target.checked ? '#f4f4f8' : '#0a0a1a');
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
  document.getElementById('toggle-colorblind').addEventListener('change', e => {{
    state.colorblindMode = e.target.checked;
    // Update pathway dot colors in sidebar
    Object.keys(DATA.pathways).forEach(k => {{
      const dot = document.querySelector(`#pathway-checkboxes .pathway-dot[style*="background"]`);
      // Update all pathway dots
    }});
    document.querySelectorAll('#pathway-checkboxes .check-row').forEach(row => {{
      const cb = row.querySelector('input[type=checkbox]');
      if (!cb) return;
      const key = cb.id.replace('pw-', '');
      const dotEl = row.querySelector('.pathway-dot');
      if (dotEl) dotEl.style.background = getPathwayColor(key);
    }});
    // Update legend pathway dots
    Object.keys(DATA.pathways).forEach(k => {{
      const dotEl = document.querySelector(`.pw-dot-${{k}}`);
      if (dotEl) dotEl.style.background = getPathwayColor(k);
    }});
    if (Graph) {{
      Graph.linkColor(lk => {{
        if (state.highlightLinks.size > 0 && !state.highlightLinks.has(lk)) return 'rgba(40,40,60,0.12)';
        if (lk.edge_sign === 'inhibitory') return 'rgba(255,68,68,0.62)';
        return getPathwayColor(lk.pathway) || '#666';
      }});
      Graph.linkDirectionalParticleColor(lk => getPathwayColor(lk.pathway) || '#888');
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
  }});
  document.getElementById('btn-screenshot').addEventListener('click', doScreenshot);
  document.getElementById('btn-export-csv').addEventListener('click', exportCSV);
  document.getElementById('btn-export-edge-csv').addEventListener('click', exportEdgeCSV);

  // Search
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
      (n.label && n.label.toLowerCase().includes(q))
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
      Graph.cameraPosition(
        {{ x: x + 200, y: y + 80, z: z + 300 }},
        {{ x, y, z }}, 800
      );
      pulseNodes(matches);
    }}
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
  }});
}}

// ============================================================
// Camera
// ============================================================
function resetCamera() {{
  if (Graph) Graph.cameraPosition({{ x: 0, y: 0, z: 1800 }}, {{ x: 0, y: 0, z: 0 }}, 1000);
}}

// ============================================================
// Keyboard
// ============================================================
document.addEventListener('keydown', e => {{
  if (e.target.tagName === 'INPUT') return;
  const key = e.key.toUpperCase();
  if (key === 'R') {{ resetCamera(); return; }}
  if (key === 'L') {{
    const cb = document.getElementById('toggle-labels');
    cb.checked = !cb.checked;
    state.showLabels = cb.checked;
    if (Graph) {{ Graph.nodeThreeObject(buildNodeObject); Graph.refresh(); }}
    return;
  }}
  if (key === 'H') {{
    const cb = document.getElementById('filter-hubs');
    cb.checked = !cb.checked;
    state.filterHubs = cb.checked; updateGraph(); return;
  }}
  if (key === ' ') {{
    e.preventDefault();
    if (!Graph) return;
    if (state.paused) {{ Graph.resumeAnimation(); Graph.d3ReheatSimulation(); }}
    else {{ Graph.pauseAnimation(); }}
    state.paused = !state.paused; return;
  }}
  const digit = parseInt(key);
  if (digit >= 1 && digit <= 7 && Graph) {{
    const pwKey = PATHWAY_KEYS[digit - 1];
    if (pwKey !== undefined) {{
      const z = LAYER_Z[pwKey] || 0;
      Graph.cameraPosition({{ x: 600, y: 200, z: z + 700 }}, {{ x: 0, y: 0, z }}, 800);
    }}
  }}
}});

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

  // Composite legend in bottom-right corner
  const legX = W - 220, legY = H - 40;
  const lineH = 18, pad = 10;
  const pwKeys = Object.keys(DATA.pathways);
  const rows = pwKeys.length + 5; // pathway rows + header + lfc bar rows + node rows
  const legH = rows * lineH + pad * 2 + 30;
  const legW = 200;
  const lx = W - legW - 12, ly = H - legH - 12;

  ctx.fillStyle = 'rgba(10,10,30,0.88)';
  ctx.beginPath();
  ctx.roundRect ? ctx.roundRect(lx, ly, legW, legH, 6) : ctx.rect(lx, ly, legW, legH);
  ctx.fill();

  ctx.font = 'bold 10px sans-serif';
  ctx.fillStyle = '#7080c0';
  ctx.textAlign = 'left';
  let cy2 = ly + pad + 12;
  ctx.fillText('LEGEND', lx + pad, cy2); cy2 += lineH;

  // LFC bar
  ctx.font = '9px sans-serif';
  ctx.fillStyle = '#7080c0';
  ctx.fillText('log2FC', lx + pad, cy2); cy2 += 4;
  const barGrad = ctx.createLinearGradient(lx + pad, 0, lx + pad + 120, 0);
  barGrad.addColorStop(0, '#2166ac'); barGrad.addColorStop(0.5, '#f7f7f7'); barGrad.addColorStop(1, '#b2182b');
  ctx.fillStyle = barGrad;
  ctx.fillRect(lx + pad, cy2, 120, 8); cy2 += 12;
  ctx.fillStyle = '#888';
  ctx.font = '8px sans-serif';
  ctx.fillText('−', lx + pad, cy2);
  ctx.textAlign = 'center';
  ctx.fillText('0', lx + pad + 60, cy2);
  ctx.textAlign = 'right';
  ctx.fillText('+', lx + pad + 120, cy2);
  ctx.textAlign = 'left';
  cy2 += lineH;

  // Node size key
  ctx.font = '9px sans-serif';
  ctx.fillStyle = '#aaa';
  ctx.fillText('● Significant  ○ Other  · Missing', lx + pad, cy2); cy2 += lineH;

  // Pathway colors
  pwKeys.forEach(k => {{
    const pdata = DATA.pathways[k];
    const color = getPathwayColor(k);
    ctx.beginPath();
    ctx.arc(lx + pad + 5, cy2 - 4, 5, 0, Math.PI * 2);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.fillStyle = '#dde';
    ctx.font = '9px sans-serif';
    ctx.fillText(pdata.display_name.slice(0, 22), lx + pad + 14, cy2);
    cy2 += lineH;
  }});

  const url = offscreen.toDataURL('image/png');
  const a = document.createElement('a');
  a.href = url; a.download = 'pathway_network_3d.png'; a.click();
}}

function exportCSV() {{
  if (!Graph) return;
  const visNodes = Graph.graphData().nodes || [];
  const header = 'id,label,pathways,modules,compartments,log2FC,padj,node_status,is_hub,matched_as';
  const rows = visNodes.map(n => [
    n.id, n.label,
    (n.pathways||[]).join(';'),
    (n.modules||[]).join(';'),
    (n.compartments||[]).join(';'),
    n.log2FC !== null && n.log2FC !== undefined ? n.log2FC : '',
    n.padj !== null && n.padj !== undefined ? n.padj : '',
    n.node_status,
    n.is_hub ? 'TRUE' : 'FALSE',
    n.matched_as || '',
  ].map(v => `"${{String(v).replace(/"/g,'""')}}"`).join(','));
  const csv = [header, ...rows].join('\n');
  const blob = new Blob([csv], {{ type: 'text/csv' }});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'pathway_network_nodes.csv'; a.click();
  URL.revokeObjectURL(url);
}}

function exportEdgeCSV() {{
  if (!Graph) return;
  const visLinks = Graph.graphData().links || [];
  const header = 'source,target,edge_class,edge_sign,edge_pathway,pathway,source_log2FC,source_padj,target_log2FC,target_padj';
  const rows = visLinks.map(lk => {{
    const srcId = typeof lk.source === 'object' ? lk.source.id : lk.source;
    const tgtId = typeof lk.target === 'object' ? lk.target.id : lk.target;
    const srcNode = nodeMap[srcId] || {{}};
    const tgtNode = nodeMap[tgtId] || {{}};
    return [
      srcId, tgtId,
      lk.edge_class || '', lk.edge_sign || 'activating',
      lk.edge_pathway || '', lk.pathway || '',
      srcNode.log2FC !== null && srcNode.log2FC !== undefined ? srcNode.log2FC : '',
      srcNode.padj !== null && srcNode.padj !== undefined ? srcNode.padj : '',
      tgtNode.log2FC !== null && tgtNode.log2FC !== undefined ? tgtNode.log2FC : '',
      tgtNode.padj !== null && tgtNode.padj !== undefined ? tgtNode.padj : '',
    ].map(v => `"${{String(v).replace(/"/g,'""')}}"`).join(',');
  }});
  const csv = [header, ...rows].join('\n');
  const blob = new Blob([csv], {{ type: 'text/csv' }});
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = 'pathway_network_edges.csv'; a.click();
  URL.revokeObjectURL(url);
}}

// ============================================================
// Isolate pathway
// ============================================================
function isolatePathway(key) {{
  Object.keys(DATA.pathways).forEach(k => {{
    const cb = document.getElementById(`pw-${{k}}`);
    if (k === key) {{
      state.visiblePathways.add(k);
      if (cb) cb.checked = true;
    }} else {{
      state.visiblePathways.delete(k);
      if (cb) cb.checked = false;
    }}
  }});
  updateGraph();
  // Fly camera to this pathway's z-layer
  const z = LAYER_Z[key] || 0;
  if (Graph) Graph.cameraPosition({{ x: 500, y: 150, z: z + 700 }}, {{ x: 0, y: 0, z }}, 800);
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
    const t = frame / totalFrames;
    const scale = 1 + 0.6 * Math.sin(t * Math.PI); // swell then return
    Graph.graphData().nodes.forEach(n => {{
      if (!ids.has(n.id) || !n.__threeObj) return;
      n.__threeObj.scale.setScalar(scale);
    }});
    if (frame < totalFrames) requestAnimationFrame(animate);
    else {{
      // Reset scale
      Graph.graphData().nodes.forEach(n => {{
        if (n.__threeObj) n.__threeObj.scale.setScalar(1);
      }});
    }}
  }}
  requestAnimationFrame(animate);
}}

// ============================================================
// Legend toggle
// ============================================================
function toggleLegend() {{
  const body = document.getElementById('legend-body');
  const icon = document.getElementById('legend-toggle-icon');
  const collapsed = body.classList.toggle('collapsed');
  icon.textContent = collapsed ? '▼' : '▲';
}}

// ============================================================
// Boot
// ============================================================
window.addEventListener('load', () => {{
  const loadMsg = document.getElementById('loading-msg');
  loadMsg.textContent = `Loading ${{DATA.nodes.length}} nodes, ${{DATA.links.length}} edges...`;

  // Check libraries loaded
  if (typeof ForceGraph3D === 'undefined') {{
    loadMsg.textContent = 'ERROR: ForceGraph3D library not loaded. Check internet connection and reload.';
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

# ---------------------------------------------------------------------------
# 7. Assemble and write HTML
# ---------------------------------------------------------------------------

PATHWAY_COLORS = {
    'cgas_sting':   '#6A5ACD',
    'apoptosis':    '#A33A3A',
    'necroptosis':  '#2F7F56',
    'inflammasome': '#FF8C00',
    'osteoclast':   '#1F4E79',
    'osteoblast':   '#1E90FF',
    'oxysterol':    '#8B6914',
}

def generate_html(pathways, nodes, links, title, p_cutoff, max_lfc, output_path):
    # Build pathways metadata for DATA
    pathways_meta = {}
    for pname, pdata in pathways.items():
        pathways_meta[pname] = {
            'display_name': pdata['display_name'],
            'color': PATHWAY_COLORS.get(pname, '#888888'),
            'validation_notes': pdata['validation_notes'],
        }

    data_obj = {
        'title': title,
        'p_cutoff': p_cutoff,
        'max_lfc': max_lfc,
        'pathways': pathways_meta,
        'nodes': nodes,
        'links': links,
    }

    data_json = json.dumps(data_obj, ensure_ascii=False, indent=None)

    pathway_colors_json = json.dumps(PATHWAY_COLORS, ensure_ascii=False)

    layer_z_json = json.dumps(PATHWAY_Z, ensure_ascii=False)

    # Ordered list of pathway keys for keyboard shortcuts 1-7
    pathway_keys = list(PATHWAY_Z.keys())
    pathway_keys_json = json.dumps(pathway_keys)

    html = HTML_TEMPLATE.format(
        title=title,
        p_cutoff=p_cutoff,
        max_lfc=round(max_lfc, 2),
        data_json=data_json,
        pathway_colors_json=pathway_colors_json,
        layer_z_json=layer_z_json,
        pathway_keys_json=pathway_keys_json,
    )

    Path(output_path).write_text(html, encoding='utf-8')
    return output_path


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------

def main():
    script_path = Path(__file__).resolve()
    script_dir = Path(__file__).resolve().parent
    default_r_script = script_dir / 'pathway_mapper_v4.R'
    default_deseq2_csv = script_dir / 'EXAMPLE_DATA_DESeq2_Free_v_Pygo_23.csv'

    parser = argparse.ArgumentParser(
        description='Build 3D pathway network HTML from R pathway definitions + DESeq2 CSV'
    )
    parser.add_argument(
        '--r-script',
        default=str(default_r_script),
        help='Path to pathway_mapper_v4.R (defaults to the copy next to this script)'
    )
    parser.add_argument(
        '--deseq2-csv',
        default=str(default_deseq2_csv),
        help='Path to DESeq2 results CSV (defaults to DESeq2_FILE_PATH_HERE.csv next to this script)'
    )
    parser.add_argument('--output', default='pathway_network_3d.html', help='Output HTML path')
    ###CHANGE THE NAMES BELOW TO MATCH THE ANALYSIS###
    parser.add_argument('--title', default='Free vs Pygo 2023— 3D Pathway Network',
                        help='Visualization title')
    parser.add_argument('--p-cutoff', type=float, default=0.05, help='Significance threshold for padj')
    args = parser.parse_args()

    def error_for_missing_path(label, supplied_path, default_path, default_var_name, cli_flag):
        supplied_path = Path(supplied_path)
        using_default = supplied_path == default_path
        edit_hint = f'edit {script_path} and update `{default_var_name}` in `main()`'
        if using_default:
            parser.error(
                f'{label} not found: {supplied_path}\n'
                f'You are using the built-in default path.\n'
                f'To fix it, either {edit_hint}, or run the script with {cli_flag} "/full/path/to/file".\n'
                f'Spyder example: %runfile "{script_path}" --args {cli_flag} "/full/path/to/file"'
            )
        parser.error(
            f'{label} not found: {supplied_path}\n'
            f'Fix the path you passed with {cli_flag}.\n'
            f'If you want to change the built-in default instead, {edit_hint}.'
        )

    if not Path(args.r_script).exists():
        error_for_missing_path('R script', args.r_script, default_r_script, 'default_r_script', '--r-script')
    if not Path(args.deseq2_csv).exists():
        error_for_missing_path(
            'DESeq2 CSV',
            args.deseq2_csv,
            default_deseq2_csv,
            'default_deseq2_csv',
            '--deseq2-csv'
        )

    print(f"[1/5] Parsing R pathway registry from: {args.r_script}")
    pathways = parse_pathway_registry(args.r_script)
    print(f"      Found {len(pathways)} pathways: {', '.join(pathways.keys())}")
    for pname, pdata in pathways.items():
        print(f"      {pname}: {len(pdata['nodes'])} nodes, {len(pdata['edges'])} edges")

    print(f"\n[2/5] Reading DESeq2 data from: {args.deseq2_csv}")
    de_table = read_deseq2(args.deseq2_csv)
    print(f"      {len(de_table)} genes loaded")

    print(f"\n[3/5] Building unified network (p-cutoff = {args.p_cutoff})")
    nodes, links = build_network(pathways, de_table, args.p_cutoff)

    # Validation summary
    status_counts = defaultdict(int)
    hub_count = 0
    for n in nodes:
        status_counts[n['node_status']] += 1
        if n['is_hub']:
            hub_count += 1

    print(f"      Total nodes: {len(nodes)}  |  Total edges: {len(links)}")
    print(f"      Hub nodes (multi-pathway): {hub_count}")
    for status, count in sorted(status_counts.items()):
        print(f"      {status}: {count}")

    print(f"\n[4/5] Known hub node validation:")
    known_hubs = [
        'TNF', 'TNFRSF1A', 'TRADD', 'CASP8', 'NFKB1', 'RELA', 'TRAF6',
        'RUNX2', 'SP7', 'TNFSF11', 'ZBP1', 'IFNB1', 'IFNAR1',
    ]
    node_by_id = {n['id']: n for n in nodes}
    for sym in known_hubs:
        if sym in node_by_id:
            n = node_by_id[sym]
            pathways_str = ', '.join(n['pathways'])
            matched = f"matched as {n['matched_as']}" if n['matched_as'] and n['matched_as'] != sym else ''
            lfc_str = f"log2FC={n['log2FC']:.3f}" if n['log2FC'] is not None else 'Missing'
            matched_note = f" [{matched}]" if matched else ''
            hub_mark = "★ HUB" if n['is_hub'] else ''
            print(f"      {sym:12s} {hub_mark:8s} {lfc_str:18s} {n['node_status']:30s} [{pathways_str}]{matched_note}")
        else:
            print(f"      {sym:12s}  !! NOT FOUND in global node list")

    print(f"\n[5/5] Synonym resolution check (STING1 → TMEM173):")
    if 'STING1' in node_by_id:
        n = node_by_id['STING1']
        print(f"      STING1 matched_as={n['matched_as']}  status={n['node_status']}  log2FC={n['log2FC']}")
    else:
        print("      STING1 NOT FOUND")

    print(f"\n[6/6] Generating HTML → {args.output}")
    max_lfc = compute_max_lfc(nodes)
    print(f"      Max |log2FC| for color scale: {max_lfc}")
    generate_html(pathways, nodes, links, args.title, args.p_cutoff, max_lfc, args.output)
    file_size_mb = Path(args.output).stat().st_size / 1_000_000
    print(f"      Done! {args.output} ({file_size_mb:.1f} MB)")
    print(f"\nOpen in browser: file://{Path(args.output).absolute()}")


if __name__ == '__main__':
    main()
