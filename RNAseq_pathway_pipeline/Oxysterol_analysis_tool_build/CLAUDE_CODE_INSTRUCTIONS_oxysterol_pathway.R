# ============================================================
# CLAUDE CODE INSTRUCTIONS
# Task: Add oxysterol pathway to pathway_mapper_v3.R
# ============================================================
#
# CONTEXT:
# pathway_mapper_v3.R uses a pathway_registry — a named list of 
# functions, each returning a config list with: name, display_name,
# nodes (tribble), edges (tribble), synonyms (named list), and 
# validation_notes (character vector). The rendering engine is 
# generic and does not need modification.
#
# YOUR TASK:
# Add the following pathway_registry entry AFTER the existing 
# pathway_registry[["osteoblast"]] block (which ends around line 513)
# and BEFORE the list_pathways() function (line 516).
#
# Do NOT modify any existing code — only INSERT this new block.
#
# ============================================================

# INSERT THIS BLOCK into pathway_mapper_v3.R:

pathway_registry[["oxysterol"]] <- function() {
  nodes <- tribble(
    ~node_id,     ~symbol_key,  ~label,            ~x,   ~y,   ~module,                    ~compartment,
    # ── Upstream: IFN input & cholesterol source ──
    "IFNB1",      "IFNB1",      "IFN\u03b2",        0.0,  4.0,  "IFN input",                "Secreted",
    "IFNAR1",     "IFNAR1",     "IFNAR1",           1.0,  4.0,  "IFN input",                "Plasma membrane",
    # ── Module 1: Enzymatic biosynthesis ──
    "CH25H",      "CH25H",      "CH25H",            2.5,  4.0,  "Biosynthesis",             "ER",
    "CYP7B1",     "CYP7B1",     "CYP7B1",           4.0,  4.0,  "Biosynthesis",             "ER",
    "CYP27A1",    "CYP27A1",    "CYP27A1",          2.5,  1.0,  "Biosynthesis",             "Mitochondria",
    "CYP11A1",    "CYP11A1",    "CYP11A1",          2.5,  2.5,  "Biosynthesis",             "Mitochondria",
    # ── Module 2A: Hedgehog / osteogenic arm ──
    "SMO",        "SMO",        "SMO",              4.0,  2.5,  "Hedgehog/Osteogenic",      "Plasma membrane",
    "PTCH1",      "PTCH1",      "PTCH1",            4.0,  1.8,  "Hedgehog/Osteogenic",      "Plasma membrane",
    "GLI1",       "GLI1",       "GLI1",             5.5,  2.5,  "Hedgehog/Osteogenic",      "Nucleus",
    "GLI2",       "GLI2",       "GLI2",             5.5,  1.8,  "Hedgehog/Osteogenic",      "Nucleus",
    "HES1",       "HES1",       "HES1",             7.0,  2.5,  "Hedgehog/Osteogenic",      "Nucleus",
    "HEY1",       "HEY1",       "HEY1",             7.0,  1.8,  "Hedgehog/Osteogenic",      "Nucleus",
    "RUNX2",      "RUNX2",      "RUNX2",            8.5,  2.5,  "OB readouts",              "Nucleus",
    "SP7",        "SP7",        "Osterix",          8.5,  1.8,  "OB readouts",              "Nucleus",
    # ── Module 2B: Immune chemotaxis arm ──
    "GPR183",     "GPR183",     "GPR183/EBI2",      5.5,  4.0,  "Immune chemotaxis",        "Plasma membrane",
    "HSD3B7",     "HSD3B7",     "HSD3B7",           5.5,  4.7,  "Immune chemotaxis",        "ER",
    # ── Module 2C: LXR / cholesterol homeostasis ──
    "NR1H3",      "NR1H3",      "LXR\u03b1",        4.0,  0.2,  "LXR/Cholesterol",          "Nucleus",
    "ABCA1",      "ABCA1",      "ABCA1",            5.5,  0.8,  "LXR/Cholesterol",          "Plasma membrane",
    "ABCG1",      "ABCG1",      "ABCG1",            5.5,  0.2,  "LXR/Cholesterol",          "Plasma membrane",
    "TNFSF11",    "TNFSF11",    "RANKL",            5.5, -0.5,  "LXR/Cholesterol",          "Secreted",
    # ── Module 2C cont: SREBP2/SCAP/inflammasome bridge ──
    "INSIG1",     "INSIG1",     "INSIG1",           4.0,  5.0,  "Cholesterol sensing",      "ER",
    "SCAP",       "SCAP",       "SCAP",             5.5,  5.5,  "Cholesterol sensing",      "ER",
    "SREBF2",     "SREBF2",     "SREBP2",           7.0,  5.5,  "Cholesterol sensing",      "Nucleus",
    "HMGCR",      "HMGCR",      "HMGCR",            8.5,  5.5,  "Cholesterol sensing",      "ER",
    # ── Module 3: Degradation ──
    "SULT2B1",    "SULT2B1",    "SULT2B1",          7.0,  4.5,  "Degradation",              "Cytosol"
  )

  edges <- tribble(
    ~from,       ~to,         ~edge_class,  ~edge_pathway,
    # IFN → CH25H (ISG induction)
    "IFNB1",     "IFNAR1",    "core",       "ifn_input",
    "IFNAR1",    "CH25H",     "core",       "ifn_input",
    # Biosynthesis: CH25H → 25-HC → CYP7B1 → 7α,25-diHC
    "CH25H",     "CYP7B1",    "core",       "biosynthesis",
    "CH25H",     "NR1H3",     "core",       "lxr_arm",
    "CH25H",     "INSIG1",    "core",       "cholesterol_sensing",
    # CYP7B1 → GPR183 ligand
    "CYP7B1",    "GPR183",    "output",     "chemotaxis_arm",
    # GPR183 gradient sculpting
    "HSD3B7",    "GPR183",    "core",       "chemotaxis_arm",
    # CYP11A1 → osteogenic oxysterols → SMO
    "CYP11A1",   "SMO",       "core",       "hedgehog_arm",
    # CYP27A1 → 27-OHC → LXR
    "CYP27A1",   "NR1H3",     "core",       "lxr_arm",
    # Hedgehog arm: SMO → GLI → HES/HEY → RUNX2/SP7
    "PTCH1",     "SMO",       "core",       "hedgehog_arm",
    "SMO",       "GLI1",      "core",       "hedgehog_arm",
    "SMO",       "GLI2",      "core",       "hedgehog_arm",
    "GLI1",      "HES1",      "core",       "hedgehog_arm",
    "GLI1",      "HEY1",      "core",       "hedgehog_arm",
    "GLI2",      "HES1",      "core",       "hedgehog_arm",
    "HES1",      "RUNX2",     "output",     "hedgehog_arm",
    "HEY1",      "RUNX2",     "output",     "hedgehog_arm",
    "RUNX2",     "SP7",       "output",     "hedgehog_arm",
    # LXR arm: NR1H3 → efflux transporters + RANKL
    "NR1H3",     "ABCA1",     "output",     "lxr_arm",
    "NR1H3",     "ABCG1",     "output",     "lxr_arm",
    "NR1H3",     "TNFSF11",   "output",     "lxr_arm",
    # Cholesterol sensing: 25-HC → INSIG1 → SCAP → SREBF2 → HMGCR
    "INSIG1",    "SCAP",      "core",       "cholesterol_sensing",
    "SCAP",      "SREBF2",    "core",       "cholesterol_sensing",
    "SREBF2",    "HMGCR",     "output",     "cholesterol_sensing",
    # Degradation
    "SULT2B1",   "CH25H",     "core",       "degradation"
  )

  synonyms <- list(
    CH25H    = c("C25H"),
    CYP7B1   = c("OAH"),
    CYP27A1  = c("CYP27", "CTX"),
    CYP11A1  = c("P450SCC", "CYP11A"),
    SMO      = c("SMOH"),
    PTCH1    = c("PTC", "PTC1"),
    GLI1     = c("GLI"),
    GPR183   = c("EBI2"),
    NR1H3    = c("LXRA", "LXR"),
    ABCA1    = c("ABC1", "TGD"),
    TNFSF11  = c("RANKL", "TRANCE", "OPGL"),
    SREBF2   = c("SREBP2"),
    HMGCR    = c("HMG-COA-R"),
    INSIG1   = c("INSIG-1"),
    HSD3B7   = c("HSD3B7"),
    SULT2B1  = c("SULT2B1B"),
    RUNX2    = c("CBFA1", "OSF2"),
    SP7      = c("OSX", "Osterix"),
    HES1     = c("HES-1"),
    HEY1     = c("HESR1", "HRT1")
  )

  validation_notes <- c(
    "CH25H is an ISG: bridges cGAS-STING/IFN module to oxysterol production. Expect co-expression with ISG signature.",
    "GPR183/EBI2: oxysterol-guided OCP chemotaxis to bone surfaces (Nevius et al. 2015 J Exp Med).",
    "Osteogenic arm: oxysterol \u2192 SMO binding is allosteric, bypassing PTCH1. Blocked by cyclopamine.",
    "27-OHC \u2192 LXR\u03b1 \u2192 RANKL: macrophage-derived oxysterol promotes osteoclastogenesis (Nelson et al. 2011).",
    "25-HC dose paradox: anti-inflammatory at nM (SREBP2 repression) vs pro-inflammatory at \u00b5M (NLRP3 activation).",
    "SULT2B1: sulfonates oxysterols to inactive forms. Low expression prolongs oxysterol signaling.",
    "HSD3B7 degrades 7\u03b1,25-diHC, shaping the GPR183 chemotactic gradient.",
    "Non-enzymatic oxysterols (7-KC, 7\u03b2-OHC) from ROS/necroptosis are metabolite targets, not gene targets.",
    "Metabolomics priorities: 25-HC, 7\u03b1,25-diHC, 7-KC, 27-OHC (Tier 1 for Free vs Fusing disc comparison).",
    "Chicken annotation caveat: CH25H is intronless \u2014 verify Ensembl ID in GRCg7b (may be LOC identifier).",
    "GPR183 annotation: verify in GRCg7b; may be listed as EBI2."
  )

  list(
    name            = "oxysterol",
    display_name    = "Oxysterol Signaling in Bone Fusion",
    nodes           = nodes,
    edges           = edges,
    synonyms        = synonyms,
    validation_notes = validation_notes
  )
}

# ============================================================
# ADDITIONAL INSTRUCTIONS FOR CLAUDE CODE:
# ============================================================
#
# 1. INSERTION POINT:
#    - Open pathway_mapper_v3.R
#    - Find the closing brace of pathway_registry[["osteoblast"]] 
#      (around line 513-514, ends with `}`)
#    - Insert the entire pathway_registry[["oxysterol"]] block above
#      BEFORE the `list_pathways <- function()` line (currently ~516)
#    - Leave a blank line before and after the new block
#
# 2. EDGE PATHWAY PALETTE (OPTIONAL BUT RECOMMENDED):
#    - In the build_pathway_plot() function, there's an 
#      `edge_pathway_palette` named vector (~line 1016-1021).
#    - Add these entries to the palette for consistent coloring:
#
#      "ifn_input"            = "#6A5ACD",
#      "biosynthesis"         = "#2E8B57",
#      "chemotaxis_arm"       = "#FF8C00",
#      "hedgehog_arm"         = "#1E90FF",
#      "lxr_arm"              = "#DC143C",
#      "cholesterol_sensing"  = "#8B6914",
#      "degradation"          = "#808080"
#
#    - Also add matching labels to `edge_pathway_labels`:
#
#      ifn_input           = "IFN input",
#      biosynthesis        = "Biosynthesis",
#      chemotaxis_arm      = "Immune chemotaxis",
#      hedgehog_arm        = "Hedgehog/Osteogenic",
#      lxr_arm             = "LXR/Resorptive",
#      cholesterol_sensing = "Cholesterol sensing",
#      degradation         = "Degradation"
#
# 3. NOTEBOOK UPDATE (pathway_mapper_v3_notebook.Rmd):
#    - In Section 4 ("Run a single pathway"), add a comment noting
#      "oxysterol" is now available
#    - In Section 5 ("Run a custom selection"), add "oxysterol" to 
#      the example MY_PATHWAYS vector
#    - In Section 6 ("Run ALL pathways"), no change needed — 
#      list_pathways() will pick it up automatically
#
# 4. VERIFICATION:
#    After insertion, source the file and run:
#      source("pathway_mapper_v3.R")
#      list_pathways()  # should now include "oxysterol"
#      cfg <- get_pathway_config("oxysterol")
#      cat("Nodes:", nrow(cfg$nodes), "\n")  # expect 25
#      cat("Edges:", nrow(cfg$edges), "\n")  # expect 25
#      cat("Synonyms:", length(cfg$synonyms), "\n")  # expect 20
#
#    Then test rendering:
#      run_pathway_maps(
#        input    = "DESeq2_Free_vs_Pygostyle_by_ENS.csv",
#        outdir   = "pathway_out",
#        title    = "Free vs Pygostyle 2023",
#        pathways = "oxysterol"
#      )
#
# 5. DO NOT MODIFY:
#    - Any existing pathway registry entries
#    - The core rendering functions (build_pathway_plot, 
#      render_pathway_map, etc.)
#    - The prepare_de_table_for_pathways function
#    - The CLI parsing code
#
# ============================================================
# BIOLOGICAL NOTES FOR REVIEW:
# ============================================================
#
# The pathway layout has three spatial tiers (top to bottom):
#   Top (y~4-5.5): IFN→CH25H→CYP7B1→GPR183 axis + INSIG/SCAP/SREBF2
#   Middle (y~1.8-2.5): CYP11A1→SMO→GLI→HES/HEY→RUNX2/SP7 (Hh/osteogenic)
#   Bottom (y~-0.5-1.0): CYP27A1→NR1H3→ABCA1/ABCG1/RANKL (LXR/resorptive)
#
# The edge from SULT2B1→CH25H represents negative regulation
# (sulfonation inactivates the 25-HC product). The rendering 
# engine currently doesn't distinguish activating vs inhibitory
# edges — this could be a future enhancement.
#
# The edge from HSD3B7→GPR183 similarly represents gradient
# sculpting (degradation of the ligand 7α,25-diHC), which 
# shapes the chemotactic field. It is not a direct activation.
#
# Shared nodes with other pathways:
#   - RUNX2, SP7: also in osteoblast pathway
#   - TNFSF11 (RANKL): also in osteoclast pathway
#   - IFNB1, IFNAR1: also in cGAS-STING pathway
# This is intentional — it shows cross-pathway convergence.
# ============================================================
