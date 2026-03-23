#!/usr/bin/env Rscript

# Pathway Mapper v3 — general-purpose pathway map generator
# Renders any signaling pathway from a DESeq2 results CSV
# (Gallus gallus, Ensembl IDs via org.Gg.eg.db)
#
# RStudio usage:
#   source("pathway_mapper_v3.R")
#   list_pathways()
#   run_pathway_maps(
#     input = "DESeq2_results.csv",
#     outdir = "pathway_out",
#     title  = "Sacral vs Free"
#   )
#
# CLI usage:
#   Rscript pathway_mapper_v3.R --input results.csv --outdir out --title "Sacral vs Free"
#   Rscript pathway_mapper_v3.R --list-pathways

suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Gg.eg.db)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(scales)
  library(svglite)
  library(grid)
  library(cowplot)
})

# ---------------------------------------------------------------------------
# Pathway Registry
# ---------------------------------------------------------------------------

pathway_registry <- list()

pathway_registry[["cgas_sting"]] <- function() {
  nodes <- tribble(
    ~node_id,   ~symbol_key, ~label,   ~x,  ~y,  ~module,           ~compartment,
    "DDX41",    "DDX41",     "DDX41",   0,   2.8, "Sensors",         "Cytosol",
    "IFI16",    "IFI16",     "IFI16",   0,   2.2, "Sensors",         "Cytosol",
    "MB21D1",   "MB21D1",    "cGAS",    0,   1.6, "Sensors",         "Cytosol",
    "ZBP1",     "ZBP1",      "ZBP1",    0,   0.8, "Sensors",         "Cytosol",
    "STING1",   "STING1",    "STING1",  1,   1.8, "STING",           "ER/Golgi",
    "TRAF3",    "TRAF3",     "TRAF3",   2,   2.7, "Trafficking",     "Cytosol",
    "TRAF6",    "TRAF6",     "TRAF6",   2,   1.1, "Trafficking",     "Cytosol",
    "TBK1",     "TBK1",      "TBK1",    2,   2.1, "Kinases",         "Cytosol",
    "IKBKB",    "IKBKB",     "IKKβ",    2,   0.4, "Kinases",         "Cytosol",
    "IKBKG",    "IKBKG",     "NEMO",    2,  -0.2, "Kinases",         "Cytosol",
    "CHUK",     "CHUK",      "IKKα",    2,  -0.8, "Kinases",         "Cytosol",
    "IRF7",     "IRF7",      "IRF7",    3,   2.9, "TFs",             "Nucleus",
    "IRF3",     "IRF3",      "IRF3",    3,   2.2, "TFs",             "Nucleus",
    "NFKB1",    "NFKB1",     "NF-κB",   3,   0.5, "TFs",             "Nucleus",
    "RELA",     "RELA",      "RelA",    3,  -0.2, "TFs",             "Nucleus",
    "IFNB1",    "IFNB1",     "IFNβ",    4,   2.3, "IFN",             "Secreted",
    "IFNAR1",   "IFNAR1",    "IFNAR1",  5,   2.8, "IFNAR/JAK-STAT",  "Plasma membrane",
    "IFNAR2",   "IFNAR2",    "IFNAR2",  5,   2.0, "IFNAR/JAK-STAT",  "Plasma membrane",
    "JAK1",     "JAK1",      "JAK1",    6,   2.8, "IFNAR/JAK-STAT",  "Cytosol",
    "TYK2",     "TYK2",      "TYK2",    6,   2.0, "IFNAR/JAK-STAT",  "Cytosol",
    "STAT1",    "STAT1",     "STAT1",   6,   3.3, "IFNAR/JAK-STAT",  "Nucleus",
    "STAT2",    "STAT2",     "STAT2",   6,   1.5, "IFNAR/JAK-STAT",  "Nucleus",
    "IRF9",     "IRF9",      "IRF9",    6,   2.3, "IFNAR/JAK-STAT",  "Nucleus",
    "CXCL10",   "CXCL10",    "CXCL10",  7,   3.4, "ISGs/Outputs",    "Secreted",
    "CCL5",     "CCL5",      "CCL5",    7,   2.8, "ISGs/Outputs",    "Secreted",
    "ISG15",    "ISG15",     "ISG15",   7,   2.4, "ISGs/Outputs",    "Cytosol",
    "MX1",      "MX1",       "MX1",     7,   2.0, "ISGs/Outputs",    "Cytosol",
    "OASL",     "OASL",      "OASL",    7,   1.6, "ISGs/Outputs",    "Cytosol",
    "IFIT5",    "IFIT5",     "IFIT5",   7,   1.2, "ISGs/Outputs",    "Cytosol",
    "RSAD2",    "RSAD2",     "RSAD2",   7,   0.8, "ISGs/Outputs",    "Cytosol",
    "IL6",      "IL6",       "IL6",     7,   0.1, "ISGs/Outputs",    "Secreted",
    "TNF",      "TNF",       "TNF",     7,  -0.5, "ISGs/Outputs",    "Secreted"
  )

  edges <- tribble(
    ~from,     ~to,       ~edge_class, ~edge_pathway,
    "DDX41",   "STING1",  "core",      "cgas_sting",
    "IFI16",   "STING1",  "core",      "cgas_sting",
    "MB21D1",  "STING1",  "core",      "cgas_sting",
    "ZBP1",    "TBK1",    "core",      "cgas_sting",
    "STING1",  "TRAF3",   "core",      "cgas_sting",
    "STING1",  "TRAF6",   "core",      "cgas_sting",
    "TRAF3",   "TBK1",    "core",      "cgas_sting",
    "TRAF6",   "IKBKB",   "core",      "cgas_sting",
    "IKBKG",   "IKBKB",   "core",      "cgas_sting",
    "STING1",  "TBK1",    "core",      "cgas_sting",
    "STING1",  "IKBKB",   "core",      "cgas_sting",
    "CHUK",    "NFKB1",   "core",      "cgas_sting",
    "IKBKB",   "NFKB1",   "core",      "cgas_sting",
    "IKBKB",   "RELA",    "core",      "cgas_sting",
    "TBK1",    "IRF3",    "core",      "cgas_sting",
    "TBK1",    "IRF7",    "core",      "cgas_sting",
    "IRF3",    "IFNB1",   "core",      "cgas_sting",
    "IRF7",    "IFNB1",   "core",      "cgas_sting",
    "IFNB1",   "IFNAR1",  "core",      "cgas_sting",
    "IFNB1",   "IFNAR2",  "core",      "cgas_sting",
    "IFNAR1",  "JAK1",    "core",      "cgas_sting",
    "IFNAR2",  "TYK2",    "core",      "cgas_sting",
    "JAK1",    "STAT1",   "core",      "cgas_sting",
    "TYK2",    "STAT2",   "core",      "cgas_sting",
    "STAT1",   "IRF9",    "core",      "cgas_sting",
    "STAT2",   "IRF9",    "core",      "cgas_sting",
    "IRF9",    "ISG15",   "output",    "cgas_sting",
    "IRF9",    "MX1",     "output",    "cgas_sting",
    "IRF9",    "OASL",    "output",    "cgas_sting",
    "IRF9",    "IFIT5",   "output",    "cgas_sting",
    "IRF9",    "RSAD2",   "output",    "cgas_sting",
    "RELA",    "IL6",     "output",    "cgas_sting",
    "RELA",    "TNF",     "output",    "cgas_sting",
    "IRF3",    "CXCL10",  "output",    "cgas_sting",
    "IRF7",    "CCL5",    "output",    "cgas_sting"
  )

  synonyms <- list(
    STING1 = c("TMEM173"),
    MB21D1 = c("CGAS"),
    ZBP1   = c("DAI")
  )

  validation_notes <- c(
    "STING axis: STING1 ER-to-Golgi trafficking by confocal.",
    "IFN output: IFNB1 / ISGs by qPCR or RNAscope.",
    "NF-kB arm: p65 nuclear translocation by IF.",
    "Chemokines: CXCL10/CCL5 by ELISA on conditioned media."
  )

  list(
    name = "cgas_sting",
    display_name = "cGAS-STING / Type I IFN Signaling",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

pathway_registry[["apoptosis"]] <- function() {
  nodes <- tribble(
    ~node_id,    ~symbol_key, ~label,    ~x,   ~y,   ~module,              ~compartment,
    "TNF",       "TNF",       "TNF",      0.0,  0.2,  "Ligand/Receptor",    NA_character_,
    "TNFRSF1A",  "TNFRSF1A",  "TNFRSF1A", 1.0,  0.2,  "Ligand/Receptor",    "Plasma membrane",
    "TRADD",     "TRADD",     "TRADD",    2.0,  0.2,  "Adaptor",            "Cytosol",
    "FADD",      "FADD",      "FADD",     3.0,  0.8,  "Caspase axis",       "Cytosol",
    "CASP8",     "CASP8",     "CASP8",    4.0,  0.8,  "Caspase axis",       "Cytosol",
    "CASP7",     "CASP7",     "CASP7",    5.0,  1.2,  "Executioners",       "Cytosol",
    "CASP3",     "CASP3",     "CASP3",    6.2,  0.8,  "Executioners",       "Cytosol",
    "BID",       "BID",       "BID",      4.8,  0.0,  "Mitochondrial axis", "Mitochondria",
    "BAX",       "BAX",       "BAX",      5.8, -0.3,  "Mitochondrial axis", "Mitochondria",
    "CYCS",      "CYCS",      "CYCS",     6.8, -0.3,  "Mitochondrial axis", "Mitochondria",
    "APAF1",     "APAF1",     "APAF1",    7.8, -0.3,  "Mitochondrial axis", "Cytosol",
    "CASP9",     "CASP9",     "CASP9",    8.8, -0.3,  "Mitochondrial axis", "Cytosol"
  )

  edges <- tribble(
    ~from,      ~to,        ~edge_class, ~edge_pathway,
    "TNF",      "TNFRSF1A", "core",      "shared_death_input",
    "TNFRSF1A", "TRADD",    "core",      "shared_death_input",
    "TRADD",    "FADD",     "core",      "apoptosis_branch",
    "FADD",     "CASP8",    "core",      "apoptosis_branch",
    "CASP8",    "CASP3",    "output",    "apoptosis_branch",
    "CASP8",    "CASP7",    "output",    "apoptosis_branch",
    "CASP8",    "BID",      "core",      "apoptosis_branch",
    "BID",      "BAX",      "core",      "apoptosis_branch",
    "BAX",      "CYCS",     "core",      "apoptosis_branch",
    "CYCS",     "APAF1",    "core",      "apoptosis_branch",
    "APAF1",    "CASP9",    "core",      "apoptosis_branch",
    "CASP9",    "CASP3",    "output",    "apoptosis_branch"
  )

  synonyms <- list(
    TNFRSF1A = c("TNFR1", "TNFRSF1"),
    CASP3    = c("CASP3A", "CASP3-like"),
    BAX      = c("BAX1")
  )

  validation_notes <- c(
    "Cleaved CASP3 / CASP8 by Western blot or IHC.",
    "TUNEL staining for DNA fragmentation.",
    "Cytochrome c release by fractionation (mitochondrial arm).",
    "Annexin V / PI flow cytometry for apoptosis vs necrosis."
  )

  list(
    name = "apoptosis",
    display_name = "Apoptosis",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

pathway_registry[["necroptosis"]] <- function() {
  nodes <- tribble(
    ~node_id,   ~symbol_key, ~label,    ~x,   ~y,   ~module,             ~compartment,
    "TNF",      "TNF",       "TNF",      0.0,  0.2,  "Ligand/Receptor",   NA_character_,
    "TNFRSF1A", "TNFRSF1A",  "TNFRSF1A", 1.0,  0.2,  "Ligand/Receptor",   "Plasma membrane",
    "TRADD",    "TRADD",     "TRADD",    2.0,  0.2,  "Adaptor",           "Cytosol",
    "RIPK1",    "RIPK1",     "RIPK1",    3.2,  0.2,  "Core necroptosis",  "Cytosol",
    "RIPK3",    "RIPK3",     "RIPK3*",   4.2,  0.2,  "Core necroptosis",  "Cytosol",
    "MLKL",     "MLKL",      "MLKL",     5.2,  0.2,  "Executioner",       "Plasma membrane",
    "ZBP1",     "ZBP1",      "ZBP1",     3.2,  1.0,  "Inputs",            "Cytosol",
    "CASP8",    "CASP8",     "CASP8",    3.2, -0.6,  "Cross-talk",        "Cytosol"
  )

  edges <- tribble(
    ~from,      ~to,        ~edge_class, ~edge_pathway,
    "TNF",      "TNFRSF1A", "core",      "shared_death_input",
    "TNFRSF1A", "TRADD",    "core",      "shared_death_input",
    "TRADD",    "RIPK1",    "core",      "necroptosis_branch",
    "RIPK1",    "RIPK3",    "core",      "necroptosis_branch",
    "RIPK3",    "MLKL",     "output",    "necroptosis_branch",
    "ZBP1",     "RIPK3",    "core",      "necroptosis_branch",
    "CASP8",    "RIPK1",    "core",      "necroptosis_branch"
  )

  synonyms <- list(
    TNFRSF1A = c("TNFR1", "TNFRSF1"),
    RIPK1    = c("RIP1"),
    RIPK3    = c("RIP3"),
    MLKL     = c("MLKL1"),
    ZBP1     = c("DAI")
  )

  validation_notes <- c(
    "*RIPK3: GRCg7b locus (GeneID 415708) is RIPK2-like by BLASTp. True RIPK3 is likely absent across Aves.",
    "pMLKL by Western blot (note: chicken MLKL may be orphaned pseudokinase).",
    "RIPK1/RIPK3 necrosome complex by co-IP.",
    "MLKL membrane translocation by fractionation or confocal.",
    "Distinguish from apoptosis: CASP8 inhibition should potentiate, not block."
  )

  list(
    name = "necroptosis",
    display_name = "Necroptosis",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

pathway_registry[["inflammasome"]] <- function() {
  nodes <- tribble(
    ~node_id,   ~symbol_key, ~label,      ~x,   ~y,   ~module,             ~compartment,
    # --- Priming arm (TLR → NF-κB) ---
    "TLR4",     "TLR4",      "TLR4",       0.0,  1.0,  "Priming",           "Plasma membrane",
    "MYD88",    "MYD88",     "MYD88",      1.0,  1.0,  "Priming",           "Cytosol",
    "IRAK4",    "IRAK4",     "IRAK4",      2.0,  1.0,  "Priming",           "Cytosol",
    "TRAF6",    "TRAF6",     "TRAF6",      3.0,  1.0,  "Priming",           "Cytosol",
    "NFKB1",    "NFKB1",     "NF-κB",      4.0,  1.0,  "Priming",           "Nucleus",
    "RELA",     "RELA",      "RelA",       4.0,  0.4,  "Priming",           "Nucleus",
    # --- Priming transcriptional targets ---
    "NLRP3_tx", "NLRP3",     "NLRP3\u2191", 5.2,  1.4,  "Priming targets",   "Cytosol",
    "IL1B_tx",  "IL1B",      "pro-IL-1\u03b2", 5.2, 0.6,  "Priming targets",   "Cytosol",
    # --- Sensor assembly ---
    "NLRP3",    "NLRP3",     "NLRP3",      3.0,  3.0,  "Sensor assembly",   "Cytosol",
    "NEK7",     "NEK7",      "NEK7",       2.0,  3.0,  "Sensor assembly",   "Cytosol",
    "TXNIP",    "TXNIP",     "TXNIP",      2.0,  3.6,  "Sensor assembly",   "Cytosol",
    "PYCARD",   "PYCARD",    "ASC",        4.0,  3.0,  "Sensor assembly",   "Cytosol",
    "AIM2",     "AIM2",      "AIM2",       3.0,  3.8,  "Sensor assembly",   "Cytosol",
    "NLRC4",    "NLRC4",     "NLRC4",      3.0,  2.4,  "Sensor assembly",   "Cytosol",
    # --- Effector caspase ---
    "CASP1",    "CASP1",     "CASP1",      5.0,  3.0,  "Effectors",         "Cytosol",
    # --- Outputs ---
    "IL1B",     "IL1B",      "IL-1\u03b2",  6.5,  3.4,  "Outputs",           "Secreted",
    "IL18",     "IL18",      "IL-18",      6.5,  2.6,  "Outputs",           "Secreted",
    "GSDMD",    "GSDMD",     "GSDMD",      6.5,  2.0,  "Outputs",           "Plasma membrane"
  )

  edges <- tribble(
    ~from,       ~to,        ~edge_class, ~edge_pathway,
    # Priming arm
    "TLR4",      "MYD88",    "core",      "priming",
    "MYD88",     "IRAK4",    "core",      "priming",
    "IRAK4",     "TRAF6",    "core",      "priming",
    "TRAF6",     "NFKB1",    "core",      "priming",
    "TRAF6",     "RELA",     "core",      "priming",
    "NFKB1",     "NLRP3_tx", "output",    "priming",
    "RELA",      "IL1B_tx",  "output",    "priming",
    # Sensor assembly
    "NEK7",      "NLRP3",    "core",      "activation",
    "TXNIP",     "NLRP3",    "core",      "activation",
    "NLRP3",     "PYCARD",   "core",      "activation",
    "AIM2",      "PYCARD",   "core",      "activation",
    "NLRC4",     "PYCARD",   "core",      "activation",
    "PYCARD",    "CASP1",    "core",      "activation",
    # Effector outputs
    "CASP1",     "IL1B",     "output",    "effector",
    "CASP1",     "IL18",     "output",    "effector",
    "CASP1",     "GSDMD",    "output",    "effector"
  )

  synonyms <- list(
    PYCARD = c("ASC", "TMS1"),
    NLRP3  = c("NALP3", "CIAS1"),
    CASP1  = c("ICE"),
    IL1B   = c("IL1BETA"),
    MYD88  = c("MYD88"),
    NLRC4  = c("IPAF"),
    GSDMD  = c("GSDMD")
  )

  validation_notes <- c(
    "CASP1 cleavage (p20 fragment) by Western blot.",
    "ASC speck formation by immunofluorescence.",
    "IL-1\u03b2 / IL-18 by ELISA on conditioned media.",
    "GSDMD N-terminal fragment by WB (note: chicken GSDMD orthology uncertain).",
    "LDH release assay for pyroptotic cell death.",
    "NLRC4 likely absent in chicken — informative if 'Missing gene'."
  )

  list(
    name = "inflammasome",
    display_name = "Inflammasome (NLRP3 Canonical)",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

pathway_registry[["osteoclast"]] <- function() {
  nodes <- tribble(
    ~node_id,      ~symbol_key,  ~label,       ~x,   ~y,   ~module,               ~compartment,
    # --- Ligand / Receptor ---
    "CSF1",        "CSF1",       "M-CSF",       0.0,  0.0,  "Ligand/Receptor",     "Secreted",
    "CSF1R",       "CSF1R",      "CSF1R",       1.0,  0.0,  "Ligand/Receptor",     "Plasma membrane",
    "TNFSF11",     "TNFSF11",    "RANKL",       0.0,  2.0,  "Ligand/Receptor",     "Secreted",
    "TNFRSF11A",   "TNFRSF11A",  "RANK",        1.0,  2.0,  "Ligand/Receptor",     "Plasma membrane",
    "TNFRSF11B",   "TNFRSF11B",  "OPG",         0.0,  3.0,  "Ligand/Receptor",     "Secreted",
    # --- Signaling cascade ---
    "TRAF6",       "TRAF6",      "TRAF6",       2.2,  2.0,  "Signaling cascade",   "Cytosol",
    "MAP3K7",      "MAP3K7",     "TAK1",        3.2,  2.0,  "Signaling cascade",   "Cytosol",
    "MAPK14",      "MAPK14",     "p38",         4.2,  2.5,  "Signaling cascade",   "Cytosol",
    "MAPK1",       "MAPK1",      "ERK",         4.2,  1.5,  "Signaling cascade",   "Cytosol",
    "NFKB1",       "NFKB1",      "NF-\u03baB",  4.2,  0.5,  "Signaling cascade",   "Nucleus",
    # --- Transcription factors ---
    "JUN",         "JUN",        "c-Jun",       5.5,  2.5,  "Transcription factors", "Nucleus",
    "FOS",         "FOS",        "c-Fos",       5.5,  1.5,  "Transcription factors", "Nucleus",
    "NFATC1",      "NFATC1",     "NFATc1",      6.5,  2.0,  "Transcription factors", "Nucleus",
    # --- Osteoclast markers ---
    "ACP5",        "ACP5",       "TRAP",        8.0,  3.0,  "OC markers",          "Cytosol",
    "CTSK",        "CTSK",       "CathK",       8.0,  2.4,  "OC markers",          "Cytosol",
    "MMP9",        "MMP9",       "MMP9",        8.0,  1.8,  "OC markers",          "Secreted",
    "CALCR",       "CALCR",      "CTR",         8.0,  1.2,  "OC markers",          "Plasma membrane",
    "DCSTAMP",     "DCSTAMP",    "DC-STAMP",    8.0,  0.6,  "OC markers",          "Plasma membrane"
  )

  edges <- tribble(
    ~from,          ~to,          ~edge_class, ~edge_pathway,
    # M-CSF survival axis
    "CSF1",         "CSF1R",      "core",      "csf1_survival",
    # RANKL-RANK core
    "TNFSF11",      "TNFRSF11A",  "core",      "rankl_rank",
    "TNFRSF11B",    "TNFSF11",    "core",      "rankl_rank",
    "TNFRSF11A",    "TRAF6",      "core",      "rankl_rank",
    "TRAF6",        "MAP3K7",     "core",      "rankl_rank",
    "MAP3K7",       "MAPK14",     "core",      "rankl_rank",
    "MAP3K7",       "MAPK1",      "core",      "rankl_rank",
    "MAP3K7",       "NFKB1",      "core",      "rankl_rank",
    # MAPK → AP-1 → NFATc1
    "MAPK14",       "JUN",        "core",      "rankl_rank",
    "MAPK1",        "FOS",        "core",      "rankl_rank",
    "JUN",          "NFATC1",     "core",      "rankl_rank",
    "FOS",          "NFATC1",     "core",      "rankl_rank",
    "NFKB1",        "NFATC1",     "core",      "rankl_rank",
    # NFATc1 → markers
    "NFATC1",       "ACP5",       "output",    "oc_output",
    "NFATC1",       "CTSK",       "output",    "oc_output",
    "NFATC1",       "MMP9",       "output",    "oc_output",
    "NFATC1",       "CALCR",      "output",    "oc_output",
    "NFATC1",       "DCSTAMP",    "output",    "oc_output"
  )

  synonyms <- list(
    TNFSF11    = c("RANKL", "TRANCE", "OPGL"),
    TNFRSF11A  = c("RANK"),
    TNFRSF11B  = c("OPG"),
    MAP3K7     = c("TAK1"),
    ACP5       = c("TRAP"),
    CTSK       = c("CATK"),
    DCSTAMP    = c("TM7SF4"),
    NFATC1     = c("NFAT2"),
    CSF1       = c("MCSF")
  )

  validation_notes <- c(
    "TRAP staining (ACP5) on tissue sections or cytospins.",
    "Cathepsin K (CTSK) activity assay or IHC.",
    "Pit/dentine resorption assay for functional osteoclast activity.",
    "NFATc1 nuclear localization by IF.",
    "OPG (TNFRSF11B) shown as inhibitor: OPG \u2192 RANKL edge is inhibitory (decoy receptor).",
    "DC-STAMP and CALCR may have limited chicken annotation."
  )

  list(
    name = "osteoclast",
    display_name = "Osteoclast Differentiation (RANKL-RANK)",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

pathway_registry[["osteoblast"]] <- function() {
  nodes <- tribble(
    ~node_id,   ~symbol_key, ~label,        ~x,   ~y,   ~module,                 ~compartment,
    # --- BMP signaling ---
    "BMP2",     "BMP2",      "BMP2",         0.0,  3.5,  "BMP signaling",         "Secreted",
    "BMP4",     "BMP4",      "BMP4",         0.0,  2.8,  "BMP signaling",         "Secreted",
    "BMP7",     "BMP7",      "BMP7",         0.0,  2.1,  "BMP signaling",         "Secreted",
    "BMPR1A",   "BMPR1A",    "BMPR1A",       1.2,  3.2,  "BMP signaling",         "Plasma membrane",
    "BMPR2",    "BMPR2",     "BMPR2",        1.2,  2.4,  "BMP signaling",         "Plasma membrane",
    "SMAD1",    "SMAD1",     "SMAD1",        2.4,  3.2,  "BMP signaling",         "Cytosol",
    "SMAD5",    "SMAD5",     "SMAD5",        2.4,  2.4,  "BMP signaling",         "Cytosol",
    "SMAD4",    "SMAD4",     "SMAD4",        3.4,  2.8,  "BMP signaling",         "Nucleus",
    # --- Wnt / β-catenin ---
    "WNT3A",    "WNT3A",     "Wnt3a",        0.0,  0.0,  "Wnt signaling",         "Secreted",
    "FZD1",     "FZD1",      "Fzd1",         1.2,  0.6,  "Wnt signaling",         "Plasma membrane",
    "LRP5",     "LRP5",      "LRP5",         1.2,  0.0,  "Wnt signaling",         "Plasma membrane",
    "LRP6",     "LRP6",      "LRP6",         1.2, -0.6,  "Wnt signaling",         "Plasma membrane",
    "DVL2",     "DVL2",      "DVL2",         2.4,  0.3,  "Wnt signaling",         "Cytosol",
    "GSK3B",    "GSK3B",     "GSK3\u03b2",   2.4, -0.5,  "Wnt signaling",         "Cytosol",
    "APC",      "APC",       "APC",          3.4, -0.5,  "Wnt signaling",         "Cytosol",
    "AXIN2",    "AXIN2",     "Axin2",        3.4, -1.1,  "Wnt signaling",         "Cytosol",
    "CTNNB1",   "CTNNB1",    "\u03b2-catenin", 3.4, 0.3,  "Wnt signaling",         "Nucleus",
    # --- Transcription factors ---
    "RUNX2",    "RUNX2",     "RUNX2",        5.0,  2.0,  "Transcription factors", "Nucleus",
    "SP7",      "SP7",       "Osterix",      5.0,  1.2,  "Transcription factors", "Nucleus",
    "ATF4",     "ATF4",      "ATF4",         5.0,  0.4,  "Transcription factors", "Nucleus",
    # --- Osteoblast markers ---
    "COL1A1",   "COL1A1",    "COL1A1",       7.0,  2.8,  "OB markers",            "ECM",
    "ALPL",     "ALPL",      "ALP",          7.0,  2.2,  "OB markers",            "Plasma membrane",
    "BGLAP",    "BGLAP",     "OCN",          7.0,  1.6,  "OB markers",            "Secreted",
    "SPP1",     "SPP1",      "OPN",          7.0,  1.0,  "OB markers",            "Secreted",
    "IBSP",     "IBSP",      "BSP",          7.0,  0.4,  "OB markers",            "ECM",
    # --- Negative regulators ---
    "SOST",     "SOST",      "Sclerostin",   0.0, -1.4,  "Negative regulators",   "Secreted",
    "DKK1",     "DKK1",      "DKK1",         0.0, -2.0,  "Negative regulators",   "Secreted"
  )

  edges <- tribble(
    ~from,      ~to,       ~edge_class, ~edge_pathway,
    # BMP axis
    "BMP2",     "BMPR1A",  "core",      "bmp_signaling",
    "BMP4",     "BMPR1A",  "core",      "bmp_signaling",
    "BMP7",     "BMPR2",   "core",      "bmp_signaling",
    "BMPR1A",   "SMAD1",   "core",      "bmp_signaling",
    "BMPR2",    "SMAD5",   "core",      "bmp_signaling",
    "SMAD1",    "SMAD4",   "core",      "bmp_signaling",
    "SMAD5",    "SMAD4",   "core",      "bmp_signaling",
    "SMAD4",    "RUNX2",   "core",      "bmp_signaling",
    # Wnt axis
    "WNT3A",    "FZD1",    "core",      "wnt_signaling",
    "WNT3A",    "LRP5",    "core",      "wnt_signaling",
    "WNT3A",    "LRP6",    "core",      "wnt_signaling",
    "FZD1",     "DVL2",    "core",      "wnt_signaling",
    "DVL2",     "GSK3B",   "core",      "wnt_signaling",
    "GSK3B",    "APC",     "core",      "wnt_signaling",
    "GSK3B",    "AXIN2",   "core",      "wnt_signaling",
    "DVL2",     "CTNNB1",  "core",      "wnt_signaling",
    "CTNNB1",   "RUNX2",   "core",      "wnt_signaling",
    # Negative regulators → Wnt
    "SOST",     "LRP5",    "core",      "wnt_inhibition",
    "DKK1",     "LRP6",    "core",      "wnt_inhibition",
    # TF cascade
    "RUNX2",    "SP7",     "core",      "ob_differentiation",
    "SP7",      "ATF4",    "core",      "ob_differentiation",
    # TF → markers
    "RUNX2",    "COL1A1",  "output",    "ob_output",
    "RUNX2",    "ALPL",    "output",    "ob_output",
    "SP7",      "BGLAP",   "output",    "ob_output",
    "SP7",      "SPP1",    "output",    "ob_output",
    "ATF4",     "BGLAP",   "output",    "ob_output",
    "ATF4",     "IBSP",    "output",    "ob_output"
  )

  synonyms <- list(
    SP7    = c("OSX", "Osterix"),
    BGLAP  = c("OCN", "osteocalcin", "BGLAP2"),
    ALPL   = c("ALP", "TNAP"),
    SPP1   = c("OPN", "osteopontin"),
    IBSP   = c("BSP", "BSP2"),
    CTNNB1 = c("BCAT", "beta-catenin"),
    SOST   = c("sclerostin"),
    DKK1   = c("DKK-1"),
    GSK3B  = c("GSK3BETA"),
    SMAD4  = c("MADH4"),
    RUNX2  = c("CBFA1", "OSF2"),
    WNT3A  = c("WNT3")
  )

  validation_notes <- c(
    "Alizarin red S staining for matrix mineralization.",
    "ALP (ALPL) activity assay — early OB marker.",
    "Osteocalcin (BGLAP) ELISA on conditioned media — late OB marker.",
    "von Kossa staining for calcium phosphate deposits.",
    "RUNX2 and Osterix nuclear localization by IF.",
    "SOST/DKK1 shown as Wnt pathway inhibitors (edges to LRP5/6).",
    "SOST and DKK1 may have limited chicken annotation — informative if 'Missing gene'."
  )

  list(
    name = "osteoblast",
    display_name = "Osteoblast Differentiation (BMP/Wnt)",
    nodes = nodes,
    edges = edges,
    synonyms = synonyms,
    validation_notes = validation_notes
  )
}

list_pathways <- function() sort(names(pathway_registry))

get_pathway_config <- function(pathway) {
  pathway <- tolower(pathway)
  if (!pathway %in% names(pathway_registry)) {
    stop(
      "Unknown pathway: ", pathway,
      "\nAvailable: ", paste(list_pathways(), collapse = ", "),
      call. = FALSE
    )
  }
  pathway_registry[[pathway]]()
}

# ---------------------------------------------------------------------------
# Core helper functions (battle-tested — do not modify internals)
# ---------------------------------------------------------------------------

prepare_de_table_for_pathways <- function(input,
                                          lfc_column_preference = c("log2FC_shrunken", "log2FC"),
                                          id_column_candidates = c("ensembl_id", "Unnamed: 0", "gene", "id")) {
  if (!file.exists(input)) {
    stop("Input file not found: ", input, call. = FALSE)
  }

  de_raw <- readr::read_csv(
    file = input,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )

  if (any(is.na(names(de_raw)))) {
    idx_na <- which(is.na(names(de_raw)))
    names(de_raw)[idx_na] <- paste0("unnamed_col_", seq_along(idx_na))
  }

  id_candidates <- unique(c(id_column_candidates, "ensembl_id", "Unnamed: 0", "...1"))
  id_col <- id_candidates[id_candidates %in% names(de_raw)][1]

  if (!is.na(id_col) && nzchar(id_col) && id_col != "ensembl_id") {
    names(de_raw)[names(de_raw) == id_col] <- "ensembl_id"
  }

  if (!"ensembl_id" %in% names(de_raw) && "" %in% names(de_raw)) {
    names(de_raw)[which(names(de_raw) == "")[1]] <- "ensembl_id"
  }

  if (!"ensembl_id" %in% names(de_raw) && ncol(de_raw) >= 1L) {
    first_col_name <- names(de_raw)[1]
    # Accept the first column as the ID column if it is unnamed/empty or looks
    # like an ID column (Ensembl IDs or gene symbols).
    if (!nzchar(first_col_name) || first_col_name %in% c("...1", "Unnamed: 0")) {
      names(de_raw)[1] <- "ensembl_id"
    } else {
      first_values <- de_raw[[1]]
      first_non_na <- as.character(first_values[which(!is.na(first_values))[1]])
      if (!is.na(first_non_na) && nzchar(first_non_na)) {
        names(de_raw)[1] <- "ensembl_id"
      }
    }
  }

  if (!"ensembl_id" %in% names(de_raw)) {
    stop(
      "Could not find ID column (Ensembl IDs or gene symbols). Checked: ",
      paste(id_column_candidates, collapse = ", "),
      call. = FALSE
    )
  }

  available_fc <- lfc_column_preference[lfc_column_preference %in% names(de_raw)]
  if (length(available_fc) >= 1L) {
    primary_fc <- available_fc[1]
    secondary_fc <- if (length(available_fc) >= 2L) available_fc[2] else NA_character_
    if (!is.na(secondary_fc)) {
      fc_column_used <- paste0(primary_fc, " (fallback to ", secondary_fc, " when NA)")
      de_raw <- de_raw %>%
        mutate(log2FC_used = ifelse(!is.na(.data[[primary_fc]]), .data[[primary_fc]], .data[[secondary_fc]]))
    } else {
      fc_column_used <- primary_fc
      de_raw <- de_raw %>% mutate(log2FC_used = .data[[primary_fc]])
    }
  } else {
    fc_column_used <- "none (all NA)"
    warning("None of the preferred log2FC columns were found: ", paste(lfc_column_preference, collapse = ", "))
    de_raw <- de_raw %>% mutate(log2FC_used = NA_real_)
  }

  if (!"padj" %in% names(de_raw)) {
    warning("padj column not found. All genes will be treated as not significant.")
    de_raw <- de_raw %>% mutate(padj = NA_real_)
  }

  de_raw <- de_raw %>%
    mutate(
      ensembl_id = as.character(ensembl_id),
      log2FC_used = suppressWarnings(as.numeric(log2FC_used)),
      padj = suppressWarnings(as.numeric(padj))
    )

  id_values   <- unique(stats::na.omit(de_raw$ensembl_id))
  is_ensembl  <- any(grepl("^ENS[A-Z]*\\d{8,}", id_values))

  if (is_ensembl) {
    annotation <- AnnotationDbi::select(
      x       = org.Gg.eg.db,
      keys    = id_values,
      columns = c("SYMBOL", "ENTREZID"),
      keytype = "ENSEMBL"
    )

    annotation_first <- annotation %>%
      as_tibble() %>%
      dplyr::rename(ensembl_id = ENSEMBL, symbol = SYMBOL, entrez_id = ENTREZID) %>%
      group_by(ensembl_id) %>%
      dplyr::slice(1) %>%
      ungroup()

    de_annot <- de_raw %>%
      left_join(annotation_first, by = "ensembl_id") %>%
      mutate(symbol = as.character(symbol))

  } else {
    # ID column contains gene symbols — rename and look up Ensembl/Entrez from symbol
    de_raw <- de_raw %>% dplyr::rename(symbol = ensembl_id)

    symbol_keys <- unique(stats::na.omit(de_raw$symbol))
    annotation <- AnnotationDbi::select(
      x       = org.Gg.eg.db,
      keys    = symbol_keys,
      columns = c("ENSEMBL", "ENTREZID"),
      keytype = "SYMBOL"
    )

    annotation_first <- annotation %>%
      as_tibble() %>%
      dplyr::rename(symbol = SYMBOL, ensembl_id = ENSEMBL, entrez_id = ENTREZID) %>%
      group_by(symbol) %>%
      dplyr::slice(1) %>%
      ungroup()

    de_annot <- de_raw %>%
      left_join(annotation_first, by = "symbol") %>%
      mutate(symbol = as.character(symbol))
  }

  de_by_symbol <- de_annot %>%
    filter(!is.na(symbol), symbol != "") %>%
    arrange(is.na(padj), padj, desc(abs(replace_na(log2FC_used, 0)))) %>%
    group_by(symbol) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(data_symbol = symbol, ensembl_id, log2FC_used, padj)

  list(
    de_raw = de_raw,
    de_annot = de_annot,
    de_by_symbol = de_by_symbol,
    fc_column_used = fc_column_used,
    input_basename = basename(input)
  )
}

resolve_symbol_match <- function(symbol_key, synonyms, available_symbols) {
  candidates <- unique(c(symbol_key, synonyms[[symbol_key]]))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0L || length(available_symbols) == 0L) {
    return(NA_character_)
  }

  available_symbols <- as.character(available_symbols)
  match_idx <- match(toupper(candidates), toupper(available_symbols), nomatch = 0L)
  match_idx <- match_idx[match_idx > 0L]
  if (length(match_idx) == 0L) {
    return(NA_character_)
  }
  available_symbols[match_idx[1]]
}

build_pathway_node_table <- function(config, de_by_symbol, p_cutoff = 0.05) {
  required_cols <- c("node_id", "symbol_key", "label", "x", "y")
  if (!all(required_cols %in% names(config$nodes))) {
    stop("Config nodes must include: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }

  nodes <- config$nodes
  if (!"module" %in% names(nodes)) {
    nodes$module <- "Pathway"
  }
  if (!"compartment" %in% names(nodes)) {
    nodes$compartment <- NA_character_
  }

  available_symbols <- de_by_symbol$data_symbol
  synonyms <- config$synonyms
  if (is.null(synonyms)) {
    synonyms <- list()
  }

  node_status_levels <- c(
    "Missing gene",
    "Present; padj NA",
    "Present; not significant",
    "Present; significant"
  )

  nodes <- nodes %>%
    rowwise() %>%
    mutate(map_symbol = resolve_symbol_match(symbol_key, synonyms, available_symbols)) %>%
    ungroup()

  nodes %>%
    left_join(de_by_symbol, by = c("map_symbol" = "data_symbol")) %>%
    mutate(
      present_in_data = !is.na(map_symbol),
      signif_flag = !is.na(padj) & padj < p_cutoff,
      node_status = case_when(
        !present_in_data ~ "Missing gene",
        is.na(padj) ~ "Present; padj NA",
        signif_flag ~ "Present; significant",
        TRUE ~ "Present; not significant"
      ),
      node_status = factor(node_status, levels = node_status_levels),
      log2FC_plot = ifelse(node_status == "Missing gene", NA_real_, log2FC_used),
      compartment = as.character(compartment),
      gene_id = ensembl_id
    ) %>%
    distinct(node_id, .keep_all = TRUE)
}

choose_compartment_label <- function(box, nodes_in_box, pad = 0.45) {
  # Always place label at top-left corner of the box, inside by a small margin.
  # This avoids the label landing in the middle of dense node clusters.
  tibble::tibble(
    label_x = box$xmin + pad * 0.6,
    label_y = box$ymax - pad * 0.3,
    hjust   = 0,
    vjust   = 1,
    pos_id  = "top_left"
  )
}

infer_compartment_boxes <- function(nodes_df,
                                    pad_x = 0.6,
                                    pad_y = 0.6,
                                    min_nodes_per_box = 2,
                                    force_singleton_compartments = NULL) {
  if (!"compartment" %in% names(nodes_df)) {
    return(tibble())
  }

  nodes_with_comp <- nodes_df %>%
    filter(!is.na(compartment), compartment != "")
  if (nrow(nodes_with_comp) == 0) {
    return(tibble())
  }

  if (is.null(force_singleton_compartments)) {
    force_singleton_compartments <- character(0)
  }

  keep_compartments <- nodes_with_comp %>%
    count(compartment, name = "n_nodes") %>%
    filter(n_nodes >= min_nodes_per_box | compartment %in% force_singleton_compartments) %>%
    pull(compartment)

  if (length(keep_compartments) == 0) {
    return(tibble())
  }

  boxes <- nodes_with_comp %>%
    filter(compartment %in% keep_compartments) %>%
    group_by(compartment) %>%
    summarise(
      xmin = min(x, na.rm = TRUE) - pad_x,
      xmax = max(x, na.rm = TRUE) + pad_x,
      ymin = min(y, na.rm = TRUE) - pad_y,
      ymax = max(y, na.rm = TRUE) + pad_y,
      .groups = "drop"
    )

  label_tbl <- purrr::map_dfr(seq_len(nrow(boxes)), function(i) {
    box <- boxes[i, , drop = FALSE]
    comp <- box$compartment[[1]]
    nodes_in_box <- nodes_with_comp %>%
      filter(compartment == comp) %>%
      dplyr::select(x, y)
    choose_compartment_label(box = box, nodes_in_box = nodes_in_box, pad = 0.45)
  }) %>%
    dplyr::select(label_x, label_y, hjust, vjust, pos_id)

  bind_cols(boxes, label_tbl)
}

compute_manual_compartment_boxes <- function(nodes_df, compartment_spec) {
  if (is.null(compartment_spec) || nrow(compartment_spec) == 0) {
    return(tibble())
  }

  if (!"compartment" %in% names(compartment_spec)) {
    stop("Manual compartment spec must include a `compartment` column.", call. = FALSE)
  }

  spec <- as_tibble(compartment_spec)
  has_explicit_bounds <- all(c("xmin", "xmax", "ymin", "ymax") %in% names(spec))
  has_node_pad_bounds <- all(c("node_ids", "pad_x", "pad_y") %in% names(spec))
  if (!has_explicit_bounds && !has_node_pad_bounds) {
    stop(
      "Manual compartment spec must include either explicit bounds (xmin/xmax/ymin/ymax) ",
      "or node_ids + pad_x/pad_y.",
      call. = FALSE
    )
  }

  if (!"label_nudge_x" %in% names(spec)) spec$label_nudge_x <- 0
  if (!"label_nudge_y" %in% names(spec)) spec$label_nudge_y <- 0
  if (!"label_pad" %in% names(spec)) spec$label_pad <- 0.45
  if (!"node_ids" %in% names(spec)) spec$node_ids <- rep(list(NULL), nrow(spec))

  boxes <- if (has_explicit_bounds) {
    spec %>%
      dplyr::transmute(
        compartment    = as.character(compartment),
        xmin           = as.numeric(xmin),
        xmax           = as.numeric(xmax),
        ymin           = as.numeric(ymin),
        ymax           = as.numeric(ymax),
        node_ids       = node_ids,
        label_nudge_x  = as.numeric(label_nudge_x),
        label_nudge_y  = as.numeric(label_nudge_y),
        label_pad      = as.numeric(label_pad)
      )
  } else {
    purrr::map_dfr(seq_len(nrow(spec)), function(i) {
      comp <- spec$compartment[[i]]
      node_ids <- spec$node_ids[[i]]
      if (is.null(node_ids)) return(tibble())
      node_ids <- as.character(unlist(node_ids, recursive = TRUE, use.names = FALSE))
      comp_nodes <- nodes_df %>% filter(node_id %in% node_ids)
      if (nrow(comp_nodes) == 0) {
        warning("Manual compartment '", comp, "' has no matching nodes in this map.")
        return(tibble())
      }
      tibble(
        compartment   = comp,
        xmin          = min(comp_nodes$x, na.rm = TRUE) - spec$pad_x[[i]],
        xmax          = max(comp_nodes$x, na.rm = TRUE) + spec$pad_x[[i]],
        ymin          = min(comp_nodes$y, na.rm = TRUE) - spec$pad_y[[i]],
        ymax          = max(comp_nodes$y, na.rm = TRUE) + spec$pad_y[[i]],
        node_ids      = list(node_ids),
        label_nudge_x = spec$label_nudge_x[[i]],
        label_nudge_y = spec$label_nudge_y[[i]],
        label_pad     = spec$label_pad[[i]]
      )
    })
  }

  label_tbl <- purrr::map_dfr(seq_len(nrow(boxes)), function(i) {
    box <- boxes[i, , drop = FALSE]
    node_ids <- box$node_ids[[1]]
    if (is.null(node_ids) || length(node_ids) == 0) {
      nodes_in_box <- nodes_df %>%
        filter(
          !is.na(compartment),
          compartment == box$compartment[[1]],
          x >= box$xmin[[1]], x <= box$xmax[[1]],
          y >= box$ymin[[1]], y <= box$ymax[[1]]
        ) %>%
        dplyr::select(x, y)
    } else {
      node_ids <- as.character(unlist(node_ids, recursive = TRUE, use.names = FALSE))
      nodes_in_box <- nodes_df %>%
        filter(node_id %in% node_ids) %>%
        dplyr::select(x, y)
    }

    lab <- choose_compartment_label(
      box = box,
      nodes_in_box = nodes_in_box,
      pad = box$label_pad[[1]]
    )
    lab$label_x <- lab$label_x + box$label_nudge_x[[1]]
    lab$label_y <- lab$label_y + box$label_nudge_y[[1]]
    lab
  }) %>%
    dplyr::select(label_x, label_y, hjust, vjust, pos_id)

  bind_cols(
    boxes %>% dplyr::select(compartment, xmin, xmax, ymin, ymax),
    label_tbl
  )
}

build_pathway_plot <- function(node_tbl,
                               edge_df,
                               title,
                               subtitle,
                               show_compartments = FALSE,
                               compartment_pad_x = 0.6,
                               compartment_pad_y = 0.6,
                               min_nodes_per_compartment = 2,
                               force_singleton_compartments = NULL,
                               compartment_boxes_spec = NULL,
                               show_edge_pathway_legend = FALSE) {
  if (!"edge_class" %in% names(edge_df)) {
    edge_df <- edge_df %>% mutate(edge_class = "core")
  }
  if (!"edge_pathway" %in% names(edge_df)) {
    edge_df <- edge_df %>% mutate(edge_pathway = "Pathway flow")
  }
  edge_df <- edge_df %>%
    mutate(
      edge_class   = as.character(edge_class),
      edge_pathway = as.character(edge_pathway),
      curvature = dplyr::case_when(
        from == "TNFRSF1A" & to == "TRADD" ~  0.25,
        from == "TNFRSF1A" & to == "RIPK1" ~ -0.30,
        from == "CASP8"    & to == "RIPK1" ~ -0.35,
        TRUE ~ 0
      ),
      is_curved = curvature != 0
    )

  edges_joined <- edge_df %>%
    left_join(node_tbl %>% dplyr::select(node_id, x, y), by = c("from" = "node_id")) %>%
    dplyr::rename(x_from = x, y_from = y) %>%
    left_join(node_tbl %>% dplyr::select(node_id, x, y), by = c("to" = "node_id")) %>%
    dplyr::rename(x_to = x, y_to = y)

  missing_xy <- edges_joined %>%
    dplyr::filter(is.na(x_from) | is.na(x_to) | is.na(y_from) | is.na(y_to))
  if (nrow(missing_xy) > 0) {
    warning(
      "Edges dropped due to missing node coords: ",
      paste0(missing_xy$from, "->", missing_xy$to, collapse = ", ")
    )
  }
  if (all(c("TNFRSF1A", "TRADD") %in% node_tbl$node_id)) {
    stopifnot(any(edges_joined$from == "TNFRSF1A" & edges_joined$to == "TRADD"))
  }

  edges_joined_valid <- edges_joined %>%
    dplyr::filter(!is.na(x_from), !is.na(y_from), !is.na(x_to), !is.na(y_to))

  edge_overlay      <- edges_joined_valid %>% dplyr::filter(is_curved)
  edges_main_joined <- edges_joined_valid %>% dplyr::filter(!is_curved)

  edge_df_main <- edges_main_joined %>%
    dplyr::select(from, to, edge_class, edge_pathway, curvature, is_curved)
  if (nrow(edge_df_main) == 0) {
    stop("No drawable edges after coordinate join.", call. = FALSE)
  }

  g <- graph_from_data_frame(
    d        = edge_df_main,
    directed = TRUE,
    vertices = node_tbl %>%
      dplyr::rename(name = node_id) %>%
      dplyr::select(name, dplyr::everything())
  )

  lim <- max(abs(node_tbl$log2FC_used), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1

  status_levels <- c(
    "Missing gene",
    "Present; padj NA",
    "Present; not significant",
    "Present; significant"
  )
  status_colours <- c(
    "Missing gene"           = "grey45",
    "Present; padj NA"       = "black",
    "Present; not significant" = "black",
    "Present; significant"   = "black"
  )
  status_linetypes <- c(
    "Missing gene"           = "22",
    "Present; padj NA"       = "13",
    "Present; not significant" = "solid",
    "Present; significant"   = "solid"
  )
  status_linewidths <- c(
    "Missing gene"           = 0.5,
    "Present; padj NA"       = 0.75,
    "Present; not significant" = 0.45,
    "Present; significant"   = 1.2
  )

  edge_width_values <- c("core" = 0.72, "output" = 0.56)
  extra_edge_classes <- setdiff(unique(edge_df_main$edge_class), names(edge_width_values))
  if (length(extra_edge_classes) > 0) {
    edge_width_values[extra_edge_classes] <- 0.72
  }

  edge_pathway_levels <- unique(edge_df_main$edge_pathway)
  use_edge_pathway_legend <- isTRUE(show_edge_pathway_legend) && length(edge_pathway_levels) > 1L
  edge_pathway_palette <- c(
    "cgas_sting"         = "black",
    "shared_death_input" = "#1F4E79",
    "apoptosis_branch"   = "#A33A3A",
    "necroptosis_branch" = "#2F7F56",
    "pathway_flow"       = "black",
    "Pathway flow"       = "black"
  )
  missing_edge_levels <- setdiff(edge_pathway_levels, names(edge_pathway_palette))
  if (length(missing_edge_levels) > 0) {
    hue_vals <- scales::hue_pal()(length(missing_edge_levels))
    names(hue_vals) <- missing_edge_levels
    edge_pathway_palette <- c(edge_pathway_palette, hue_vals)
  }
  edge_pathway_palette <- edge_pathway_palette[unique(c(edge_pathway_levels, names(edge_pathway_palette)))]
  edge_pathway_labels <- c(
    cgas_sting         = "cGAS-STING",
    shared_death_input = "Shared death input",
    apoptosis_branch   = "Apoptosis branch",
    necroptosis_branch = "Necroptosis branch",
    pathway_flow       = "Pathway flow",
    `Pathway flow`     = "Pathway flow"
  )

  if (use_edge_pathway_legend) {
    edge_df_main <- edge_df_main %>%
      mutate(edge_pathway = factor(edge_pathway, levels = edge_pathway_levels))
  }

  edge_layer_straight <- if (use_edge_pathway_legend) {
    geom_edge_link(
      aes(
        start_cap = label_rect(node1.label, padding = margin(2.5, 3.5, 2.5, 3.5, "mm")),
        end_cap   = label_rect(node2.label, padding = margin(2.5, 3.5, 2.5, 3.5, "mm")),
        edge_width   = edge_class,
        edge_colour  = edge_pathway,
        filter       = !is_curved
      ),
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      lineend = "round",
      show.legend = TRUE
    )
  } else {
    geom_edge_link(
      aes(
        start_cap = label_rect(node1.label, padding = margin(2.5, 3.5, 2.5, 3.5, "mm")),
        end_cap   = label_rect(node2.label, padding = margin(2.5, 3.5, 2.5, 3.5, "mm")),
        edge_width = edge_class,
        filter     = !is_curved
      ),
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      edge_colour = "black",
      lineend = "round",
      show.legend = FALSE
    )
  }

  edge_layer_curved <- if (use_edge_pathway_legend) {
    geom_edge_arc(
      aes(
        edge_width  = edge_class,
        edge_colour = edge_pathway,
        filter      = is_curved
      ),
      strength = 0.35,
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      lineend = "round",
      show.legend = FALSE
    )
  } else {
    geom_edge_arc(
      aes(edge_width = edge_class, filter = is_curved),
      strength = 0.35,
      arrow = arrow(length = unit(2.8, "mm"), type = "closed"),
      edge_colour = "black",
      lineend = "round",
      show.legend = FALSE
    )
  }

  compartment_boxes <- if (isTRUE(show_compartments)) {
    if (!is.null(compartment_boxes_spec) && nrow(compartment_boxes_spec) > 0) {
      compute_manual_compartment_boxes(nodes_df = node_tbl, compartment_spec = compartment_boxes_spec)
    } else {
      infer_compartment_boxes(
        nodes_df = node_tbl,
        pad_x = compartment_pad_x,
        pad_y = compartment_pad_y,
        min_nodes_per_box = min_nodes_per_compartment,
        force_singleton_compartments = force_singleton_compartments
      )
    }
  } else {
    tibble()
  }

  if (nrow(edge_overlay) > 0) {
    edge_overlay <- edge_overlay %>%
      mutate(
        dx               = x_to - x_from,
        dy               = y_to - y_from,
        dist             = sqrt(dx^2 + dy^2),
        shrink_amt       = ifelse(dist > 0, pmin(0.25, pmax((dist / 2) - 1e-6, 0)), 0),
        ux               = ifelse(dist > 0, dx / dist, 0),
        uy               = ifelse(dist > 0, dy / dist, 0),
        x_from_plot      = x_from + ux * shrink_amt,
        y_from_plot      = y_from + uy * shrink_amt,
        x_to_plot        = x_to   - ux * shrink_amt,
        y_to_plot        = y_to   - uy * shrink_amt,
        edge_colour_plot = unname(edge_pathway_palette[as.character(edge_pathway)])
      )
    if (any(is.na(edge_overlay$edge_colour_plot))) {
      edge_overlay$edge_colour_plot[is.na(edge_overlay$edge_colour_plot)] <- "black"
    }
  }

  x_candidates <- c(node_tbl$x)
  y_candidates <- c(node_tbl$y)
  if (nrow(compartment_boxes) > 0) {
    x_candidates <- c(x_candidates, compartment_boxes$xmin, compartment_boxes$xmax)
    y_candidates <- c(y_candidates, compartment_boxes$ymin, compartment_boxes$ymax)
  }
  x_lim <- c(min(x_candidates, na.rm = TRUE) - 1.20, max(x_candidates, na.rm = TRUE) + 0.80)
  y_lim <- c(min(y_candidates, na.rm = TRUE) - 0.80, max(y_candidates, na.rm = TRUE) + 0.80)

  ggraph(g, layout = "manual", x = x, y = y) +
    {
      if (nrow(compartment_boxes) > 0) {
        geom_rect(
          data = compartment_boxes,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill      = "grey95",
          alpha     = 0.10,
          color     = "grey85",
          linewidth = 0.3
        )
      } else {
        NULL
      }
    } +
    {
      if (nrow(compartment_boxes) > 0) {
        geom_text(
          data = compartment_boxes,
          aes(x = label_x, y = label_y, label = compartment, hjust = hjust, vjust = vjust),
          inherit.aes = FALSE,
          color    = "grey40",
          size     = 3.0,
          fontface = "italic",
          alpha    = 0.85
        )
      } else {
        NULL
      }
    } +
    edge_layer_straight +
    edge_layer_curved +
    scale_edge_width_manual(values = edge_width_values, guide = "none") +
    geom_node_label(
      aes(
        label     = label,
        fill      = log2FC_plot,
        colour    = node_status,
        linetype  = node_status,
        linewidth = node_status
      ),
      size          = 3.0,
      family        = "sans",
      fontface      = "bold",
      text.colour   = "black",
      label.padding = unit(0.22, "lines"),
      label.r       = unit(0.08, "lines"),
      lineheight    = 0.95
    ) +
    {
      if (nrow(edge_overlay) > 0) {
        purrr::map(seq_len(nrow(edge_overlay)), function(i) {
          row <- edge_overlay[i, , drop = FALSE]
          geom_curve(
            data = row,
            aes(x = x_from_plot, y = y_from_plot, xend = x_to_plot, yend = y_to_plot),
            inherit.aes = FALSE,
            curvature   = row$curvature,
            linewidth   = 0.9,
            lineend     = "round",
            colour      = row$edge_colour_plot,
            arrow       = grid::arrow(type = "closed", length = grid::unit(0.20, "cm")),
            show.legend = FALSE
          )
        })
      } else {
        NULL
      }
    } +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-lim, lim),
      oob      = squish,
      na.value = "grey90",
      name     = "log2FC"
    ) +
    scale_colour_manual(
      values = status_colours,
      breaks = status_levels,
      drop   = FALSE,
      name   = "Node status"
    ) +
    scale_linetype_manual(
      values = status_linetypes,
      breaks = status_levels,
      drop   = FALSE,
      name   = "Node status"
    ) +
    scale_linewidth_manual(
      values = status_linewidths,
      breaks = status_levels,
      drop   = FALSE,
      name   = "Node status"
    ) +
    guides(
      fill   = guide_colorbar(order = 1),
      colour = guide_legend(
        order = 2,
        override.aes = list(
          fill      = c("grey90", "grey85", "grey85", "grey85"),
          linetype  = unname(status_linetypes[status_levels]),
          linewidth = unname(status_linewidths[status_levels]),
          colour    = unname(status_colours[status_levels])
        )
      ),
      linetype  = "none",
      linewidth = "none"
    ) +
    {
      if (use_edge_pathway_legend) {
        scale_edge_colour_manual(
          values = edge_pathway_palette,
          breaks = edge_pathway_levels,
          drop   = FALSE,
          labels = ifelse(
            edge_pathway_levels %in% names(edge_pathway_labels),
            unname(edge_pathway_labels[edge_pathway_levels]),
            tools::toTitleCase(gsub("_", " ", edge_pathway_levels))
          ),
          name   = "Edge pathway"
        )
      } else {
        scale_edge_colour_manual(values = c("Pathway flow" = "black"), guide = "none")
      }
    } +
    {
      if (use_edge_pathway_legend) {
        guides(edge_colour = guide_legend(order = 3))
      } else {
        guides(edge_colour = "none")
      }
    } +
    labs(
      title    = title,
      subtitle = subtitle,
      caption  = "Missing gene = gray dashed outline; padj NA = dotted outline; significant = thick outline."
    ) +
    coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE, clip = "off") +
    theme_void(base_size = 12) +
    theme(
      plot.background       = element_rect(fill = "white", colour = NA),
      panel.background      = element_rect(fill = "white", colour = NA),
      legend.background     = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      legend.position       = "right",
      plot.title            = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle         = element_text(size = 10, hjust = 0.5),
      plot.caption          = element_text(size = 9, colour = "grey30")
    )
}

write_pathway_outputs <- function(prefix, outdir, plot_obj, node_tbl, de_annot) {
  png_file        <- file.path(outdir, paste0(prefix, "_map.png"))
  nodes_file      <- file.path(outdir, paste0(prefix, "_nodes.tsv"))
  missing_file    <- file.path(outdir, paste0(prefix, "_missing.tsv"))
  all_results_file <- file.path(outdir, paste0(prefix, "_all_results.csv"))

  ggsave(
    filename = png_file,
    plot     = plot_obj,
    width    = 14,
    height   = 7.5,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )

  node_tbl %>%
    dplyr::transmute(
      gene_id        = gene_id,
      symbol         = symbol_key,
      label          = label,
      log2FC_used    = log2FC_used,
      padj           = padj,
      present_in_data = present_in_data,
      node_status    = node_status,
      signif_flag    = signif_flag
    ) %>%
    readr::write_tsv(nodes_file, na = "NA")

  node_tbl %>%
    filter(node_status == "Missing gene") %>%
    dplyr::transmute(symbol = symbol_key, label = label, module = module) %>%
    readr::write_tsv(missing_file, na = "NA")

  de_annot %>%
    dplyr::select(ensembl_id, symbol, entrez_id, log2FC_used, padj, dplyr::everything()) %>%
    readr::write_csv(all_results_file, na = "NA")

  list(
    png         = png_file,
    nodes       = nodes_file,
    missing     = missing_file,
    all_results = all_results_file
  )
}

make_interpretation_panel <- function(nodes_annot,
                                      edges_cfg = NULL,
                                      title = "Interpretation notes",
                                      top_n = 3,
                                      p_cutoff = 0.05,
                                      validation_notes = NULL) {
  boundary_text <- paste(
    "- DESeq2 overlay shows transcriptional change, not protein activation.",
    "- Edges are canonical pathway wiring, not inferred causality in this dataset.",
    "- Pathway assignment requires orthogonal validation in tissue.",
    sep = "\n"
  )

  if (!is.null(validation_notes) && length(validation_notes) > 0) {
    validation_text <- paste(
      "Suggested validation:",
      paste0("- ", validation_notes, collapse = "\n"),
      sep = "\n"
    )
  } else {
    validation_text <- "No pathway-specific validation notes provided."
  }

  top_text <- "Top changing genes by module: no mapped genes with finite log2FC."
  req_cols <- c("module", "label", "log2FC_used", "present_in_data", "padj")
  if (all(req_cols %in% names(nodes_annot))) {
    top_tbl <- nodes_annot %>%
      filter(present_in_data, !is.na(log2FC_used), !is.na(module), module != "") %>%
      group_by(module) %>%
      slice_max(order_by = abs(log2FC_used), n = top_n, with_ties = FALSE) %>%
      ungroup()

    if (nrow(top_tbl) > 0) {
      module_lines <- top_tbl %>%
        mutate(
          sig_tag  = ifelse(!is.na(padj) & padj < p_cutoff, "*", ""),
          gene_txt = paste0(label, "(", format(round(log2FC_used, 2), nsmall = 2), sig_tag, ")")
        ) %>%
        group_by(module) %>%
        summarise(
          line = paste0(dplyr::first(module), ": ", paste(gene_txt, collapse = ", ")),
          .groups = "drop"
        ) %>%
        pull(line)

      top_text <- paste(
        "Top changing genes by module (abs log2FC; * = padj < 0.05):",
        paste0("- ", module_lines, collapse = "\n"),
        sep = "\n"
      )
    }
  }

  ggplot() +
    annotate("text", x = 0, y = 1.00, label = title,           hjust = 0, vjust = 1, size = 5.0, fontface = "bold") +
    annotate("text", x = 0, y = 0.88, label = boundary_text,   hjust = 0, vjust = 1, size = 3.25, lineheight = 1.15) +
    annotate("text", x = 0, y = 0.58, label = validation_text, hjust = 0, vjust = 1, size = 3.20, lineheight = 1.15) +
    annotate("text", x = 0, y = 0.22, label = top_text,        hjust = 0, vjust = 1, size = 3.00, lineheight = 1.10) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin      = margin(12, 8, 12, 12)
    )
}

save_multipanel_figure <- function(p_left,
                                   p_right,
                                   outfile,
                                   width      = 16,
                                   height     = 8,
                                   dpi        = 300,
                                   rel_widths = c(3.2, 1.4)) {
  grDevices::png(
    filename = outfile,
    width    = width,
    height   = height,
    units    = "in",
    res      = dpi,
    bg       = "white"
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  layout <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(rel_widths, "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  print(p_left,  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::popViewport()
  invisible(outfile)
}

render_pathway_map <- function(config,
                               prepared_de,
                               outdir,
                               title,
                               p_cutoff = 0.05,
                               show_compartments = FALSE,
                               compartment_pad_x = 0.6,
                               compartment_pad_y = 0.6,
                               min_nodes_per_compartment = 2,
                               force_singleton_compartments = NULL,
                               return_plot = FALSE) {
  node_tbl  <- build_pathway_node_table(config, prepared_de$de_by_symbol, p_cutoff = p_cutoff)
  matched_n <- sum(node_tbl$present_in_data, na.rm = TRUE)
  total_n   <- nrow(node_tbl)
  subtitle_txt <- paste0(prepared_de$input_basename, " | matched: ", matched_n, "/", total_n, " pathway genes")
  plot_title   <- paste0(title, " - ", config$display_name)

  plot_obj <- build_pathway_plot(
    node_tbl                    = node_tbl,
    edge_df                     = config$edges,
    title                       = plot_title,
    subtitle                    = subtitle_txt,
    show_compartments           = show_compartments,
    compartment_pad_x           = compartment_pad_x,
    compartment_pad_y           = compartment_pad_y,
    min_nodes_per_compartment   = min_nodes_per_compartment,
    force_singleton_compartments = force_singleton_compartments,
    compartment_boxes_spec      = if ("compartment_boxes" %in% names(config)) config$compartment_boxes else NULL
  )

  if (isTRUE(return_plot)) {
    return(list(
      plot     = plot_obj,
      node_tbl = node_tbl,
      subtitle = subtitle_txt,
      title    = plot_title
    ))
  }

  write_pathway_outputs(
    prefix  = config$name,
    outdir  = outdir,
    plot_obj = plot_obj,
    node_tbl = node_tbl,
    de_annot = prepared_de$de_annot
  )
}

# ---------------------------------------------------------------------------
# run_pathway_maps() — top-level orchestrator
# ---------------------------------------------------------------------------

run_pathway_maps <- function(input,
                             outdir,
                             title,
                             pathways = list_pathways(),
                             p_cutoff = 0.05,
                             lfc_column_preference = c("log2FC_shrunken", "log2FC"),
                             id_column_candidates  = c("ensembl_id", "Unnamed: 0", "gene", "id"),
                             stitch_panels = TRUE) {
  if (missing(input) || !nzchar(input)) {
    stop("`input` is required.", call. = FALSE)
  }
  if (missing(outdir) || !nzchar(outdir)) {
    stop("`outdir` is required.", call. = FALSE)
  }
  if (missing(title) || !nzchar(title)) {
    title <- "Pathway maps"
  }

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  prepared_de <- prepare_de_table_for_pathways(
    input                 = input,
    lfc_column_preference = lfc_column_preference,
    id_column_candidates  = id_column_candidates
  )

  valid     <- list_pathways()
  requested <- unique(tolower(pathways))
  invalid   <- setdiff(requested, valid)
  if (length(invalid) > 0) {
    stop(
      "Unknown pathway(s): ", paste(invalid, collapse = ", "),
      "\nAvailable: ", paste(valid, collapse = ", "),
      call. = FALSE
    )
  }

  out          <- list()
  plot_objects <- list()
  all_node_tbls <- list()

  for (pathway_name in requested) {
    cfg <- get_pathway_config(pathway_name)

    result <- render_pathway_map(
      config                       = cfg,
      prepared_de                  = prepared_de,
      outdir                       = outdir,
      title                        = title,
      p_cutoff                     = p_cutoff,
      show_compartments            = TRUE,
      min_nodes_per_compartment    = 2,
      force_singleton_compartments = c("Plasma membrane", "ER/Golgi"),
      return_plot                  = TRUE
    )

    write_pathway_outputs(
      prefix   = cfg$name,
      outdir   = outdir,
      plot_obj = result$plot,
      node_tbl = result$node_tbl,
      de_annot = prepared_de$de_annot
    )

    p_interpret <- make_interpretation_panel(
      nodes_annot      = result$node_tbl,
      edges_cfg        = cfg$edges,
      title            = paste0(cfg$display_name, " \u2014 Interpretation"),
      validation_notes = cfg$validation_notes
    )
    panel_png <- file.path(outdir, paste0(cfg$name, "_panel.png"))
    save_multipanel_figure(
      p_left  = result$plot,
      p_right = p_interpret,
      outfile = panel_png
    )

    all_node_tbls[[pathway_name]] <- result$node_tbl %>%
      dplyr::mutate(pathway = pathway_name, .before = 1)

    out[[pathway_name]]            <- result
    out[[pathway_name]]$panel_png  <- panel_png
    plot_objects[[pathway_name]]   <- result$plot
  }

  if (length(all_node_tbls) > 0) {
    combined_nodes_file <- file.path(outdir, "all_pathways_nodes.csv")
    dplyr::bind_rows(all_node_tbls) %>%
      dplyr::transmute(
        pathway,
        gene_id         = gene_id,
        symbol          = symbol_key,
        label           = label,
        module          = module,
        log2FC_used     = log2FC_used,
        padj            = padj,
        present_in_data = present_in_data,
        node_status     = node_status,
        signif_flag     = signif_flag
      ) %>%
      readr::write_csv(combined_nodes_file, na = "NA")
    out[["all_pathways_nodes_csv"]] <- combined_nodes_file
  }

  if (isTRUE(stitch_panels) && length(plot_objects) > 1) {
    overview_png <- file.path(outdir, "all_pathways_overview.png")
    n_panels <- length(plot_objects)
    nc       <- min(n_panels, 2)
    nr       <- ceiling(n_panels / nc)
    combined <- cowplot::plot_grid(
      plotlist = plot_objects,
      ncol     = nc,
      labels   = "AUTO"
    )
    cowplot::ggsave2(
      filename  = overview_png,
      plot      = combined,
      width     = 14 * nc,
      height    = 7.5 * nr,
      dpi       = 200,
      bg        = "white",
      limitsize = FALSE
    )
    out[["overview_png"]] <- overview_png
  }

  invisible(out)
}

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

parse_args <- function(args) {
  out <- list(
    input     = NULL,
    outdir    = "pathway_out",
    title     = "Pathway Map",
    pathways  = NULL,
    help      = FALSE,
    list_pathways = FALSE
  )

  if (length(args) == 0L) return(out)

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      out$help <- TRUE
      return(out)
    }
    if (key == "--list-pathways") {
      out$list_pathways <- TRUE
      return(out)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (key == "--input") {
      out$input <- value
    } else if (key == "--outdir") {
      out$outdir <- value
    } else if (key == "--title") {
      out$title <- value
    } else if (key == "--pathways") {
      out$pathways <- value
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    i <- i + 2L
  }
  out
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) == 0L || any(args %in% c("-h", "--help"))) {
    cat(
      "Usage:\n",
      "  Rscript pathway_mapper_v3.R --input <csv> --outdir <dir> [options]\n\n",
      "Options:\n",
      "  --input <csv>       DESeq2 results CSV (Ensembl IDs)\n",
      "  --outdir <dir>      Output directory\n",
      "  --pathways <list>   Comma-separated pathway names (default: all)\n",
      "  --title <string>    Plot title prefix\n",
      "  --list-pathways     Print available pathways and exit\n",
      "  -h, --help          Show this help\n",
      sep = ""
    )
    return(invisible(NULL))
  }

  if ("--list-pathways" %in% args) {
    cat("Available pathways:\n")
    cat(paste0("  - ", list_pathways(), collapse = "\n"), "\n")
    return(invisible(NULL))
  }

  parsed <- parse_args(args)

  pathways <- if (!is.null(parsed$pathways)) {
    strsplit(parsed$pathways, ",")[[1]]
  } else {
    list_pathways()
  }

  if (is.null(parsed$input) || !nzchar(parsed$input)) {
    stop("--input is required.", call. = FALSE)
  }

  run_pathway_maps(
    input    = parsed$input,
    outdir   = parsed$outdir,
    title    = parsed$title,
    pathways = pathways
  )
}

# Execute only when run as a script (not when sourced in RStudio).
if (sys.nframe() == 0L) {
  main()
}
