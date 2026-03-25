cGAS-STING + death-module plotting package

Included files:
- cgast_sting_map.R (latest script)
- cgast_sting_map_notebook_v6.Rmd (R Notebook)
- DESeq2_GRCg7b_Sacral_vs_Free_by_ENS.csv (example input)

Quick run (R console):
source("cgast_sting_map.R")
run_cgast_sting_map(
  input = "DESeq2_GRCg7b_Sacral_vs_Free_by_ENS.csv",
  outdir = "outputs",
  title = "Sacral vs Free (cGAS-STING map)"
)
extra <- run_pathway_maps(
  input = "DESeq2_GRCg7b_Sacral_vs_Free_by_ENS.csv",
  outdir = "outputs",
  title = "Sacral vs Free (cGAS-STING map)",
  pathways = c("apoptosis", "necroptosis"),
  include_combined = TRUE
)

Expected outputs in ./outputs:
- cgast_sting_map_SVG.svg
- cgast_sting_map_PNG.png
- apoptosis_map.png
- necroptosis_map.png
- combined_death_map.png
- combined_death_multipanel.png
- and node/missing/all_results tables
