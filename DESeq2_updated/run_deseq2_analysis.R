options(stringsAsFactors = FALSE)

load_driver_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing required package '", pkg, "'. Install it before running the pipeline.",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_driver_package("yaml")

source("R/io.R")
source("R/cleaning.R")
source("R/validation.R")
source("R/deseq_model.R")
source("R/qc_plots.R")
source("R/contrasts.R")
source("R/gsea.R")
source("R/exports.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else "config/analysis_config.yaml"

config <- load_config(config_path)
config <- load_required_packages(config)
dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- load_inputs(config)
cleaned <- clean_inputs(inputs, config)
validate_inputs(cleaned, config)
fit <- fit_deseq_pipeline(cleaned, config)
make_qc_outputs(fit, config)
contrast_results <- run_requested_contrasts(fit, cleaned, config)
write_all_outputs(fit, cleaned, contrast_results, config)
print_run_summary(fit, cleaned, contrast_results, config)
