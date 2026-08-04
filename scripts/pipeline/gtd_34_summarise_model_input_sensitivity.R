#!/usr/bin/env Rscript

# Consolidate the primary and alternative phylogeographic-exposure input
# definitions after each has been propagated through a full model refit.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage: gtd_34_summarise_model_input_sensitivity.R",
      "RESULTS_ROOT OUTPUT_DIR"
    )
  )
}

results_root <- normalizePath(args[[1L]], mustWork = TRUE)
outdir <- args[[2L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

analyses <- data.table(
  analysis_id = c(
    "primary_lower_0_5",
    "threshold_0_7",
    "threshold_0_9",
    "time_midpoint",
    "time_interval_uniform",
    "alternative_root"
  ),
  analysis_label = c(
    "Primary: threshold 0.5, lower-bound time",
    "Transition threshold 0.7",
    "Transition threshold 0.9",
    "Earliest-sample interval midpoint",
    "Earliest-sample interval-uniform",
    "Alternative root"
  ),
  result_directory = c(
    "model_main",
    file.path("model_input_sensitivity", c(
      "threshold_0_7",
      "threshold_0_9",
      "time_midpoint",
      "time_interval_uniform",
      "alternative_root"
    ))
  )
)

read_result <- function(spec, filename) {
  path <- file.path(results_root, spec$result_directory, filename)
  if (!file.exists(path)) {
    stop("Missing sensitivity result: ", path)
  }
  frame <- fread(path)
  frame[
    ,
    `:=`(
      analysis_id = spec$analysis_id,
      analysis_label = spec$analysis_label
    )
  ]
  frame
}

growth <- rbindlist(lapply(seq_len(nrow(analyses)), function(index) {
  spec <- analyses[index]
  read_result(spec, "lineage_relative_transmission.tsv")[
    lineage == "L1_02.07"
  ]
}), fill = TRUE)
setcolorder(
  growth,
  c(
    "analysis_id", "analysis_label", "lineage",
    "mean", "median", "lower_95", "upper_95"
  )
)
fwrite(
  growth,
  file.path(outdir, "l10207_input_sensitivity.tsv"),
  sep = "\t"
)

pairwise <- rbindlist(lapply(seq_len(nrow(analyses)), function(index) {
  read_result(analyses[index], "l10207_pairwise_growth.tsv")
}), fill = TRUE)
setcolorder(
  pairwise,
  c(
    "analysis_id", "analysis_label", "numerator", "denominator",
    "mean", "median", "lower_95", "upper_95",
    "posterior_probability_above_one"
  )
)
fwrite(
  pairwise,
  file.path(outdir, "l10207_pairwise_input_sensitivity.tsv"),
  sep = "\t"
)

scenario_keys <- c(
  "no_new_exposure_case_difference_fraction",
  "l10207_growth_scenario_difference_fraction"
)
scenarios <- rbindlist(lapply(seq_len(nrow(analyses)), function(index) {
  read_result(analyses[index], "counterfactual_summary.tsv")[
    scenario %chin% scenario_keys
  ]
}), fill = TRUE)
scenarios <- scenarios[
  ,
  .(
    analysis_id, analysis_label, country_iso3, scenario,
    mean, median, lower_95, upper_95
  )
]
setcolorder(
  scenarios,
  c(
    "analysis_id", "analysis_label", "country_iso3", "scenario",
    "mean", "median", "lower_95", "upper_95"
  )
)
fwrite(
  scenarios,
  file.path(outdir, "scenario_input_sensitivity.tsv"),
  sep = "\t"
)

diagnostics <- rbindlist(lapply(seq_len(nrow(analyses)), function(index) {
  spec <- analyses[index]
  path <- file.path(
    results_root,
    spec$result_directory,
    "sampling_diagnostics.json"
  )
  payload <- fromJSON(path, simplifyVector = TRUE)
  data.table(
    analysis_id = spec$analysis_id,
    analysis_label = spec$analysis_label,
    chains = payload$chains,
    iterations_per_chain = payload$iterations_per_chain,
    post_warmup_draws = payload$post_warmup_draws,
    divergent_transitions = payload$divergent_transitions,
    maximum_treedepth_hits = payload$maximum_treedepth_hits,
    maximum_rhat = payload$maximum_rhat,
    minimum_effective_sample_size = if (
      !is.null(payload$minimum_neff)
    ) {
      payload$minimum_neff
    } else {
      payload$minimum_bulk_neff
    }
  )
}), fill = TRUE)
fwrite(
  diagnostics,
  file.path(outdir, "model_input_sensitivity_diagnostics.tsv"),
  sep = "\t"
)

fwrite(
  analyses,
  file.path(outdir, "model_input_sensitivity_manifest.tsv"),
  sep = "\t"
)
