#!/usr/bin/env Rscript

# Refit the complete joint model to an alternative frozen input definition.
# The compact posterior output is sufficient for the standard model summaries
# and conditional scenarios, while avoiding storage of the full rstan object.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7L) {
  stop(
    paste(
      "Usage: gtd_32_fit_input_sensitivity.R COMPILED_MODEL_RDS DATA_JSON",
      "ANALYSIS_LABEL OUTDIR CHAINS ITER SEED"
    )
  )
}

compiled_model_file <- normalizePath(args[[1L]], mustWork = TRUE)
data_file <- normalizePath(args[[2L]], mustWork = TRUE)
analysis_label <- args[[3L]]
outdir <- args[[4L]]
chains <- as.integer(args[[5L]])
iter <- as.integer(args[[6L]])
seed <- as.integer(args[[7L]])
lineages <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
model_data <- fromJSON(data_file, simplifyVector = TRUE)
model <- readRDS(compiled_model_file)
options(mc.cores = min(chains, parallel::detectCores()))
rstan_options(auto_write = TRUE)

fit <- sampling(
  model,
  data = model_data,
  chains = chains,
  iter = iter,
  warmup = floor(iter / 2),
  seed = seed,
  refresh = max(1L, floor(iter / 20)),
  control = list(adapt_delta = 0.95, max_treedepth = 13)
)

draws <- extract(
  fit,
  pars = c(
    "r_coef", "log_theta", "import_scale", "density_feedback",
    "reporting_jump", "lineage_relative_transmission", "q",
    "mu_cases", "cases_rep"
  ),
  permuted = TRUE
)
saveRDS(
  draws,
  file.path(outdir, "posterior_outputs.rds"),
  compress = "xz"
)

qsum <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    lower_95 = unname(quantile(x, 0.025)),
    upper_95 = unname(quantile(x, 0.975))
  )
}

lineage_summary <- rbindlist(lapply(seq_along(lineages), function(l) {
  values <- draws$lineage_relative_transmission[, l]
  data.table(
    analysis = analysis_label,
    lineage = lineages[[l]],
    t(qsum(values)),
    posterior_probability_above_one = mean(values > 1)
  )
}))
fwrite(
  lineage_summary,
  file.path(outdir, "lineage_relative_transmission.tsv"),
  sep = "\t"
)

l10207_idx <- match("L1_02.07", lineages)
pairwise <- rbindlist(lapply(
  setdiff(seq_along(lineages), l10207_idx),
  function(l) {
    ratio <- draws$lineage_relative_transmission[, l10207_idx] /
      draws$lineage_relative_transmission[, l]
    data.table(
      analysis = analysis_label,
      numerator = "L1_02.07",
      denominator = lineages[[l]],
      t(qsum(ratio)),
      posterior_probability_above_one = mean(ratio > 1)
    )
  }
))
fwrite(pairwise, file.path(outdir, "l10207_pairwise_growth.tsv"), sep = "\t")

sampler <- get_sampler_params(fit, inc_warmup = FALSE)
fit_summary <- summary(
  fit,
  pars = c(
    "log_theta_raw", "import_scale", "density_feedback",
    "reporting_jump", "phi_cases", "sigma_project",
    "genome_concentration"
  )
)$summary
diagnostics <- list(
  analysis = analysis_label,
  data_file = data_file,
  compiled_model_file = compiled_model_file,
  chains = chains,
  iterations_per_chain = iter,
  post_warmup_draws = chains * (iter - floor(iter / 2)),
  seed = seed,
  divergent_transitions = sum(vapply(
    sampler,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  )),
  maximum_treedepth_hits = sum(vapply(
    sampler,
    function(x) sum(x[, "treedepth__"] >= 13),
    numeric(1)
  )),
  maximum_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
  minimum_neff = min(fit_summary[, "n_eff"], na.rm = TRUE)
)
write_json(
  diagnostics,
  file.path(outdir, "sampling_diagnostics.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
