#!/usr/bin/env Rscript

# Fit the frozen joint model to a pre-specified country subset. This sensitivity
# analysis tests whether shared lineage-growth estimates depend on any one
# country and reports single-country estimates without treating them as
# independent replications of the primary model.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) {
  stop(
    paste(
      "Usage: gtd_29_fit_country_subset.R COMPILED_MODEL_RDS DATA_JSON",
      "KEEP_COUNTRIES_CSV OUTDIR CHAINS ITER [SEED]"
    )
  )
}

compiled_model_file <- normalizePath(args[[1L]], mustWork = TRUE)
data_file <- normalizePath(args[[2L]], mustWork = TRUE)
keep_names <- strsplit(args[[3L]], ",", fixed = TRUE)[[1L]]
outdir <- args[[4L]]
chains <- as.integer(args[[5L]])
iter <- as.integer(args[[6L]])
seed <- if (length(args) >= 7L) as.integer(args[[7L]]) else 20260811L

all_countries <- c("AUS", "CHN", "JPN")
lineages <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
if (!length(keep_names) || any(!keep_names %chin% all_countries)) {
  stop("KEEP_COUNTRIES_CSV must contain AUS, CHN, and/or JPN")
}
keep_names <- all_countries[all_countries %chin% keep_names]
keep_old <- match(keep_names, all_countries)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
model_data <- fromJSON(data_file, simplifyVector = TRUE)

obs_keep <- which(model_data$obs_country %in% keep_old)
old_project <- model_data$obs_project[obs_keep]
used_projects <- sort(unique(old_project))
project_map <- setNames(seq_along(used_projects), used_projects)
country_map <- setNames(seq_along(keep_old), keep_old)

subset_data <- list(
  C = length(keep_old),
  T = model_data$T,
  L = model_data$L,
  P = length(used_projects),
  J = length(obs_keep),
  K = model_data$K,
  use_project_effects = model_data$use_project_effects,
  cases = model_data$cases[keep_old, , drop = FALSE],
  B = model_data$B,
  reporting_change = model_data$reporting_change[keep_old, , drop = FALSE],
  initial_alpha = model_data$initial_alpha[keep_old, , drop = FALSE],
  import_exposure = model_data$import_exposure[keep_old, , , drop = FALSE],
  obs_country = unname(as.integer(country_map[as.character(
    model_data$obs_country[obs_keep]
  )])),
  obs_month = model_data$obs_month[obs_keep],
  obs_project = unname(as.integer(project_map[as.character(old_project)])),
  y_genome = model_data$y_genome[obs_keep, , drop = FALSE]
)

write_json(
  subset_data,
  file.path(outdir, "joint_model_data.json"),
  auto_unbox = TRUE,
  digits = NA
)

model <- readRDS(compiled_model_file)
options(mc.cores = min(chains, parallel::detectCores()))
rstan_options(auto_write = TRUE)

fit <- sampling(
  model,
  data = subset_data,
  chains = chains,
  iter = iter,
  warmup = floor(iter / 2),
  seed = seed,
  refresh = max(1L, floor(iter / 20)),
  control = list(adapt_delta = 0.95, max_treedepth = 13)
)

theta <- rstan::extract(
  fit,
  pars = "lineage_relative_transmission",
  permuted = TRUE
)$lineage_relative_transmission

summaries <- rbindlist(lapply(seq_along(lineages), function(l) {
  values <- theta[, l]
  data.table(
    country_subset = paste(keep_names, collapse = "+"),
    n_countries = length(keep_names),
    lineage = lineages[[l]],
    mean = mean(values),
    median = median(values),
    lower_95 = unname(quantile(values, 0.025)),
    upper_95 = unname(quantile(values, 0.975)),
    posterior_probability_above_one = mean(values > 1)
  )
}))
fwrite(
  summaries,
  file.path(outdir, "lineage_relative_transmission.tsv"),
  sep = "\t"
)
saveRDS(
  list(
    country_subset = keep_names,
    lineage_relative_transmission = theta
  ),
  file.path(outdir, "lineage_relative_transmission_draws.rds"),
  compress = "xz"
)

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
  country_subset = keep_names,
  chains = chains,
  iterations_per_chain = iter,
  post_warmup_draws = chains * (iter - floor(iter / 2)),
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
