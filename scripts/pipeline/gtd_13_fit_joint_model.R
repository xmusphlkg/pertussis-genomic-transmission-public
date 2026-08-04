#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(rstan)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: gtd_13_fit_joint_model.R STAN DATA_JSON OUTDIR [CHAINS] [ITER] [SEED]")
}
stan_file <- normalizePath(args[[1L]], mustWork = TRUE)
data_file <- normalizePath(args[[2L]], mustWork = TRUE)
outdir <- args[[3L]]
chains <- if (length(args) >= 4L) as.integer(args[[4L]]) else 4L
iter <- if (length(args) >= 5L) as.integer(args[[5L]]) else 2000L
seed <- if (length(args) >= 6L) as.integer(args[[6L]]) else 20260725L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

options(mc.cores = min(chains, parallel::detectCores()))
rstan_options(auto_write = TRUE)
stan_data <- fromJSON(data_file, simplifyVector = TRUE)

message("Compiling joint transmission-sampling model")
model <- stan_model(file = stan_file, model_name = "gtd_joint_transmission_sampling")
saveRDS(model, file.path(outdir, "compiled_model.rds"))

message(sprintf("Sampling %d chains, %d iterations per chain", chains, iter))
fit <- sampling(
  model,
  data = stan_data,
  chains = chains,
  iter = iter,
  warmup = floor(iter / 2),
  seed = seed,
  refresh = max(1L, floor(iter / 20)),
  control = list(adapt_delta = 0.95, max_treedepth = 13)
)
saveRDS(fit, file.path(outdir, "joint_model_fit.rds"), compress = "xz")

key_pars <- c(
  "lineage_relative_transmission", "post_import_fraction", "import_scale",
  "density_feedback",
  "reporting_jump", "phi_cases", "sigma_project", "genome_concentration"
)
sm <- as.data.table(summary(fit, pars = key_pars, probs = c(0.025, 0.5, 0.975))$summary,
                    keep.rownames = "parameter")
fwrite(sm, file.path(outdir, "posterior_key_parameters.tsv"), sep = "\t")

all_sm <- summary(
  fit,
  pars = c("r_coef", "log_theta_raw", "import_scale", "density_feedback",
           "q0", "project_raw",
           "sigma_project", "reporting_jump", "phi_cases",
           "genome_concentration"),
  probs = c(0.025, 0.5, 0.975)
)$summary
all_dt <- as.data.table(all_sm, keep.rownames = "parameter")
fwrite(all_dt, file.path(outdir, "posterior_sampling_parameters.tsv"), sep = "\t")

sampler <- get_sampler_params(fit, inc_warmup = FALSE)
n_divergent <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
n_max_depth <- sum(vapply(
  sampler, function(x) sum(x[, "treedepth__"] >= 13), numeric(1)
))
diagnostics <- list(
  chains = chains,
  iterations_per_chain = iter,
  post_warmup_draws = chains * (iter - floor(iter / 2)),
  seed = seed,
  divergent_transitions = n_divergent,
  maximum_treedepth_hits = n_max_depth,
  maximum_rhat = max(all_sm[, "Rhat"], na.rm = TRUE),
  minimum_bulk_neff = min(all_sm[, "n_eff"], na.rm = TRUE)
)
write_json(
  diagnostics,
  file.path(outdir, "sampling_diagnostics.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

draws <- extract(
  fit,
  pars = c(
    "r_coef", "log_theta", "import_scale", "density_feedback",
    "reporting_jump", "lineage_relative_transmission", "post_import_fraction", "q",
    "project_effect", "mu_cases", "cases_rep"
  ),
  permuted = TRUE
)
saveRDS(draws, file.path(outdir, "posterior_outputs.rds"), compress = "xz")
