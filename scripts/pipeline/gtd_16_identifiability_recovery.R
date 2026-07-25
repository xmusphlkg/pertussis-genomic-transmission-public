#!/usr/bin/env Rscript

# Realistic parameter-recovery experiment for the joint model. Each replicate
# draws a truth from the fitted posterior, simulates cases and genomic strata
# under the fitted design, and refits the original model.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop(paste(
    "Usage: gtd_16_identifiability_recovery.R",
    "COMPILED_MODEL_RDS FIT_RDS DATA_JSON OUTDIR [N_REP] [CHAINS] [ITER] [SEED]"
  ))
}
model_file <- args[[1L]]
fit_file <- args[[2L]]
data_file <- args[[3L]]
outdir <- args[[4L]]
n_rep <- if (length(args) >= 5L) as.integer(args[[5L]]) else 6L
chains <- if (length(args) >= 6L) as.integer(args[[6L]]) else 2L
iter <- if (length(args) >= 7L) as.integer(args[[7L]]) else 1200L
seed <- if (length(args) >= 8L) as.integer(args[[8L]]) else 20260801L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

model <- readRDS(model_file)
fit <- readRDS(fit_file)
base <- fromJSON(data_file, simplifyVector = TRUE)
draws <- rstan::extract(
  fit,
  pars = c(
    "r_coef", "log_theta_raw", "import_scale", "density_feedback", "q0",
    "project_raw", "sigma_project", "reporting_jump", "phi_cases",
    "genome_concentration"
  ),
  permuted = TRUE
)

softmax <- function(x) {
  z <- exp(x - max(x))
  z / sum(z)
}

simulate_data <- function(draw_id, replicate_seed) {
  set.seed(replicate_seed)
  C <- base$C
  T <- base$T
  L <- base$L
  log_theta <- draws$log_theta_raw[draw_id, ] -
    mean(draws$log_theta_raw[draw_id, ])
  theta <- exp(log_theta)
  cases <- matrix(0L, nrow = C, ncol = T)
  cases[, 1] <- base$cases[, 1]
  q <- array(0, dim = c(C, T, L))
  q[, 1, ] <- draws$q0[draw_id, , ]

  for (c in seq_len(C)) {
    for (t in 2:T) {
      log_r <- sum(base$B[t, ] * draws$r_coef[draw_id, c, ]) -
        draws$density_feedback[draw_id, c] *
        log1p(cases[c, t - 1] / 1000)
      local <- exp(log_r + log_theta) * q[c, t - 1, ] *
        (cases[c, t - 1] + 0.5)
      imported <- draws$import_scale[draw_id, c] *
        base$import_exposure[c, t, ]
      component <- local + imported + 1e-9
      q[c, t, ] <- component / sum(component)
      mu <- sum(component) * exp(
        draws$reporting_jump[draw_id, c] *
        (base$reporting_change[c, t] - base$reporting_change[c, t - 1])
      )
      cases[c, t] <- rnbinom(
        1, mu = mu, size = draws$phi_cases[draw_id, c]
      )
    }
  }

  y_genome <- matrix(0L, nrow = base$J, ncol = L)
  for (j in seq_len(base$J)) {
    c <- base$obs_country[j]
    t <- base$obs_month[j]
    p <- base$obs_project[j]
    project <- base$use_project_effects *
      draws$sigma_project[draw_id] *
      (draws$project_raw[draw_id, p, ] -
         mean(draws$project_raw[draw_id, p, ]))
    prob <- softmax(log(q[c, t, ]) + project)
    alpha <- draws$genome_concentration[draw_id] * prob
    dm_prob <- rgamma(L, shape = alpha, rate = 1)
    dm_prob <- dm_prob / sum(dm_prob)
    n <- sum(base$y_genome[j, ])
    y_genome[j, ] <- as.vector(rmultinom(1, size = n, prob = dm_prob))
  }

  simulated <- base
  simulated$cases <- unname(cases)
  simulated$y_genome <- unname(y_genome)
  truth <- list(
    theta = theta,
    import_scale = draws$import_scale[draw_id, ]
  )
  list(data = simulated, truth = truth)
}

set.seed(seed)
available_draws <- dim(draws$r_coef)[1]
draw_ids <- sample(seq_len(available_draws), n_rep, replace = FALSE)

fit_one <- function(rep_id) {
  rep_seed <- seed + rep_id * 100L
  sim <- simulate_data(draw_ids[rep_id], rep_seed)
  rep_dir <- file.path(outdir, sprintf("replicate_%02d", rep_id))
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(
    sim$data,
    file.path(rep_dir, "simulated_data.json"),
    auto_unbox = TRUE,
    digits = NA
  )
  write_json(
    sim$truth,
    file.path(rep_dir, "truth.json"),
    auto_unbox = TRUE,
    digits = NA,
    pretty = TRUE
  )
  rep_fit <- sampling(
    model,
    data = sim$data,
    chains = chains,
    cores = chains,
    iter = iter,
    warmup = floor(iter / 2),
    seed = rep_seed,
    refresh = 0,
    control = list(adapt_delta = 0.95, max_treedepth = 13)
  )
  pars <- c("lineage_relative_transmission", "import_scale")
  sm <- as.data.table(
    summary(rep_fit, pars = pars, probs = c(0.025, 0.5, 0.975))$summary,
    keep.rownames = "parameter"
  )
  sm[, replicate := rep_id]
  sm[, truth := c(sim$truth$theta, sim$truth$import_scale)]
  sm[, covered_95 := truth >= `2.5%` & truth <= `97.5%`]
  sm[, log_error_median := log(`50%` / truth)]
  sampler <- get_sampler_params(rep_fit, inc_warmup = FALSE)
  diagnostics <- data.table(
    replicate = rep_id,
    divergent_transitions = sum(vapply(
      sampler, function(x) sum(x[, "divergent__"]), numeric(1)
    )),
    maximum_treedepth_hits = sum(vapply(
      sampler, function(x) sum(x[, "treedepth__"] >= 13), numeric(1)
    )),
    maximum_rhat = max(sm$Rhat, na.rm = TRUE),
    minimum_neff = min(sm$n_eff, na.rm = TRUE),
    true_highest_lineage = which.max(sim$truth$theta),
    recovered_highest_lineage = which.max(sm[1:base$L, `50%`])
  )
  fwrite(sm, file.path(rep_dir, "recovery.tsv"), sep = "\t")
  fwrite(diagnostics, file.path(rep_dir, "diagnostics.tsv"), sep = "\t")
  list(recovery = sm, diagnostics = diagnostics)
}

options(mc.cores = min(n_rep * chains, parallel::detectCores()))
results <- parallel::mclapply(
  seq_len(n_rep),
  fit_one,
  mc.cores = min(3L, n_rep),
  mc.preschedule = FALSE
)
if (any(vapply(results, inherits, logical(1), what = "try-error"))) {
  stop("At least one recovery replicate failed")
}

recovery <- rbindlist(lapply(results, `[[`, "recovery"))
diagnostics <- rbindlist(lapply(results, `[[`, "diagnostics"))
recovery[, parameter_type := fifelse(
  grepl("^lineage_relative_transmission", parameter), "lineage_growth",
  "import_scale"
)]
summary_table <- recovery[, .(
  n_parameters = .N,
  coverage_95 = mean(covered_95),
  median_absolute_log_error = median(abs(log_error_median)),
  correlation_truth_posterior_median = cor(truth, `50%`)
), by = parameter_type]

gate <- list(
  thresholds_defined_before_recovery_results = list(
    minimum_95_interval_coverage = 0.80,
    maximum_median_absolute_log_error_lineage_growth = 0.20,
    maximum_median_absolute_log_error_import_scale = 0.50,
    minimum_truth_median_correlation = 0.70,
    minimum_highest_lineage_rank_recovery = 0.75
  ),
  observed = list(
    lineage_growth_coverage = summary_table[
      parameter_type == "lineage_growth", coverage_95
    ],
    import_scale_coverage = summary_table[
      parameter_type == "import_scale", coverage_95
    ],
    lineage_growth_median_absolute_log_error = summary_table[
      parameter_type == "lineage_growth", median_absolute_log_error
    ],
    import_scale_median_absolute_log_error = summary_table[
      parameter_type == "import_scale", median_absolute_log_error
    ],
    lineage_growth_truth_median_correlation = summary_table[
      parameter_type == "lineage_growth", correlation_truth_posterior_median
    ],
    import_scale_truth_median_correlation = summary_table[
      parameter_type == "import_scale", correlation_truth_posterior_median
    ],
    highest_lineage_rank_recovery = mean(
      diagnostics$true_highest_lineage == diagnostics$recovered_highest_lineage
    )
  )
)
o <- gate$observed
t <- gate$thresholds_defined_before_recovery_results
gate$pass <- (
  o$lineage_growth_coverage >= t$minimum_95_interval_coverage &&
  o$import_scale_coverage >= t$minimum_95_interval_coverage &&
  o$lineage_growth_median_absolute_log_error <=
    t$maximum_median_absolute_log_error_lineage_growth &&
  o$import_scale_median_absolute_log_error <=
    t$maximum_median_absolute_log_error_import_scale &&
  o$lineage_growth_truth_median_correlation >=
    t$minimum_truth_median_correlation &&
  o$import_scale_truth_median_correlation >=
    t$minimum_truth_median_correlation &&
  o$highest_lineage_rank_recovery >= t$minimum_highest_lineage_rank_recovery &&
  sum(diagnostics$divergent_transitions) == 0
)

fwrite(recovery, file.path(outdir, "all_recovery.tsv"), sep = "\t")
fwrite(diagnostics, file.path(outdir, "all_diagnostics.tsv"), sep = "\t")
fwrite(summary_table, file.path(outdir, "recovery_summary.tsv"), sep = "\t")
write_json(gate, file.path(outdir, "identifiability_gate.json"),
           pretty = TRUE, auto_unbox = TRUE)
message(sprintf("Identifiability gate pass: %s", gate$pass))
