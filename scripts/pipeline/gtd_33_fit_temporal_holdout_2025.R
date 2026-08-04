#!/usr/bin/env Rscript

# Fit the frozen joint model through December 2024, then recursively forecast
# the 12 months of 2025 without using 2025 cases, genomes, or exposure in the
# fit. Projection uses rows 73-84 of the original pre-frozen spline basis.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L) {
  stop(
    paste(
      "Usage: gtd_33_fit_temporal_holdout_2025.R COMPILED_MODEL_RDS",
      "TRAINING_DATA_JSON FULL_DATA_JSON PREFIT_AUDIT_JSON OUTDIR",
      "POSTERIOR_CACHE_RDS [CHAINS] [ITER] [SEED]"
    )
  )
}

compiled_model_file <- normalizePath(args[[1L]], mustWork = TRUE)
training_data_file <- normalizePath(args[[2L]], mustWork = TRUE)
full_data_file <- normalizePath(args[[3L]], mustWork = TRUE)
prefit_audit_file <- normalizePath(args[[4L]], mustWork = TRUE)
outdir <- args[[5L]]
posterior_cache_file <- args[[6L]]
chains <- if (length(args) >= 7L) as.integer(args[[7L]]) else 4L
iter <- if (length(args) >= 8L) as.integer(args[[8L]]) else 2000L
seed <- if (length(args) >= 9L) as.integer(args[[9L]]) else 20260821L
warmup <- floor(iter / 2)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(
  dirname(posterior_cache_file),
  recursive = TRUE,
  showWarnings = FALSE
)

countries <- c("AUS", "CHN", "JPN")
lineages <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
training_data <- fromJSON(training_data_file, simplifyVector = TRUE)
if (training_data$T != 72L ||
    training_data$C != length(countries) ||
    training_data$L != length(lineages)) {
  stop("The temporal-holdout training dimensions are not the frozen 2019-2024 design")
}
if (any(training_data$obs_month > training_data$T)) {
  stop("A post-cutoff genome observation leaked into the training data")
}

training_md5 <- unname(tools::md5sum(training_data_file))
compiled_model_md5 <- unname(tools::md5sum(compiled_model_file))
expected_draws <- chains * (iter - warmup)

summarise_fit <- function(fit) {
  sampler <- get_sampler_params(fit, inc_warmup = FALSE)
  monitored <- c(
    "r_coef", "log_theta_raw", "import_scale", "density_feedback",
    "q0", "project_raw", "sigma_project", "reporting_jump",
    "phi_cases", "genome_concentration"
  )
  fit_summary <- summary(fit, pars = monitored)$summary
  list(
    chains = chains,
    iterations_per_chain = iter,
    warmup_per_chain = warmup,
    post_warmup_draws = expected_draws,
    seed = seed,
    divergent_transitions = sum(vapply(
      sampler,
      function(x) sum(x[, "divergent__"]),
      numeric(1L)
    )),
    maximum_treedepth_hits = sum(vapply(
      sampler,
      function(x) sum(x[, "treedepth__"] >= 13),
      numeric(1L)
    )),
    maximum_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
    minimum_neff = min(fit_summary[, "n_eff"], na.rm = TRUE)
  )
}

if (file.exists(posterior_cache_file)) {
  message("Reading validated temporal-holdout posterior cache")
  cache <- readRDS(posterior_cache_file)
  if (!identical(cache$training_md5, training_md5) ||
      !identical(cache$compiled_model_md5, compiled_model_md5)) {
    stop("The posterior cache does not match the current model and training data")
  }
  if (cache$fit_config$chains != chains ||
      cache$fit_config$iter != iter ||
      cache$fit_config$seed != seed) {
    stop("The posterior cache was produced with a different fit configuration")
  }
} else {
  message("Reading compiled Stan model")
  model <- readRDS(compiled_model_file)
  options(mc.cores = min(chains, parallel::detectCores()))
  rstan_options(auto_write = TRUE)

  message(sprintf(
    "Sampling leakage-free training fit: %d chains x %d iterations",
    chains,
    iter
  ))
  fit <- sampling(
    model,
    data = training_data,
    chains = chains,
    iter = iter,
    warmup = warmup,
    seed = seed,
    refresh = max(1L, floor(iter / 20)),
    control = list(adapt_delta = 0.95, max_treedepth = 13)
  )
  diagnostics <- summarise_fit(fit)
  key_summary <- as.data.table(
    summary(
      fit,
      pars = c(
        "lineage_relative_transmission", "import_scale",
        "density_feedback", "phi_cases", "sigma_project",
        "genome_concentration"
      ),
      probs = c(0.025, 0.5, 0.975)
    )$summary,
    keep.rownames = "parameter"
  )
  posterior_draws <- extract(
    fit,
    pars = c(
      "r_coef", "log_theta", "import_scale", "density_feedback",
      "reporting_jump", "phi_cases", "q",
      "lineage_relative_transmission"
    ),
    permuted = TRUE
  )
  cache <- list(
    cache_version = 1L,
    training_md5 = training_md5,
    compiled_model_md5 = compiled_model_md5,
    fit_config = list(chains = chains, iter = iter, seed = seed),
    diagnostics = diagnostics,
    key_summary = key_summary,
    draws = posterior_draws
  )
  message("Saving compact posterior cache")
  saveRDS(cache, posterior_cache_file, compress = "xz")
  rm(fit, model, posterior_draws)
  invisible(gc())
}

draws <- cache$draws
n_draws <- dim(draws$q)[1L]
if (n_draws != expected_draws) {
  stop("The posterior draw count does not match the requested fit")
}
fwrite(
  as.data.table(cache$key_summary),
  file.path(outdir, "temporal_holdout_2025_posterior_parameters.tsv"),
  sep = "\t"
)
write_json(
  cache$diagnostics,
  file.path(outdir, "temporal_holdout_2025_sampling_diagnostics.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

# The full data, including observed holdout outcomes, are intentionally loaded
# only after fitting (or after validating a cache fitted to the truncated JSON).
full_data <- fromJSON(full_data_file, simplifyVector = TRUE)
if (full_data$T != 84L ||
    full_data$C != training_data$C ||
    full_data$L != training_data$L ||
    full_data$K != training_data$K) {
  stop("Full and training model dimensions are incompatible")
}
if (!isTRUE(all.equal(
  training_data$cases,
  full_data$cases[, seq_len(training_data$T), drop = FALSE],
  tolerance = 0
))) {
  stop("Training cases are not the exact through-2024 prefix")
}
if (!isTRUE(all.equal(
  training_data$B,
  full_data$B[seq_len(training_data$T), , drop = FALSE],
  tolerance = 0
))) {
  stop("Training spline design is not the exact pre-frozen prefix")
}

all_months <- seq(
  as.Date("2019-01-01"),
  by = "month",
  length.out = full_data$T
)
forecast_indices <- seq.int(training_data$T + 1L, full_data$T)
forecast_months <- all_months[forecast_indices]
horizon <- length(forecast_indices)
if (horizon != 12L ||
    format(forecast_months[[1L]], "%Y-%m") != "2025-01" ||
    format(forecast_months[[horizon]], "%Y-%m") != "2025-12") {
  stop("The requested holdout is not January-December 2025")
}

draw_country_matrix <- function(x, country_index) {
  out <- x[, country_index, , drop = FALSE]
  dim(out) <- c(dim(x)[1L], dim(x)[3L])
  out
}

draw_lineage_matrix <- function(x, country_index, time_index) {
  out <- x[, country_index, time_index, , drop = FALSE]
  dim(out) <- c(dim(x)[1L], dim(x)[4L])
  out
}

theta <- exp(draws$log_theta)
expected_cases <- array(
  NA_real_,
  dim = c(n_draws, full_data$C, horizon)
)
predictive_cases <- array(
  NA_real_,
  dim = c(n_draws, full_data$C, horizon)
)
predictive_conditional_mean <- array(
  NA_real_,
  dim = c(n_draws, full_data$C, horizon)
)

set.seed(seed + 1L)
for (country_index in seq_len(full_data$C)) {
  r_coef <- draw_country_matrix(draws$r_coef, country_index)
  q_expected <- draw_lineage_matrix(
    draws$q,
    country_index,
    training_data$T
  )
  q_predictive <- q_expected
  previous_expected <- rep(
    training_data$cases[country_index, training_data$T],
    n_draws
  )
  previous_predictive <- previous_expected

  for (horizon_index in seq_len(horizon)) {
    time_index <- forecast_indices[[horizon_index]]
    spline_term <- as.vector(r_coef %*% full_data$B[time_index, ])
    reporting_delta <-
      full_data$reporting_change[country_index, time_index] -
      full_data$reporting_change[country_index, time_index - 1L]
    reporting_factor <- exp(
      draws$reporting_jump[, country_index] * reporting_delta
    )

    log_r_expected <- spline_term -
      draws$density_feedback[, country_index] *
        log1p(previous_expected / 1000)
    local_expected <- sweep(
      theta * q_expected,
      1L,
      exp(log_r_expected) * (previous_expected + 0.5),
      "*"
    )
    # Exposure is intentionally and explicitly zero for every 2025 month.
    component_expected <- local_expected + 1e-9
    total_expected <- rowSums(component_expected)
    current_expected <- total_expected * reporting_factor
    q_expected <- component_expected / total_expected

    log_r_predictive <- spline_term -
      draws$density_feedback[, country_index] *
        log1p(previous_predictive / 1000)
    local_predictive <- sweep(
      theta * q_predictive,
      1L,
      exp(log_r_predictive) * (previous_predictive + 0.5),
      "*"
    )
    component_predictive <- local_predictive + 1e-9
    total_predictive <- rowSums(component_predictive)
    current_predictive_mean <- total_predictive * reporting_factor
    current_predictive <- rnbinom(
      n_draws,
      size = draws$phi_cases[, country_index],
      mu = current_predictive_mean
    )
    q_predictive <- component_predictive / total_predictive

    if (any(!is.finite(current_expected)) ||
        any(!is.finite(current_predictive_mean)) ||
        any(!is.finite(current_predictive))) {
      stop("A non-finite 2025 forecast was generated")
    }
    expected_cases[, country_index, horizon_index] <- current_expected
    predictive_conditional_mean[, country_index, horizon_index] <-
      current_predictive_mean
    predictive_cases[, country_index, horizon_index] <- current_predictive
    previous_expected <- current_expected
    previous_predictive <- current_predictive
  }
}

# Seasonal-naive comparator:
#   point forecast = the observed count in the same month of 2024;
#   predictive uncertainty = bootstrap of median-centred, one-year log1p
#   differences from all 2020-2024 training months within that country.
seasonal_naive_cases <- array(
  NA_real_,
  dim = c(n_draws, full_data$C, horizon)
)
seasonal_naive_point <- matrix(
  NA_real_,
  nrow = full_data$C,
  ncol = horizon
)
baseline_calibration <- list()
set.seed(seed + 2L)
for (country_index in seq_len(full_data$C)) {
  training_cases <- training_data$cases[country_index, ]
  residuals <- log1p(training_cases[13:training_data$T]) -
    log1p(training_cases[1:(training_data$T - 12L)])
  residual_median <- median(residuals)
  centred_residuals <- residuals - residual_median
  residual_sample <- matrix(
    sample(
      centred_residuals,
      n_draws * horizon,
      replace = TRUE
    ),
    nrow = n_draws,
    ncol = horizon
  )
  same_month_2024 <- training_cases[
    (training_data$T - horizon + 1L):training_data$T
  ]
  seasonal_naive_cases[, country_index, ] <- pmax(
    0,
    exp(sweep(
      residual_sample,
      2L,
      log1p(same_month_2024),
      "+"
    )) - 1
  )
  seasonal_naive_point[country_index, ] <- same_month_2024
  baseline_calibration[[country_index]] <- data.table(
    country_iso3 = countries[[country_index]],
    residual_definition =
      "log1p(case_t)-log1p(case_t_minus_12), median-centred",
    residual_months = length(residuals),
    raw_residual_median = residual_median,
    raw_residual_mad = mad(residuals),
    bootstrap_draws = n_draws
  )
}
fwrite(
  rbindlist(baseline_calibration),
  file.path(outdir, "temporal_holdout_2025_seasonal_naive_calibration.tsv"),
  sep = "\t"
)

quantile_probabilities <- c(0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.975)
quantile_columns <- c(
  "lower_95", "lower_90", "lower_80", "lower_50", "median",
  "upper_50", "upper_80", "upper_90", "upper_95"
)
forecast_methods <- c(
  "joint_model_expected",
  "joint_model_negative_binomial",
  "seasonal_naive"
)
uncertainty_scope <- c(
  joint_model_expected = "posterior expected-case uncertainty",
  joint_model_negative_binomial =
    "recursive negative-binomial posterior predictive uncertainty",
  seasonal_naive =
    "empirical one-year log-difference predictive uncertainty"
)
forecast_arrays <- list(
  joint_model_expected = expected_cases,
  joint_model_negative_binomial = predictive_cases,
  seasonal_naive = seasonal_naive_cases
)

monthly_rows <- list()
row_index <- 0L
for (method in forecast_methods) {
  values <- forecast_arrays[[method]]
  for (country_index in seq_len(full_data$C)) {
    for (horizon_index in seq_len(horizon)) {
      forecast_quantiles <- quantile(
        values[, country_index, horizon_index],
        quantile_probabilities,
        names = FALSE
      )
      point_forecast <- if (method == "seasonal_naive") {
        seasonal_naive_point[country_index, horizon_index]
      } else {
        forecast_quantiles[[5L]]
      }
      row_index <- row_index + 1L
      monthly_rows[[row_index]] <- data.table(
        forecast_method = method,
        uncertainty_scope = uncertainty_scope[[method]],
        country_iso3 = countries[[country_index]],
        forecast_month = format(
          forecast_months[[horizon_index]],
          "%Y-%m"
        ),
        horizon_month = horizon_index,
        observed_cases =
          full_data$cases[country_index, forecast_indices[[horizon_index]]],
        point_forecast = point_forecast,
        lower_95 = forecast_quantiles[[1L]],
        lower_90 = forecast_quantiles[[2L]],
        lower_80 = forecast_quantiles[[3L]],
        lower_50 = forecast_quantiles[[4L]],
        median = forecast_quantiles[[5L]],
        upper_50 = forecast_quantiles[[6L]],
        upper_80 = forecast_quantiles[[7L]],
        upper_90 = forecast_quantiles[[8L]],
        upper_95 = forecast_quantiles[[9L]]
      )
    }
  }
}
monthly <- rbindlist(monthly_rows)
fwrite(
  monthly,
  file.path(outdir, "temporal_holdout_2025_monthly_forecasts.tsv"),
  sep = "\t"
)

interval_score <- function(observed, lower, upper, alpha) {
  (upper - lower) +
    (2 / alpha) * (lower - observed) * (observed < lower) +
    (2 / alpha) * (observed - upper) * (observed > upper)
}

wis_from_row <- function(row) {
  alphas <- c(0.50, 0.20, 0.10, 0.05)
  lowers <- c(row$lower_50, row$lower_80, row$lower_90, row$lower_95)
  uppers <- c(row$upper_50, row$upper_80, row$upper_90, row$upper_95)
  interval_scores <- vapply(
    seq_along(alphas),
    function(index) {
      interval_score(
        row$observed_cases,
        lowers[[index]],
        uppers[[index]],
        alphas[[index]]
      )
    },
    numeric(1L)
  )
  (
    0.5 * abs(row$observed_cases - row$point_forecast) +
      sum(alphas / 2 * interval_scores)
  ) / (length(alphas) + 0.5)
}

monthly[, wis := vapply(
  seq_len(.N),
  function(index) wis_from_row(monthly[index]),
  numeric(1L)
)]
metrics_country <- monthly[, .(
  n_months = .N,
  log_rmse = sqrt(mean(
    (log1p(observed_cases) - log1p(point_forecast))^2
  )),
  mean_absolute_error = mean(abs(observed_cases - point_forecast)),
  coverage_95 = mean(
    observed_cases >= lower_95 & observed_cases <= upper_95
  ),
  median_95_interval_width = median(upper_95 - lower_95),
  mean_wis = mean(wis)
), by = .(forecast_method, country_iso3)]
metrics_all <- monthly[, .(
  n_months = .N,
  log_rmse = sqrt(mean(
    (log1p(observed_cases) - log1p(point_forecast))^2
  )),
  mean_absolute_error = mean(abs(observed_cases - point_forecast)),
  coverage_95 = mean(
    observed_cases >= lower_95 & observed_cases <= upper_95
  ),
  median_95_interval_width = median(upper_95 - lower_95),
  mean_wis = mean(wis)
), by = forecast_method]
metrics_all[, country_iso3 := "ALL"]
metrics <- rbindlist(list(metrics_country, metrics_all), use.names = TRUE)
baseline_metrics <- metrics[
  forecast_method == "seasonal_naive",
  .(
    country_iso3,
    seasonal_naive_log_rmse = log_rmse,
    seasonal_naive_mean_wis = mean_wis
  )
]
metrics[
  baseline_metrics,
  on = "country_iso3",
  `:=`(
    log_rmse_ratio_to_seasonal_naive =
      log_rmse / i.seasonal_naive_log_rmse,
    mean_wis_ratio_to_seasonal_naive =
      mean_wis / i.seasonal_naive_mean_wis
  )
]
metrics[, wis_interval_levels := "50%;80%;90%;95%"]
setorder(metrics, country_iso3, forecast_method)
fwrite(
  metrics,
  file.path(outdir, "temporal_holdout_2025_metrics.tsv"),
  sep = "\t"
)

prefit_audit <- fromJSON(prefit_audit_file, simplifyVector = FALSE)
prefit_audit$postfit_prediction_audit <- list(
  posterior_cache = normalizePath(posterior_cache_file, mustWork = TRUE),
  posterior_draws = n_draws,
  fit_chains = chains,
  fit_iterations_per_chain = iter,
  fit_seed = seed,
  full_outcome_data_loaded_only_after_fit_or_cache_validation = TRUE,
  prediction_anchor_month = "2024-12",
  prediction_anchor_case_source = (
    "Observed December 2024 cases contained in the training JSON"
  ),
  prediction_anchor_lineage_source = (
    "Posterior q in December 2024 from the truncated training fit"
  ),
  full_prefrozen_spline_projection_rows = "73-84",
  spline_basis_rebuilt_after_truncation = FALSE,
  forecast_exposure_forced_zero = TRUE,
  forecast_exposure_weight_used = setNames(
    as.list(rep(0, length(countries))),
    countries
  ),
  import_scale_used_in_forecast = FALSE,
  holdout_cases_used_only_for_scoring = TRUE,
  holdout_genomes_used_for_fit = FALSE,
  recursive_expected_forecast = (
    "Each month uses the previous month's posterior expected case path"
  ),
  recursive_predictive_forecast = paste(
    "Each month uses the previous simulated count and draws a new count",
    "from negative-binomial-2"
  ),
  seasonal_naive_definition = paste(
    "Point forecast is the same calendar month of 2024; intervals bootstrap",
    "median-centred 12-month log1p differences from 2020-2024"
  ),
  wis_definition = (
    "Weighted interval score using central 50%, 80%, 90%, and 95% intervals"
  )
)
write_json(
  prefit_audit,
  file.path(outdir, "temporal_holdout_2025_leakage_audit.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA,
  null = "null"
)
