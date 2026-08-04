#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: gtd_15_summarise_joint_model.R POSTERIOR_RDS DATA_JSON OBS_TSV EVENTS_TSV OUTDIR")
}
draws <- readRDS(args[[1L]])
model_data <- fromJSON(args[[2L]], simplifyVector = TRUE)
observations <- fread(args[[3L]])
events <- fread(args[[4L]])
outdir <- args[[5L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

countries <- c("AUS", "CHN", "JPN")
lineages <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
months <- seq(as.Date("2019-01-01"), as.Date("2025-12-01"), by = "month")
post_idx <- which(months >= as.Date("2023-01-01"))
anchor_idx <- which(months == as.Date("2022-12-01"))

qsum <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    lower_95 = unname(quantile(x, 0.025)),
    upper_95 = unname(quantile(x, 0.975))
  )
}

rr <- rbindlist(lapply(seq_along(lineages), function(l) {
  z <- qsum(draws$lineage_relative_transmission[, l])
  data.table(lineage = lineages[[l]], t(z))
}))
fwrite(rr, file.path(outdir, "lineage_relative_transmission.tsv"), sep = "\t")

l10207_idx <- match("L1_02.07", lineages)
pairwise_growth <- rbindlist(lapply(
  setdiff(seq_along(lineages), l10207_idx),
  function(l) {
    ratio <- draws$lineage_relative_transmission[, l10207_idx] /
      draws$lineage_relative_transmission[, l]
    data.table(
      numerator = "L1_02.07",
      denominator = lineages[[l]],
      t(qsum(ratio)),
      posterior_probability_above_one = mean(ratio > 1)
    )
  }
))
fwrite(
  pairwise_growth,
  file.path(outdir, "l10207_pairwise_growth.tsv"),
  sep = "\t"
)

required_scenario_draws <- c(
  "r_coef", "log_theta", "import_scale", "density_feedback",
  "reporting_jump", "q"
)
if (!all(required_scenario_draws %in% names(draws))) {
  stop(
    "Posterior output is missing parameters required for anchored scenarios: ",
    paste(setdiff(required_scenario_draws, names(draws)), collapse = ", ")
  )
}

# Conditional scenarios are anchored to the fitted lineage state and observed
# reported-case count in December 2022. From January 2023 onward, each trajectory
# evolves recursively. The no-new-exposure scenario changes only post-2022
# import exposure; it does not also reweight the initial lineage state by
# phylogeographic persistence support.
n_draws <- dim(draws$q)[1L]
n_countries <- length(countries)
n_months <- length(months)
n_lineages <- length(lineages)
scenario_baseline <- array(NA_real_, c(n_draws, n_countries, n_months))
scenario_no_new_exposure <- array(NA_real_, c(n_draws, n_countries, n_months))
scenario_equal_lineage <- array(NA_real_, c(n_draws, n_countries, n_months))
scenario_l10207_reference <- array(
  NA_real_,
  c(n_draws, n_countries, n_months)
)

for (c in seq_along(countries)) {
  historical_cases <- matrix(
    rep(model_data$cases[c, seq_len(anchor_idx)], each = n_draws),
    nrow = n_draws
  )
  scenario_baseline[, c, seq_len(anchor_idx)] <- historical_cases
  scenario_no_new_exposure[, c, seq_len(anchor_idx)] <- historical_cases
  scenario_equal_lineage[, c, seq_len(anchor_idx)] <- historical_cases
  scenario_l10207_reference[, c, seq_len(anchor_idx)] <- historical_cases

  q_base <- draws$q[, c, anchor_idx, ]
  q_no_exposure <- q_base
  q_equal <- q_base
  q_l10207_reference <- q_base

  for (t in post_idx) {
    spline_term <- as.vector(draws$r_coef[, c, ] %*% model_data$B[t, ])
    previous_base <- scenario_baseline[, c, t - 1L]
    previous_no_exposure <- scenario_no_new_exposure[, c, t - 1L]
    previous_equal <- scenario_equal_lineage[, c, t - 1L]
    previous_l10207_reference <-
      scenario_l10207_reference[, c, t - 1L]
    log_r_base <- spline_term -
      draws$density_feedback[, c] * log1p(previous_base / 1000)
    log_r_no_exposure <- spline_term -
      draws$density_feedback[, c] * log1p(previous_no_exposure / 1000)
    log_r_equal <- spline_term -
      draws$density_feedback[, c] * log1p(previous_equal / 1000)
    log_r_l10207_reference <- spline_term -
      draws$density_feedback[, c] *
        log1p(previous_l10207_reference / 1000)

    comp_base <- matrix(0, nrow = n_draws, ncol = n_lineages)
    comp_no_exposure <- matrix(0, nrow = n_draws, ncol = n_lineages)
    comp_equal <- matrix(0, nrow = n_draws, ncol = n_lineages)
    comp_l10207_reference <- matrix(
      0,
      nrow = n_draws,
      ncol = n_lineages
    )
    for (l in seq_along(lineages)) {
      theta <- exp(draws$log_theta[, l])
      theta_l10207_reference <- if (l == l10207_idx) 1 else theta
      imported <- draws$import_scale[, c] * model_data$import_exposure[c, t, l]
      comp_base[, l] <- exp(log_r_base) * theta * q_base[, l] *
        (previous_base + 0.5) + imported + 1e-9
      comp_no_exposure[, l] <- exp(log_r_no_exposure) * theta *
        q_no_exposure[, l] * (previous_no_exposure + 0.5) + 1e-9
      comp_equal[, l] <- exp(log_r_equal) * q_equal[, l] *
        (previous_equal + 0.5) + imported + 1e-9
      comp_l10207_reference[, l] <- exp(log_r_l10207_reference) *
        theta_l10207_reference * q_l10207_reference[, l] *
        (previous_l10207_reference + 0.5) + imported + 1e-9
    }

    total_base <- rowSums(comp_base)
    total_no_exposure <- rowSums(comp_no_exposure)
    total_equal <- rowSums(comp_equal)
    total_l10207_reference <- rowSums(comp_l10207_reference)
    reporting_factor <- exp(
      draws$reporting_jump[, c] *
        (model_data$reporting_change[c, t] -
           model_data$reporting_change[c, t - 1L])
    )
    scenario_baseline[, c, t] <- total_base * reporting_factor
    scenario_no_new_exposure[, c, t] <- total_no_exposure * reporting_factor
    scenario_equal_lineage[, c, t] <- total_equal * reporting_factor
    scenario_l10207_reference[, c, t] <-
      total_l10207_reference * reporting_factor
    q_base <- comp_base / total_base
    q_no_exposure <- comp_no_exposure / total_no_exposure
    q_equal <- comp_equal / total_equal
    q_l10207_reference <-
      comp_l10207_reference / total_l10207_reference
  }
}

conditional_scenarios <- rbindlist(lapply(seq_along(countries), function(c) {
  base <- rowSums(scenario_baseline[, c, post_idx, drop = FALSE])
  no_exposure <- rowSums(
    scenario_no_new_exposure[, c, post_idx, drop = FALSE]
  )
  equal <- rowSums(scenario_equal_lineage[, c, post_idx, drop = FALSE])
  l10207_reference <- rowSums(
    scenario_l10207_reference[, c, post_idx, drop = FALSE]
  )
  no_exposure_difference <- 1 - no_exposure / base
  lineage_difference <- 1 - equal / base
  l10207_difference <- 1 - l10207_reference / base
  rbind(
    data.table(country_iso3 = countries[[c]], scenario = "baseline",
               t(qsum(base)), change_fraction_median = NA_real_),
    data.table(country_iso3 = countries[[c]], scenario = "no_new_exposure",
               t(qsum(no_exposure)),
               change_fraction_median = median(no_exposure / base - 1)),
    data.table(country_iso3 = countries[[c]], scenario = "equal_lineage_growth",
               t(qsum(equal)),
               change_fraction_median = median(equal / base - 1)),
    data.table(
      country_iso3 = countries[[c]],
      scenario = "l10207_reference_growth",
      t(qsum(l10207_reference)),
      change_fraction_median = median(l10207_reference / base - 1)
    ),
    data.table(
      country_iso3 = countries[[c]],
      scenario = "no_new_exposure_case_difference_fraction",
      t(qsum(no_exposure_difference)),
      change_fraction_median = NA_real_
    ),
    data.table(
      country_iso3 = countries[[c]],
      scenario = "lineage_growth_scenario_difference_fraction",
      t(qsum(lineage_difference)),
      change_fraction_median = NA_real_
    ),
    data.table(
      country_iso3 = countries[[c]],
      scenario = "l10207_growth_scenario_difference_fraction",
      t(qsum(l10207_difference)),
      change_fraction_median = NA_real_
    )
  )
}))
fwrite(
  conditional_scenarios,
  file.path(outdir, "counterfactual_summary.tsv"),
  sep = "\t",
  na = "NA"
)

# Project-adjusted modelled lineage shares are weighted by the observed national
# case series, while raw shares are the unadjusted project archive proportions.
case_array <- model_data$cases
share_rows <- list()
for (c in seq_along(countries)) {
  raw_country <- observations[country_iso3 == countries[[c]] & month_id %in% post_idx]
  count_cols <- grep("^n_", names(raw_country), value = TRUE)
  raw_counts <- colSums(as.matrix(raw_country[, ..count_cols]))
  raw_share <- raw_counts / sum(raw_counts)
  for (l in seq_along(lineages)) {
    corrected <- rowSums(
      sweep(draws$q[, c, post_idx, l, drop = FALSE][, 1, , 1],
            2, case_array[c, post_idx], "*")
    ) / sum(case_array[c, post_idx])
    z <- qsum(corrected)
    share_rows[[length(share_rows) + 1L]] <- data.table(
      country_iso3 = countries[[c]],
      lineage = lineages[[l]],
      raw_public_tree_share = raw_share[[l]],
      t(z)
    )
  }
}
shares <- rbindlist(share_rows)
fwrite(shares, file.path(outdir, "sampling_corrected_lineage_shares.tsv"), sep = "\t")

# Monthly modelled lineage shares for trajectory panels. These are posterior
# summaries of q and are distinct from raw monthly public-genome proportions.
monthly_share_rows <- list()
for (c in seq_along(countries)) {
  for (t in seq_along(months)) {
    for (l in seq_along(lineages)) {
      z <- qsum(draws$q[, c, t, l])
      monthly_share_rows[[length(monthly_share_rows) + 1L]] <- data.table(
        country_iso3 = countries[[c]],
        model_month = months[[t]],
        lineage = lineages[[l]],
        t(z)
      )
    }
  }
}
monthly_shares <- rbindlist(monthly_share_rows)
fwrite(
  monthly_shares,
  file.path(outdir, "monthly_sampling_corrected_lineage_shares.tsv"),
  sep = "\t"
)

# Monthly source data for model fit and anchored conditional-scenario figures.
monthly <- list()
for (c in seq_along(countries)) {
  for (t in seq_along(months)) {
    base <- qsum(scenario_baseline[, c, t])
    no_exposure <- qsum(scenario_no_new_exposure[, c, t])
    equal <- qsum(scenario_equal_lineage[, c, t])
    l10207_reference <- qsum(scenario_l10207_reference[, c, t])
    mu <- qsum(draws$mu_cases[, c, t])
    monthly[[length(monthly) + 1L]] <- data.table(
      country_iso3 = countries[[c]],
      model_month = months[[t]],
      observed_cases = case_array[c, t],
      fitted_median = mu[["median"]],
      fitted_lower_95 = mu[["lower_95"]],
      fitted_upper_95 = mu[["upper_95"]],
      baseline_median = base[["median"]],
      baseline_lower_95 = base[["lower_95"]],
      baseline_upper_95 = base[["upper_95"]],
      no_new_exposure_median = no_exposure[["median"]],
      no_new_exposure_lower_95 = no_exposure[["lower_95"]],
      no_new_exposure_upper_95 = no_exposure[["upper_95"]],
      equal_lineage_median = equal[["median"]],
      equal_lineage_lower_95 = equal[["lower_95"]],
      equal_lineage_upper_95 = equal[["upper_95"]],
      l10207_reference_median = l10207_reference[["median"]],
      l10207_reference_lower_95 = l10207_reference[["lower_95"]],
      l10207_reference_upper_95 = l10207_reference[["upper_95"]]
    )
  }
}
monthly <- rbindlist(monthly)
fwrite(monthly, file.path(outdir, "monthly_model_and_counterfactuals.tsv"), sep = "\t")

fit_metrics <- monthly[, .(
  log_rmse = sqrt(mean((log1p(observed_cases) - log1p(fitted_median))^2)),
  log_correlation = cor(log1p(observed_cases), log1p(fitted_median)),
  expected_interval_coverage = mean(
    observed_cases >= fitted_lower_95 & observed_cases <= fitted_upper_95
  )
), by = country_iso3]
if (!is.null(draws$cases_rep)) {
  fit_metrics[, posterior_predictive_coverage := vapply(
    seq_along(countries),
    function(c) {
      lo <- apply(draws$cases_rep[, c, , drop = FALSE], 3, quantile, 0.025)
      hi <- apply(draws$cases_rep[, c, , drop = FALSE], 3, quantile, 0.975)
      mean(case_array[c, ] >= lo & case_array[c, ] <= hi)
    },
    numeric(1)
  )]
}
fwrite(fit_metrics, file.path(outdir, "posterior_predictive_metrics.tsv"), sep = "\t")

# Beta-binomial descriptive intervals quantify how often eligible sampled-tree
# edges meet either the conservative minimum-span definition or the
# interval-compatible maximum-span definition. They are not individual-level
# transmission probabilities.
if (!"robust_success" %in% names(events)) {
  events[, robust_success := successful_sampled_cluster]
}
if (!"interval_compatible_success" %in% names(events)) {
  events[, interval_compatible_success := successful_sampled_cluster]
}
cluster_definition_columns <- c(
  robust_minimum_span = "robust_success",
  interval_compatible_maximum_span = "interval_compatible_success"
)
event_summary <- rbindlist(lapply(
  names(cluster_definition_columns),
  function(cluster_definition_name) {
    success_column <- cluster_definition_columns[[cluster_definition_name]]
    events[high_support_post_reseeding == TRUE, .(
      n_high_support_events = .N,
      n_successful_sampled_clusters = sum(get(success_column))
    ), by = destination_country][
      ,
      cluster_definition := cluster_definition_name
    ]
  }
), use.names = TRUE)
event_summary[, `:=`(
  posterior_mean_success_probability =
    (n_successful_sampled_clusters + 1) / (n_high_support_events + 2),
  lower_95 = qbeta(0.025, n_successful_sampled_clusters + 1,
                   n_high_support_events - n_successful_sampled_clusters + 1),
  upper_95 = qbeta(0.975, n_successful_sampled_clusters + 1,
                   n_high_support_events - n_successful_sampled_clusters + 1)
)]
fwrite(event_summary, file.path(outdir, "successful_sampled_cluster_probability.tsv"),
       sep = "\t")
