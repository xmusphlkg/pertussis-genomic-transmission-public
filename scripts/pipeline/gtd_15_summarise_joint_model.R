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

counterfactual <- rbindlist(lapply(seq_along(countries), function(c) {
  base <- rowSums(draws$expected_cases_baseline[, c, post_idx, drop = FALSE])
  no_import <- rowSums(draws$expected_cases_no_import[, c, post_idx, drop = FALSE])
  equal <- rowSums(draws$expected_cases_equal_lineage[, c, post_idx, drop = FALSE])
  # Report a bounded conditional fractional reduction, not an attributable fraction.
  # Recursive nonlinear counterfactuals can otherwise give tiny negative
  # Monte Carlo values when the no-import trajectory slightly exceeds baseline.
  import_attr <- pmax(0, pmin(1, 1 - no_import / base))
  lineage_attr <- 1 - equal / base
  rbind(
    data.table(country_iso3 = countries[[c]], scenario = "baseline",
               t(qsum(base)), change_fraction_median = NA_real_),
    data.table(country_iso3 = countries[[c]], scenario = "no_import",
               t(qsum(no_import)),
               change_fraction_median = median(no_import / base - 1)),
    data.table(country_iso3 = countries[[c]], scenario = "equal_lineage_transmission",
               t(qsum(equal)),
               change_fraction_median = median(equal / base - 1)),
    data.table(country_iso3 = countries[[c]], scenario = "no_new_introduction_case_reduction_fraction",
               t(qsum(import_attr)), change_fraction_median = NA_real_),
    data.table(country_iso3 = countries[[c]], scenario = "lineage_difference_effect_fraction",
               t(qsum(lineage_attr)), change_fraction_median = NA_real_)
  )
}))
fwrite(counterfactual, file.path(outdir, "counterfactual_summary.tsv"), sep = "\t")

# Sampling-corrected lineage shares are weighted by the observed national case
# series, while raw shares are the unadjusted project archive proportions.
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

# Monthly latent lineage shares for trajectory panels. These are posterior
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

# Monthly source data for model fit and counterfactual figures.
monthly <- list()
for (c in seq_along(countries)) {
  for (t in seq_along(months)) {
    base <- qsum(draws$expected_cases_baseline[, c, t])
    no_imp <- qsum(draws$expected_cases_no_import[, c, t])
    equal <- qsum(draws$expected_cases_equal_lineage[, c, t])
    mu <- qsum(draws$mu_cases[, c, t])
    monthly[[length(monthly) + 1L]] <- data.table(
      country_iso3 = countries[[c]],
      model_month = months[[t]],
      observed_cases = case_array[c, t],
      fitted_median = mu[["median"]],
      fitted_lower_95 = mu[["lower_95"]],
      fitted_upper_95 = mu[["upper_95"]],
      baseline_median = base[["median"]],
      no_import_median = no_imp[["median"]],
      equal_lineage_median = equal[["median"]]
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

# A beta-binomial descriptive interval quantifies the probability that a
# high-support introduction becomes a sampled cluster of >=5 genomes spanning
# >=6 months. This is not called a transmission probability in individuals.
event_summary <- events[high_support_post_reseeding == TRUE, .(
  n_high_support_events = .N,
  n_successful_sampled_clusters = sum(successful_sampled_cluster)
), by = destination_country]
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
