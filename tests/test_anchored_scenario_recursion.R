source("scripts/lib/gtd_anchored_scenarios.R")

n_draws <- 3L
n_countries <- 1L
n_months <- 3L
n_lineages <- 2L
n_basis <- 1L
months <- as.Date(c("2022-12-01", "2023-01-01", "2023-02-01"))

model_data <- list(
  C = n_countries,
  T = n_months,
  L = n_lineages,
  K = n_basis,
  cases = matrix(c(10, 20, 30), nrow = 1L),
  B = matrix(c(0.2, 0.4, 0.6), ncol = 1L),
  reporting_change = matrix(c(0, 1, 1), nrow = 1L),
  import_exposure = array(
    c(0, 0, 2, 3, 1, 4),
    dim = c(n_countries, n_months, n_lineages)
  )
)

draws <- list(
  r_coef = array(c(0.1, -0.2, 0.3), dim = c(n_draws, 1L, 1L)),
  log_theta = log(matrix(
    c(0.8, 1.5, 1.0, 1.8, 1.2, 2.0),
    nrow = n_draws,
    byrow = TRUE
  )),
  import_scale = matrix(c(0.5, 1.0, 1.5), ncol = 1L),
  density_feedback = matrix(c(0.1, 0.2, 0.3), ncol = 1L),
  reporting_jump = matrix(c(-0.2, 0.0, 0.2), ncol = 1L),
  q = array(NA_real_, dim = c(n_draws, 1L, n_months, n_lineages)),
  mu_cases = array(NA_real_, dim = c(n_draws, 1L, n_months))
)
for (time_index in seq_len(n_months)) {
  draws$q[, 1L, time_index, ] <- matrix(
    rep(c(0.4, 0.6), each = n_draws),
    nrow = n_draws
  )
}
draws$mu_cases[, 1L, 1L] <- model_data$cases[1L, 1L] + 1e-6

# Independent transcription of the Stan transformed-parameter expression for
# the first month after the anchor.
previous_cases <- rep(model_data$cases[1L, 1L], n_draws)
previous_q <- draws$q[, 1L, 1L, ]
log_r <- draws$r_coef[, 1L, 1L] * model_data$B[2L, 1L] -
  draws$density_feedback[, 1L] * log1p(previous_cases / 1000)
manual_component <- matrix(NA_real_, nrow = n_draws, ncol = n_lineages)
for (lineage_index in seq_len(n_lineages)) {
  manual_component[, lineage_index] <-
    exp(log_r + draws$log_theta[, lineage_index]) *
      previous_q[, lineage_index] * (previous_cases + 0.5) +
    draws$import_scale[, 1L] *
      model_data$import_exposure[1L, 2L, lineage_index] +
    1e-9
}
draws$mu_cases[, 1L, 2L] <- rowSums(manual_component) *
  exp(draws$reporting_jump[, 1L])
draws$mu_cases[, 1L, 3L] <- 1

baseline_step <- gtd_scenario_step(
  previous_cases,
  previous_q,
  2L,
  1L,
  draws,
  model_data
)
stopifnot(isTRUE(all.equal(
  baseline_step$cases,
  draws$mu_cases[, 1L, 2L],
  tolerance = 1e-12
)))

# Removing exposure changes only the import component at an otherwise
# identical one-step state.
no_exposure_step <- gtd_scenario_step(
  previous_cases,
  previous_q,
  2L,
  1L,
  draws,
  model_data,
  remove_import = TRUE
)
stopifnot(isTRUE(all.equal(
  no_exposure_step$local_component,
  baseline_step$local_component,
  tolerance = 0
)))
stopifnot(all(no_exposure_step$import_component == 0))

# Standardising L1_02.07 (represented by column 2 here) changes only that
# lineage's local multiplier; exposure and every other lineage are untouched.
neutral_step <- gtd_scenario_step(
  previous_cases,
  previous_q,
  2L,
  1L,
  draws,
  model_data,
  neutral_lineage_index = 2L
)
stopifnot(isTRUE(all.equal(
  neutral_step$import_component,
  baseline_step$import_component,
  tolerance = 0
)))
stopifnot(isTRUE(all.equal(
  neutral_step$local_component[, 1L],
  baseline_step$local_component[, 1L],
  tolerance = 0
)))
stopifnot(isTRUE(all.equal(
  neutral_step$local_component[, 2L],
  baseline_step$local_component[, 2L] /
    exp(draws$log_theta[, 2L]),
  tolerance = 1e-12
)))

# With zero exposure in a country, the no-exposure recursion is draw-wise
# identical to baseline throughout.
zero_exposure_data <- model_data
zero_exposure_data$import_exposure[] <- 0
zero_exposure_draws <- draws
first_zero_step <- gtd_scenario_step(
  previous_cases,
  previous_q,
  2L,
  1L,
  zero_exposure_draws,
  zero_exposure_data
)
zero_exposure_draws$mu_cases[, 1L, 2L] <- first_zero_step$cases
zero_exposure_run <- gtd_run_anchored_scenarios(
  zero_exposure_draws,
  zero_exposure_data,
  months,
  1L,
  2L
)
stopifnot(identical(
  zero_exposure_run$cases[["baseline"]],
  zero_exposure_run$cases[["no_post_2022_exposure"]]
))

# Cumulative proportions are paired within draws. They are not a ratio of
# marginal summaries, and no clipping is applied.
baseline_cases <- rbind(c(1, 2, 8), c(1, 8, 12))
scenario_cases <- rbind(c(1, 1, 4), c(1, 6, 9))
contrast <- gtd_drawwise_cumulative_difference(
  baseline_cases,
  scenario_cases,
  2:3
)
stopifnot(identical(contrast$difference_fraction, c(0.5, 0.25)))
stopifnot(
  median(contrast$difference_fraction) !=
    1 - median(contrast$scenario_cumulative) /
      median(contrast$baseline_cumulative)
)

unclipped <- gtd_drawwise_cumulative_difference(
  baseline_cases,
  baseline_cases * 2,
  2:3
)
stopifnot(identical(unclipped$difference_fraction, c(-1, -1)))
