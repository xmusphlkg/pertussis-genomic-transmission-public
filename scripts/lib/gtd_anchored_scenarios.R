# Shared deterministic recursion for conditional scenarios anchored to an
# observed monthly case count and the fitted lineage state at that month.
#
# The functions in this file deliberately return unbounded draw-wise contrasts.
# In particular, a negative difference fraction is retained when a scenario has
# more cumulative cases than its paired baseline draw.

gtd_scenario_required_draws <- c(
  "r_coef", "log_theta", "import_scale", "density_feedback",
  "reporting_jump", "q", "mu_cases"
)

gtd_validate_scenario_inputs <- function(draws, model_data, months) {
  missing_draws <- setdiff(gtd_scenario_required_draws, names(draws))
  if (length(missing_draws) > 0L) {
    stop(
      "Posterior output is missing scenario parameters: ",
      paste(missing_draws, collapse = ", ")
    )
  }

  required_data <- c(
    "C", "T", "L", "K", "cases", "B", "reporting_change",
    "import_exposure"
  )
  missing_data <- setdiff(required_data, names(model_data))
  if (length(missing_data) > 0L) {
    stop(
      "Model data are missing scenario inputs: ",
      paste(missing_data, collapse = ", ")
    )
  }

  n_draws <- dim(draws$q)[1L]
  expected_q <- c(n_draws, model_data$C, model_data$T, model_data$L)
  if (!identical(as.integer(dim(draws$q)), as.integer(expected_q))) {
    stop("draws$q dimensions do not match the frozen model data")
  }
  if (length(months) != model_data$T) {
    stop("The month vector length does not match model_data$T")
  }
  if (!identical(
    as.integer(dim(draws$mu_cases)),
    as.integer(c(n_draws, model_data$C, model_data$T))
  )) {
    stop("draws$mu_cases dimensions do not match the frozen model data")
  }

  invisible(TRUE)
}

gtd_draw_lineage_matrix <- function(x, country_index, time_index) {
  out <- x[, country_index, time_index, , drop = FALSE]
  dim(out) <- c(dim(x)[1L], dim(x)[4L])
  out
}

gtd_draw_country_matrix <- function(x, country_index) {
  out <- x[, country_index, , drop = FALSE]
  dim(out) <- c(dim(x)[1L], dim(x)[3L])
  out
}

gtd_scenario_step <- function(
    previous_cases,
    previous_q,
    time_index,
    country_index,
    draws,
    model_data,
    remove_import = FALSE,
    neutral_lineage_index = integer()) {
  n_draws <- length(previous_cases)
  n_lineages <- model_data$L
  if (!identical(dim(previous_q), c(n_draws, n_lineages))) {
    stop("previous_q must be an n_draws by n_lineages matrix")
  }
  if (time_index < 2L || time_index > model_data$T) {
    stop("time_index must identify a month after the first model month")
  }
  if (length(neutral_lineage_index) > 1L ||
      any(neutral_lineage_index < 1L | neutral_lineage_index > n_lineages)) {
    stop("neutral_lineage_index must be empty or one valid lineage index")
  }

  r_coef <- gtd_draw_country_matrix(draws$r_coef, country_index)
  spline_term <- as.vector(r_coef %*% model_data$B[time_index, ])
  log_r <- spline_term -
    draws$density_feedback[, country_index] *
      log1p(previous_cases / 1000)

  theta <- exp(draws$log_theta)
  if (length(neutral_lineage_index) == 1L) {
    theta[, neutral_lineage_index] <- 1
  }
  local_multiplier <- exp(log_r) * (previous_cases + 0.5)
  local_component <- sweep(
    theta * previous_q,
    1L,
    local_multiplier,
    "*"
  )

  exposure <- model_data$import_exposure[
    country_index,
    time_index,
    ,
    drop = TRUE
  ]
  if (remove_import) {
    import_component <- matrix(0, nrow = n_draws, ncol = n_lineages)
  } else {
    import_component <- outer(
      draws$import_scale[, country_index],
      exposure,
      "*"
    )
  }

  component <- local_component + import_component + 1e-9
  total_component <- rowSums(component)
  reporting_delta <-
    model_data$reporting_change[country_index, time_index] -
    model_data$reporting_change[country_index, time_index - 1L]
  reporting_factor <- exp(
    draws$reporting_jump[, country_index] * reporting_delta
  )

  list(
    cases = total_component * reporting_factor,
    q = component / total_component,
    local_component = local_component,
    import_component = import_component,
    component = component,
    reporting_factor = reporting_factor
  )
}

gtd_run_anchored_scenarios <- function(
    draws,
    model_data,
    months,
    anchor_index,
    neutral_lineage_index) {
  gtd_validate_scenario_inputs(draws, model_data, months)
  if (anchor_index < 1L || anchor_index >= model_data$T) {
    stop("anchor_index must precede the final model month")
  }

  n_draws <- dim(draws$q)[1L]
  scenario_names <- c(
    "baseline",
    "no_post_2022_exposure",
    "l10207_neutral_growth"
  )
  cases <- setNames(
    lapply(
      scenario_names,
      function(x) {
        array(
          NA_real_,
          dim = c(n_draws, model_data$C, model_data$T)
        )
      }
    ),
    scenario_names
  )

  for (country_index in seq_len(model_data$C)) {
    anchor_cases <- rep(
      model_data$cases[country_index, anchor_index],
      n_draws
    )
    q_anchor <- gtd_draw_lineage_matrix(
      draws$q,
      country_index,
      anchor_index
    )
    q_state <- setNames(
      lapply(scenario_names, function(x) q_anchor),
      scenario_names
    )
    for (scenario in scenario_names) {
      cases[[scenario]][, country_index, anchor_index] <- anchor_cases
    }

    for (time_index in seq.int(anchor_index + 1L, model_data$T)) {
      baseline_step <- gtd_scenario_step(
        cases[["baseline"]][, country_index, time_index - 1L],
        q_state[["baseline"]],
        time_index,
        country_index,
        draws,
        model_data
      )
      no_exposure_step <- gtd_scenario_step(
        cases[["no_post_2022_exposure"]][
          ,
          country_index,
          time_index - 1L
        ],
        q_state[["no_post_2022_exposure"]],
        time_index,
        country_index,
        draws,
        model_data,
        remove_import = months[[time_index]] >= as.Date("2023-01-01")
      )
      neutral_step <- gtd_scenario_step(
        cases[["l10207_neutral_growth"]][
          ,
          country_index,
          time_index - 1L
        ],
        q_state[["l10207_neutral_growth"]],
        time_index,
        country_index,
        draws,
        model_data,
        neutral_lineage_index = neutral_lineage_index
      )

      cases[["baseline"]][, country_index, time_index] <-
        baseline_step$cases
      cases[["no_post_2022_exposure"]][, country_index, time_index] <-
        no_exposure_step$cases
      cases[["l10207_neutral_growth"]][, country_index, time_index] <-
        neutral_step$cases
      q_state[["baseline"]] <- baseline_step$q
      q_state[["no_post_2022_exposure"]] <- no_exposure_step$q
      q_state[["l10207_neutral_growth"]] <- neutral_step$q
    }
  }

  first_month_index <- anchor_index + 1L
  first_month_max_abs_error <- vapply(
    seq_len(model_data$C),
    function(country_index) {
      max(abs(
        cases[["baseline"]][, country_index, first_month_index] -
          draws$mu_cases[, country_index, first_month_index]
      ))
    },
    numeric(1L)
  )

  list(
    cases = cases,
    anchor_index = anchor_index,
    first_month_index = first_month_index,
    first_month_max_abs_error = first_month_max_abs_error
  )
}

gtd_drawwise_cumulative_difference <- function(
    baseline_cases,
    scenario_cases,
    evaluation_indices) {
  if (!identical(dim(baseline_cases), dim(scenario_cases))) {
    stop("baseline_cases and scenario_cases must have identical dimensions")
  }
  if (length(dim(baseline_cases)) != 2L) {
    stop("Case inputs must be n_draws by n_months matrices")
  }
  baseline_cumulative <- rowSums(
    baseline_cases[, evaluation_indices, drop = FALSE]
  )
  scenario_cumulative <- rowSums(
    scenario_cases[, evaluation_indices, drop = FALSE]
  )
  if (any(!is.finite(baseline_cumulative)) ||
      any(baseline_cumulative <= 0)) {
    stop("Every cumulative baseline draw must be finite and positive")
  }
  case_difference <- baseline_cumulative - scenario_cumulative

  list(
    baseline_cumulative = baseline_cumulative,
    scenario_cumulative = scenario_cumulative,
    case_difference = case_difference,
    difference_fraction = case_difference / baseline_cumulative
  )
}

gtd_scenario_quantiles <- function(x) {
  c(
    mean = mean(x),
    median = median(x),
    lower_95 = unname(quantile(x, 0.025)),
    upper_95 = unname(quantile(x, 0.975))
  )
}
