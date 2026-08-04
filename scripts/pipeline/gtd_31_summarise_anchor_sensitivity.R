#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: gtd_31_summarise_anchor_sensitivity.R",
      "MAIN_POSTERIOR_RDS NO_PROJECT_POSTERIOR_RDS DATA_JSON OUTDIR"
    )
  )
}

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the script location")
}
script_file <- normalizePath(sub("^--file=", "", script_arg))
source(file.path(dirname(script_file), "..", "lib", "gtd_anchored_scenarios.R"))

posterior_paths <- c(main = args[[1L]], no_project = args[[2L]])
model_data <- fromJSON(args[[3L]], simplifyVector = TRUE)
outdir <- args[[4L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

countries <- c("AUS", "CHN", "JPN")
lineages <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
months <- seq(as.Date("2019-01-01"), as.Date("2025-12-01"), by = "month")
anchor_months <- as.Date(c(
  "2022-06-01", "2022-12-01", "2023-06-01", "2023-12-01"
))
evaluation_indices <- which(
  months >= as.Date("2024-01-01") &
    months <= as.Date("2025-12-01")
)
l10207_index <- match("L1_02.07", lineages)

if (model_data$C != length(countries) ||
    model_data$L != length(lineages) ||
    model_data$T != length(months)) {
  stop("Frozen country, lineage, or month dimensions have changed")
}

summary_rows <- list()
row_index <- 0L
for (model_specification in names(posterior_paths)) {
  message("Reading ", model_specification, " posterior draws")
  draws <- readRDS(posterior_paths[[model_specification]])
  gtd_validate_scenario_inputs(draws, model_data, months)
  n_draws <- dim(draws$q)[1L]

  for (anchor_index in match(anchor_months, months)) {
    if (is.na(anchor_index)) {
      stop("At least one anchor month is outside the frozen model period")
    }
    anchor_month <- months[[anchor_index]]
    scenario_run <- gtd_run_anchored_scenarios(
      draws,
      model_data,
      months,
      anchor_index,
      l10207_index
    )

    for (country_index in seq_along(countries)) {
      baseline <- scenario_run$cases[["baseline"]][
        ,
        country_index,
        ,
        drop = FALSE
      ]
      dim(baseline) <- c(n_draws, model_data$T)

      for (scenario in c(
        "no_post_2022_exposure",
        "l10207_neutral_growth"
      )) {
        scenario_cases <- scenario_run$cases[[scenario]][
          ,
          country_index,
          ,
          drop = FALSE
        ]
        dim(scenario_cases) <- c(n_draws, model_data$T)
        contrast <- gtd_drawwise_cumulative_difference(
          baseline,
          scenario_cases,
          evaluation_indices
        )
        case_summary <- gtd_scenario_quantiles(contrast$case_difference)
        fraction_summary <- gtd_scenario_quantiles(
          contrast$difference_fraction
        )
        scenario_cumulative_summary <- gtd_scenario_quantiles(
          contrast$scenario_cumulative
        )
        baseline_cumulative_summary <- gtd_scenario_quantiles(
          contrast$baseline_cumulative
        )

        intervention_start <- if (
          scenario == "no_post_2022_exposure"
        ) {
          max(
            months[[anchor_index + 1L]],
            as.Date("2023-01-01")
          )
        } else {
          months[[anchor_index + 1L]]
        }
        row_index <- row_index + 1L
        summary_rows[[row_index]] <- data.table(
          model_specification = model_specification,
          anchor_month = format(anchor_month, "%Y-%m"),
          intervention_start_month = format(
            intervention_start,
            "%Y-%m"
          ),
          evaluation_start_month = "2024-01",
          evaluation_end_month = "2025-12",
          country_iso3 = countries[[country_index]],
          scenario = scenario,
          n_posterior_draws = n_draws,
          baseline_cumulative_cases_median =
            baseline_cumulative_summary[["median"]],
          scenario_cumulative_cases_median =
            scenario_cumulative_summary[["median"]],
          cumulative_case_difference_mean = case_summary[["mean"]],
          cumulative_case_difference_median = case_summary[["median"]],
          cumulative_case_difference_lower_95 =
            case_summary[["lower_95"]],
          cumulative_case_difference_upper_95 =
            case_summary[["upper_95"]],
          cumulative_difference_fraction_mean =
            fraction_summary[["mean"]],
          cumulative_difference_fraction_median =
            fraction_summary[["median"]],
          cumulative_difference_fraction_lower_95 =
            fraction_summary[["lower_95"]],
          cumulative_difference_fraction_upper_95 =
            fraction_summary[["upper_95"]],
          posterior_probability_difference_above_zero =
            mean(contrast$difference_fraction > 0),
          scenario_equals_baseline_all_draws =
            all(contrast$case_difference == 0),
          first_month_baseline_max_abs_error_from_stan_mu =
            scenario_run$first_month_max_abs_error[[country_index]]
        )
      }
    }
  }

  rm(draws)
  invisible(gc())
}

summary_table <- rbindlist(summary_rows)
setorder(
  summary_table,
  model_specification,
  anchor_month,
  country_iso3,
  scenario
)
fwrite(
  summary_table,
  file.path(outdir, "anchor_scenario_sensitivity.tsv"),
  sep = "\t"
)

metadata <- list(
  analysis = "anchored conditional-scenario sensitivity",
  anchors = format(anchor_months, "%Y-%m"),
  evaluation_period = c("2024-01", "2025-12"),
  contrast_definition = paste(
    "(baseline cumulative expected reported cases - scenario cumulative",
    "expected reported cases) / baseline cumulative expected reported cases;",
    "computed within each posterior draw before summarisation and not clipped"
  ),
  no_post_2022_exposure = paste(
    "Import exposure is set to zero from 2023-01 or the first month after",
    "the anchor, whichever is later; all other terms are unchanged."
  ),
  l10207_neutral_growth = paste(
    "Only the L1_02.07 lineage relative net-growth multiplier is set to 1",
    "from the first month after the anchor; exposure and all other lineage",
    "multipliers are retained."
  ),
  baseline_check = paste(
    "The first recursively projected baseline month is compared draw-wise",
    "with the corresponding Stan transformed-parameter mu_cases."
  )
)
write_json(
  metadata,
  file.path(outdir, "anchor_scenario_sensitivity_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
