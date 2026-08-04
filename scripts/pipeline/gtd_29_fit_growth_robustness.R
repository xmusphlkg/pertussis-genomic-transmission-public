#!/usr/bin/env Rscript

# Formal robustness refits for the shared relative lineage-growth result.
#
# The script fits:
#   1. each model country separately with and without project effects; and
#   2. the project-adjusted model after omitting each country in turn; and
#   3. the full project-adjusted model after omitting the dominant genomic
#      project from each country in turn.
#
# Only compact posterior summaries and diagnostics are written. Large stanfit
# objects and posterior arrays remain in memory and are not deposited.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(paste(
    "Usage: gtd_29_fit_growth_robustness.R",
    "STAN_OR_COMPILED_RDS DATA_JSON OBS_TSV PROJECT_INDEX_TSV OUTDIR",
    paste(
      "[CHAINS] [ITER] [SEED] [PARALLEL_FITS] [ANALYSIS_REGEX]",
      "[ADAPT_DELTA] [MAX_TREEDEPTH] [INIT_MODE]"
    )
  ))
}

model_path <- normalizePath(args[[1L]], mustWork = TRUE)
data_path <- normalizePath(args[[2L]], mustWork = TRUE)
observation_path <- normalizePath(args[[3L]], mustWork = TRUE)
project_index_path <- normalizePath(args[[4L]], mustWork = TRUE)
outdir <- args[[5L]]
chains <- if (length(args) >= 6L) as.integer(args[[6L]]) else 4L
iter <- if (length(args) >= 7L) as.integer(args[[7L]]) else 2000L
seed <- if (length(args) >= 8L) as.integer(args[[8L]]) else 20261300L
parallel_fits <- if (length(args) >= 9L) as.integer(args[[9L]]) else 3L
analysis_regex <- if (length(args) >= 10L) args[[10L]] else ".*"
adapt_delta <- if (length(args) >= 11L) as.numeric(args[[11L]]) else 0.95
max_treedepth <- if (length(args) >= 12L) as.integer(args[[12L]]) else 13L
init_mode <- if (length(args) >= 13L) args[[13L]] else "random"
warmup <- floor(iter / 2)

if (chains < 2L || iter < 400L || warmup >= iter) {
  stop("At least two chains and 400 iterations are required")
}
if (parallel_fits < 1L) {
  stop("PARALLEL_FITS must be at least one")
}
if (
  !is.finite(adapt_delta) ||
    adapt_delta <= 0 ||
    adapt_delta >= 1 ||
    max_treedepth < 10L
) {
  stop("Invalid ADAPT_DELTA or MAX_TREEDEPTH")
}
if (
  !init_mode %in%
    c("random", "zero", "prior_centered", "mode_audit")
) {
  stop(
    "INIT_MODE must be random, zero, prior_centered, or mode_audit"
  )
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rstan_options(auto_write = FALSE)

message("Loading frozen model inputs")
base <- fromJSON(data_path, simplifyVector = TRUE)
observations <- fread(observation_path)
project_index <- fread(project_index_path)

required_observation_columns <- c(
  "country_iso3", "country_id", "project_id", "project_id_numeric"
)
if (!all(required_observation_columns %in% names(observations))) {
  stop(
    "Observation table is missing columns: ",
    paste(
      setdiff(required_observation_columns, names(observations)),
      collapse = ", "
    )
  )
}
count_columns <- grep("^n_", names(observations), value = TRUE)
if (length(count_columns) != base$L) {
  stop("Observation lineage columns do not match L in the frozen JSON")
}
lineages <- sub("^n_", "", count_columns)
target_lineage <- "L1_02.07"
target_index <- match(target_lineage, lineages)
if (is.na(target_index)) {
  stop("The target lineage L1_02.07 is absent from the observation table")
}
if (nrow(observations) != base$J) {
  stop("Observation row count does not match J in the frozen JSON")
}
if (!identical(
  unname(as.matrix(observations[, ..count_columns])),
  unname(base$y_genome)
)) {
  stop("Observation lineage counts do not match y_genome in the frozen JSON")
}
if (!identical(
  as.integer(observations$country_id),
  as.integer(base$obs_country)
)) {
  stop("Observation country indices do not match the frozen JSON")
}
if (!identical(
  as.integer(observations$project_id_numeric),
  as.integer(base$obs_project)
)) {
  stop("Observation project indices do not match the frozen JSON")
}

country_lookup <- unique(
  observations[, .(country_id = as.integer(country_id), country_iso3)]
)
setorder(country_lookup, country_id)
if (!identical(country_lookup$country_id, seq_len(base$C))) {
  stop("Country indices are not complete and consecutive")
}
countries <- country_lookup$country_iso3

project_lookup <- unique(
  project_index[
    ,
    .(
      project_id_numeric = as.integer(project_id_numeric),
      project_id
    )
  ]
)
setorder(project_lookup, project_id_numeric)
if (!all(unique(base$obs_project) %in% project_lookup$project_id_numeric)) {
  stop("The project index does not cover all projects in the frozen JSON")
}

observations[, n_genomes := rowSums(.SD), .SDcols = count_columns]
dominant_projects <- observations[
  ,
  .(n_genomes = sum(n_genomes)),
  by = .(country_iso3, project_id, project_id_numeric)
][
  order(country_iso3, -n_genomes, project_id)
][
  ,
  .SD[1L],
  by = country_iso3
]
if (!identical(sort(dominant_projects$country_iso3), sort(countries))) {
  stop("A dominant project could not be resolved for every model country")
}

if (grepl("\\.stan$", model_path, ignore.case = TRUE)) {
  message("Compiling the current Stan source")
  model <- stan_model(
    file = model_path,
    model_name = "gtd_growth_robustness"
  )
  model_kind <- "stan_source"
} else {
  message("Loading the supplied compiled Stan model")
  model <- readRDS(model_path)
  if (!inherits(model, "stanmodel")) {
    stop("The first argument is neither a Stan source file nor a stanmodel RDS")
  }
  model_kind <- "compiled_rds"
}

analysis_specs <- list()
for (country_id in seq_len(base$C)) {
  country <- countries[[country_id]]
  analysis_specs[[length(analysis_specs) + 1L]] <- list(
    analysis_id = sprintf("country_%s_project_adjusted", country),
    analysis_type = "country_only",
    keep_countries = country_id,
    use_project_effects = 1L,
    omitted_country = NA_character_,
    omitted_project_id_numeric = NA_integer_,
    omitted_project_id = NA_character_,
    omitted_project_country = NA_character_,
    n_omitted_genomes = 0L
  )
  analysis_specs[[length(analysis_specs) + 1L]] <- list(
    analysis_id = sprintf("country_%s_no_project", country),
    analysis_type = "country_only",
    keep_countries = country_id,
    use_project_effects = 0L,
    omitted_country = NA_character_,
    omitted_project_id_numeric = NA_integer_,
    omitted_project_id = NA_character_,
    omitted_project_country = NA_character_,
    n_omitted_genomes = 0L
  )
}
for (country_id in seq_len(base$C)) {
  country <- countries[[country_id]]
  n_country_genomes <- as.integer(
    observations[country_iso3 == country, sum(n_genomes)]
  )
  analysis_specs[[length(analysis_specs) + 1L]] <- list(
    analysis_id = sprintf("omit_country_%s", country),
    analysis_type = "country_omission",
    keep_countries = setdiff(seq_len(base$C), country_id),
    use_project_effects = 1L,
    omitted_country = country,
    omitted_project_id_numeric = NA_integer_,
    omitted_project_id = NA_character_,
    omitted_project_country = NA_character_,
    n_omitted_genomes = n_country_genomes
  )
  analysis_specs[[length(analysis_specs) + 1L]] <- list(
    analysis_id = sprintf("omit_country_%s_no_project", country),
    analysis_type = "country_omission",
    keep_countries = setdiff(seq_len(base$C), country_id),
    use_project_effects = 0L,
    omitted_country = country,
    omitted_project_id_numeric = NA_integer_,
    omitted_project_id = NA_character_,
    omitted_project_country = NA_character_,
    n_omitted_genomes = n_country_genomes
  )
}
for (row_id in seq_len(nrow(dominant_projects))) {
  row <- dominant_projects[row_id]
  analysis_specs[[length(analysis_specs) + 1L]] <- list(
    analysis_id = sprintf(
      "omit_dominant_project_%s",
      row$country_iso3
    ),
    analysis_type = "dominant_project_omission",
    keep_countries = seq_len(base$C),
    use_project_effects = 1L,
    omitted_country = NA_character_,
    omitted_project_id_numeric = as.integer(row$project_id_numeric),
    omitted_project_id = row$project_id,
    omitted_project_country = row$country_iso3,
    n_omitted_genomes = as.integer(row$n_genomes)
  )
}
analysis_specs <- Filter(
  function(spec) grepl(analysis_regex, spec$analysis_id),
  analysis_specs
)
if (!length(analysis_specs)) {
  stop("ANALYSIS_REGEX did not match any robustness analysis")
}

subset_model_data <- function(spec) {
  keep_observations <- base$obs_country %in% spec$keep_countries
  if (!is.na(spec$omitted_project_id_numeric)) {
    keep_observations <- keep_observations &
      base$obs_project != spec$omitted_project_id_numeric
  }
  keep_rows <- which(keep_observations)
  used_projects <- sort(unique(base$obs_project[keep_rows]))
  project_remap <- setNames(seq_along(used_projects), used_projects)
  country_remap <- setNames(
    seq_along(spec$keep_countries),
    spec$keep_countries
  )

  if (!length(keep_rows) || !length(used_projects)) {
    stop("A robustness specification has no genomic observations")
  }

  list(
    C = length(spec$keep_countries),
    T = base$T,
    L = base$L,
    P = length(used_projects),
    J = length(keep_rows),
    K = base$K,
    use_project_effects = spec$use_project_effects,
    cases = base$cases[spec$keep_countries, , drop = FALSE],
    B = base$B,
    reporting_change = base$reporting_change[
      spec$keep_countries,
      ,
      drop = FALSE
    ],
    initial_alpha = base$initial_alpha[
      spec$keep_countries,
      ,
      drop = FALSE
    ],
    import_exposure = base$import_exposure[
      spec$keep_countries,
      ,
      ,
      drop = FALSE
    ],
    obs_country = as.integer(
      country_remap[as.character(base$obs_country[keep_rows])]
    ),
    obs_month = as.integer(base$obs_month[keep_rows]),
    obs_project = as.integer(
      project_remap[as.character(base$obs_project[keep_rows])]
    ),
    y_genome = base$y_genome[keep_rows, , drop = FALSE]
  )
}

make_initial_values <- function(stan_data, fit_seed) {
  if (init_mode == "random") {
    return("random")
  }
  if (init_mode == "zero") {
    return(0)
  }
  lapply(seq_len(chains), function(chain_id) {
    set.seed(fit_seed + chain_id)
    q0 <- stan_data$initial_alpha /
      rowSums(stan_data$initial_alpha)
    initial_values <- list(
      r_coef = matrix(
        rnorm(stan_data$C * stan_data$K, 0, 0.02),
        nrow = stan_data$C,
        ncol = stan_data$K
      ),
      log_theta_raw = rnorm(stan_data$L, 0, 0.02),
      import_scale = exp(
        rnorm(stan_data$C, log(100), 0.05)
      ),
      density_feedback = pmax(
        0.01,
        rnorm(stan_data$C, 0.5, 0.03)
      ),
      q0 = q0,
      project_raw = matrix(
        rnorm(stan_data$P * stan_data$L, 0, 0.05),
        nrow = stan_data$P,
        ncol = stan_data$L
      ),
      sigma_project = exp(rnorm(1L, log(0.5), 0.05)),
      reporting_jump = rnorm(stan_data$C, 0, 0.02),
      phi_cases = exp(rnorm(stan_data$C, log(10), 0.05)),
      genome_concentration = exp(rnorm(1L, log(20), 0.05))
    )
    if (init_mode == "mode_audit" && chain_id == 2L) {
      # Deliberately initialise one chain near the empirically observed
      # low-growth/high-project-heterogeneity mode. Convergence is accepted
      # only if this chain joins the other dispersed chains.
      initial_values$log_theta_raw[] <- 0
      initial_values$log_theta_raw[4L] <- -0.65
      initial_values$import_scale[] <- 100
      initial_values$import_scale[1L] <- 6
      initial_values$project_raw <- matrix(
        rnorm(stan_data$P * stan_data$L, 0, 1),
        nrow = stan_data$P,
        ncol = stan_data$L
      )
      initial_values$sigma_project <- 4.8
      initial_values$genome_concentration <- 3.5
    } else if (init_mode == "mode_audit" && chain_id == chains) {
      # A second overdispersed chain explores the broad prior region without
      # being targeted to either identified mode.
      initial_values$r_coef <- matrix(
        rnorm(stan_data$C * stan_data$K, 0, 0.4),
        nrow = stan_data$C,
        ncol = stan_data$K
      )
      initial_values$log_theta_raw <- rnorm(stan_data$L, 0, 0.2)
      initial_values$import_scale <- exp(
        rnorm(stan_data$C, log(100), 0.8)
      )
      initial_values$project_raw <- matrix(
        rnorm(stan_data$P * stan_data$L, 0, 1),
        nrow = stan_data$P,
        ncol = stan_data$L
      )
      initial_values$sigma_project <- exp(rnorm(1L, log(0.8), 0.5))
      initial_values$genome_concentration <- exp(
        rnorm(1L, log(20), 0.5)
      )
    }
    initial_values
  })
}

summarise_fit <- function(spec, fit, stan_data, fit_seed, elapsed_seconds) {
  rr_draws <- rstan::extract(
    fit,
    pars = "lineage_relative_transmission",
    permuted = TRUE
  )$lineage_relative_transmission

  rr_summary <- rbindlist(lapply(seq_along(lineages), function(lineage_id) {
    values <- rr_draws[, lineage_id]
    data.table(
      analysis_id = spec$analysis_id,
      lineage = lineages[[lineage_id]],
      mean = mean(values),
      median = median(values),
      lower_95 = unname(quantile(values, 0.025)),
      upper_95 = unname(quantile(values, 0.975)),
      probability_above_reference = mean(values > 1)
    )
  }))

  target_draws <- rr_draws[, target_index]
  other_indices <- setdiff(seq_along(lineages), target_index)
  probability_highest <- mean(
    target_draws > apply(
      rr_draws[, other_indices, drop = FALSE],
      1L,
      max
    )
  )
  target_summary <- rr_summary[lineage == target_lineage]
  target_summary[
    ,
    `:=`(
      analysis_type = spec$analysis_type,
      country_subset = paste(
        countries[spec$keep_countries],
        collapse = "+"
      ),
      observation_specification = if (
        spec$use_project_effects == 1L
      ) {
        "project_adjusted"
      } else {
        "no_project"
      },
      omitted_country = spec$omitted_country,
      omitted_project_country = spec$omitted_project_country,
      omitted_project_id = spec$omitted_project_id,
      n_countries = stan_data$C,
      n_genome_strata = stan_data$J,
      n_genomes = sum(stan_data$y_genome),
      n_omitted_genomes = spec$n_omitted_genomes,
      probability_highest_lineage = probability_highest
    )
  ]
  setcolorder(
    target_summary,
    c(
      "analysis_id", "analysis_type", "country_subset",
      "observation_specification", "omitted_country",
      "omitted_project_country",
      "omitted_project_id", "n_countries", "n_genome_strata",
      "n_genomes", "n_omitted_genomes", "lineage", "mean", "median",
      "lower_95", "upper_95", "probability_above_reference",
      "probability_highest_lineage"
    )
  )

  pairwise <- rbindlist(lapply(other_indices, function(other_id) {
    data.table(
      analysis_id = spec$analysis_id,
      target_lineage = target_lineage,
      comparator_lineage = lineages[[other_id]],
      probability_target_greater = mean(
        target_draws > rr_draws[, other_id]
      ),
      median_log_ratio = median(
        log(target_draws / rr_draws[, other_id])
      ),
      lower_95_log_ratio = unname(quantile(
        log(target_draws / rr_draws[, other_id]),
        0.025
      )),
      upper_95_log_ratio = unname(quantile(
        log(target_draws / rr_draws[, other_id]),
        0.975
      ))
    )
  }))

  diagnostic_pars <- c(
    "r_coef", "log_theta_raw", "import_scale", "density_feedback",
    "q0", "reporting_jump", "phi_cases", "genome_concentration"
  )
  if (spec$use_project_effects == 1L) {
    diagnostic_pars <- c(
      diagnostic_pars,
      "project_raw",
      "sigma_project"
    )
  }
  parameter_summary <- summary(
    fit,
    pars = diagnostic_pars,
    probs = c(0.025, 0.5, 0.975)
  )$summary
  target_parameter_summary <- summary(
    fit,
    pars = "lineage_relative_transmission",
    probs = c(0.025, 0.5, 0.975)
  )$summary[target_index, , drop = FALSE]
  sampler <- get_sampler_params(fit, inc_warmup = FALSE)
  diagnostics <- data.table(
    analysis_id = spec$analysis_id,
    analysis_type = spec$analysis_type,
    chains = chains,
    iterations_per_chain = iter,
    warmup_per_chain = warmup,
    post_warmup_draws = chains * (iter - warmup),
    seed = fit_seed,
    elapsed_seconds = elapsed_seconds,
    divergent_transitions = sum(vapply(
      sampler,
      function(x) sum(x[, "divergent__"]),
      numeric(1)
    )),
    maximum_treedepth_hits = sum(vapply(
      sampler,
      function(x) sum(x[, "treedepth__"] >= max_treedepth),
      numeric(1)
    )),
    maximum_rhat = max(parameter_summary[, "Rhat"], na.rm = TRUE),
    minimum_effective_sample_size = min(
      parameter_summary[, "n_eff"],
      na.rm = TRUE
    ),
    target_lineage_rhat = target_parameter_summary[1L, "Rhat"],
    target_lineage_effective_sample_size =
      target_parameter_summary[1L, "n_eff"]
  )
  diagnostics[
    ,
    diagnostic_pass := (
      maximum_rhat < 1.01 &
        minimum_effective_sample_size >= 400 &
        target_lineage_rhat < 1.01 &
        target_lineage_effective_sample_size >= 400 &
        divergent_transitions == 0 &
        maximum_treedepth_hits == 0
    )
  ]
  target_summary[
    ,
    diagnostic_pass := diagnostics$diagnostic_pass
  ]
  setcolorder(
    target_summary,
    c(
      setdiff(names(target_summary), "diagnostic_pass"),
      "diagnostic_pass"
    )
  )

  list(
    all_lineages = rr_summary,
    target = target_summary,
    pairwise = pairwise,
    diagnostics = diagnostics
  )
}

fit_one <- function(spec_id) {
  spec <- analysis_specs[[spec_id]]
  fit_seed <- seed + spec_id * 100L
  stan_data <- subset_model_data(spec)
  initial_values <- make_initial_values(stan_data, fit_seed)
  message(sprintf(
    "Fitting %s (%d countries, %d genomic strata, %d genomes)",
    spec$analysis_id,
    stan_data$C,
    stan_data$J,
    sum(stan_data$y_genome)
  ))
  timing <- system.time({
    fit <- sampling(
      model,
      data = stan_data,
      chains = chains,
      cores = chains,
      iter = iter,
      warmup = warmup,
      seed = fit_seed,
      init = initial_values,
      refresh = 0,
      control = list(
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth
      )
    )
  })
  result <- summarise_fit(
    spec,
    fit,
    stan_data,
    fit_seed,
    unname(timing[["elapsed"]])
  )
  rm(fit)
  invisible(gc())
  result
}

run_started <- Sys.time()
worker <- function(spec_id) {
  tryCatch(
    fit_one(spec_id),
    error = function(e) {
      structure(
        list(
          analysis_id = analysis_specs[[spec_id]]$analysis_id,
          message = conditionMessage(e)
        ),
        class = "growth_robustness_fit_error"
      )
    }
  )
}

if (.Platform$OS.type == "unix" && parallel_fits > 1L) {
  results <- parallel::mclapply(
    seq_along(analysis_specs),
    worker,
    mc.cores = min(parallel_fits, length(analysis_specs)),
    mc.preschedule = FALSE
  )
} else {
  results <- lapply(seq_along(analysis_specs), worker)
}

failed <- vapply(
  results,
  inherits,
  logical(1),
  what = "growth_robustness_fit_error"
)
if (any(failed)) {
  failure_table <- rbindlist(lapply(results[failed], as.data.table))
  fwrite(
    failure_table,
    file.path(outdir, "fit_failures.tsv"),
    sep = "\t"
  )
  stop(
    "At least one growth-robustness fit failed; see fit_failures.tsv"
  )
}

all_lineages <- rbindlist(lapply(results, `[[`, "all_lineages"))
target_summary <- rbindlist(lapply(results, `[[`, "target"))
pairwise <- rbindlist(lapply(results, `[[`, "pairwise"))
diagnostics <- rbindlist(lapply(results, `[[`, "diagnostics"))

fwrite(
  all_lineages,
  file.path(outdir, "lineage_relative_growth_by_analysis.tsv"),
  sep = "\t"
)
fwrite(
  target_summary,
  file.path(outdir, "l1_02_07_growth_robustness.tsv"),
  sep = "\t"
)
fwrite(
  pairwise,
  file.path(outdir, "l1_02_07_pairwise_probabilities.tsv"),
  sep = "\t"
)
fwrite(
  diagnostics,
  file.path(outdir, "fit_diagnostics.tsv"),
  sep = "\t"
)

run_finished <- Sys.time()
run_configuration <- data.table(
  model_path = model_path,
  model_kind = model_kind,
  model_sha256 = digest(model_path, algo = "sha256", file = TRUE),
  data_path = data_path,
  data_sha256 = digest(data_path, algo = "sha256", file = TRUE),
  observation_path = observation_path,
  observation_sha256 = digest(
    observation_path,
    algo = "sha256",
    file = TRUE
  ),
  project_index_path = project_index_path,
  project_index_sha256 = digest(
    project_index_path,
    algo = "sha256",
    file = TRUE
  ),
  chains = chains,
  iterations_per_chain = iter,
  warmup_per_chain = warmup,
  base_seed = seed,
  parallel_fits = parallel_fits,
  analysis_regex = analysis_regex,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth,
  init_mode = init_mode,
  n_analyses = length(analysis_specs),
  started_at = format(run_started, "%Y-%m-%dT%H:%M:%S%z"),
  finished_at = format(run_finished, "%Y-%m-%dT%H:%M:%S%z"),
  elapsed_minutes = as.numeric(
    difftime(run_finished, run_started, units = "mins")
  )
)
fwrite(
  run_configuration,
  file.path(outdir, "run_configuration.tsv"),
  sep = "\t"
)

message(sprintf(
  "Completed %d formal growth-robustness fits in %.1f minutes",
  length(analysis_specs),
  run_configuration$elapsed_minutes
))
