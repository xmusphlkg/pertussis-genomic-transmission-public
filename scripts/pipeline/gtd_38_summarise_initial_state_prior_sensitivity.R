#!/usr/bin/env Rscript

# Compare the strict pre-2019 historical initial-state prior with a weak
# symmetric Dirichlet prior using compact posterior outputs from complete
# model refits.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "Usage: gtd_38_summarise_initial_state_prior_sensitivity.R",
      "STRICT_PRIMARY_POSTERIOR_RDS STRICT_PRIMARY_DIAGNOSTICS_JSON",
      "SYMMETRIC_POSTERIOR_RDS SYMMETRIC_DIAGNOSTICS_JSON OUTDIR"
    )
  )
}

lineages <- c(
  "L1_01.02",
  "L1_02.05",
  "L1_02.06",
  "L1_02.07",
  "Other"
)
target_lineage <- "L1_02.07"

inputs <- data.table(
  analysis_id = c("historical_pre2019", "symmetric_dirichlet"),
  analysis_label = c(
    "Strict primary: historical genomes before 2019",
    "Sensitivity: symmetric Dirichlet(0.5)"
  ),
  initial_state_prior = c(
    "country_specific_historical_genomes_strictly_before_2019",
    "symmetric_dirichlet_0.5_each_lineage"
  ),
  posterior_rds = c(args[[1L]], args[[3L]]),
  diagnostics_json = c(args[[2L]], args[[4L]])
)

inputs[
  ,
  `:=`(
    posterior_rds = vapply(
      posterior_rds,
      normalizePath,
      character(1L),
      mustWork = TRUE
    ),
    diagnostics_json = vapply(
      diagnostics_json,
      normalizePath,
      character(1L),
      mustWork = TRUE
    )
  )
]

outdir <- args[[5L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)

scalar_number <- function(payload, field, source) {
  value <- payload[[field]]
  if (
    is.null(value) ||
      length(value) != 1L ||
      !is.numeric(value) ||
      !is.finite(value)
  ) {
    stop(
      "Diagnostic field '", field,
      "' must be one finite number in ", source
    )
  }
  as.numeric(value)
}

read_diagnostics <- function(path) {
  payload <- fromJSON(path, simplifyVector = TRUE)
  if (!is.list(payload)) {
    stop("Diagnostics JSON must contain an object: ", path)
  }

  minimum_neff_field <- if (!is.null(payload$minimum_neff)) {
    "minimum_neff"
  } else if (!is.null(payload$minimum_bulk_neff)) {
    "minimum_bulk_neff"
  } else {
    stop(
      "Diagnostics JSON must contain minimum_neff or minimum_bulk_neff: ",
      path
    )
  }

  values <- list(
    chains = scalar_number(payload, "chains", path),
    iterations_per_chain = scalar_number(
      payload,
      "iterations_per_chain",
      path
    ),
    post_warmup_draws = scalar_number(payload, "post_warmup_draws", path),
    seed = scalar_number(payload, "seed", path),
    divergent_transitions = scalar_number(
      payload,
      "divergent_transitions",
      path
    ),
    maximum_treedepth_hits = scalar_number(
      payload,
      "maximum_treedepth_hits",
      path
    ),
    maximum_rhat = scalar_number(payload, "maximum_rhat", path),
    minimum_effective_sample_size = scalar_number(
      payload,
      minimum_neff_field,
      path
    ),
    minimum_effective_sample_size_field = minimum_neff_field
  )

  integer_fields <- c(
    "chains",
    "iterations_per_chain",
    "post_warmup_draws",
    "seed",
    "divergent_transitions",
    "maximum_treedepth_hits"
  )
  for (field in integer_fields) {
    if (values[[field]] != floor(values[[field]])) {
      stop("Diagnostic field '", field, "' must be integer-valued in ", path)
    }
    values[[field]] <- as.integer(values[[field]])
  }
  if (
    values$chains <= 0L ||
      values$iterations_per_chain <= 0L ||
      values$post_warmup_draws <= 0L ||
      values$divergent_transitions < 0L ||
      values$maximum_treedepth_hits < 0L ||
      values$maximum_rhat <= 0 ||
      values$minimum_effective_sample_size <= 0
  ) {
    stop("Diagnostics contain an out-of-range value: ", path)
  }
  values
}

read_lineage_draws <- function(path, expected_draws) {
  posterior <- readRDS(path)
  if (
    !is.list(posterior) ||
      !"lineage_relative_transmission" %in% names(posterior)
  ) {
    stop(
      "Posterior RDS is missing lineage_relative_transmission: ",
      path
    )
  }

  draws <- posterior$lineage_relative_transmission
  if (!is.matrix(draws) || length(dim(draws)) != 2L) {
    stop(
      "lineage_relative_transmission must be a two-dimensional matrix: ",
      path
    )
  }
  if (ncol(draws) != length(lineages)) {
    stop(
      "lineage_relative_transmission must have exactly 5 columns; found ",
      ncol(draws), " in ", path
    )
  }
  if (nrow(draws) != expected_draws) {
    stop(
      "Posterior draw count (", nrow(draws),
      ") does not match post_warmup_draws (", expected_draws,
      ") in diagnostics for ", path
    )
  }
  if (
    !is.numeric(draws) ||
      any(!is.finite(draws)) ||
      any(draws <= 0)
  ) {
    stop(
      "lineage_relative_transmission draws must be finite and positive: ",
      path
    )
  }
  colnames(draws) <- lineages
  draws
}

summarise_draws <- function(values) {
  list(
    mean = mean(values),
    median = median(values),
    lower_95 = unname(quantile(values, 0.025)),
    upper_95 = unname(quantile(values, 0.975)),
    posterior_probability_above_one = mean(values > 1)
  )
}

diagnostic_payloads <- lapply(inputs$diagnostics_json, read_diagnostics)
posterior_draws <- lapply(seq_len(nrow(inputs)), function(index) {
  read_lineage_draws(
    inputs$posterior_rds[[index]],
    diagnostic_payloads[[index]]$post_warmup_draws
  )
})

lineage_summary <- rbindlist(lapply(seq_len(nrow(inputs)), function(index) {
  rbindlist(lapply(seq_along(lineages), function(lineage_index) {
    cbind(
      data.table(
        analysis_id = inputs$analysis_id[[index]],
        analysis_label = inputs$analysis_label[[index]],
        lineage = lineages[[lineage_index]]
      ),
      as.data.table(
        summarise_draws(posterior_draws[[index]][, lineage_index])
      )
    )
  }))
}))
fwrite(
  lineage_summary,
  file.path(
    outdir,
    "lineage_relative_transmission_by_initial_state_prior.tsv"
  ),
  sep = "\t"
)

target_index <- match(target_lineage, lineages)
comparator_indices <- setdiff(seq_along(lineages), target_index)
pairwise_summary <- rbindlist(lapply(seq_len(nrow(inputs)), function(index) {
  rbindlist(lapply(comparator_indices, function(comparator_index) {
    ratio <- posterior_draws[[index]][, target_index] /
      posterior_draws[[index]][, comparator_index]
    cbind(
      data.table(
        analysis_id = inputs$analysis_id[[index]],
        analysis_label = inputs$analysis_label[[index]],
        numerator = target_lineage,
        denominator = lineages[[comparator_index]]
      ),
      as.data.table(summarise_draws(ratio))
    )
  }))
}))
fwrite(
  pairwise_summary,
  file.path(outdir, "l10207_pairwise_growth_by_initial_state_prior.tsv"),
  sep = "\t"
)

diagnostics <- rbindlist(lapply(seq_len(nrow(inputs)), function(index) {
  cbind(
    data.table(
      analysis_id = inputs$analysis_id[[index]],
      analysis_label = inputs$analysis_label[[index]]
    ),
    as.data.table(diagnostic_payloads[[index]])
  )
}))
fwrite(
  diagnostics,
  file.path(outdir, "initial_state_prior_sampling_diagnostics.tsv"),
  sep = "\t"
)

manifest <- copy(inputs)
manifest[
  ,
  `:=`(
    n_posterior_draws = vapply(posterior_draws, nrow, integer(1L)),
    n_lineages = length(lineages),
    lineage_order = paste(lineages, collapse = ";"),
    target_lineage = target_lineage
  )
]
fwrite(
  manifest,
  file.path(outdir, "initial_state_prior_sensitivity_manifest.tsv"),
  sep = "\t"
)
