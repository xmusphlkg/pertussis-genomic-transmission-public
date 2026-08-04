#!/usr/bin/env Rscript

# Summarise the existing simulation-refit experiment by estimand domain without
# replacing or relaxing the original joint all-parameter identifiability gate.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: gtd_39_summarise_identifiability_domains.R",
    "ALL_RECOVERY_TSV ALL_DIAGNOSTICS_TSV IDENTIFIABILITY_GATE_JSON OUTDIR"
  ))
}

recovery_file <- args[[1L]]
diagnostics_file <- args[[2L]]
gate_file <- args[[3L]]
outdir <- args[[4L]]

for (path in c(recovery_file, diagnostics_file, gate_file)) {
  if (!file.exists(path)) {
    stop(sprintf("Required input does not exist: %s", path))
  }
}

recovery <- fread(recovery_file)
diagnostics <- fread(diagnostics_file)
gate <- fromJSON(gate_file, simplifyVector = TRUE)

required_recovery_columns <- c(
  "parameter", "50%", "Rhat", "n_eff", "replicate", "truth",
  "covered_95", "log_error_median", "parameter_type"
)
required_diagnostic_columns <- c(
  "replicate", "divergent_transitions", "maximum_treedepth_hits",
  "maximum_rhat", "minimum_neff", "true_highest_lineage",
  "recovered_highest_lineage"
)
missing_recovery <- setdiff(required_recovery_columns, names(recovery))
missing_diagnostics <- setdiff(required_diagnostic_columns, names(diagnostics))
if (length(missing_recovery) > 0L) {
  stop(sprintf(
    "Recovery input is missing required columns: %s",
    paste(missing_recovery, collapse = ", ")
  ))
}
if (length(missing_diagnostics) > 0L) {
  stop(sprintf(
    "Diagnostics input is missing required columns: %s",
    paste(missing_diagnostics, collapse = ", ")
  ))
}
if (nrow(recovery) == 0L || nrow(diagnostics) == 0L) {
  stop("Recovery and diagnostics inputs must both contain at least one row")
}
if (anyDuplicated(recovery[, .(replicate, parameter)]) > 0L) {
  stop("Recovery input has duplicate replicate-parameter rows")
}
if (anyDuplicated(diagnostics$replicate) > 0L) {
  stop("Diagnostics input has duplicate replicate rows")
}
if (!setequal(unique(recovery$replicate), diagnostics$replicate)) {
  stop("Recovery and diagnostics inputs do not contain the same replicates")
}

lineage_parameter_map <- c(
  "lineage_relative_transmission[1]" = "L1_01.02",
  "lineage_relative_transmission[2]" = "L1_02.05",
  "lineage_relative_transmission[3]" = "L1_02.06",
  "lineage_relative_transmission[4]" = "L1_02.07",
  "lineage_relative_transmission[5]" = "Other"
)
lineage_parameters <- names(lineage_parameter_map)
target_parameter <- "lineage_relative_transmission[4]"
target_analysis_identifier <- unname(lineage_parameter_map[[target_parameter]])

observed_lineage_parameters <- unique(
  recovery[parameter_type == "lineage_growth", parameter]
)
if (!setequal(observed_lineage_parameters, lineage_parameters)) {
  stop(paste(
    "Lineage-growth rows must contain the frozen five-parameter ordering:",
    paste(lineage_parameters, collapse = ", ")
  ))
}
if (any(
  recovery[parameter %chin% lineage_parameters, parameter_type] !=
    "lineage_growth"
)) {
  stop("At least one frozen lineage parameter has an inconsistent parameter_type")
}
if (any(
  recovery[grepl("^import_scale\\[", parameter), parameter_type] !=
    "import_scale"
)) {
  stop("At least one import-scale parameter has an inconsistent parameter_type")
}

numeric_recovery_columns <- c(
  "50%", "Rhat", "n_eff", "truth", "log_error_median"
)
numeric_diagnostic_columns <- c(
  "divergent_transitions", "maximum_treedepth_hits", "maximum_rhat",
  "minimum_neff", "true_highest_lineage", "recovered_highest_lineage"
)
if (any(!is.finite(as.matrix(recovery[, ..numeric_recovery_columns])))) {
  stop("Recovery input contains a non-finite required numeric value")
}
if (any(!is.finite(as.matrix(diagnostics[, ..numeric_diagnostic_columns])))) {
  stop("Diagnostics input contains a non-finite required numeric value")
}
if (!is.logical(recovery$covered_95) || anyNA(recovery$covered_95)) {
  stop("covered_95 must be a complete logical column")
}

threshold_names <- c(
  "minimum_95_interval_coverage",
  "maximum_median_absolute_log_error_lineage_growth",
  "maximum_median_absolute_log_error_import_scale",
  "minimum_truth_median_correlation",
  "minimum_highest_lineage_rank_recovery"
)
thresholds <- gate$thresholds_defined_before_recovery_results
if (is.null(thresholds) || !all(threshold_names %in% names(thresholds))) {
  stop("The existing gate does not contain all prespecified thresholds")
}
for (name in threshold_names) {
  value <- thresholds[[name]]
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value)) {
    stop(sprintf("Gate threshold is not one finite numeric scalar: %s", name))
  }
}

if (!is.null(gate$joint_all_parameters_pass)) {
  joint_all_parameters_pass <- gate$joint_all_parameters_pass
  joint_gate_source_field <- "joint_all_parameters_pass"
} else {
  joint_all_parameters_pass <- gate$pass
  joint_gate_source_field <- "pass"
}
if (
  length(joint_all_parameters_pass) != 1L ||
    !is.logical(joint_all_parameters_pass) ||
    is.na(joint_all_parameters_pass)
) {
  stop("The existing joint gate result must be one non-missing JSON boolean")
}

summarise_recovery <- function(x) {
  correlation <- cor(x$truth, x[["50%"]])
  if (!is.finite(correlation)) {
    stop("A requested recovery-domain correlation is not finite")
  }
  list(
    n_replicates = uniqueN(x$replicate),
    n_parameters = uniqueN(x$parameter),
    n_parameter_replicate_pairs = nrow(x),
    n_covered_95 = sum(x$covered_95),
    coverage_95 = mean(x$covered_95),
    median_absolute_log_error = median(abs(x$log_error_median)),
    truth_posterior_median_correlation = correlation
  )
}

lineage_recovery <- recovery[
  parameter_type == "lineage_growth" & parameter %chin% lineage_parameters
]
target_recovery <- recovery[parameter == target_parameter]
import_recovery <- recovery[parameter_type == "import_scale"]
if (nrow(target_recovery) == 0L || nrow(import_recovery) == 0L) {
  stop("Target-lineage and import-scale recovery rows must both be present")
}

lineage_metrics <- summarise_recovery(lineage_recovery)
target_metrics <- summarise_recovery(target_recovery)
import_metrics <- summarise_recovery(import_recovery)

sampler_metrics <- list(
  n_replicates = nrow(diagnostics),
  divergent_transitions = sum(diagnostics$divergent_transitions),
  maximum_treedepth_hits = sum(diagnostics$maximum_treedepth_hits),
  maximum_rhat = max(diagnostics$maximum_rhat),
  minimum_neff = min(diagnostics$minimum_neff)
)
sampler_metrics$zero_divergence_constraint_pass <-
  sampler_metrics$divergent_transitions == 0

rank_correct <- diagnostics$true_highest_lineage ==
  diagnostics$recovered_highest_lineage
rank_metrics <- list(
  n_replicates = nrow(diagnostics),
  n_correct = sum(rank_correct),
  recovery_fraction = mean(rank_correct)
)

# The JSON written by the original recovery script is rounded for display.
# This tolerance verifies that the supplied TSVs and gate belong together; it
# is not an inferential threshold and is never used in a domain pass decision.
gate_rounding_tolerance <- 5.1e-5
assert_matches_gate <- function(calculated, gate_name) {
  gate_value <- gate$observed[[gate_name]]
  if (
    length(gate_value) != 1L ||
      !is.numeric(gate_value) ||
      !is.finite(gate_value) ||
      abs(calculated - gate_value) > gate_rounding_tolerance
  ) {
    stop(sprintf(
      "Calculated metric does not match existing gate field '%s'",
      gate_name
    ))
  }
}
assert_matches_gate(
  lineage_metrics$coverage_95,
  "lineage_growth_coverage"
)
assert_matches_gate(
  lineage_metrics$median_absolute_log_error,
  "lineage_growth_median_absolute_log_error"
)
assert_matches_gate(
  lineage_metrics$truth_posterior_median_correlation,
  "lineage_growth_truth_median_correlation"
)
assert_matches_gate(
  import_metrics$coverage_95,
  "import_scale_coverage"
)
assert_matches_gate(
  import_metrics$median_absolute_log_error,
  "import_scale_median_absolute_log_error"
)
assert_matches_gate(
  import_metrics$truth_posterior_median_correlation,
  "import_scale_truth_median_correlation"
)
assert_matches_gate(
  rank_metrics$recovery_fraction,
  "highest_lineage_rank_recovery"
)

lineage_growth_domain_pass <- (
  lineage_metrics$coverage_95 >=
    thresholds$minimum_95_interval_coverage &&
  lineage_metrics$median_absolute_log_error <=
    thresholds$maximum_median_absolute_log_error_lineage_growth &&
  lineage_metrics$truth_posterior_median_correlation >=
    thresholds$minimum_truth_median_correlation &&
  rank_metrics$recovery_fraction >=
    thresholds$minimum_highest_lineage_rank_recovery &&
  sampler_metrics$zero_divergence_constraint_pass
)
import_scale_domain_pass <- (
  import_metrics$coverage_95 >=
    thresholds$minimum_95_interval_coverage &&
  import_metrics$median_absolute_log_error <=
    thresholds$maximum_median_absolute_log_error_import_scale &&
  import_metrics$truth_posterior_median_correlation >=
    thresholds$minimum_truth_median_correlation &&
  sampler_metrics$zero_divergence_constraint_pass
)

common_summary <- list(
  divergent_transitions = sampler_metrics$divergent_transitions,
  maximum_treedepth_hits = sampler_metrics$maximum_treedepth_hits,
  maximum_rhat = sampler_metrics$maximum_rhat,
  minimum_neff = sampler_metrics$minimum_neff
)

summary_rows <- rbindlist(list(
  as.data.table(c(
    list(
    estimand_group = "all_lineage_growth",
    analysis_identifier = "all_frozen_model_lineages",
    recovery_status = "threshold_evaluated"
    ),
    lineage_metrics,
    list(
    highest_lineage_rank_n_correct = rank_metrics$n_correct,
    highest_lineage_rank_n_replicates = rank_metrics$n_replicates,
    highest_lineage_rank_recovery = rank_metrics$recovery_fraction
    ),
    common_summary,
    list(
    minimum_95_interval_coverage_threshold =
      thresholds$minimum_95_interval_coverage,
    maximum_median_absolute_log_error_threshold =
      thresholds$maximum_median_absolute_log_error_lineage_growth,
    minimum_truth_median_correlation_threshold =
      thresholds$minimum_truth_median_correlation,
    minimum_highest_lineage_rank_recovery_threshold =
      thresholds$minimum_highest_lineage_rank_recovery,
    zero_divergence_constraint_pass =
      sampler_metrics$zero_divergence_constraint_pass,
    domain_pass = lineage_growth_domain_pass
    )
  )),
  as.data.table(c(
    list(
    estimand_group = "target_lineage_growth",
    analysis_identifier = target_analysis_identifier,
    recovery_status = sprintf(
      "descriptive_only_n_%d_replicates",
      target_metrics$n_replicates
    )
    ),
    target_metrics,
    list(
    highest_lineage_rank_n_correct = NA_integer_,
    highest_lineage_rank_n_replicates = NA_integer_,
    highest_lineage_rank_recovery = NA_real_
    ),
    common_summary,
    list(
    minimum_95_interval_coverage_threshold = NA_real_,
    maximum_median_absolute_log_error_threshold = NA_real_,
    minimum_truth_median_correlation_threshold = NA_real_,
    minimum_highest_lineage_rank_recovery_threshold = NA_real_,
    zero_divergence_constraint_pass =
      sampler_metrics$zero_divergence_constraint_pass,
    domain_pass = NA
    )
  )),
  as.data.table(c(
    list(
    estimand_group = "import_scale",
    analysis_identifier = "all_country_import_scale_parameters",
    recovery_status = "threshold_evaluated"
    ),
    import_metrics,
    list(
    highest_lineage_rank_n_correct = NA_integer_,
    highest_lineage_rank_n_replicates = NA_integer_,
    highest_lineage_rank_recovery = NA_real_
    ),
    common_summary,
    list(
    minimum_95_interval_coverage_threshold =
      thresholds$minimum_95_interval_coverage,
    maximum_median_absolute_log_error_threshold =
      thresholds$maximum_median_absolute_log_error_import_scale,
    minimum_truth_median_correlation_threshold =
      thresholds$minimum_truth_median_correlation,
    minimum_highest_lineage_rank_recovery_threshold = NA_real_,
    zero_divergence_constraint_pass =
      sampler_metrics$zero_divergence_constraint_pass,
    domain_pass = import_scale_domain_pass
    )
  ))
), use.names = TRUE, fill = TRUE)

interpretation <- list(
  schema_version = "1.0",
  joint_all_parameters_pass = joint_all_parameters_pass,
  lineage_growth_domain_pass = lineage_growth_domain_pass,
  import_scale_domain_pass = import_scale_domain_pass,
  target_specific_inference_status = sprintf(
    "descriptive_only_n_%d_replicates",
    target_metrics$n_replicates
  ),
  domains = list(
    all_lineage_growth = c(
      lineage_metrics,
      list(
        highest_lineage_rank_recovery = rank_metrics$recovery_fraction,
        thresholds_copied_from_existing_gate = list(
          minimum_95_interval_coverage =
            thresholds$minimum_95_interval_coverage,
          maximum_median_absolute_log_error =
            thresholds$maximum_median_absolute_log_error_lineage_growth,
          minimum_truth_median_correlation =
            thresholds$minimum_truth_median_correlation,
          minimum_highest_lineage_rank_recovery =
            thresholds$minimum_highest_lineage_rank_recovery
        ),
        zero_divergence_constraint_pass =
          sampler_metrics$zero_divergence_constraint_pass,
        domain_pass = lineage_growth_domain_pass
      )
    ),
    target_lineage_growth = c(
      list(
        analysis_identifier = target_analysis_identifier,
        model_parameter = target_parameter
      ),
      target_metrics,
      list(
        inference_status = sprintf(
          "descriptive_only_n_%d_replicates",
          target_metrics$n_replicates
        ),
        threshold_evaluated = FALSE,
        reason = paste(
          "Target-specific recovery is reported descriptively because only",
          sprintf("%d recovery replicates were run.", target_metrics$n_replicates)
        )
      )
    ),
    import_scale = c(
      import_metrics,
      list(
        thresholds_copied_from_existing_gate = list(
          minimum_95_interval_coverage =
            thresholds$minimum_95_interval_coverage,
          maximum_median_absolute_log_error =
            thresholds$maximum_median_absolute_log_error_import_scale,
          minimum_truth_median_correlation =
            thresholds$minimum_truth_median_correlation
        ),
        zero_divergence_constraint_pass =
          sampler_metrics$zero_divergence_constraint_pass,
        domain_pass = import_scale_domain_pass
      )
    )
  ),
  highest_lineage_rank_recovery = rank_metrics,
  sampler_diagnostics = c(
    sampler_metrics,
    list(
      gating_note = paste(
        "Only the zero-divergence constraint belonged to the original joint",
        "gate. Maximum-treedepth hits, R-hat, and effective sample size are",
        "reported descriptively without newly introduced cutoffs."
      )
    )
  ),
  interpretation = list(
    joint_gate = paste(
      "The original joint all-parameter gate is retained exactly and is not",
      "replaced by either domain-specific result."
    ),
    lineage_growth = paste(
      "The all-lineage growth domain is evaluated with only the lineage",
      "thresholds, highest-lineage rank threshold, and zero-divergence",
      "constraint prespecified in the existing gate."
    ),
    import_scale = paste(
      "Within this simulation-refit design, import-scale parameters are",
      "structurally less recoverable than shared lineage growth. Scenario",
      "evidence that depends on import scale must remain conditional."
    )
  ),
  provenance = list(
    input_files = list(
      all_recovery = list(
        basename = basename(recovery_file),
        md5 = unname(tools::md5sum(recovery_file))
      ),
      all_diagnostics = list(
        basename = basename(diagnostics_file),
        md5 = unname(tools::md5sum(diagnostics_file))
      ),
      identifiability_gate = list(
        basename = basename(gate_file),
        md5 = unname(tools::md5sum(gate_file))
      )
    ),
    joint_gate_source_field = joint_gate_source_field,
    joint_gate_value_copied_not_recomputed = TRUE,
    domain_pass_does_not_replace_joint_gate = TRUE,
    threshold_source =
      "thresholds_defined_before_recovery_results in the supplied gate",
    new_inferential_thresholds_added = FALSE,
    diagnostic_constraint_reused =
      "zero total divergent transitions",
    diagnostic_metrics_reported_without_new_cutoffs = c(
      "maximum_treedepth_hits", "maximum_rhat", "minimum_neff"
    ),
    frozen_lineage_parameter_map = as.list(lineage_parameter_map),
    gate_input_consistency_tolerance = gate_rounding_tolerance,
    gate_input_consistency_tolerance_role =
      "display-rounding check only; never used for inference or domain pass"
  )
)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  summary_rows,
  file.path(outdir, "estimand_recovery_summary.tsv"),
  sep = "\t",
  na = ""
)
write_json(
  interpretation,
  file.path(outdir, "identifiability_domain_interpretation.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  digits = NA
)

message(sprintf(
  paste(
    "Joint gate=%s; lineage-growth domain=%s;",
    "import-scale domain=%s"
  ),
  joint_all_parameters_pass,
  lineage_growth_domain_pass,
  import_scale_domain_pass
))
