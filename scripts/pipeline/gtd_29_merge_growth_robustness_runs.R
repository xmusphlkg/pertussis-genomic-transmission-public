#!/usr/bin/env Rscript

# Merge the three prespecified growth-robustness run components:
#   1. country-only, project-adjusted country omissions, and dominant-project
#      omissions;
#   2. no-project country omissions; and
#   3. the deliberately overdispersed omit-Japan mode audit.
#
# Every input is validated before any formal output is replaced. This prevents
# accidental de-duplication on analysis_id alone from collapsing the lineage
# and pairwise tables.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    paste(
      "Usage: gtd_29_merge_growth_robustness_runs.R",
      "BASE_RUN_DIR NO_PROJECT_LOO_DIR OMIT_JPN_AUDIT_DIR OUTDIR"
    )
  )
}

run_dirs <- setNames(
  lapply(args[1:3], normalizePath, mustWork = TRUE),
  c("base_refits", "no_project_country_omissions", "omit_japan_mode_audit")
)
outdir <- args[[4L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

expected_analysis_ids <- c(
  "country_AUS_project_adjusted",
  "country_AUS_no_project",
  "country_CHN_project_adjusted",
  "country_CHN_no_project",
  "country_JPN_project_adjusted",
  "country_JPN_no_project",
  "omit_country_AUS",
  "omit_country_CHN",
  "omit_country_JPN",
  "omit_country_AUS_no_project",
  "omit_country_CHN_no_project",
  "omit_country_JPN_no_project",
  "omit_dominant_project_AUS",
  "omit_dominant_project_CHN",
  "omit_dominant_project_JPN"
)
expected_lineages <- c(
  "L1_01.02",
  "L1_02.05",
  "L1_02.06",
  "L1_02.07",
  "Other"
)
expected_comparators <- setdiff(expected_lineages, "L1_02.07")

table_specs <- list(
  lineage_relative_growth_by_analysis.tsv = list(
    keys = c("analysis_id", "lineage"),
    rows_per_analysis = length(expected_lineages)
  ),
  l1_02_07_growth_robustness.tsv = list(
    keys = "analysis_id",
    rows_per_analysis = 1L
  ),
  l1_02_07_pairwise_probabilities.tsv = list(
    keys = c("analysis_id", "comparator_lineage"),
    rows_per_analysis = length(expected_comparators)
  ),
  fit_diagnostics.tsv = list(
    keys = "analysis_id",
    rows_per_analysis = 1L
  )
)

read_component_table <- function(filename) {
  tables <- lapply(names(run_dirs), function(component) {
    path <- file.path(run_dirs[[component]], filename)
    if (!file.exists(path)) {
      stop("Missing component table: ", path)
    }
    value <- fread(path)
    value[, run_component := component]
    value
  })
  common_names <- lapply(tables, function(value) {
    setdiff(names(value), "run_component")
  })
  if (!all(vapply(
    common_names[-1L],
    identical,
    logical(1),
    common_names[[1L]]
  ))) {
    stop("Component schemas differ for ", filename)
  }
  rbindlist(tables, use.names = TRUE)
}

merged <- lapply(names(table_specs), read_component_table)
names(merged) <- names(table_specs)

for (filename in names(merged)) {
  value <- merged[[filename]]
  spec <- table_specs[[filename]]
  if (anyDuplicated(value[, ..spec$keys])) {
    stop("Duplicate composite keys in merged ", filename)
  }
  if (!setequal(value$analysis_id, expected_analysis_ids)) {
    stop(
      "Unexpected analysis IDs in ",
      filename,
      ": ",
      paste(setdiff(value$analysis_id, expected_analysis_ids), collapse = ", ")
    )
  }
  counts <- value[, .N, by = analysis_id]
  if (
    nrow(counts) != length(expected_analysis_ids) ||
      any(counts$N != spec$rows_per_analysis)
  ) {
    stop("Unexpected per-analysis row counts in ", filename)
  }
}

lineage_table <- merged[["lineage_relative_growth_by_analysis.tsv"]]
if (!all(
  lineage_table[
    ,
    setequal(lineage, expected_lineages),
    by = analysis_id
  ]$V1
)) {
  stop("At least one analysis has an incomplete lineage set")
}

pairwise_table <- merged[["l1_02_07_pairwise_probabilities.tsv"]]
if (
  !"target_lineage" %in% names(pairwise_table) ||
    !all(pairwise_table$target_lineage == "L1_02.07") ||
    !all(
      pairwise_table[
        ,
        setequal(comparator_lineage, expected_comparators),
        by = analysis_id
      ]$V1
    )
) {
  stop("At least one analysis has an invalid pairwise-comparator set")
}

target_table <- merged[["l1_02_07_growth_robustness.tsv"]]
diagnostic_table <- merged[["fit_diagnostics.tsv"]]
if (
  !"diagnostic_pass" %in% names(target_table) ||
    !"diagnostic_pass" %in% names(diagnostic_table) ||
    !all(target_table$diagnostic_pass) ||
    !all(diagnostic_table$diagnostic_pass)
) {
  stop("A formal growth-robustness result failed its diagnostic gate")
}

audit_target <- target_table[analysis_id == "omit_country_JPN"]
audit_diagnostic <- diagnostic_table[analysis_id == "omit_country_JPN"]
if (
  nrow(audit_target) != 1L ||
    nrow(audit_diagnostic) != 1L ||
    audit_target$lower_95 <= 1 ||
    audit_diagnostic$chains != 4L ||
    audit_diagnostic$iterations_per_chain != 4000L ||
    audit_diagnostic$warmup_per_chain != 2000L ||
    audit_diagnostic$divergent_transitions != 0L ||
    audit_diagnostic$maximum_treedepth_hits != 0L ||
    audit_diagnostic$maximum_rhat >= 1.01
) {
  stop("The omit-Japan mode audit does not satisfy the formal acceptance gate")
}

run_configs <- rbindlist(
  lapply(names(run_dirs), function(component) {
    path <- file.path(run_dirs[[component]], "run_configuration.tsv")
    if (!file.exists(path)) {
      stop("Missing component run configuration: ", path)
    }
    value <- fread(path)
    value[, run_component := component]
    value
  }),
  use.names = TRUE,
  fill = TRUE
)
audit_config <- run_configs[run_component == "omit_japan_mode_audit"]
if (
  nrow(audit_config) != 1L ||
    audit_config$chains != 4L ||
    audit_config$iterations_per_chain != 4000L ||
    audit_config$warmup_per_chain != 2000L ||
    audit_config$adapt_delta != 0.99 ||
    audit_config$max_treedepth != 15L ||
    audit_config$init_mode != "mode_audit"
) {
  stop("The omit-Japan run configuration is not the prespecified mode audit")
}

analysis_order <- setNames(seq_along(expected_analysis_ids), expected_analysis_ids)
for (filename in names(merged)) {
  value <- copy(merged[[filename]])
  value[, analysis_order__ := unname(analysis_order[analysis_id])]
  extra_order <- intersect(
    c("lineage", "comparator_lineage"),
    names(value)
  )
  setorderv(value, c("analysis_order__", extra_order))
  value[, c("analysis_order__", "run_component") := NULL]
  merged[[filename]] <- value
}

# All validation above completes before formal outputs are written.
for (filename in names(merged)) {
  fwrite(merged[[filename]], file.path(outdir, filename), sep = "\t")
}
fwrite(
  run_configs,
  file.path(outdir, "run_configuration_merged_components.tsv"),
  sep = "\t"
)
fwrite(
  run_configs,
  file.path(outdir, "run_configuration.tsv"),
  sep = "\t"
)

message(
  "Merged ",
  length(expected_analysis_ids),
  " analyses: ",
  nrow(lineage_table),
  " lineage rows, ",
  nrow(pairwise_table),
  " pairwise rows."
)
