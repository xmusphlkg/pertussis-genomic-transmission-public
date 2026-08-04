#!/usr/bin/env Rscript

# Design-weighted L1_02.07 shares under the deterministic focal-tree sampling
# caps. The weights address the study's cohort subsampling only; they do not
# address public-archive submission, sequence availability, or post-selection
# QC mechanisms.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop(paste(
    "Usage: gtd_30_build_selection_cap_weighting.R",
    "INVENTORY_TSV SELECTION_TSV UNIFORM_QC_TSV LINEAGE_TSV OUT_TSV",
    "[COUNTRIES_COMMA_SEPARATED]"
  ))
}

inventory_path <- normalizePath(args[[1L]], mustWork = TRUE)
selection_path <- normalizePath(args[[2L]], mustWork = TRUE)
uniform_qc_path <- normalizePath(args[[3L]], mustWork = TRUE)
lineage_path <- normalizePath(args[[4L]], mustWork = TRUE)
output_path <- args[[5L]]
model_countries <- if (length(args) >= 6L) {
  strsplit(args[[6L]], ",", fixed = TRUE)[[1L]]
} else {
  c("AUS", "CHN", "JPN")
}
model_countries <- trimws(model_countries)

inventory <- fread(inventory_path)
selection <- fread(selection_path)
uniform_qc <- fread(uniform_qc_path)
lineages <- fread(lineage_path)

required_inventory_columns <- c(
  "genome_record_id", "country_iso3", "year", "date_lower",
  "date_resolution", "project_id", "preliminary_lineage_id"
)
required_selection_columns <- c(
  required_inventory_columns,
  "selection_reason"
)
required_lineage_columns <- c(
  "tree_sample_id", "tree_role", "country_iso3", "epidemic_period",
  "primary_model_lineage_id"
)
if (!all(required_inventory_columns %in% names(inventory))) {
  stop("The inventory table lacks required sampling-stratum fields")
}
if (!all(required_selection_columns %in% names(selection))) {
  stop("The selection table lacks required sampling-stratum fields")
}
if (!all(c("tree_sample_id", "genome_record_id") %in% names(uniform_qc))) {
  stop("The uniform-QC table cannot map tree tips to genome records")
}
if (!all(required_lineage_columns %in% names(lineages))) {
  stop("The lineage table lacks required formal-lineage fields")
}
if (anyDuplicated(inventory$genome_record_id)) {
  stop("Inventory genome_record_id values are not unique")
}
if (anyDuplicated(selection$genome_record_id)) {
  stop("Selection genome_record_id values are not unique")
}

rescue_countries <- unique(
  selection[
    grepl(
      "^historical_all_available",
      selection_reason
    ),
    country_iso3
  ]
)

five_year_bin <- function(year) {
  start <- (as.integer(year) %/% 5L) * 5L
  fifelse(
    is.na(start),
    "unknown",
    sprintf("%d_%d", start, start + 4L)
  )
}

modern_time_bin <- function(
  year,
  date_lower,
  date_resolution
) {
  precise <- date_resolution %chin% c("day", "month", "quarter")
  lower <- as.IDate(date_lower)
  fifelse(
    precise & !is.na(lower),
    substr(as.character(lower), 1L, 7L),
    fifelse(is.na(year), "unknown", as.character(as.integer(year)))
  )
}

add_sampling_stratum <- function(frame) {
  frame <- copy(frame)
  frame[, year_numeric := as.integer(year)]
  frame[
    ,
    lineage_stratum := fifelse(
      is.na(preliminary_lineage_id) | preliminary_lineage_id == "",
      "UNKNOWN_LINEAGE",
      preliminary_lineage_id
    )
  ]
  frame[
    ,
    time_stratum := fifelse(
      year_numeric <= 2019L,
      fifelse(
        country_iso3 %chin% rescue_countries,
        paste0("all_available_", genome_record_id),
        five_year_bin(year_numeric)
      ),
      modern_time_bin(
        year_numeric,
        date_lower,
        date_resolution
      )
    )
  ]
  frame[
    ,
    sampling_stratum := paste(
      country_iso3,
      time_stratum,
      fifelse(is.na(project_id), "", project_id),
      lineage_stratum,
      sep = "|"
    )
  ]
  frame
}

inventory <- add_sampling_stratum(
  inventory[country_iso3 %chin% model_countries]
)
selection <- add_sampling_stratum(
  selection[country_iso3 %chin% model_countries]
)

inventory_counts <- inventory[
  ,
  .(n_inventory = .N),
  by = sampling_stratum
]
selection_counts <- selection[
  ,
  .(n_initially_selected = .N),
  by = sampling_stratum
]
stratum_counts <- merge(
  inventory_counts,
  selection_counts,
  by = "sampling_stratum",
  all.x = TRUE
)
stratum_counts[
  is.na(n_initially_selected),
  n_initially_selected := 0L
]

tip_map <- unique(
  uniform_qc[, .(tree_sample_id, genome_record_id)]
)
if (anyDuplicated(tip_map$tree_sample_id)) {
  stop("Uniform-QC tree_sample_id values are not unique")
}
formal <- lineages[
  tree_role == "focal" &
    country_iso3 %chin% model_countries &
    epidemic_period %chin% c("prepandemic", "resurgence")
]
formal <- merge(
  formal,
  tip_map,
  by = "tree_sample_id",
  all.x = TRUE
)
if (formal[is.na(genome_record_id) | genome_record_id == "", .N]) {
  stop("Some formal focal tips cannot be mapped to genome_record_id")
}
formal <- merge(
  formal,
  selection[, .(genome_record_id, sampling_stratum)],
  by = "genome_record_id",
  all.x = TRUE
)
if (formal[is.na(sampling_stratum), .N]) {
  stop("Some formal focal tips are absent from the frozen selection table")
}
formal <- merge(
  formal,
  stratum_counts,
  by = "sampling_stratum",
  all.x = TRUE
)
if (
  formal[
    is.na(n_inventory) |
      is.na(n_initially_selected) |
      n_initially_selected <= 0L,
    .N
  ]
) {
  stop("Invalid deterministic-sampling counts for at least one formal tip")
}

formal[
  ,
  selection_cap_weight := n_inventory / n_initially_selected
]
formal[
  ,
  is_l1_02_07 := primary_model_lineage_id == "L1_02.07"
]

weighted <- formal[
  ,
  {
    total_weight <- sum(selection_cap_weight)
    weighted_positive <- sum(
      selection_cap_weight * is_l1_02_07
    )
    .(
      n_final_tree_tips = .N,
      n_l1_02_07 = sum(is_l1_02_07),
      unweighted_share = mean(is_l1_02_07),
      weighted_l1_02_07_total = weighted_positive,
      weighted_total = total_weight,
      selection_cap_weighted_share =
        weighted_positive / total_weight,
      kish_effective_n =
        total_weight^2 / sum(selection_cap_weight^2),
      minimum_weight = min(selection_cap_weight),
      maximum_weight = max(selection_cap_weight)
    )
  },
  by = .(country_iso3, epidemic_period)
]

# Quantify how much of the public candidate inventory lies in strata that had
# at least one initially selected genome. This is a design-support diagnostic,
# not a correction for post-selection QC.
inventory[
  ,
  epidemic_period := fifelse(
    year_numeric <= 2019L,
    "prepandemic",
    fifelse(year_numeric >= 2023L, "resurgence", "pandemic")
  )
]
represented_strata <- stratum_counts[
  n_initially_selected > 0L,
  sampling_stratum
]
support <- inventory[
  epidemic_period %chin% c("prepandemic", "resurgence"),
  .(
    n_inventory_records = .N,
    n_inventory_records_in_selected_strata =
      sum(sampling_stratum %chin% represented_strata),
    inventory_stratum_support_fraction =
      mean(sampling_stratum %chin% represented_strata)
  ),
  by = .(country_iso3, epidemic_period)
]
weighted <- merge(
  weighted,
  support,
  by = c("country_iso3", "epidemic_period"),
  all.x = TRUE
)
setorder(weighted, country_iso3, epidemic_period)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
fwrite(weighted, output_path, sep = "\t")
message(sprintf(
  "Wrote deterministic selection-cap weighting for %d country-period cells",
  nrow(weighted)
))
