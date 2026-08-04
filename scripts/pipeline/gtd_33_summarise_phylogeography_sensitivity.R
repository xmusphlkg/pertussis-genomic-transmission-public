#!/usr/bin/env Rscript

# Consolidate deterministic tree/root sensitivities and sampling jackknives
# into compact country-score and country-rank tables.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    paste(
      "Usage: gtd_33_summarise_phylogeography_sensitivity.R",
      "RESULTS_ROOT OUTPUT_DIR"
    )
  )
}

results_root <- normalizePath(args[[1L]], mustWork = TRUE)
outdir <- args[[2L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

focal_countries <- c("AUS", "BEL", "CHN", "FRA", "JPN")
score_semantics <- paste(
  "Non-complementary heuristic ancestry scores conditional on the sampled",
  "tree, geographic-state reconstruction, and retained reference set"
)

analyses <- data.table(
  analysis_id = c(
    "primary_ml_midpoint",
    "primary_ml_outgroup",
    "bootstrap_consensus_midpoint",
    "bionj_midpoint",
    "historical_reference_keep80",
    "global_background_keep80"
  ),
  analysis_class = c(
    "primary",
    "rooting",
    "topology",
    "topology",
    "reference_sampling",
    "background_sampling"
  ),
  analysis_label = c(
    "Primary ML tree, midpoint root",
    "Primary ML tree, outgroup root",
    "Bootstrap consensus tree, midpoint root",
    "BIONJ distance tree, midpoint root",
    "Historical focal reference, 80% retained",
    "Global background, 80% retained"
  ),
  tree_basis = c(
    "maximum_likelihood",
    "maximum_likelihood",
    "bootstrap_consensus",
    "bionj_distance",
    "maximum_likelihood_pruned",
    "maximum_likelihood_pruned"
  ),
  rooting = c(
    "midpoint",
    "outgroup_mrca",
    "midpoint",
    "midpoint",
    "fixed_primary_root",
    "fixed_primary_root"
  ),
  topology_file_or_pattern = c(
    "results/phylogeny/gtd_gate_rescue4_depth_core_snp.treefile",
    paste0(
      "results/phylogeography_alternative_root/",
      "gtd_primary_final_outgroup_rooted.tree"
    ),
    "results/phylogeny/gtd_gate_rescue4_depth_core_snp.contree",
    paste0(
      "results/phylogeography_bionj/",
      "gtd_gate_rescue4_depth_core_snp.bionj"
    ),
    paste0(
      "results/phylogeography_reference_sensitivity/",
      "replicate_*/historical_jackknife.tree"
    ),
    paste0(
      "results/phylogeography_background_sensitivity/",
      "replicate_*/historical_jackknife.tree"
    )
  ),
  sampling_perturbation = c(
    "none",
    "none",
    "none",
    "none",
    "stratified_20_percent_historical_focal_deletion",
    "stratified_20_percent_global_background_deletion"
  ),
  result_type = c("point", "point", "point", "point", "jackknife", "jackknife"),
  source_file = c(
    "results/phylogeography/events_thr0_5/country_phylogeographic_summary.tsv",
    paste0(
      "results/phylogeography_alternative_root/events_thr0_5/",
      "country_phylogeographic_summary.tsv"
    ),
    paste0(
      "results/phylogeography_consensus/events_thr0_5/",
      "country_phylogeographic_summary.tsv"
    ),
    paste0(
      "results/phylogeography_bionj/events_thr0_5/",
      "country_phylogeographic_summary.tsv"
    ),
    paste0(
      "results/phylogeography_reference_sensitivity/",
      "country_score_jackknife_summary.tsv"
    ),
    paste0(
      "results/phylogeography_background_sensitivity/",
      "country_score_jackknife_summary.tsv"
    )
  ),
  replicate_source_file = c(
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    paste0(
      "results/phylogeography_reference_sensitivity/",
      "country_score_replicates.tsv"
    ),
    paste0(
      "results/phylogeography_background_sensitivity/",
      "country_score_replicates.tsv"
    )
  )
)

resolve_source <- function(relative_path) {
  sub("^results/", paste0(results_root, "/"), relative_path)
}

required_files <- unique(
  c(
    analyses$source_file,
    analyses[!is.na(replicate_source_file), replicate_source_file]
  )
)
missing_files <- required_files[
  !file.exists(vapply(required_files, resolve_source, character(1)))
]
if (length(missing_files)) {
  stop(
    paste(
      "Missing required phylogeography sensitivity inputs:",
      paste(missing_files, collapse = ", ")
    )
  )
}

read_point_scores <- function(spec) {
  frame <- fread(resolve_source(spec$source_file))
  if (!setequal(frame$country_iso3, focal_countries)) {
    stop(sprintf("Country set mismatch in %s", spec$source_file))
  }
  frame[
    ,
    .(
      analysis_id = spec$analysis_id,
      analysis_class = spec$analysis_class,
      analysis_label = spec$analysis_label,
      country_iso3,
      n_resurgence_tips,
      n_replicates = 1L,
      central_statistic = "point_estimate",
      reseeding_central = mean_post_import_support,
      reseeding_lower_95 = mean_post_import_support,
      reseeding_upper_95 = mean_post_import_support,
      local_central = mean_local_persistence_support,
      local_lower_95 = mean_local_persistence_support,
      local_upper_95 = mean_local_persistence_support,
      source_file = spec$source_file
    )
  ]
}

read_jackknife_scores <- function(spec) {
  frame <- fread(resolve_source(spec$source_file))
  replicates <- fread(resolve_source(spec$replicate_source_file))[status == "ok"]
  if (!setequal(frame$country_iso3, focal_countries)) {
    stop(sprintf("Country set mismatch in %s", spec$source_file))
  }
  n_tips <- unique(
    replicates[
      ,
      .(country_iso3, n_resurgence_tips)
    ]
  )
  if (n_tips[, .N, by = country_iso3][N != 1L, .N]) {
    stop(sprintf("Resurgence-tip counts vary in %s", spec$replicate_source_file))
  }
  frame <- merge(frame, n_tips, by = "country_iso3", all.x = TRUE)
  frame[
    ,
    .(
      analysis_id = spec$analysis_id,
      analysis_class = spec$analysis_class,
      analysis_label = spec$analysis_label,
      country_iso3,
      n_resurgence_tips,
      n_replicates,
      central_statistic = "jackknife_median",
      reseeding_central = reseeding_median,
      reseeding_lower_95,
      reseeding_upper_95,
      local_central = local_median,
      local_lower_95,
      local_upper_95,
      source_file = spec$source_file
    )
  ]
}

score_rows <- lapply(
  seq_len(nrow(analyses)),
  function(index) {
    spec <- analyses[index]
    if (spec$result_type == "point") {
      read_point_scores(spec)
    } else {
      read_jackknife_scores(spec)
    }
  }
)
score_comparison <- rbindlist(score_rows, use.names = TRUE)
primary_scores <- score_comparison[
  analysis_id == "primary_ml_midpoint",
  .(
    country_iso3,
    primary_reseeding = reseeding_central,
    primary_local = local_central
  )
]
score_comparison <- merge(
  score_comparison,
  primary_scores,
  by = "country_iso3",
  all.x = TRUE
)
score_comparison[
  ,
  `:=`(
    reseeding_delta_from_primary = reseeding_central - primary_reseeding,
    local_delta_from_primary = local_central - primary_local,
    primary_reseeding = NULL,
    primary_local = NULL
  )
]
setcolorder(
  score_comparison,
  c(
    "analysis_id",
    "analysis_class",
    "analysis_label",
    "country_iso3",
    "n_resurgence_tips",
    "n_replicates",
    "central_statistic",
    "reseeding_central",
    "reseeding_lower_95",
    "reseeding_upper_95",
    "reseeding_delta_from_primary",
    "local_central",
    "local_lower_95",
    "local_upper_95",
    "local_delta_from_primary",
    "source_file"
  )
)
score_comparison[
  ,
  `:=`(
    analysis_order = match(analysis_id, analyses$analysis_id),
    country_order = match(country_iso3, focal_countries)
  )
]
setorder(score_comparison, analysis_order, country_order)
score_comparison[, c("analysis_order", "country_order") := NULL]
fwrite(
  score_comparison,
  file.path(outdir, "country_ancestry_score_comparison.tsv"),
  sep = "\t"
)

rank_point <- function(spec) {
  frame <- fread(resolve_source(spec$source_file))
  frame[
    ,
    `:=`(
      local_rank = frank(
        -mean_local_persistence_support,
        ties.method = "min"
      ),
      reseeding_rank = frank(
        -mean_post_import_support,
        ties.method = "min"
      )
    )
  ]
  frame[
    ,
    .(
      analysis_id = spec$analysis_id,
      analysis_class = spec$analysis_class,
      analysis_label = spec$analysis_label,
      country_iso3,
      n_replicates = 1L,
      local_rank_median = local_rank,
      local_rank_lower_95 = local_rank,
      local_rank_upper_95 = local_rank,
      probability_highest_local = as.numeric(local_rank == 1L),
      reseeding_rank_median = reseeding_rank,
      reseeding_rank_lower_95 = reseeding_rank,
      reseeding_rank_upper_95 = reseeding_rank,
      probability_highest_reseeding = as.numeric(reseeding_rank == 1L)
    )
  ]
}

rank_jackknife <- function(spec) {
  frame <- fread(resolve_source(spec$replicate_source_file))[status == "ok"]
  frame[
    ,
    `:=`(
      local_rank = frank(
        -mean_local_persistence_support,
        ties.method = "min"
      ),
      reseeding_rank = frank(
        -mean_post_import_support,
        ties.method = "min"
      )
    ),
    by = replicate
  ]
  frame[
    ,
    .(
      analysis_id = spec$analysis_id,
      analysis_class = spec$analysis_class,
      analysis_label = spec$analysis_label,
      n_replicates = uniqueN(replicate),
      local_rank_median = median(local_rank),
      local_rank_lower_95 = quantile(local_rank, 0.025),
      local_rank_upper_95 = quantile(local_rank, 0.975),
      probability_highest_local = mean(local_rank == 1L),
      reseeding_rank_median = median(reseeding_rank),
      reseeding_rank_lower_95 = quantile(reseeding_rank, 0.025),
      reseeding_rank_upper_95 = quantile(reseeding_rank, 0.975),
      probability_highest_reseeding = mean(reseeding_rank == 1L)
    ),
    by = country_iso3
  ]
}

rank_rows <- lapply(
  seq_len(nrow(analyses)),
  function(index) {
    spec <- analyses[index]
    if (spec$result_type == "point") {
      rank_point(spec)
    } else {
      rank_jackknife(spec)
    }
  }
)
rank_stability <- rbindlist(rank_rows, use.names = TRUE)
rank_stability[
  ,
  `:=`(
    analysis_order = match(analysis_id, analyses$analysis_id),
    country_order = match(country_iso3, focal_countries)
  )
]
setorder(rank_stability, analysis_order, country_order)
rank_stability[, c("analysis_order", "country_order") := NULL]
fwrite(
  rank_stability,
  file.path(outdir, "country_ancestry_rank_stability.tsv"),
  sep = "\t"
)

manifest <- analyses[
  ,
  .(
    analysis_id,
    analysis_class,
    analysis_label,
    tree_basis,
    rooting,
    topology_file_or_pattern,
    sampling_perturbation,
    result_type,
    n_replicates = ifelse(
      result_type == "point",
      1L,
      vapply(
        replicate_source_file,
        function(path) {
          if (is.na(path)) {
            return(NA_integer_)
          }
          uniqueN(fread(resolve_source(path))[status == "ok", replicate])
        },
        integer(1)
      )
    ),
    score_semantics,
    source_file,
    replicate_source_file
  )
]
fwrite(
  manifest,
  file.path(outdir, "phylogeography_sensitivity_manifest.tsv"),
  sep = "\t",
  na = ""
)

message(
  sprintf(
    "Wrote %d country-score rows and %d rank-stability rows to %s",
    nrow(score_comparison),
    nrow(rank_stability),
    normalizePath(outdir)
  )
)
