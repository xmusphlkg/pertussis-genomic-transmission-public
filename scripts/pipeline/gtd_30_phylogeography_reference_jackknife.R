#!/usr/bin/env Rscript

# Reconstruct geographic ancestry after deterministic stratified deletion of
# historical focal tips. All resurgence and global-background tips are retained.
# This tests whether country-level ancestry scores depend on a small subset of
# the historical focal reference.

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9L) {
  stop(
    paste(
      "Usage: gtd_30_phylogeography_reference_jackknife.R TREE METADATA",
      "TIP_GEOGRAPHY LINEAGES PASTML_BIN PYTHON_BIN EVENT_SCRIPT OUTDIR",
      "N_REPLICATES [HISTORICAL_KEEP_FRACTION] [CORES]",
      "[BACKGROUND_KEEP_FRACTION]"
    )
  )
}

tree_file <- normalizePath(args[[1L]], mustWork = TRUE)
metadata_file <- normalizePath(args[[2L]], mustWork = TRUE)
tip_geography_file <- normalizePath(args[[3L]], mustWork = TRUE)
lineages_file <- normalizePath(args[[4L]], mustWork = TRUE)
pastml_bin <- normalizePath(args[[5L]], mustWork = TRUE)
python_bin <- normalizePath(args[[6L]], mustWork = TRUE)
event_script <- normalizePath(args[[7L]], mustWork = TRUE)
outdir <- args[[8L]]
n_replicates <- as.integer(args[[9L]])
keep_fraction <- if (length(args) >= 10L) as.numeric(args[[10L]]) else 0.80
cores <- if (length(args) >= 11L) as.integer(args[[11L]]) else 8L
background_keep_fraction <- if (length(args) >= 12L) {
  as.numeric(args[[12L]])
} else {
  1
}

if (
  !is.finite(keep_fraction) ||
    keep_fraction <= 0 ||
    keep_fraction > 1
) {
  stop("HISTORICAL_KEEP_FRACTION must be greater than 0 and at most 1")
}
if (
  !is.finite(background_keep_fraction) ||
    background_keep_fraction <= 0 ||
    background_keep_fraction > 1
) {
  stop("BACKGROUND_KEEP_FRACTION must be greater than 0 and at most 1")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
tree <- read.tree(tree_file)
meta <- fread(metadata_file)
geo <- fread(tip_geography_file)
focal_countries <- c("AUS", "BEL", "CHN", "FRA", "JPN")

if (!setequal(tree$tip.label, meta$tree_sample_id)) {
  stop("Tree tips and metadata identifiers differ")
}

historical <- meta[
  tree_role == "focal" &
    country_iso3 %chin% focal_countries &
    epidemic_period %chin% c("prepandemic", "pandemic")
]
background <- copy(meta[tree_role == "global_background"])
background[
  is.na(preliminary_lineage_id) | !nzchar(preliminary_lineage_id),
  preliminary_lineage_id := "Unknown"
]
background[
  is.na(epidemic_period) | !nzchar(epidemic_period),
  epidemic_period := "Unknown"
]
background[
  is.na(continent) | !nzchar(continent),
  continent := "Unknown"
]

run_replicate <- function(replicate_id) {
  set.seed(20260830L + replicate_id)
  retained_historical <- if (keep_fraction < 1) {
    historical[
      ,
      {
        retain_n <- max(1L, ceiling(.N * keep_fraction))
        .SD[sample.int(.N, retain_n)]
      },
      by = .(country_iso3, epidemic_period)
    ]$tree_sample_id
  } else {
    historical$tree_sample_id
  }
  retained_background <- if (background_keep_fraction < 1) {
    background[
      ,
      {
        retain_n <- max(1L, ceiling(.N * background_keep_fraction))
        .SD[sample.int(.N, retain_n)]
      },
      by = .(continent, epidemic_period, preliminary_lineage_id)
    ]$tree_sample_id
  } else {
    background$tree_sample_id
  }
  remove_historical <- setdiff(
    historical$tree_sample_id,
    retained_historical
  )
  remove_background <- setdiff(
    background$tree_sample_id,
    retained_background
  )
  remove_tips <- union(remove_historical, remove_background)
  pruned <- drop.tip(tree, remove_tips, trim.internal = TRUE)
  pruned$node.label <- sprintf(
    "JK%03d_NODE_%04d",
    replicate_id,
    seq_len(pruned$Nnode)
  )

  rep_dir <- file.path(outdir, sprintf("replicate_%03d", replicate_id))
  pastml_dir <- file.path(rep_dir, "pastml")
  event_dir <- file.path(rep_dir, "events")
  dir.create(pastml_dir, recursive = TRUE, showWarnings = FALSE)
  tree_out <- file.path(rep_dir, "historical_jackknife.tree")
  geo_out <- file.path(rep_dir, "tip_geography.tsv")
  metadata_out <- file.path(rep_dir, "tree_tip_metadata.tsv")
  write.tree(pruned, tree_out)
  fwrite(
    geo[tree_sample_id %chin% pruned$tip.label],
    geo_out,
    sep = "\t"
  )
  fwrite(
    meta[tree_sample_id %chin% pruned$tip.label],
    metadata_out,
    sep = "\t"
  )

  pastml_status <- system2(
    pastml_bin,
    c(
      "--tree", tree_out,
      "--data", geo_out,
      "--columns", "geo_state",
      "--prediction_method", "MPPA",
      "--model", "F81",
      "--out_data", file.path(pastml_dir, "ancestral_states.tsv"),
      "--work_dir", pastml_dir,
      "--threads", "1"
    ),
    stdout = file.path(rep_dir, "pastml.stdout.log"),
    stderr = file.path(rep_dir, "pastml.stderr.log")
  )
  if (pastml_status != 0L) {
    return(data.table(replicate = replicate_id, status = "pastml_failed"))
  }

  marginals <- file.path(
    pastml_dir,
    "marginal_probabilities.character_geo_state.model_F81.tab"
  )
  event_status <- system2(
    python_bin,
    c(
      event_script,
      tree_out,
      marginals,
      metadata_out,
      lineages_file,
      event_dir,
      "--transition-threshold", "0.5"
    ),
    stdout = file.path(rep_dir, "events.stdout.log"),
    stderr = file.path(rep_dir, "events.stderr.log")
  )
  if (event_status != 0L) {
    return(data.table(replicate = replicate_id, status = "events_failed"))
  }

  summary_file <- file.path(event_dir, "country_phylogeographic_summary.tsv")
  result <- fread(summary_file)
  result[, `:=`(
    replicate = replicate_id,
    status = "ok",
    retained_tree_tips = length(pruned$tip.label),
    removed_historical_tips = length(remove_historical),
    removed_background_tips = length(remove_background)
  )]
  result
}

results <- rbindlist(
  parallel::mclapply(
    seq_len(n_replicates),
    run_replicate,
    mc.cores = min(cores, n_replicates)
  ),
  fill = TRUE
)
fwrite(results, file.path(outdir, "country_score_replicates.tsv"), sep = "\t")

successful <- results[status == "ok"]
if (uniqueN(successful$replicate) != n_replicates) {
  warning(
    sprintf(
      "%d of %d historical-reference replicates completed",
      uniqueN(successful$replicate),
      n_replicates
    )
  )
}
if (!nrow(successful)) {
  stop("No historical-reference replicate completed")
}

score_summary <- successful[
  ,
  .(
    n_replicates = .N,
    reseeding_median = median(mean_post_import_support),
    reseeding_lower_95 = quantile(mean_post_import_support, 0.025),
    reseeding_upper_95 = quantile(mean_post_import_support, 0.975),
    local_median = median(mean_local_persistence_support),
    local_lower_95 = quantile(mean_local_persistence_support, 0.025),
    local_upper_95 = quantile(mean_local_persistence_support, 0.975)
  ),
  by = country_iso3
]
fwrite(
  score_summary,
  file.path(outdir, "country_score_jackknife_summary.tsv"),
  sep = "\t"
)

rank_summary <- successful[
  ,
  .(
    highest_local_country = country_iso3[
      which.max(mean_local_persistence_support)
    ],
    highest_reseeding_country = country_iso3[
      which.max(mean_post_import_support)
    ]
  ),
  by = replicate
][
  ,
  .N,
  by = .(highest_local_country, highest_reseeding_country)
][order(-N)]
rank_summary[, fraction := N / sum(N)]
fwrite(
  rank_summary,
  file.path(outdir, "country_rank_stability.tsv"),
  sep = "\t"
)
