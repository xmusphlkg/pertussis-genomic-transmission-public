#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(rhierbaps)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: gtd_09_define_lineages.R ALIGNMENT METADATA OUTDIR [N_POPS] [SEED]")
}

alignment_path <- normalizePath(args[[1L]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], mustWork = TRUE)
outdir <- args[[3L]]
n_pops <- if (length(args) >= 4L) as.integer(args[[4L]]) else 20L
seed <- if (length(args) >= 5L) as.integer(args[[5L]]) else 20260725L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
message("Loading SNP alignment")
snp_matrix <- load_fasta(alignment_path, keep.singletons = FALSE)
message(sprintf(
  "Running hierBAPS on %d genomes and %d parsimony-informative SNPs",
  nrow(snp_matrix), ncol(snp_matrix)
))
fit <- hierBAPS(
  snp_matrix,
  max.depth = 2,
  n.pops = n_pops,
  quiet = FALSE,
  n.extra.rounds = 1,
  assignment.probs = TRUE,
  n.cores = 4
)

assignments <- as.data.table(fit$partition.df)
id_col <- names(assignments)[[1L]]
setnames(assignments, id_col, "tree_sample_id")

meta <- fread(metadata_path, na.strings = c("", "NA"))
if (!setequal(assignments$tree_sample_id, meta$tree_sample_id)) {
  stop("hierBAPS identifiers do not match the frozen tree metadata")
}
assignments <- merge(
  assignments,
  meta[, .(
    tree_sample_id, tree_role, country_iso3, date_lower, date_upper,
    date_resolution, year, month, epidemic_period, project_id,
    preliminary_lineage_id
  )],
  by = "tree_sample_id",
  all.x = TRUE,
  sort = FALSE
)

level_cols <- grep("^level", names(assignments), value = TRUE)
if (length(level_cols) < 1L) {
  stop("hierBAPS did not return a level-1 population assignment")
}
assignments[, model_lineage_id := sprintf("L1_%02d", as.integer(get(level_cols[[1L]])))]
if (length(level_cols) >= 2L) {
  assignments[, model_sublineage_id := sprintf(
    "%s.%02d", model_lineage_id, as.integer(get(level_cols[[2L]]))
  )]
} else {
  assignments[, model_sublineage_id := model_lineage_id]
}
assignments[, lineage_definition_method :=
  "hierBAPS_v1.1.4_parsimony_informative_core_SNPs_pre_effect_freeze"]
assignments[, lineage_definition_seed := seed]
assignments[, lineage_definition_n_pops := n_pops]

summary_l1 <- assignments[, .(
  n_genomes = .N,
  n_focal = sum(tree_role == "focal", na.rm = TRUE),
  n_countries = uniqueN(country_iso3, na.rm = TRUE),
  countries = paste(sort(unique(na.omit(country_iso3))), collapse = ";"),
  n_periods = uniqueN(epidemic_period[!is.na(epidemic_period) & epidemic_period != "unknown"]),
  periods = paste(sort(unique(na.omit(epidemic_period[epidemic_period != "unknown"]))), collapse = ";"),
  min_year = suppressWarnings(min(year, na.rm = TRUE)),
  max_year = suppressWarnings(max(year, na.rm = TRUE))
), by = model_lineage_id]
summary_l1[!is.finite(min_year), min_year := NA_real_]
summary_l1[!is.finite(max_year), max_year := NA_real_]
summary_l1[, primary_model_eligible :=
  n_genomes >= 20L & (n_countries >= 2L | n_periods >= 2L)]
summary_l1[, exclusion_reason := fifelse(
  primary_model_eligible, "",
  "Fewer than 20 genomes or insufficient cross-country/cross-period coverage"
)]

fwrite(
  assignments,
  file.path(outdir, sprintf("hierbaps_assignments_npops%d_seed%d.tsv", n_pops, seed)),
  sep = "\t", na = ""
)
fwrite(
  summary_l1[order(-n_genomes)],
  file.path(outdir, sprintf("hierbaps_lineage_summary_npops%d_seed%d.tsv", n_pops, seed)),
  sep = "\t", na = ""
)
saveRDS(
  fit,
  file.path(outdir, sprintf("hierbaps_fit_npops%d_seed%d.rds", n_pops, seed))
)

report <- c(
  "# Pre-effect transmission-lineage freeze",
  "",
  sprintf("- Alignment: `%s`", alignment_path),
  sprintf("- Genomes: %d", nrow(snp_matrix)),
  sprintf("- Parsimony-informative core SNPs: %d", ncol(snp_matrix)),
  sprintf("- hierBAPS maximum populations: %d", n_pops),
  sprintf("- Random seed: %d", seed),
  "- The clustering was defined without using case counts or post-clustering growth estimates.",
  "- Level-1 populations are candidate model lineages; level-2 populations are sensitivity sublineages.",
  "- A primary lineage requires at least 20 genomes and cross-country or cross-period coverage.",
  "",
  "## Level-1 populations",
  "",
  paste(capture.output(print(summary_l1[order(-n_genomes)])), collapse = "\n")
)
writeLines(
  report,
  file.path(outdir, sprintf("LINEAGE_FREEZE_npops%d_seed%d.md", n_pops, seed))
)
