#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
  library(phangorn)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Usage: gtd_10_prepare_phylogeography.R TREE METADATA OUTDIR [OUTGROUP_TIPS_CSV]")
}

tree_path <- normalizePath(args[[1L]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2L]], mustWork = TRUE)
outdir <- args[[3L]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tree <- read.tree(tree_path)
meta <- fread(metadata_path, na.strings = c("", "NA"))
if (!setequal(tree$tip.label, meta$tree_sample_id)) {
  stop("Tree tips do not match the frozen metadata")
}

# The ML topology is unrooted. Midpoint rooting is primary; a prespecified
# divergent outgroup clade can be supplied for root-placement sensitivity.
if (length(args) == 4L) {
  outgroup_tips <- strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]
  if (!all(outgroup_tips %in% tree$tip.label)) {
    stop("One or more requested outgroup tips are absent from the tree")
  }
  outgroup_node <- getMRCA(tree, outgroup_tips)
  rooted <- root(tree, node = outgroup_node, resolve.root = TRUE)
  rooting_label <- paste0("outgroup_mrca:", paste(outgroup_tips, collapse = ";"))
  tree_filename <- "gtd_primary_final_outgroup_rooted.tree"
} else {
  rooted <- midpoint(tree)
  rooting_label <- "midpoint"
  tree_filename <- "gtd_primary_final_midpoint_rooted.tree"
}
rooted$node.label <- sprintf("NODE_%04d", seq_len(rooted$Nnode))
write.tree(rooted, file.path(outdir, tree_filename))

focal <- c("AUS", "BEL", "CHN", "FRA", "JPN")
meta[, geo_state := fifelse(
  country_iso3 %in% focal,
  country_iso3,
  fifelse(
    continent %chin% c("Africa", "Asia", "Europe", "North America",
                       "Oceania", "South America"),
    gsub(" ", "_", paste0(continent, "_other")),
    "Unknown_other"
  )
)]

geo <- meta[, .(tree_sample_id, geo_state)]
geo[, tree_order := match(tree_sample_id, rooted$tip.label)]
setorder(geo, tree_order)
geo[, tree_order := NULL]
fwrite(geo, file.path(outdir, "tip_geography.tsv"), sep = "\t", na = "")

counts <- geo[, .N, by = geo_state][order(-N)]
fwrite(counts, file.path(outdir, "geography_state_counts.tsv"), sep = "\t")

report <- c(
  "# Phylogeography preparation",
  "",
  sprintf("- Tips: %d", Ntip(rooted)),
  sprintf("- Internal nodes: %d", rooted$Nnode),
  sprintf("- Rooting: %s (orientation only; no exact tMRCA claims).", rooting_label),
  "- Focal countries are separate states; non-focal tips are grouped by continent.",
  "- The global background was already country-time-lineage balanced before tree construction.",
  "",
  "## State counts",
  "",
  paste(capture.output(print(counts)), collapse = "\n")
)
writeLines(report, file.path(outdir, "PHYLOGEOGRAPHY_PREPARATION.md"))
