#!/usr/bin/env Rscript

# Crosswalk the study-defined rhierBAPS population L1_02.07 to published
# MT28 sublineage terminology without changing any frozen lineage assignment.
#
# The source-published `Sublineages` field was imported from the 8,117-genome
# global compilation (Zhang et al., Journal of Infection 2026;92:106718).
# This script tests label concordance and the MRCA span of MT28-labelled tips.
# It does not infer missing MLVA types and does not use MT28 labels to redefine
# the lineage or refit any model.

options(stringsAsFactors = FALSE, width = 180)

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (is.na(hit)) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", flag, call. = FALSE)
  }
  args[[hit + 1L]]
}

repo_root <- normalizePath(get_arg("--repo-root", "."), mustWork = TRUE)
output_dir <- get_arg(
  "--output-dir",
  file.path(repo_root, "results", "lineage_nomenclature")
)

paths <- c(
  assignments = file.path(
    repo_root,
    "results",
    "lineages",
    "primary_finalized",
    "model_lineage_assignments.tsv"
  ),
  tree_metadata = file.path(
    repo_root,
    "results",
    "phylogeny",
    "tree_tip_metadata.tsv"
  ),
  focal_records = file.path(
    repo_root,
    "data",
    "derived",
    "transmission_genome_records.tsv"
  ),
  background_records = file.path(
    repo_root,
    "data",
    "derived",
    "global_tree_background_records.tsv"
  ),
  tree = file.path(
    repo_root,
    "results",
    "phylogeny",
    "gtd_gate_rescue4_depth_core_snp.treefile"
  )
)

if (any(!file.exists(paths))) {
  stop(
    "Missing required inputs:\n",
    paste(paths[!file.exists(paths)], collapse = "\n"),
    call. = FALSE
  )
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required.", call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The jsonlite package is required.", call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(
    path,
    check.names = FALSE,
    na.strings = c("", "NA"),
    quote = "",
    comment.char = ""
  )
}

write_tsv <- function(frame, filename) {
  write.table(
    frame,
    file.path(output_dir, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

pct <- function(numerator, denominator) {
  ifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
}

assignments <- read_tsv(paths[["assignments"]])
tree_metadata <- read_tsv(paths[["tree_metadata"]])
focal_records <- read_tsv(paths[["focal_records"]])
background_records <- read_tsv(paths[["background_records"]])

stopifnot(
  !anyDuplicated(assignments$tree_sample_id),
  !anyDuplicated(tree_metadata$tree_sample_id),
  !anyDuplicated(focal_records$genome_record_id)
)

tip_metadata <- tree_metadata[
  ,
  c(
    "tree_sample_id",
    "genome_record_id",
    "sample_id",
    "biosample_accession",
    "assembly_accession",
    "country_iso3",
    "epidemic_period"
  )
]
crosswalk <- merge(
  assignments[
    ,
    c(
      "tree_sample_id",
      "tree_role",
      "model_sublineage_id",
      "primary_model_lineage_id"
    )
  ],
  tip_metadata,
  by = "tree_sample_id",
  all.x = TRUE,
  sort = FALSE
)

annotation_fields <- c(
  "published_branch",
  "published_lineage",
  "published_sublineage"
)
for (field in annotation_fields) {
  crosswalk[[field]] <- NA_character_
}
crosswalk$published_annotation_source <- NA_character_

focal_index <- match(
  crosswalk$genome_record_id,
  focal_records$genome_record_id
)
use_focal <- !is.na(focal_index)
for (field in annotation_fields) {
  crosswalk[[field]][use_focal] <-
    focal_records[[field]][focal_index[use_focal]]
}
crosswalk$published_annotation_source[use_focal] <-
  "transmission_genome_records"

# One background BioSample has two run-level records, so background matching
# uses the exact genome-record/sample pair rather than genome_record_id alone.
background_key <- paste(
  background_records$genome_record_id,
  background_records$sample_id,
  sep = "|||"
)
tip_background_key <- paste(
  crosswalk$genome_record_id,
  crosswalk$sample_id,
  sep = "|||"
)
background_index <- match(tip_background_key, background_key)
use_background <- !is.na(background_index)
for (field in annotation_fields) {
  crosswalk[[field]][use_background] <-
    background_records[[field]][background_index[use_background]]
}
crosswalk$published_annotation_source[use_background] <-
  "global_tree_background_records"

if (any(is.na(crosswalk$published_annotation_source))) {
  stop(
    "One or more tree tips could not be linked to a source annotation row.",
    call. = FALSE
  )
}

crosswalk$is_l10207 <- crosswalk$model_sublineage_id == "L1_02.07"
crosswalk$published_mt28_label <- ifelse(
  !is.na(crosswalk$published_sublineage) &
    grepl("MT28", crosswalk$published_sublineage, fixed = TRUE),
  crosswalk$published_sublineage,
  NA_character_
)
crosswalk$has_published_mt28_label <- !is.na(
  crosswalk$published_mt28_label
)

tree <- ape::read.tree(paths[["tree"]])
if (!setequal(tree$tip.label, crosswalk$tree_sample_id)) {
  stop("Tree tips and lineage crosswalk identifiers differ.", call. = FALSE)
}

mt28_tips <- crosswalk$tree_sample_id[
  crosswalk$has_published_mt28_label
]
target_tips <- crosswalk$tree_sample_id[crosswalk$is_l10207]
if (length(mt28_tips) < 2L) {
  stop("At least two published MT28-labelled tips are required.", call. = FALSE)
}

mt28_mrca_node <- ape::getMRCA(tree, mt28_tips)
mt28_mrca_tips <- ape::extract.clade(
  tree,
  node = mt28_mrca_node
)$tip.label
crosswalk$inside_published_mt28_mrca <- crosswalk$tree_sample_id %in%
  mt28_mrca_tips

n_tree <- nrow(crosswalk)
n_target <- sum(crosswalk$is_l10207)
n_target_focal <- sum(
  crosswalk$is_l10207 & crosswalk$tree_role == "focal"
)
n_target_background <- sum(
  crosswalk$is_l10207 & crosswalk$tree_role == "global_background"
)
n_mt28 <- sum(crosswalk$has_published_mt28_label)
n_mt28_in_target <- sum(
  crosswalk$has_published_mt28_label & crosswalk$is_l10207
)
n_mt28_outside_target <- sum(
  crosswalk$has_published_mt28_label & !crosswalk$is_l10207
)
n_target_mt28_unlabelled <- sum(
  crosswalk$is_l10207 & !crosswalk$has_published_mt28_label
)
n_mrca <- length(mt28_mrca_tips)
n_target_in_mrca <- sum(target_tips %in% mt28_mrca_tips)
n_non_target_in_mrca <- sum(!mt28_mrca_tips %in% target_tips)
n_target_outside_mrca <- sum(!target_tips %in% mt28_mrca_tips)
target_monophyletic <- isTRUE(ape::is.monophyletic(tree, target_tips))
mrca_exact_match <- setequal(mt28_mrca_tips, target_tips)

summary_table <- data.frame(
  metric = c(
    "tree_tips",
    "l10207_tree_tips",
    "l10207_focal_tips",
    "l10207_background_tips",
    "published_mt28_labelled_tips",
    "published_mt28_labels_in_l10207",
    "published_mt28_labels_outside_l10207",
    "l10207_published_mt28_annotation_coverage",
    "published_mt28_label_positive_predictive_value",
    "published_mt28_mrca_descendant_tips",
    "published_mt28_mrca_l10207_recovery",
    "published_mt28_mrca_l10207_precision",
    "l10207_target_monophyletic",
    "published_mt28_mrca_exactly_matches_l10207"
  ),
  numerator = c(
    n_tree,
    n_target,
    n_target_focal,
    n_target_background,
    n_mt28,
    n_mt28_in_target,
    n_mt28_outside_target,
    n_mt28_in_target,
    n_mt28_in_target,
    n_mrca,
    n_target_in_mrca,
    n_target_in_mrca,
    as.integer(target_monophyletic),
    as.integer(mrca_exact_match)
  ),
  denominator = c(
    n_tree,
    n_tree,
    n_target,
    n_target,
    n_tree,
    n_mt28,
    n_mt28,
    n_target,
    n_mt28,
    n_tree,
    n_target,
    n_mrca,
    1L,
    1L
  ),
  percent = c(
    100,
    pct(n_target, n_tree),
    pct(n_target_focal, n_target),
    pct(n_target_background, n_target),
    pct(n_mt28, n_tree),
    pct(n_mt28_in_target, n_mt28),
    pct(n_mt28_outside_target, n_mt28),
    pct(n_mt28_in_target, n_target),
    pct(n_mt28_in_target, n_mt28),
    pct(n_mrca, n_tree),
    pct(n_target_in_mrca, n_target),
    pct(n_target_in_mrca, n_mrca),
    100 * as.integer(target_monophyletic),
    100 * as.integer(mrca_exact_match)
  ),
  interpretation = c(
    "Frozen core-SNP tree denominator.",
    "Study-defined rhierBAPS level-2 population.",
    "Focal members used in genomic observation and case-linked models.",
    "Global-background members used only to stabilize phylogenetic context.",
    "Tips with an MT28 category in the source-published Sublineages field.",
    "All source-published MT28-labelled tips mapped to L1_02.07.",
    "No source-published MT28-labelled tip mapped outside L1_02.07.",
    "Published MT28 labels directly covered a subset of L1_02.07.",
    "A published MT28 label predicted L1_02.07 membership without an observed exception.",
    "Descendant span of the MRCA of all source-published MT28-labelled tips.",
    "The published-MT28 MRCA recovered every L1_02.07 member.",
    "The published-MT28 MRCA contained no non-L1_02.07 tip.",
    "L1_02.07 was monophyletic in the frozen core-SNP tree.",
    "The published-MT28 MRCA descendant set equalled L1_02.07."
  ),
  limitation = c(
    "Conditional on the frozen public-tree cohort.",
    "The identifier is analysis-specific and is retained in machine-readable outputs.",
    "Focal public genomes are not a national prevalence sample.",
    "Background tips are excluded from national observation denominators.",
    "Published labels are unavailable for many tree tips.",
    "Concordance is a nomenclature crosswalk, not an independent expansion test.",
    "The source-published field may be missing even when a tip belongs to the same clade.",
    paste0(
      n_target_mt28_unlabelled,
      " of ",
      n_target,
      " L1_02.07 tips lacked a source-published MT28 sublineage label."
    ),
    "This does not imply that every unlabelled L1_02.07 tip was directly MLVA-typed as MT28.",
    "MRCA correspondence supports a clade-level display name, not imputation of per-tip MLVA type.",
    "Recovery is evaluated in the same frozen tree used to define L1_02.07.",
    "Precision is evaluated in the same frozen tree used to define L1_02.07.",
    "Monophyly does not establish a mechanism for expansion.",
    "The professional alias should remain linked to the analysis identifier in Methods and reproducibility files."
  ),
  stringsAsFactors = FALSE
)

label_levels <- sort(unique(na.omit(crosswalk$published_mt28_label)))
label_counts <- do.call(
  rbind,
  lapply(
    label_levels,
    function(label) {
      keep <- crosswalk$published_mt28_label == label &
        !is.na(crosswalk$published_mt28_label)
      data.frame(
        published_mt28_label = label,
        n_tree_tips = sum(keep),
        n_focal = sum(keep & crosswalk$tree_role == "focal"),
        n_global_background = sum(
          keep & crosswalk$tree_role == "global_background"
        ),
        n_l10207 = sum(keep & crosswalk$is_l10207),
        n_outside_l10207 = sum(keep & !crosswalk$is_l10207),
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(label_counts) <- NULL

crosswalk <- crosswalk[
  order(
    !crosswalk$is_l10207,
    !crosswalk$has_published_mt28_label,
    crosswalk$tree_role,
    crosswalk$tree_sample_id
  ),
  c(
    "tree_sample_id",
    "genome_record_id",
    "tree_role",
    "country_iso3",
    "epidemic_period",
    "model_sublineage_id",
    "primary_model_lineage_id",
    "published_branch",
    "published_lineage",
    "published_sublineage",
    "published_mt28_label",
    "has_published_mt28_label",
    "is_l10207",
    "inside_published_mt28_mrca",
    "published_annotation_source"
  )
]
rownames(crosswalk) <- NULL

validation <- list(
  status = if (
    n_tree == 989L &&
      n_target == 288L &&
      n_target_focal == 271L &&
      n_target_background == 17L &&
      n_mt28 == 99L &&
      n_mt28_in_target == 99L &&
      n_mt28_outside_target == 0L &&
      n_mrca == 288L &&
      n_non_target_in_mrca == 0L &&
      n_target_outside_mrca == 0L &&
      target_monophyletic &&
      mrca_exact_match
  ) {
    "PASS"
  } else {
    "FAIL"
  },
  analysis_identifier = "L1_02.07",
  recommended_display_name = "MT28-associated genomic lineage",
  published_annotation_source = paste(
    "Zhang et al. Evolutionary dynamics and global spread of",
    "macrolide-resistant Bordetella pertussis during the post-pandemic",
    "pertussis resurgence. Journal of Infection. 2026;92:106718."
  ),
  published_annotation_source_doi = "10.1016/j.jinf.2026.106718",
  n_tree_tips = n_tree,
  n_l10207_tree_tips = n_target,
  n_l10207_focal_tips = n_target_focal,
  n_l10207_global_background_tips = n_target_background,
  n_published_mt28_labelled_tips = n_mt28,
  n_published_mt28_labels_in_l10207 = n_mt28_in_target,
  n_published_mt28_labels_outside_l10207 = n_mt28_outside_target,
  n_l10207_without_published_mt28_label = n_target_mt28_unlabelled,
  published_mt28_mrca_node = mt28_mrca_node,
  n_published_mt28_mrca_descendants = n_mrca,
  n_non_l10207_inside_published_mt28_mrca = n_non_target_in_mrca,
  n_l10207_outside_published_mt28_mrca = n_target_outside_mrca,
  l10207_monophyletic = target_monophyletic,
  published_mt28_mrca_exactly_matches_l10207 = mrca_exact_match,
  interpretation = paste(
    "The source-published MT28 labels and their MRCA support the reader-facing",
    "alias 'MT28-associated genomic lineage' for the study-defined cluster.",
    "Machine-readable identifiers remain unchanged."
  ),
  limitation = paste(
    "The crosswalk does not impute MLVA type or resistance status to unlabelled",
    "tips and is not an independent test of lineage expansion."
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(summary_table, "mt28_lineage_nomenclature_summary.tsv")
write_tsv(label_counts, "published_mt28_label_counts.tsv")
write_tsv(crosswalk, "mt28_lineage_tip_crosswalk.tsv")
jsonlite::write_json(
  validation,
  file.path(output_dir, "mt28_lineage_nomenclature_validation.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
cat(
  "MT28 lineage nomenclature validation: ",
  validation$status,
  "\n",
  sep = ""
)
