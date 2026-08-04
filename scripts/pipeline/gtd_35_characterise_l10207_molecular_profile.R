#!/usr/bin/env Rscript

# Descriptive molecular characterisation of the frozen L1_02.07 population.
#
# This analysis deliberately remains conditional on the focal public-tree
# sample. It does not estimate national antigen, PRN, or resistance prevalence.
# Molecular fields were excluded from the upstream lineage definition.

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
  file.path(repo_root, "results", "molecular_characterisation")
)

resolve_archive_path <- function() {
  explicit <- get_arg("--archive", "")
  if (nzchar(explicit)) {
    return(normalizePath(explicit, mustWork = TRUE))
  }

  registry_path <- file.path(repo_root, "data", "derived", "source_registry.tsv")
  registry <- read.delim(
    registry_path,
    check.names = FALSE,
    na.strings = c("", "NA"),
    quote = ""
  )
  archive_rel <- registry$path[registry$source_id == "archive"]
  if (length(archive_rel) != 1L || is.na(archive_rel)) {
    stop("Could not identify the frozen archive in source_registry.tsv", call. = FALSE)
  }

  data_roots <- unique(c(
    Sys.getenv("PERTUSSIS_DATA_ROOT", unset = ""),
    repo_root,
    dirname(repo_root),
    dirname(dirname(repo_root))
  ))
  data_roots <- data_roots[nzchar(data_roots)]
  candidates <- file.path(data_roots, archive_rel)
  hits <- candidates[file.exists(candidates)]
  if (!length(hits)) {
    stop(
      "Frozen archive not found. Pass --archive /absolute/path/frozen_archive_isolates.tsv",
      call. = FALSE
    )
  }
  normalizePath(hits[[1L]], mustWork = TRUE)
}

paths <- c(
  assignments = file.path(
    repo_root,
    "results",
    "lineages",
    "primary_finalized",
    "model_lineage_assignments.tsv"
  ),
  tree_metadata = file.path(repo_root, "results", "phylogeny", "tree_tip_metadata.tsv"),
  source_records = file.path(
    repo_root,
    "data",
    "derived",
    "transmission_genome_records.tsv"
  ),
  frozen_archive = resolve_archive_path()
)

if (any(!file.exists(paths))) {
  stop(
    "Missing required inputs:\n",
    paste(paths[!file.exists(paths)], collapse = "\n"),
    call. = FALSE
  )
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

write_tsv <- function(frame, filename, na_value = "") {
  write.table(
    frame,
    file.path(output_dir, filename),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = na_value
  )
}

pct <- function(numerator, denominator) {
  ifelse(denominator > 0, 100 * numerator / denominator, NA_real_)
}

assignments <- read_tsv(paths[["assignments"]])
tree_metadata <- read_tsv(paths[["tree_metadata"]])
source_records <- read_tsv(paths[["source_records"]])
frozen <- read_tsv(paths[["frozen_archive"]])

stopifnot(
  !anyDuplicated(assignments$tree_sample_id),
  !anyDuplicated(tree_metadata$tree_sample_id),
  !anyDuplicated(source_records$genome_record_id),
  !anyDuplicated(frozen$sample_id_canonical),
  !anyDuplicated(frozen$biosample_accession[!is.na(frozen$biosample_accession)]),
  !anyDuplicated(frozen$assembly_accession[!is.na(frozen$assembly_accession)])
)

tree_identifiers <- tree_metadata[
  ,
  c(
    "tree_sample_id",
    "genome_record_id",
    "sample_id",
    "biosample_accession",
    "assembly_accession",
    "run_accession"
  )
]

source_annotations <- source_records[
  ,
  c(
    "genome_record_id",
    "preliminary_lineage_id",
    "published_sublineage",
    "ptxP_label",
    "fim3_label",
    "marker_23s_status"
  )
]
names(source_annotations)[names(source_annotations) == "preliminary_lineage_id"] <-
  "preliminary_lineage_id_source"
names(source_annotations)[names(source_annotations) == "marker_23s_status"] <-
  "marker_23s_status_source"

analysis <- merge(
  assignments,
  tree_identifiers,
  by = "tree_sample_id",
  all.x = TRUE,
  sort = FALSE
)
analysis <- merge(
  analysis,
  source_annotations,
  by = "genome_record_id",
  all.x = TRUE,
  sort = FALSE
)

# Frozen annotations are matched hierarchically using stable public identifiers.
frozen_by_sample <- match(analysis$sample_id, frozen$sample_id_canonical)
frozen_by_biosample <- match(analysis$biosample_accession, frozen$biosample_accession)
frozen_by_assembly <- match(analysis$assembly_accession, frozen$assembly_accession)

identifier_conflict <- vapply(
  seq_len(nrow(analysis)),
  function(index) {
    candidates <- unique(na.omit(c(
      frozen_by_sample[[index]],
      frozen_by_biosample[[index]],
      frozen_by_assembly[[index]]
    )))
    length(candidates) > 1L
  },
  logical(1)
)
if (any(identifier_conflict)) {
  stop("Stable identifiers mapped to conflicting frozen archive rows", call. = FALSE)
}

frozen_index <- frozen_by_sample
match_method <- ifelse(!is.na(frozen_index), "sample_id", NA_character_)
use_biosample <- is.na(frozen_index) & !is.na(frozen_by_biosample)
frozen_index[use_biosample] <- frozen_by_biosample[use_biosample]
match_method[use_biosample] <- "biosample_accession"
use_assembly <- is.na(frozen_index) & !is.na(frozen_by_assembly)
frozen_index[use_assembly] <- frozen_by_assembly[use_assembly]
match_method[use_assembly] <- "assembly_accession"

analysis$frozen_match_method <- match_method
analysis$frozen_archive_match <- !is.na(frozen_index)
analysis$frozen_prn_outcome <- frozen$frozen_outcome_status[frozen_index]
analysis$frozen_prn_mechanism <- frozen$frozen_prn_mechanism_call[frozen_index]
analysis$frozen_23S_status <- frozen$marker_23s_status[frozen_index]
analysis$frozen_prn_primary_eligible <-
  frozen$frozen_primary_outcome_eligible[frozen_index] == "true"

analysis <- analysis[analysis$tree_role == "focal", , drop = FALSE]
if (nrow(analysis) != 774L) {
  stop("Expected 774 focal tree tips; observed ", nrow(analysis), call. = FALSE)
}
if (any(is.na(analysis$genome_record_id))) {
  stop("One or more focal lineage assignments lack tree metadata", call. = FALSE)
}
if (any(!analysis$genome_record_id %in% source_records$genome_record_id)) {
  stop("One or more focal tips could not be recovered in source records", call. = FALSE)
}

analysis$lineage_group <- ifelse(
  analysis$primary_model_lineage_id == "L1_02.07",
  "L1_02.07",
  "Other_model_lineages"
)
analysis$source_antigen_profile <- ifelse(
  !is.na(analysis$ptxP_label) & !is.na(analysis$fim3_label),
  paste(analysis$ptxP_label, analysis$fim3_label, sep = "/"),
  NA_character_
)
analysis$published_MT_label <- ifelse(
  !is.na(analysis$published_sublineage) &
    grepl("MT28", analysis$published_sublineage, fixed = TRUE),
  analysis$published_sublineage,
  NA_character_
)

period_order <- c("prepandemic", "pandemic", "resurgence")
analysis$epidemic_period <- factor(
  analysis$epidemic_period,
  levels = period_order
)
analysis <- analysis[
  order(
    analysis$country_iso3,
    analysis$epidemic_period,
    analysis$lineage_group,
    analysis$tree_sample_id
  ),
  ,
  drop = FALSE
]

# Tree counts make the denominator for every molecular summary explicit.
count_by_group <- function(frame) {
  data.frame(
    n_tree_tips = nrow(frame),
    n_l10207 = sum(frame$lineage_group == "L1_02.07"),
    n_other_lineages = sum(frame$lineage_group == "Other_model_lineages"),
    l10207_percent = pct(
      sum(frame$lineage_group == "L1_02.07"),
      nrow(frame)
    )
  )
}

tree_share <- do.call(
  rbind,
  lapply(
    split(
      analysis,
      interaction(
        analysis$country_iso3,
        analysis$epidemic_period,
        drop = TRUE,
        lex.order = TRUE
      )
    ),
    function(frame) {
      cbind(
        country_iso3 = frame$country_iso3[[1L]],
        epidemic_period = as.character(frame$epidemic_period[[1L]]),
        count_by_group(frame)
      )
    }
  )
)
rownames(tree_share) <- NULL

# Annotation coverage is always reported against all tree tips in the
# country-period-lineage cell; non-interpretable PRN calls are not silently
# pooled with interpretable calls.
coverage_for_group <- function(frame) {
  n_total <- nrow(frame)
  n_frozen <- sum(frame$frozen_archive_match)
  n_prn_interpretable <- sum(
    frame$frozen_prn_outcome %in% c("intact", "disrupted"),
    na.rm = TRUE
  )
  n_prn_non_interpretable <- sum(
    frame$frozen_prn_outcome == "non_interpretable",
    na.rm = TRUE
  )
  data.frame(
    n_tree_tips = n_total,
    n_source_ptxP_annotated = sum(!is.na(frame$ptxP_label)),
    pct_source_ptxP_annotated = pct(sum(!is.na(frame$ptxP_label)), n_total),
    n_source_fim3_annotated = sum(!is.na(frame$fim3_label)),
    pct_source_fim3_annotated = pct(sum(!is.na(frame$fim3_label)), n_total),
    n_source_antigen_profile_annotated = sum(!is.na(frame$source_antigen_profile)),
    pct_source_antigen_profile_annotated = pct(
      sum(!is.na(frame$source_antigen_profile)),
      n_total
    ),
    n_source_23S_raw_annotated = sum(!is.na(frame$marker_23s_status_source)),
    pct_source_23S_raw_annotated = pct(
      sum(!is.na(frame$marker_23s_status_source)),
      n_total
    ),
    n_frozen_archive_matched = n_frozen,
    pct_frozen_archive_matched = pct(n_frozen, n_total),
    n_frozen_23S_annotated = sum(!is.na(frame$frozen_23S_status)),
    pct_frozen_23S_annotated = pct(
      sum(!is.na(frame$frozen_23S_status)),
      n_total
    ),
    n_frozen_PRN_interpretable = n_prn_interpretable,
    pct_frozen_PRN_interpretable = pct(n_prn_interpretable, n_total),
    n_frozen_PRN_non_interpretable = n_prn_non_interpretable,
    pct_frozen_PRN_non_interpretable = pct(
      n_prn_non_interpretable,
      n_total
    ),
    n_frozen_PRN_not_matched = n_total - n_frozen,
    n_published_MT28_annotated = sum(!is.na(frame$published_MT_label)),
    pct_published_MT28_annotated = pct(
      sum(!is.na(frame$published_MT_label)),
      n_total
    )
  )
}

coverage_groups <- split(
  analysis,
  interaction(
    analysis$country_iso3,
    analysis$epidemic_period,
    analysis$lineage_group,
    drop = TRUE,
    lex.order = TRUE
  )
)
coverage <- do.call(
  rbind,
  lapply(
    coverage_groups,
    function(frame) {
      cbind(
        country_iso3 = frame$country_iso3[[1L]],
        epidemic_period = as.character(frame$epidemic_period[[1L]]),
        lineage_group = frame$lineage_group[[1L]],
        coverage_for_group(frame)
      )
    }
  )
)
rownames(coverage) <- NULL

# Long-form exact category counts preserve source labels instead of imposing an
# unsupported cross-study harmonisation of 23S allele, genotype, and phenotype.
feature_vectors <- list(
  source_ptxP = analysis$ptxP_label,
  source_fim3 = analysis$fim3_label,
  source_antigen_profile = analysis$source_antigen_profile,
  source_23S_raw = analysis$marker_23s_status_source,
  frozen_23S_status = analysis$frozen_23S_status,
  frozen_PRN_outcome = analysis$frozen_prn_outcome,
  frozen_PRN_mechanism = analysis$frozen_prn_mechanism,
  published_MT_label = analysis$published_MT_label
)

profile_rows <- list()
row_index <- 0L
profile_keys <- unique(
  analysis[, c("country_iso3", "epidemic_period", "lineage_group")]
)
for (key_index in seq_len(nrow(profile_keys))) {
  key <- profile_keys[key_index, , drop = FALSE]
  keep <- analysis$country_iso3 == key$country_iso3 &
    analysis$epidemic_period == key$epidemic_period &
    analysis$lineage_group == key$lineage_group
  n_total <- sum(keep)

  for (feature_name in names(feature_vectors)) {
    values <- feature_vectors[[feature_name]][keep]
    values_with_missing <- ifelse(is.na(values), "MISSING", values)
    counts <- sort(table(values_with_missing), decreasing = TRUE)
    n_annotated <- sum(!is.na(values))
    n_interpretable <- n_annotated
    if (feature_name == "frozen_PRN_outcome") {
      n_interpretable <- sum(values %in% c("intact", "disrupted"), na.rm = TRUE)
    } else if (feature_name == "frozen_PRN_mechanism") {
      n_interpretable <- sum(
        values %in% c(
          "intact",
          "coding_disrupted_is481",
          "coding_disrupted_inversion_or_rearrangement",
          "coding_disrupted_other"
        ),
        na.rm = TRUE
      )
    }

    for (category in names(counts)) {
      row_index <- row_index + 1L
      n_category <- as.integer(counts[[category]])
      profile_rows[[row_index]] <- data.frame(
        country_iso3 = key$country_iso3,
        epidemic_period = as.character(key$epidemic_period),
        lineage_group = key$lineage_group,
        feature = feature_name,
        category = category,
        n_category = n_category,
        n_feature_annotated = n_annotated,
        n_feature_interpretable = n_interpretable,
        n_feature_missing = n_total - n_annotated,
        n_lineage_tree_total = n_total,
        category_percent_of_annotated = ifelse(
          category == "MISSING",
          NA_real_,
          pct(n_category, n_annotated)
        ),
        category_percent_of_tree_total = pct(n_category, n_total)
      )
    }
  }
}
profile_counts <- do.call(rbind, profile_rows)

# Project-level annotation availability exposes study-confounded missingness.
project_groups <- split(
  analysis,
  interaction(
    analysis$country_iso3,
    analysis$epidemic_period,
    analysis$project_id,
    drop = TRUE,
    lex.order = TRUE
  )
)
project_coverage <- do.call(
  rbind,
  lapply(
    project_groups,
    function(frame) {
      data.frame(
        country_iso3 = frame$country_iso3[[1L]],
        epidemic_period = as.character(frame$epidemic_period[[1L]]),
        project_id = frame$project_id[[1L]],
        n_tree_tips = nrow(frame),
        n_l10207 = sum(frame$lineage_group == "L1_02.07"),
        pct_l10207 = pct(
          sum(frame$lineage_group == "L1_02.07"),
          nrow(frame)
        ),
        n_source_antigen_profile_annotated = sum(
          !is.na(frame$source_antigen_profile)
        ),
        pct_source_antigen_profile_annotated = pct(
          sum(!is.na(frame$source_antigen_profile)),
          nrow(frame)
        ),
        n_frozen_23S_annotated = sum(!is.na(frame$frozen_23S_status)),
        pct_frozen_23S_annotated = pct(
          sum(!is.na(frame$frozen_23S_status)),
          nrow(frame)
        ),
        n_frozen_PRN_interpretable = sum(
          frame$frozen_prn_outcome %in% c("intact", "disrupted"),
          na.rm = TRUE
        ),
        pct_frozen_PRN_interpretable = pct(
          sum(
            frame$frozen_prn_outcome %in% c("intact", "disrupted"),
            na.rm = TRUE
          ),
          nrow(frame)
        )
      )
    }
  )
)
rownames(project_coverage) <- NULL

# Period-shift checks ask whether the L1_02.07 share changes from the
# prepandemic to resurgence periods within an observed molecular stratum.
# Fisher estimates are emitted only when each period contributes >=5 records
# and both lineage margins are non-zero.
shift_strata <- list(
  c("source_ptxP", "ptxP3"),
  c("source_ptxP", "ptxP1"),
  c("source_fim3", "fim3-1"),
  c("source_fim3", "fim3-2"),
  c("source_antigen_profile", "ptxP3/fim3-1"),
  c("source_antigen_profile", "ptxP1/fim3-1"),
  c("source_antigen_profile", "ptxP3/fim3-2"),
  c("frozen_23S_status", "MR_A2047G"),
  c("frozen_23S_status", "MS"),
  c("frozen_PRN_outcome", "intact"),
  c("frozen_PRN_outcome", "disrupted"),
  c("published_MT_label", "MS-MT28"),
  c("published_MT_label", "MR-MT28-PG1")
)

shift_source <- list(
  source_ptxP = analysis$ptxP_label,
  source_fim3 = analysis$fim3_label,
  source_antigen_profile = analysis$source_antigen_profile,
  frozen_23S_status = analysis$frozen_23S_status,
  frozen_PRN_outcome = analysis$frozen_prn_outcome,
  published_MT_label = analysis$published_MT_label
)

fisher_period_shift <- function(pre_l1, pre_other, res_l1, res_other) {
  contingency <- matrix(
    c(res_l1, res_other, pre_l1, pre_other),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(
      period = c("resurgence", "prepandemic"),
      lineage = c("L1_02.07", "Other")
    )
  )
  fit <- fisher.test(contingency)
  c(
    odds_ratio = unname(fit$estimate),
    conf_low = fit$conf.int[[1L]],
    conf_high = fit$conf.int[[2L]],
    p_value = fit$p.value
  )
}

shift_rows <- list()
shift_index <- 0L
countries <- sort(unique(analysis$country_iso3))
for (country in countries) {
  for (stratum in shift_strata) {
    feature_name <- stratum[[1L]]
    category <- stratum[[2L]]
    values <- shift_source[[feature_name]]
    in_stratum <- analysis$country_iso3 == country &
      !is.na(values) &
      values == category &
      analysis$epidemic_period %in% c("prepandemic", "resurgence")
    frame <- analysis[in_stratum, , drop = FALSE]

    pre_l1 <- sum(
      frame$epidemic_period == "prepandemic" &
        frame$lineage_group == "L1_02.07"
    )
    pre_other <- sum(
      frame$epidemic_period == "prepandemic" &
        frame$lineage_group == "Other_model_lineages"
    )
    res_l1 <- sum(
      frame$epidemic_period == "resurgence" &
        frame$lineage_group == "L1_02.07"
    )
    res_other <- sum(
      frame$epidemic_period == "resurgence" &
        frame$lineage_group == "Other_model_lineages"
    )
    pre_total <- pre_l1 + pre_other
    res_total <- res_l1 + res_other

    pre_projects <- unique(frame$project_id[frame$epidemic_period == "prepandemic"])
    res_projects <- unique(frame$project_id[frame$epidemic_period == "resurgence"])
    shared_projects <- intersect(pre_projects, res_projects)
    shared_min5 <- sum(vapply(
      shared_projects,
      function(project) {
        project_frame <- frame[frame$project_id == project, , drop = FALSE]
        all(table(
          factor(
            project_frame$epidemic_period,
            levels = c("prepandemic", "resurgence")
          )
        ) >= 5L)
      },
      logical(1)
    ))

    status <- "fisher_exact_prespecified"
    estimate <- c(
      odds_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    )
    if (pre_total == 0L || res_total == 0L) {
      status <- "not_estimable_period_absent"
    } else if (pre_total < 5L || res_total < 5L) {
      status <- "descriptive_only_sparse_period"
    } else if ((pre_l1 + res_l1) == 0L || (pre_other + res_other) == 0L) {
      status <- "not_estimable_zero_lineage_margin"
    } else {
      estimate <- fisher_period_shift(pre_l1, pre_other, res_l1, res_other)
    }

    project_overlap_note <- if (length(shared_projects) == 0L) {
      "no project represented in both periods"
    } else if (shared_min5 == 0L) {
      "shared project exists, but none has >=5 annotated records in each period"
    } else {
      paste0(
        shared_min5,
        " shared project(s) have >=5 annotated records in each period"
      )
    }

    shift_index <- shift_index + 1L
    shift_rows[[shift_index]] <- data.frame(
      country_iso3 = country,
      feature = feature_name,
      molecular_stratum = category,
      prepandemic_l10207 = pre_l1,
      prepandemic_other = pre_other,
      prepandemic_total = pre_total,
      resurgence_l10207 = res_l1,
      resurgence_other = res_other,
      resurgence_total = res_total,
      prepandemic_l10207_percent = pct(pre_l1, pre_total),
      resurgence_l10207_percent = pct(res_l1, res_total),
      percentage_point_change = pct(res_l1, res_total) - pct(pre_l1, pre_total),
      fisher_odds_ratio_resurgence_vs_prepandemic = estimate[["odds_ratio"]],
      fisher_conf_low = estimate[["conf_low"]],
      fisher_conf_high = estimate[["conf_high"]],
      fisher_p_value = estimate[["p_value"]],
      analysis_status = status,
      n_prepandemic_projects = length(pre_projects),
      n_resurgence_projects = length(res_projects),
      n_projects_in_both_periods = length(shared_projects),
      n_shared_projects_min5_each_period = shared_min5,
      project_overlap_note = project_overlap_note
    )
  }
}
period_shift <- do.call(rbind, shift_rows)
period_shift$fisher_p_bh_exploratory <- NA_real_
valid_shift <- period_shift$analysis_status == "fisher_exact_prespecified"
period_shift$fisher_p_bh_exploratory[valid_shift] <- p.adjust(
  period_shift$fisher_p_value[valid_shift],
  method = "BH"
)

# Within-country/period comparisons describe how observed features are
# distributed between L1_02.07 and all other frozen model lineages.
contrast_definitions <- list(
  list(
    contrast = "ptxP3_vs_other_ptxP",
    values = analysis$ptxP_label,
    exposed = "ptxP3",
    reference = setdiff(unique(na.omit(analysis$ptxP_label)), "ptxP3")
  ),
  list(
    contrast = "fim3-1_vs_other_fim3",
    values = analysis$fim3_label,
    exposed = "fim3-1",
    reference = setdiff(unique(na.omit(analysis$fim3_label)), "fim3-1")
  ),
  list(
    contrast = "MR_A2047G_vs_MS",
    values = analysis$frozen_23S_status,
    exposed = "MR_A2047G",
    reference = "MS"
  ),
  list(
    contrast = "PRN_disrupted_vs_intact",
    values = analysis$frozen_prn_outcome,
    exposed = "disrupted",
    reference = "intact"
  )
)

contrast_rows <- list()
contrast_index <- 0L
for (country in countries) {
  for (period in period_order) {
    for (definition in contrast_definitions) {
      values <- definition$values
      in_cell <- analysis$country_iso3 == country &
        analysis$epidemic_period == period
      in_contrast <- in_cell &
        !is.na(values) &
        values %in% c(definition$exposed, definition$reference)
      frame <- analysis[in_contrast, , drop = FALSE]
      frame_values <- values[in_contrast]

      l1_exposed <- sum(
        frame$lineage_group == "L1_02.07" &
          frame_values == definition$exposed
      )
      l1_reference <- sum(
        frame$lineage_group == "L1_02.07" &
          frame_values %in% definition$reference
      )
      other_exposed <- sum(
        frame$lineage_group == "Other_model_lineages" &
          frame_values == definition$exposed
      )
      other_reference <- sum(
        frame$lineage_group == "Other_model_lineages" &
          frame_values %in% definition$reference
      )

      l1_total <- l1_exposed + l1_reference
      other_total <- other_exposed + other_reference
      exposed_total <- l1_exposed + other_exposed
      reference_total <- l1_reference + other_reference
      status <- "fisher_exact_prespecified"
      estimate <- c(
        odds_ratio = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = NA_real_
      )
      if (l1_total == 0L || other_total == 0L) {
        status <- "not_estimable_lineage_group_absent"
      } else if (l1_total < 5L || other_total < 5L) {
        status <- "descriptive_only_sparse_lineage_group"
      } else if (exposed_total == 0L || reference_total == 0L) {
        status <- "not_estimable_zero_feature_margin"
      } else {
        contingency <- matrix(
          c(l1_exposed, l1_reference, other_exposed, other_reference),
          nrow = 2L,
          byrow = TRUE,
          dimnames = list(
            lineage = c("L1_02.07", "Other"),
            feature = c("exposed", "reference")
          )
        )
        fit <- fisher.test(contingency)
        estimate <- c(
          odds_ratio = unname(fit$estimate),
          conf_low = fit$conf.int[[1L]],
          conf_high = fit$conf.int[[2L]],
          p_value = fit$p.value
        )
      }

      contrast_index <- contrast_index + 1L
      contrast_rows[[contrast_index]] <- data.frame(
        country_iso3 = country,
        epidemic_period = period,
        contrast = definition$contrast,
        exposed_category = definition$exposed,
        reference_categories = paste(definition$reference, collapse = ";"),
        l10207_exposed = l1_exposed,
        l10207_reference = l1_reference,
        other_exposed = other_exposed,
        other_reference = other_reference,
        l10207_exposed_percent = pct(l1_exposed, l1_total),
        other_exposed_percent = pct(other_exposed, other_total),
        fisher_odds_ratio_l10207_vs_other = estimate[["odds_ratio"]],
        fisher_conf_low = estimate[["conf_low"]],
        fisher_conf_high = estimate[["conf_high"]],
        fisher_p_value = estimate[["p_value"]],
        analysis_status = status
      )
    }
  }
}
feature_contrasts <- do.call(rbind, contrast_rows)
feature_contrasts$fisher_p_bh_exploratory <- NA_real_
valid_contrast <- feature_contrasts$analysis_status == "fisher_exact_prespecified"
feature_contrasts$fisher_p_bh_exploratory[valid_contrast] <- p.adjust(
  feature_contrasts$fisher_p_value[valid_contrast],
  method = "BH"
)

get_shift <- function(country, feature, stratum) {
  hit <- period_shift[
    period_shift$country_iso3 == country &
      period_shift$feature == feature &
      period_shift$molecular_stratum == stratum,
    ,
    drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop("Expected one period-shift result for ", country, "/", feature, "/", stratum)
  }
  hit
}

l1 <- analysis[analysis$lineage_group == "L1_02.07", , drop = FALSE]
l1_ptxp <- l1[!is.na(l1$ptxP_label), , drop = FALSE]
l1_antigen <- l1[!is.na(l1$source_antigen_profile), , drop = FALSE]
l1_23s <- l1[!is.na(l1$frozen_23S_status), , drop = FALSE]
l1_prn <- l1[
  l1$frozen_prn_outcome %in% c("intact", "disrupted"),
  ,
  drop = FALSE
]

key_rows <- list(
  data.frame(
    finding_id = "focal_tree_l10207",
    scope = "all focal tree tips",
    numerator = nrow(l1),
    denominator = nrow(analysis),
    percent = pct(nrow(l1), nrow(analysis)),
    estimate = paste0(nrow(l1), "/", nrow(analysis)),
    interpretation = "Frozen tree membership of L1_02.07",
    limitation = "Conditional on selected public-tree genomes"
  ),
  data.frame(
    finding_id = "l10207_ptxP3_among_typed",
    scope = "L1_02.07 with source ptxP annotation",
    numerator = sum(l1_ptxp$ptxP_label == "ptxP3"),
    denominator = nrow(l1_ptxp),
    percent = pct(sum(l1_ptxp$ptxP_label == "ptxP3"), nrow(l1_ptxp)),
    estimate = paste0(
      sum(l1_ptxp$ptxP_label == "ptxP3"),
      "/",
      nrow(l1_ptxp)
    ),
    interpretation = "Every ptxP-typed L1_02.07 tip was ptxP3",
    limitation = paste0(
      nrow(l1) - nrow(l1_ptxp),
      " of ",
      nrow(l1),
      " L1_02.07 tips lacked ptxP annotation"
    )
  ),
  data.frame(
    finding_id = "l10207_ptxP3_fim3_1_among_jointly_typed",
    scope = "L1_02.07 with both source antigen annotations",
    numerator = sum(l1_antigen$source_antigen_profile == "ptxP3/fim3-1"),
    denominator = nrow(l1_antigen),
    percent = pct(
      sum(l1_antigen$source_antigen_profile == "ptxP3/fim3-1"),
      nrow(l1_antigen)
    ),
    estimate = paste0(
      sum(l1_antigen$source_antigen_profile == "ptxP3/fim3-1"),
      "/",
      nrow(l1_antigen)
    ),
    interpretation = "Every jointly typed L1_02.07 tip was ptxP3/fim3-1",
    limitation = paste0(
      nrow(l1) - nrow(l1_antigen),
      " of ",
      nrow(l1),
      " L1_02.07 tips lacked a complete ptxP/fim3 annotation"
    )
  ),
  data.frame(
    finding_id = "l10207_frozen_23S_MR_A2047G",
    scope = "L1_02.07 with frozen standardised 23S status",
    numerator = sum(l1_23s$frozen_23S_status == "MR_A2047G"),
    denominator = nrow(l1_23s),
    percent = pct(
      sum(l1_23s$frozen_23S_status == "MR_A2047G"),
      nrow(l1_23s)
    ),
    estimate = paste0(
      sum(l1_23s$frozen_23S_status == "MR_A2047G"),
      " MR_A2047G; ",
      sum(l1_23s$frozen_23S_status == "MS"),
      " MS"
    ),
    interpretation = "Observed L1_02.07 tips included both 23S resistance-marker categories",
    limitation = paste0(
      nrow(l1) - nrow(l1_23s),
      " of ",
      nrow(l1),
      " L1_02.07 tips lacked frozen 23S status"
    )
  ),
  data.frame(
    finding_id = "l10207_frozen_PRN_disrupted",
    scope = "L1_02.07 with interpretable frozen PRN structure",
    numerator = sum(l1_prn$frozen_prn_outcome == "disrupted"),
    denominator = nrow(l1_prn),
    percent = pct(
      sum(l1_prn$frozen_prn_outcome == "disrupted"),
      nrow(l1_prn)
    ),
    estimate = paste0(
      sum(l1_prn$frozen_prn_outcome == "disrupted"),
      " disrupted; ",
      sum(l1_prn$frozen_prn_outcome == "intact"),
      " intact"
    ),
    interpretation = "L1_02.07 was not restricted to one frozen PRN structural state",
    limitation = paste0(
      nrow(l1) - nrow(l1_prn),
      " of ",
      nrow(l1),
      " L1_02.07 tips lacked an interpretable frozen PRN endpoint"
    )
  ),
  data.frame(
    finding_id = "l10207_frozen_PRN_disruption_mechanisms",
    scope = "structurally disrupted L1_02.07 tips",
    numerator = sum(l1_prn$frozen_prn_outcome == "disrupted"),
    denominator = sum(l1_prn$frozen_prn_outcome == "disrupted"),
    percent = 100,
    estimate = paste0(
      sum(
        l1$frozen_prn_mechanism == "coding_disrupted_is481",
        na.rm = TRUE
      ),
      " IS481; ",
      sum(
        l1$frozen_prn_mechanism ==
          "coding_disrupted_inversion_or_rearrangement",
        na.rm = TRUE
      ),
      " inversion/rearrangement; ",
      sum(
        l1$frozen_prn_mechanism == "coding_disrupted_other",
        na.rm = TRUE
      ),
      " other"
    ),
    interpretation = "Observed structurally disrupted L1_02.07 tips represented multiple PRN mechanisms",
    limitation = "Mechanism counts are conditional on frozen interpretable genomic calls, not protein phenotype"
  ),
  data.frame(
    finding_id = "published_MT28_labels_in_l10207",
    scope = "focal tips with explicit published MT28 labels",
    numerator = sum(
      !is.na(analysis$published_MT_label) &
        analysis$lineage_group == "L1_02.07"
    ),
    denominator = sum(!is.na(analysis$published_MT_label)),
    percent = pct(
      sum(
        !is.na(analysis$published_MT_label) &
          analysis$lineage_group == "L1_02.07"
      ),
      sum(!is.na(analysis$published_MT_label))
    ),
    estimate = paste0(
      sum(l1$published_MT_label == "MS-MT28", na.rm = TRUE),
      " MS-MT28; ",
      sum(l1$published_MT_label == "MR-MT28-others", na.rm = TRUE),
      " MR-MT28-others; ",
      sum(l1$published_MT_label == "MR-MT28-PG1", na.rm = TRUE),
      " MR-MT28-PG1; ",
      sum(l1$published_MT_label == "MR-MT28-PG2", na.rm = TRUE),
      " MR-MT28-PG2"
    ),
    interpretation = "All focal tips carrying a source-published MT28 sublineage label mapped to L1_02.07",
    limitation = paste0(
      "Source-published MT28 sublineage labels were available for ",
      sum(!is.na(l1$published_MT_label)),
      " of ",
      nrow(l1),
      " focal L1_02.07 tips; absence was retained as missing and does not rule out clade membership"
    )
  )
)

shift_finding <- function(
  finding_id,
  country,
  feature,
  stratum,
  interpretation,
  limitation
) {
  row <- get_shift(country, feature, stratum)
  data.frame(
    finding_id = finding_id,
    scope = paste(country, feature, stratum, sep = " / "),
    numerator = row$resurgence_l10207,
    denominator = row$resurgence_total,
    percent = row$resurgence_l10207_percent,
    estimate = paste0(
      "pre ",
      row$prepandemic_l10207,
      "/",
      row$prepandemic_total,
      "; resurgence ",
      row$resurgence_l10207,
      "/",
      row$resurgence_total,
      "; Fisher OR ",
      signif(row$fisher_odds_ratio_resurgence_vs_prepandemic, 5),
      ", p=",
      format(row$fisher_p_value, digits = 4, scientific = TRUE)
    ),
    interpretation = interpretation,
    limitation = paste0(limitation, "; ", row$project_overlap_note)
  )
}

key_rows <- c(
  key_rows,
  list(
    shift_finding(
      "CHN_shift_within_MR_A2047G",
      "CHN",
      "frozen_23S_status",
      "MR_A2047G",
      "The Chinese L1_02.07 period shift persisted within MR_A2047G-annotated tips",
      "Within-annotated-sample check; resurgence 23S coverage was incomplete"
    ),
    shift_finding(
      "CHN_shift_within_MS",
      "CHN",
      "frozen_23S_status",
      "MS",
      "The Chinese L1_02.07 period shift also persisted within MS-annotated tips",
      "Within-annotated-sample check; resurgence 23S coverage was incomplete"
    ),
    shift_finding(
      "JPN_shift_within_MS",
      "JPN",
      "frozen_23S_status",
      "MS",
      "The Japanese L1_02.07 period shift was present among MS-annotated tips",
      paste0(
        "Only ",
        sum(
          analysis$country_iso3 == "JPN" &
            analysis$epidemic_period == "resurgence" &
            !is.na(analysis$frozen_23S_status)
        ),
        " of ",
        sum(
          analysis$country_iso3 == "JPN" &
            analysis$epidemic_period == "resurgence"
        ),
        " Japanese resurgence tips had frozen 23S status"
      )
    ),
    shift_finding(
      "CHN_shift_within_PRN_intact",
      "CHN",
      "frozen_PRN_outcome",
      "intact",
      "The Chinese L1_02.07 period shift persisted among frozen-PRN-intact tips",
      "Only interpretable frozen PRN calls were included"
    ),
    shift_finding(
      "JPN_shift_within_PRN_intact",
      "JPN",
      "frozen_PRN_outcome",
      "intact",
      "The Japanese L1_02.07 period shift persisted among frozen-PRN-intact tips",
      "Only interpretable frozen PRN calls were included"
    ),
    shift_finding(
      "CHN_shift_within_ptxP3_fim3_1",
      "CHN",
      "source_antigen_profile",
      "ptxP3/fim3-1",
      "The Chinese L1_02.07 period shift persisted within the observed ptxP3/fim3-1 stratum",
      "Only 46 of 153 Chinese resurgence tips had both antigen annotations"
    )
  )
)
key_findings <- do.call(rbind, key_rows)

# Machine-readable stop rules explain why this layer is descriptive.
limitations <- data.frame(
  issue_id = c(
    "public_tree_conditioning",
    "molecular_missingness",
    "project_confounding",
    "heterogeneous_23S_semantics",
    "PRN_structure_not_protein",
    "no_causal_mechanism_test",
    "multiple_exploratory_tests",
    "published_sublineage_crosswalk"
  ),
  assessment = c(
    "All estimates condition on 774 selected focal public-tree genomes.",
    "Molecular fields are incomplete, especially for recent China and Japan; every table retains full tree denominators and explicit missing counts.",
    "Annotation completeness is project-specific. For the China/Japan strata supporting the L1_02.07 expansion, no shared project had at least five annotated records in both prepandemic and resurgence periods.",
    "Source 23S labels mix allele names, marker calls, and reported phenotype categories; raw source labels are not harmonised. Primary resistance-stratum checks use only frozen MS and MR_A2047G labels.",
    "Frozen PRN calls describe genomic structural state, not measured pertactin protein expression.",
    "Lineages were defined without molecular markers, but descriptive associations do not identify a causal fitness, immune-escape, or resistance mechanism.",
    "Fisher tests are exploratory. Raw and Benjamini-Hochberg adjusted p values are supplied without treating them as confirmatory.",
    "The MT28 sublineage field was imported from the source-published 8,117-genome compilation and was unavailable for many tips."
  ),
  consequence = c(
    "Do not describe percentages as national or circulating-population prevalence.",
    "Interpret within-stratum changes as robustness checks among annotated records, not missing-at-random estimates.",
    "Project-independent molecular period effects are not identifiable from these data.",
    "Do not collapse 23S_rRNA_1, MS, Susceptible, or mutation-detected labels into a single phenotype.",
    "Use 'structurally disrupted PRN', not 'PRN-negative', unless protein phenotype is measured.",
    "Use the layer to show molecular heterogeneity and rule out confinement to one observed stratum.",
    "Focus on effect sizes, denominators, and consistency rather than thresholded significance.",
    "Use MT28-associated as a clade-level nomenclature crosswalk; do not impute per-tip MLVA type or treat the crosswalk as an independent expansion test."
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv(tree_share, "l10207_tree_share_by_country_period.tsv")
write_tsv(coverage, "molecular_annotation_coverage.tsv")
write_tsv(profile_counts, "molecular_profile_counts.tsv")
write_tsv(project_coverage, "annotation_coverage_by_project.tsv")
write_tsv(
  period_shift,
  "l10207_period_shift_within_molecular_strata.tsv",
  na_value = "NA"
)
write_tsv(feature_contrasts, "l10207_feature_contrasts.tsv")
write_tsv(key_findings, "key_findings.tsv")
write_tsv(limitations, "interpretation_limits.tsv")

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (!length(output)) {
    return(NA_character_)
  }
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

input_manifest <- data.frame(
  input_id = names(paths),
  path = unname(paths),
  sha256 = vapply(paths, sha256, character(1)),
  role = c(
    "frozen phylogenetic lineage membership",
    "stable tree-tip identifiers",
    "focal source antigen and raw 23S annotations",
    "frozen PRN structural and standardised 23S annotations"
  )
)
write_tsv(input_manifest, "input_manifest.tsv")

validation <- list(
  status = "PASS",
  scope = "focal public-tree genomes only",
  n_tree_tips = nrow(assignments),
  n_focal_tree_tips = nrow(analysis),
  n_global_background_excluded = sum(assignments$tree_role == "global_background"),
  n_focal_l10207 = nrow(l1),
  n_focal_source_record_matches = sum(
    analysis$genome_record_id %in% source_records$genome_record_id
  ),
  n_focal_frozen_archive_matches = sum(analysis$frozen_archive_match),
  n_identifier_conflicts = sum(identifier_conflict),
  n_l10207_ptxP_annotated = nrow(l1_ptxp),
  n_l10207_joint_antigen_annotated = nrow(l1_antigen),
  n_l10207_frozen_23S_annotated = nrow(l1_23s),
  n_l10207_frozen_PRN_interpretable = nrow(l1_prn),
  n_valid_period_shift_fisher_tests = sum(valid_shift),
  n_valid_feature_contrast_fisher_tests = sum(valid_contrast),
  project_independent_molecular_period_effect_identifiable_in_CHN_or_JPN = FALSE,
  reason_project_independent_effect_not_identifiable =
    "For the China/Japan expansion strata, no shared project had >=5 annotated records in both periods.",
  analysis_date = as.character(Sys.Date())
)
jsonlite::write_json(
  validation,
  file.path(output_dir, "molecular_characterisation_validation.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

message(
  "Molecular characterisation complete: ",
  nrow(analysis),
  " focal tips; ",
  nrow(l1),
  " L1_02.07; ",
  sum(analysis$frozen_archive_match),
  " frozen-archive matches."
)
