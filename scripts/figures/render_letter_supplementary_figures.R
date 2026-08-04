#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(ragg)
  library(scales)
  library(systemfonts)
  library(tidyr)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(repo_dir, "figures", "letter")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

available_fonts <- unique(systemfonts::system_fonts()$family)
font_family <- if ("Arial" %in% available_fonts) {
  "Arial"
} else if ("Liberation Sans" %in% available_fonts) {
  "Liberation Sans"
} else if ("DejaVu Sans" %in% available_fonts) {
  "DejaVu Sans"
} else {
  "sans"
}

country_order <- c("AUS", "BEL", "CHN", "FRA", "JPN")
country_labels <- c(
  AUS = "Australia", BEL = "Belgium", CHN = "China",
  FRA = "France", JPN = "Japan"
)
country_colors <- c(
  AUS = "#3B78A8", BEL = "#C69232", CHN = "#C75245",
  FRA = "#7C65A8", JPN = "#2F8F83"
)
lineage_order <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
lineage_display_map <- readr::read_tsv(
  file.path(repo_dir, "data", "derived", "model_lineage_display_map.tsv"),
  show_col_types = FALSE
)
lineage_display_labels <- setNames(
  lineage_display_map$short_label,
  lineage_display_map$internal_id
)

read_tsv_repo <- function(...) {
  readr::read_tsv(
    file.path(repo_dir, ...),
    show_col_types = FALSE,
    progress = FALSE
  )
}

theme_letter_supplement <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_letter_supplement())

publication_annotation <- function() {
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 9,
        face = "bold",
        family = font_family,
        hjust = 0
      ),
      plot.tag.position = "topleft",
      plot.tag.location = "margin"
    )
  )
}

save_letter_supplement_png <- function(
    plot,
    stem,
    width_mm = 183,
    height_mm = 120,
    dpi = 600) {
  ragg::agg_png(
    file.path(output_dir, paste0(stem, ".png")),
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = dpi,
    background = "white"
  )
  print(plot)
  dev.off()
}

# Supplementary Figure S1: independent genomic validation of the frozen lineages
pairwise <- read_tsv_repo(
  "results", "cgmlst", "cgmlst_core_snp_pairwise_comparison.tsv"
)
set.seed(20260725)
if (nrow(pairwise) > 30000) {
  pairwise <- dplyr::slice_sample(pairwise, n = 30000)
}

s1a <- ggplot(pairwise, aes(.data[["Distance"]], cgmlst_allelic_mismatches)) +
  geom_point(size = 0.35, alpha = 0.10, colour = "#3B78A8") +
  geom_smooth(
    method = "lm", formula = y ~ x, se = FALSE,
    linewidth = 0.55, colour = "#C75245"
  ) +
  labs(x = "Core-SNP distance", y = "cgMLST allelic mismatches")

nearest_neighbours <- read_tsv_repo(
  "results", "cgmlst", "cgmlst_nearest_neighbour_lineage_concordance.tsv"
) %>%
  filter(primary_model_lineage_id %in% lineage_order) %>%
  group_by(primary_model_lineage_id) %>%
  summarise(
    nearest = mean(nearest_cgmlst_same_lineage),
    five = mean(five_neighbour_same_lineage_fraction),
    .groups = "drop"
  ) %>%
  pivot_longer(c(nearest, five), names_to = "metric", values_to = "agreement") %>%
  mutate(
    lineage = factor(primary_model_lineage_id, levels = lineage_order),
    metric = recode(
      metric,
      nearest = "Nearest neighbour",
      five = "Five nearest neighbours"
    )
  )

s1b <- ggplot(nearest_neighbours, aes(agreement, lineage, shape = metric)) +
  geom_point(size = 2.1, colour = "#3B78A8") +
  scale_x_continuous(limits = c(0.90, 1), labels = label_percent(accuracy = 1)) +
  scale_y_discrete(labels = lineage_display_labels) +
  scale_shape_manual(
    values = c("Nearest neighbour" = 16, "Five nearest neighbours" = 17),
    name = NULL
  ) +
  labs(x = "Same-lineage agreement", y = NULL)

partition <- read_tsv_repo(
  "results", "lineages", "sensitivity", "lineage_partition_sensitivity.tsv"
) %>%
  filter(primary_cluster %in% lineage_order) %>%
  mutate(primary_cluster = factor(primary_cluster, levels = lineage_order))

s1c <- ggplot(partition, aes(jaccard, primary_cluster)) +
  geom_segment(
    aes(x = 0.90, xend = jaccard, yend = primary_cluster),
    linewidth = 0.55, colour = "#B8B8B8"
  ) +
  geom_point(size = 2.1, colour = "#C75245") +
  scale_x_continuous(limits = c(0.90, 1), labels = label_percent(accuracy = 1)) +
  scale_y_discrete(labels = lineage_display_labels) +
  labs(x = "Jaccard similarity", y = NULL)

figure_s1 <- (s1a | (s1b / s1c)) +
  plot_layout(widths = c(1.25, 1), guides = "keep", tag_level = "new") +
  publication_annotation()
save_letter_supplement_png(
  figure_s1,
  "Supplementary_Figure_S1_genomic_validation",
  height_mm = 105
)

# Supplementary Figure S2: fitted country-level posterior predictive checks
monthly <- read_tsv_repo(
  "figures", "source_data", "figure4abc_monthly_counterfactuals.tsv"
) %>%
  mutate(model_month = as.Date(model_month)) %>%
  filter(model_month >= as.Date("2022-01-01"))

posterior_predictive_panel <- function(country) {
  dat <- filter(monthly, country_iso3 == country)
  ggplot(dat, aes(model_month)) +
    geom_ribbon(
      aes(ymin = fitted_lower_95, ymax = fitted_upper_95),
      fill = country_colors[[country]], alpha = 0.18
    ) +
    geom_line(
      aes(y = fitted_median),
      colour = country_colors[[country]], linewidth = 0.7
    ) +
    geom_point(aes(y = observed_cases), size = 0.8, colour = "#222222") +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = NULL, y = "Monthly reported cases")
}

posterior_predictive_metrics <- read_tsv_repo(
  "results", "model_main", "posterior_predictive_metrics.tsv"
) %>%
  select(country_iso3, log_correlation, posterior_predictive_coverage) %>%
  pivot_longer(-country_iso3, names_to = "metric", values_to = "value") %>%
  mutate(
    country = factor(
      country_iso3,
      levels = c("AUS", "CHN", "JPN"),
      labels = country_labels[c("AUS", "CHN", "JPN")]
    ),
    metric = recode(
      metric,
      log_correlation = "Observed-fitted correlation",
      posterior_predictive_coverage = "95% predictive coverage"
    )
  )

s2d <- ggplot(
  posterior_predictive_metrics,
  aes(value, country, shape = metric)
) +
  geom_point(size = 2.3, colour = "#3B78A8") +
  scale_x_continuous(limits = c(0.90, 1), labels = label_percent(accuracy = 1)) +
  scale_shape_manual(
    values = c(
      "Observed-fitted correlation" = 16,
      "95% predictive coverage" = 17
    ),
    name = NULL
  ) +
  labs(x = "Model performance", y = NULL)

figure_s2 <- (
  posterior_predictive_panel("AUS") |
    posterior_predictive_panel("CHN")
) / (
  posterior_predictive_panel("JPN") | s2d
) +
  plot_layout(guides = "collect", tag_level = "new") +
  publication_annotation()
save_letter_supplement_png(
  figure_s2,
  "Supplementary_Figure_S2_posterior_predictive_checks",
  height_mm = 140
)

# Supplementary Figure S3: sampled-ancestry ranks and growth-input sensitivity
analysis_labels <- c(
  primary_ml_midpoint = "Primary",
  primary_ml_outgroup = "Outgroup root",
  bootstrap_consensus_midpoint = "Consensus",
  bionj_midpoint = "BIONJ",
  historical_reference_keep80 = "80% historical",
  global_background_keep80 = "80% background"
)
ancestry_ranks <- read_tsv_repo(
  "results", "phylogeography_sensitivity_summary",
  "country_ancestry_rank_stability.tsv"
) %>%
  select(analysis_id, country_iso3, local_rank_median, reseeding_rank_median) %>%
  pivot_longer(
    c(local_rank_median, reseeding_rank_median),
    names_to = "score_type", values_to = "rank"
  ) %>%
  mutate(
    analysis = factor(
      analysis_id,
      levels = names(analysis_labels),
      labels = analysis_labels
    ),
    country = factor(
      country_iso3,
      levels = rev(country_order),
      labels = rev(country_labels[country_order])
    ),
    score_type = recode(
      score_type,
      local_rank_median = "Sampled historical-local\nancestry-support rank",
      reseeding_rank_median = "Reseeding-compatible sampled\nancestry-support rank"
    )
  )

s3a <- ggplot(ancestry_ranks, aes(analysis, country, fill = rank)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = number(rank, accuracy = 0.1)), size = 2.0) +
  facet_wrap(~score_type, ncol = 1) +
  scale_fill_gradient(
    low = "#C75245", high = "#DCE6EB",
    limits = c(1, 5), breaks = 1:5,
    name = "Country rank\n(1 = highest)"
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

input_growth_labels <- c(
  primary_lower_0_5 = "Primary",
  threshold_0_7 = "Threshold 0.7",
  threshold_0_9 = "Threshold 0.9",
  time_midpoint = "Date midpoint",
  time_interval_uniform = "Date-uniform",
  alternative_root = "Alternative root"
)
input_growth <- read_tsv_repo(
  "results", "model_input_sensitivity_summary", "l10207_input_sensitivity.tsv"
) %>%
  mutate(
    analysis = factor(
      analysis_id,
      levels = rev(names(input_growth_labels)),
      labels = rev(input_growth_labels)
    )
  )

s3b <- ggplot(input_growth, aes(median, analysis)) +
  geom_vline(
    xintercept = 1, linetype = "22",
    linewidth = 0.4, colour = "#777777"
  ) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y", width = 0,
    linewidth = 0.55, colour = "#C75245"
  ) +
  geom_point(
    shape = 21, fill = "#C75245", colour = "#333333", size = 2.0
  ) +
  scale_x_continuous(
    limits = c(0.99, 1.16),
    breaks = c(1, 1.05, 1.10, 1.15)
  ) +
  labs(
    x = "MT28-associated relative net-growth multiplier (95% CrI)",
    y = NULL
  )

figure_s3 <- (s3a | s3b) +
  plot_layout(widths = c(1.25, 1), guides = "keep", tag_level = "new") +
  publication_annotation()
save_letter_supplement_png(
  figure_s3,
  "Supplementary_Figure_S3_ancestry_input_sensitivity",
  height_mm = 115
)

readr::write_tsv(
  tibble::tribble(
    ~figure, ~stem, ~width_mm, ~height_mm, ~panels,
    "Supplementary Figure S1", "Supplementary_Figure_S1_genomic_validation", 183, 105, "A-C",
    "Supplementary Figure S2", "Supplementary_Figure_S2_posterior_predictive_checks", 183, 140, "A-D",
    "Supplementary Figure S3", "Supplementary_Figure_S3_ancestry_input_sensitivity", 183, 115, "A-B"
  ),
  file.path(output_dir, "SUPPLEMENTARY_RENDER_MANIFEST.tsv")
)

readr::write_tsv(
  tibble::tribble(
    ~figure, ~source_file,
    "Supplementary Figure S1", "results/cgmlst/cgmlst_core_snp_pairwise_comparison.tsv",
    "Supplementary Figure S1", "results/cgmlst/cgmlst_nearest_neighbour_lineage_concordance.tsv",
    "Supplementary Figure S1", "results/lineages/sensitivity/lineage_partition_sensitivity.tsv",
    "Supplementary Figure S2", "figures/source_data/figure4abc_monthly_counterfactuals.tsv",
    "Supplementary Figure S2", "results/model_main/posterior_predictive_metrics.tsv",
    "Supplementary Figure S3", "results/phylogeography_sensitivity_summary/country_ancestry_rank_stability.tsv",
    "Supplementary Figure S3", "results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv"
  ),
  file.path(output_dir, "LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv")
)
