#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(repo_dir, "figures", "supplementary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "sans"
country_order <- c("AUS", "BEL", "CHN", "FRA", "JPN")
country_labels <- c(AUS = "Australia", BEL = "Belgium", CHN = "China", FRA = "France", JPN = "Japan")
country_colors <- c(AUS = "#3B78A8", BEL = "#C69232", CHN = "#C75245", FRA = "#7C65A8", JPN = "#2F8F83")
lineage_order <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
lineage_colors <- c(L1_01.02 = "#3B78A8", L1_02.05 = "#7C65A8", L1_02.06 = "#2F8F83", L1_02.07 = "#C75245", Other = "#777777")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.tag = element_text(size = 9, face = "bold"),
      plot.tag.position = c(0.01, 0.99),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}
theme_set(theme_pub())

tag_panel <- function(plot, label) {
  plot + labs(tag = label) +
    theme(plot.tag.location = "panel", plot.tag.position = c(0.01, 0.99))
}

save_png <- function(plot, stem, width_mm = 183, height_mm = 145) {
  ragg::agg_png(
    file.path(output_dir, paste0(stem, ".png")),
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 600,
    background = "white"
  )
  print(plot)
  dev.off()
}

read_tsv_repo <- function(...) {
  read_tsv(file.path(repo_dir, ...), show_col_types = FALSE, progress = FALSE)
}

# Figure S1: cohort, QC, temporal coverage, and clock signal
qc <- read_tsv_repo("results", "phylogeny", "uniform_sequence_qc.tsv") %>%
  filter(tree_role == "focal", country_iso3 %in% country_order) %>%
  mutate(
    country = factor(country_iso3, levels = country_order, labels = country_labels[country_order]),
    qc_result = if_else(tree_include_after_uniform_qc, "Passed uniform QC", "Excluded")
  )

s1a <- ggplot(qc, aes(country, fill = qc_result)) +
  geom_bar(position = "stack", width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("Passed uniform QC" = "#3B78A8", "Excluded" = "#D3D3D3"), name = NULL) +
  labs(x = NULL, y = "Selected focal genomes") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

missingness <- read_tsv_repo("results", "phylogeny", "alignment_missingness_final.tsv")
missing_col <- intersect(c("missing_fraction", "alignment_missingness", "missingness"), names(missingness))[[1]]
id_col <- intersect(c("tree_sample_id", "sample_id", "name"), names(missingness))[[1]]
tree_meta <- read_tsv_repo("results", "phylogeny", "tree_tip_metadata.tsv")
missingness <- missingness %>%
  rename(tree_sample_id = all_of(id_col), missing_fraction = all_of(missing_col)) %>%
  left_join(tree_meta %>% select(tree_sample_id, tree_role), by = "tree_sample_id") %>%
  mutate(tree_role = recode(tree_role, focal = "Focal", background = "Background"))

s1b <- ggplot(missingness, aes(missing_fraction, fill = tree_role)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.68, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(Focal = "#3B78A8", Background = "#AFAFAF"), name = NULL) +
  scale_x_continuous(labels = label_percent(accuracy = 1)) +
  labs(x = "Final-alignment missingness", y = "Genomes")

focal_period <- tree_meta %>%
  filter(tree_role == "focal", country_iso3 %in% country_order) %>%
  mutate(
    country = factor(country_iso3, levels = country_order, labels = country_labels[country_order]),
    epidemic_period = factor(
      epidemic_period,
      levels = c("prepandemic", "pandemic", "resurgence"),
      labels = c("Pre-pandemic", "Pandemic", "Resurgence")
    )
  )

s1c <- ggplot(focal_period, aes(country, fill = epidemic_period)) +
  geom_bar(width = 0.68, colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("Pre-pandemic" = "#808080", "Pandemic" = "#9DBDD5", "Resurgence" = "#DF7D6D"),
    name = "Sampling period"
  ) +
  labs(x = NULL, y = "Final-tree focal genomes") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

clock_data <- tree_meta %>%
  filter(is.finite(decimal_date), is.finite(root_to_tip)) %>%
  mutate(role = if_else(tree_role == "focal", "Focal", "Background"))

s1d <- ggplot(clock_data, aes(decimal_date, root_to_tip, colour = role)) +
  geom_point(size = 0.65, alpha = 0.55) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "#222222") +
  scale_colour_manual(values = c(Focal = "#3B78A8", Background = "#B0B0B0"), name = NULL) +
  labs(x = "Sampling year", y = "Root-to-tip distance")

figure_s1 <- (tag_panel(s1a, "A") | tag_panel(s1b, "B")) /
  (tag_panel(s1c, "C") | tag_panel(s1d, "D")) +
  plot_layout(guides = "collect")
save_png(figure_s1, "Figure_S1_cohort_qc_and_temporal_signal")

# Figure S2: independent genomic representations and sensitivity
pairwise <- read_tsv_repo("results", "cgmlst", "cgmlst_core_snp_pairwise_comparison.tsv")
set.seed(20260725)
if (nrow(pairwise) > 30000) pairwise <- slice_sample(pairwise, n = 30000)

s2a <- ggplot(pairwise, aes(.data[["Distance"]], cgmlst_allelic_mismatches)) +
  geom_point(size = 0.35, alpha = 0.10, colour = "#3B78A8") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "#C75245") +
  labs(x = "Core-SNP distance", y = "cgMLST allelic mismatches")

nn <- read_tsv_repo("results", "cgmlst", "cgmlst_nearest_neighbour_lineage_concordance.tsv") %>%
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
    metric = recode(metric, nearest = "Nearest neighbour", five = "Five nearest neighbours")
  )

s2b <- ggplot(nn, aes(agreement, lineage, shape = metric)) +
  geom_point(size = 2.1, colour = "#3B78A8") +
  scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
  scale_shape_manual(values = c("Nearest neighbour" = 16, "Five nearest neighbours" = 17), name = NULL) +
  labs(x = "Same-lineage agreement", y = NULL)

partition <- read_tsv_repo("results", "lineages", "sensitivity", "lineage_partition_sensitivity.tsv") %>%
  filter(primary_cluster %in% lineage_order) %>%
  mutate(primary_cluster = factor(primary_cluster, levels = lineage_order))

s2c <- ggplot(partition, aes(jaccard, primary_cluster)) +
  geom_segment(aes(x = 0.9, xend = jaccard, yend = primary_cluster), linewidth = 0.55, colour = "#B8B8B8") +
  geom_point(size = 2.1, colour = "#C75245") +
  scale_x_continuous(limits = c(0.9, 1), labels = label_percent(accuracy = 1)) +
  labs(x = "Jaccard similarity", y = NULL)

threshold <- read_tsv_repo("figures", "source_data", "figure3d_threshold_sensitivity.tsv") %>%
  select(
    transition_threshold,
    `Reseeding edges` = n_high_support_post_reseeding_edges,
    `Successful sampled clusters` = n_successful_sampled_clusters
  ) %>%
  pivot_longer(-transition_threshold, names_to = "quantity", values_to = "count")

s2d <- ggplot(threshold, aes(transition_threshold, count, colour = quantity, group = quantity)) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Reseeding edges" = "#3B78A8", "Successful sampled clusters" = "#C75245"), name = NULL) +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
  labs(x = "Transition-support threshold", y = "Count")

figure_s2 <- (tag_panel(s2a, "A") | tag_panel(s2b, "B")) /
  (tag_panel(s2c, "C") | tag_panel(s2d, "D")) +
  plot_layout(guides = "collect")
save_png(figure_s2, "Figure_S2_genomic_concordance_and_sensitivity")

# Figure S3: country-level posterior predictive checks
monthly <- read_tsv_repo("figures", "source_data", "figure5abc_monthly_counterfactuals.tsv") %>%
  mutate(model_month = as.Date(model_month)) %>%
  filter(model_month >= as.Date("2022-01-01"))

ppc_panel <- function(country, label) {
  dat <- filter(monthly, country_iso3 == country)
  panel <- ggplot(dat, aes(model_month)) +
    geom_ribbon(aes(ymin = fitted_lower_95, ymax = fitted_upper_95), fill = country_colors[[country]], alpha = 0.18) +
    geom_line(aes(y = fitted_median), colour = country_colors[[country]], linewidth = 0.7) +
    geom_point(aes(y = observed_cases), size = 0.8, colour = "#222222") +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = NULL, y = "Monthly cases")
  tag_panel(panel, label)
}

ppc_metrics <- read_tsv_repo("results", "model_main", "posterior_predictive_metrics.tsv") %>%
  select(country_iso3, log_correlation, posterior_predictive_coverage) %>%
  pivot_longer(-country_iso3, names_to = "metric", values_to = "value") %>%
  mutate(
    country = factor(country_iso3, levels = c("AUS", "CHN", "JPN"), labels = country_labels[c("AUS", "CHN", "JPN")]),
    metric = recode(metric, log_correlation = "Observed–fitted correlation", posterior_predictive_coverage = "95% predictive coverage")
  )

s3d <- ggplot(ppc_metrics, aes(value, country, shape = metric)) +
  geom_point(size = 2.3, colour = "#3B78A8") +
  scale_x_continuous(limits = c(0.9, 1), labels = label_percent(accuracy = 1)) +
  scale_shape_manual(values = c("Observed–fitted correlation" = 16, "95% predictive coverage" = 17), name = NULL) +
  labs(x = "Model performance", y = NULL)

figure_s3 <- (ppc_panel("AUS", "A") | ppc_panel("CHN", "B")) /
  (ppc_panel("JPN", "C") | tag_panel(s3d, "D")) +
  plot_layout(guides = "collect")
save_png(figure_s3, "Figure_S3_posterior_predictive_checks")

# Figure S4: sampling, counterfactual, and recovery sensitivity
growth <- read_tsv_repo("figures", "source_data", "figure4b_lineage_growth_sensitivity.tsv") %>%
  select(
    lineage,
    median_main, lower_95_main, upper_95_main,
    median_no_project, lower_95_no_project, upper_95_no_project
  ) %>%
  pivot_longer(
    -lineage,
    names_to = c(".value", "model"),
    names_pattern = "(median|lower_95|upper_95)_(main|no_project)"
  ) %>%
  mutate(
    lineage = factor(lineage, levels = lineage_order),
    model = recode(model, main = "Project-adjusted", no_project = "No project effects")
  )

s4a <- ggplot(growth, aes(median, lineage, colour = model)) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.4, colour = "#777777") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), position = position_dodge(width = 0.45), height = 0) +
  geom_point(position = position_dodge(width = 0.45), size = 1.8) +
  scale_colour_manual(values = c("Project-adjusted" = "#3B78A8", "No project effects" = "#C75245"), name = NULL) +
  labs(x = "Relative net-growth multiplier", y = NULL) +
  guides(colour = guide_legend(nrow = 1))

shares <- read_tsv_repo("figures", "source_data", "figure4d_raw_vs_corrected_shares.tsv") %>%
  mutate(
    country = factor(country_iso3, levels = c("AUS", "CHN", "JPN"), labels = country_labels[c("AUS", "CHN", "JPN")]),
    lineage = factor(lineage, levels = lineage_order)
  )

s4b <- ggplot(shares, aes(raw_public_tree_share, median, colour = lineage)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.4, colour = "#777777") +
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95), width = 0, linewidth = 0.35) +
  geom_point(size = 1.5) +
  facet_wrap(~country, nrow = 1) +
  scale_colour_manual(values = lineage_colors, name = "Lineage") +
  scale_x_continuous(labels = label_percent()) +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "Raw public-tree share", y = "Sampling-corrected share") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

counterfactual <- read_tsv_repo("figures", "source_data", "figure5d_counterfactual_summary.tsv") %>%
  filter(scenario %in% c("no_new_introduction_case_reduction_fraction", "lineage_difference_effect_fraction")) %>%
  select(
    country_iso3, scenario,
    median_main, lower_95_main, upper_95_main,
    median_no_project, lower_95_no_project, upper_95_no_project
  ) %>%
  pivot_longer(
    -c(country_iso3, scenario),
    names_to = c(".value", "model"),
    names_pattern = "(median|lower_95|upper_95)_(main|no_project)"
  ) %>%
  mutate(
    country = factor(country_iso3, levels = c("AUS", "CHN", "JPN"), labels = country_labels[c("AUS", "CHN", "JPN")]),
    scenario = recode(
      scenario,
      no_new_introduction_case_reduction_fraction = "No new introduction",
      lineage_difference_effect_fraction = "Equal lineage growth"
    ),
    model = recode(model, main = "Project-adjusted", no_project = "No project effects")
  )

s4c <- ggplot(counterfactual, aes(median, country, colour = model)) +
  geom_vline(xintercept = 0, linetype = "22", linewidth = 0.4, colour = "#777777") +
  geom_errorbarh(aes(xmin = lower_95, xmax = upper_95), position = position_dodge(width = 0.45), height = 0) +
  geom_point(position = position_dodge(width = 0.45), size = 1.7) +
  facet_wrap(~scenario, scales = "free_x") +
  scale_colour_manual(values = c("Project-adjusted" = "#3B78A8", "No project effects" = "#C75245"), name = NULL) +
  scale_x_continuous(labels = label_percent()) +
  labs(x = "Proportional case reduction", y = NULL) +
  guides(colour = "none")

recovery <- read_tsv_repo("results", "model_main", "all_recovery.tsv") %>%
  mutate(parameter_type = recode(parameter_type, lineage_growth = "Lineage growth", import_scale = "Import scale"))

s4d <- ggplot(recovery, aes(truth, `50%`, colour = parameter_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.4, colour = "#777777") +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0, linewidth = 0.3, alpha = 0.55) +
  geom_point(size = 1.2, alpha = 0.8) +
  facet_wrap(~parameter_type, scales = "free") +
  scale_colour_manual(values = c("Lineage growth" = "#3B78A8", "Import scale" = "#C75245"), guide = "none") +
  labs(x = "Simulated truth", y = "Posterior median")

figure_s4 <- (tag_panel(s4a, "A") | tag_panel(s4b, "B")) /
  (tag_panel(s4c, "C") | tag_panel(s4d, "D")) +
  plot_layout(guides = "keep")
save_png(figure_s4, "Figure_S4_sampling_counterfactual_and_recovery_sensitivity", height_mm = 155)

write_tsv(
  tibble(
    figure = paste0("Figure S", 1:4),
    file = c(
      "Figure_S1_cohort_qc_and_temporal_signal.png",
      "Figure_S2_genomic_concordance_and_sensitivity.png",
      "Figure_S3_posterior_predictive_checks.png",
      "Figure_S4_sampling_counterfactual_and_recovery_sensitivity.png"
    ),
    format = "PNG-600dpi",
    backend = "R",
    rendered_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
  ),
  file.path(output_dir, "RENDER_MANIFEST.tsv")
)

message("Rendered four R multi-panel supplementary figures to: ", output_dir)
