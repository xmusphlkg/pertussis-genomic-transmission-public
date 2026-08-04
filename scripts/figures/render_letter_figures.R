#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggtree)
  library(patchwork)
  library(phangorn)
  library(readr)
  library(ragg)
  library(scales)
  library(systemfonts)
  library(tibble)
  library(tidyr)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
source_dir <- file.path(repo_dir, "figures", "source_data")
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

country_order <- c("AUS", "CHN", "JPN")
country_order_full <- c("AUS", "BEL", "CHN", "FRA", "JPN")
country_labels <- c(
  AUS = "Australia", BEL = "Belgium", CHN = "China",
  FRA = "France", JPN = "Japan", Background = "Background"
)
country_colors <- c(
  AUS = "#3B78A8", BEL = "#C69232", CHN = "#C75245",
  FRA = "#7C65A8", JPN = "#2F8F83", Background = "#D4D4D4"
)
period_colors <- c(
  `Pre-pandemic` = "#808080",
  Pandemic = "#9DBDD5",
  Resurgence = "#DF7D6D"
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
lineage_colors <- c(
  L1_01.02 = "#3B78A8",
  L1_02.05 = "#7C65A8",
  L1_02.06 = "#2F8F83",
  L1_02.07 = "#C75245",
  Other = "#777777",
  Background = "#D0D0D0"
)

read_source <- function(filename) {
  readr::read_tsv(
    file.path(source_dir, filename),
    show_col_types = FALSE,
    progress = FALSE
  )
}

theme_letter <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#222222"),
      axis.ticks.length = unit(1.1, "mm"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.45, colour = "#222222"),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.55),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(4.3, "mm"),
      legend.spacing.x = unit(1.0, "mm"),
      legend.box.spacing = unit(0.7, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size - 0.1, face = "bold"),
      plot.title = element_text(size = base_size + 0.35, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.35, colour = "#555555"),
      plot.caption = element_text(size = base_size - 0.8, colour = "#555555"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

theme_set(theme_letter())

save_letter_figure <- function(plot, stem, width_mm = 183, height_mm = 150) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::cairo_pdf(
    file.path(output_dir, paste0(stem, ".pdf")),
    width = width_in,
    height = height_in,
    family = font_family,
    bg = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    file.path(output_dir, paste0(stem, ".png")),
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 300,
    background = "white"
  )
  print(plot)
  dev.off()
}

publication_annotation <- function() {
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 8,
        face = "bold",
        family = font_family,
        margin = margin(0, 1.4, 1.3, 0),
        hjust = 0
      ),
      plot.tag.position = "topleft",
      plot.tag.location = "margin"
    )
  )
}

fmt_interval <- function(median, lower, upper, digits = 3) {
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
    median,
    lower,
    upper
  )
}

# ---------------------------------------------------------------------------
# Letter Figure 1: phylogenetic definition and observation context
# ---------------------------------------------------------------------------

tree_manifest <- read_source("figure2_tree_manifest.tsv")
tree <- phangorn::midpoint(ape::read.tree(tree_manifest$tree_file[[1]]))
tree_meta <- read_source("figure2a_tree_tip_metadata.tsv") %>%
  select(tree_sample_id, tree_role, country_iso3, date_lower, date_upper)
tip_ancestry <- read_source("figure2b_tip_ancestry_support.tsv") %>%
  select(tree_sample_id, primary_model_lineage_id)
tree_meta <- tree_meta %>%
  left_join(tip_ancestry, by = "tree_sample_id") %>%
  transmute(
    label = tree_sample_id,
    country_display = if_else(tree_role == "focal", country_iso3, "Background"),
    lineage_display = if_else(
      tree_role == "focal",
      coalesce(primary_model_lineage_id, "Other"),
      "Background"
    ),
    sampling_year = {
      lower_day <- as.numeric(as.Date(date_lower))
      upper_day <- as.numeric(as.Date(date_upper))
      midpoint <- as.Date((lower_day + upper_day) / 2, origin = "1970-01-01")
      as.numeric(format(midpoint, "%Y")) +
        (as.numeric(format(midpoint, "%j")) - 1) / 365.25
    }
  )

p_tree_base <- ggtree(tree, linewidth = 0.10, colour = "#A7A7A7", ladderize = TRUE) %<+%
  tree_meta
target_tips <- p_tree_base$data %>%
  filter(isTip, lineage_display == "L1_02.07")
tip_strip_wide <- p_tree_base$data %>%
  filter(isTip) %>%
  transmute(
    y,
    Lineage = factor(
      coalesce(as.character(lineage_display), "Background"),
      levels = c(lineage_order, "Background")
    ),
    Sampling_year = sampling_year
  )
tip_y_limits <- c(0.5, max(tip_strip_wide$y) + 0.5)
sampling_year_limits <- range(tip_strip_wide$Sampling_year, na.rm = TRUE)
sampling_year_palette <- c(
  "#9A133DFF", "#B93961FF", "#D8527CFF", "#F28AAAFF",
  "#F9B4C9FF", "#F9E0E8FF", "#FFFFFFFF", "#EAF3FFFF",
  "#C5DAF6FF", "#A1C2EDFF", "#6996E3FF", "#4060C8FF",
  "#1A318BFF"
)
sampling_year_breaks <- breaks_pretty(n = 4)(sampling_year_limits)

p_tree_base$data <- p_tree_base$data %>%
  mutate(
    high_support = if_else(
      !isTip & grepl("/", label, fixed = TRUE),
      vapply(
        strsplit(if_else(is.na(label), "", label), "/", fixed = TRUE),
        function(z) {
          if (length(z) != 2) return(FALSE)
          values <- suppressWarnings(as.numeric(z))
          all(is.finite(values)) && all(values >= 95)
        },
        logical(1)
      ),
      FALSE
    )
  )
target_mrca <- ape::getMRCA(tree, target_tips$label)
tree_x_max <- max(p_tree_base$data$x, na.rm = TRUE)
target_tip_x_max <- max(target_tips$x, na.rm = TRUE)

p1a_tree <- p_tree_base +
  geom_hilight(
    node = target_mrca,
    fill = "#C75245",
    alpha = 0.060,
    extend = 0.00012
  ) +
  geom_point(
    data = p_tree_base$data %>% filter(!isTip, high_support),
    aes(x = x, y = y, shape = "SH-aLRT and UFBoot ≥95%"),
    inherit.aes = FALSE,
    size = 0.24,
    colour = "#333333",
    alpha = 0.65
  ) +
  geom_tippoint(
    aes(colour = lineage_display),
    shape = 22,
    size = 0.1,
    alpha = 0,
    show.legend = TRUE
  ) +
  geom_tippoint(
    aes(fill = country_display),
    shape = 21, size = 0.50, colour = "#555555",
    stroke = 0.14, alpha = 0.88, show.legend = TRUE
  ) +
  geom_tippoint(
    data = target_tips,
    shape = 21, size = 0.82, fill = "#C75245",
    colour = "#8F2E2A", stroke = 0.22, alpha = 0.95,
    inherit.aes = FALSE,
    aes(x = x, y = y)
  ) +
  annotate(
    "text",
    x = target_tip_x_max + 0.00006,
    y = median(target_tips$y),
    label = "MT28-associated\n288-tip clade",
    hjust = 0,
    size = 1.78,
    lineheight = 0.88,
    colour = "#8F2E2A",
    family = font_family
  ) +
  scale_fill_manual(
    values = country_colors,
    breaks = c(country_order_full, "Background"),
    labels = country_labels[c(country_order_full, "Background")],
    name = "Country",
    guide = guide_legend(
      order = 1,
      nrow = 3,
      byrow = TRUE,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        shape = 21, size = 1.85, colour = "#555555",
        stroke = 0.22, alpha = 1
      )
    )
  ) +
  scale_colour_manual(
    values = lineage_colors,
    breaks = c(lineage_order, "Background"),
    labels = c(lineage_display_labels, Background = "Background"),
    name = "Frozen lineage",
    guide = guide_legend(
      order = 2,
      nrow = 3,
      byrow = TRUE,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        shape = 22,
        size = 2.25,
        fill = unname(lineage_colors[c(lineage_order, "Background")]),
        colour = unname(lineage_colors[c(lineage_order, "Background")]),
        stroke = 0.18,
        alpha = 1
      )
    )
  ) +
  scale_shape_manual(
    values = c(`SH-aLRT and UFBoot ≥95%` = 16),
    name = "Internal-node support",
    guide = guide_legend(
      order = 3,
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0,
      override.aes = list(size = 1.25, colour = "#333333", alpha = 0.9)
    )
  ) +
  geom_treescale(x = 0, y = 11, width = 0.001, fontsize = 1.9, linesize = 0.30) +
  coord_cartesian(xlim = c(0, tree_x_max * 1.06), clip = "off") +
  theme_tree2() +
  labs(x = "Substitutions per site") +
  theme(
    text = element_text(family = font_family, size = 6.0),
    axis.title.x = element_text(size = 5.2, colour = "#333333", margin = margin(t = 1.2)),
    axis.text.x = element_text(size = 4.8, colour = "#333333"),
    axis.ticks.x = element_line(linewidth = 0.22, colour = "#555555"),
    axis.line.x = element_line(linewidth = 0.25, colour = "#555555"),
    legend.position = c(0.018, 0.985),
    legend.justification = c(0, 1),
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.title = element_text(size = 4.7, face = "bold"),
    legend.text = element_text(size = 4.15),
    legend.key.height = unit(1.9, "mm"),
    legend.key.width = unit(2.4, "mm"),
    legend.spacing.x = unit(0.45, "mm"),
    legend.spacing.y = unit(0.10, "mm"),
    legend.background = element_blank(),
    legend.box.background = element_rect(
      fill = scales::alpha("white", 0.94),
      colour = "#D2D2D2",
      linewidth = 0.22
    ),
    legend.box.margin = margin(0.7, 1.1, 0.7, 1.1),
    plot.margin = margin(3, 8, 3, 3)
  )

p1a_lineage_strip <- ggplot(tip_strip_wide, aes(1, y, fill = Lineage)) +
  geom_tile(width = 0.82, height = 1.0) +
  scale_fill_manual(values = lineage_colors, guide = "none", drop = FALSE) +
  scale_x_continuous(
    position = "top", breaks = 1, labels = "Lineage",
    expand = expansion(add = 0.12)
  ) +
  scale_y_continuous(limits = tip_y_limits, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = font_family, base_size = 6.0) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0.5,
      size = 4.5, colour = "#333333"
    ),
    plot.margin = margin(3, 1, 3, 1)
  )

p1a_time_strip <- ggplot(tip_strip_wide, aes(1, y, fill = Sampling_year)) +
  geom_tile(width = 0.82, height = 1.0) +
  scale_fill_gradientn(
    colours = sampling_year_palette,
    limits = sampling_year_limits,
    breaks = sampling_year_breaks,
    labels = label_number(accuracy = 1),
    na.value = "#B3B3B3",
    guide = "none"
  ) +
  scale_x_continuous(
    position = "top", breaks = 1, labels = "Year",
    expand = expansion(add = 0.12)
  ) +
  scale_y_continuous(limits = tip_y_limits, expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = font_family, base_size = 6.0) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0.5,
      size = 4.5, colour = "#333333"
    ),
    plot.margin = margin(3, 1, 3, 1)
  )

p1a <- wrap_elements(
  full = p1a_tree + p1a_lineage_strip + p1a_time_strip +
    plot_layout(widths = c(1, 0.055, 0.075))
)

case_data <- read_source("figure1a_cases.tsv") %>%
  mutate(model_month = as.Date(model_month)) %>%
  filter(
    country_iso3 %in% country_order,
    model_month >= as.Date("2019-01-01"),
    model_month <= as.Date("2025-12-31"),
    case_data_available
  ) %>%
  mutate(country = factor(country_labels[country_iso3], levels = country_labels[country_order]))

p1c <- ggplot(case_data, aes(model_month, cases, colour = country_iso3)) +
  annotate(
    "rect",
    xmin = as.Date("2020-03-01"),
    xmax = as.Date("2022-12-31"),
    ymin = -Inf,
    ymax = Inf,
    fill = "#DCE9F2",
    alpha = 0.55
  ) +
  geom_line(linewidth = 0.45, lineend = "round") +
  facet_wrap(~country, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = country_colors, guide = "none") +
  scale_x_date(
    breaks = as.Date(c("2019-01-01", "2022-01-01", "2025-01-01")),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Date", y = "Reported cases per month")

annual_genomes <- read_source("figure1b_annual_genomes.tsv") %>%
  filter(country_iso3 %in% country_order, sampling_year >= 2015) %>%
  mutate(
    country = factor(country_labels[country_iso3], levels = country_labels[country_order]),
    epidemic_period = factor(
      recode(
        epidemic_period,
        prepandemic = "Pre-pandemic",
        pandemic = "Pandemic",
        resurgence = "Resurgence"
      ),
      levels = names(period_colors)
    )
  )

p1d <- ggplot(annual_genomes, aes(sampling_year, n_sampled_genomes, fill = epidemic_period)) +
  geom_col(width = 0.82, colour = NA) +
  facet_wrap(~country, nrow = 1, scales = "free_y") +
  scale_fill_manual(
    values = period_colors,
    name = "Sampling period",
    guide = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0
    )
  ) +
  scale_x_continuous(breaks = c(2015, 2020, 2025), expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Sampling year", y = "Focal genomes") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.4),
    legend.key.width = unit(3.2, "mm"),
    legend.key.height = unit(2.2, "mm"),
    legend.spacing.x = unit(0.7, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  )

ct_curve <- read_source("figure4e_australia_ct_curve.tsv")
p1f <- ggplot(ct_curve, aes(ct, success_probability)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#9CBED5", alpha = 0.38) +
  geom_line(colour = "#3B78A8", linewidth = 0.75) +
  annotate(
    "text",
    x = 30.0,
    y = 0.82,
    label = "OR 0.684\n95% CI: 0.62-0.76",
    hjust = 0.5,
    size = 1.95,
    colour = "#333333",
    family = font_family,
    lineheight = 0.92
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = label_percent()) +
  labs(x = "PCR cycle threshold (Ct)", y = "Complete-profile probability")

selection_weighted <- read_source("figure3d_selection_cap_weighted_l10207_shares.tsv") %>%
  filter(country_iso3 %in% country_order) %>%
  select(country_iso3, epidemic_period, unweighted_share, selection_cap_weighted_share) %>%
  pivot_longer(
    c(unweighted_share, selection_cap_weighted_share),
    names_to = "estimate",
    values_to = "share"
  ) %>%
  mutate(
    country = factor(country_labels[country_iso3], levels = country_labels[country_order]),
    period = factor(
      recode(epidemic_period, prepandemic = "Pre-pandemic", resurgence = "Resurgence"),
      levels = c("Pre-pandemic", "Resurgence")
    ),
    estimate = recode(
      estimate,
      unweighted_share = "Final tree",
      selection_cap_weighted_share = "Cap weighted"
    )
  )

p1e <- ggplot(
  selection_weighted,
  aes(period, share, colour = estimate, shape = estimate, group = estimate)
) +
  geom_line(linewidth = 0.52) +
  geom_point(size = 1.9, stroke = 0.45) +
  facet_wrap(~country, nrow = 1) +
  scale_colour_manual(values = c(`Final tree` = "#3B78A8", `Cap weighted` = "#C75245"), name = NULL) +
  scale_shape_manual(values = c(`Final tree` = 16, `Cap weighted` = 1), name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = label_percent()) +
  labs(x = NULL, y = "MT28-associated share") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom"
  )

figure1 <- p1a + p1c + p1d + p1e + p1f +
  plot_layout(
    design = "
AB
AC
AD
AE
",
    widths = c(1.35, 1.0),
    heights = c(1, 1, 1, 1),
    guides = "keep",
    tag_level = "new"
  ) +
  publication_annotation()

save_letter_figure(figure1, "Figure_1_observation_structure", 183, 190)

# ---------------------------------------------------------------------------
# Letter Figure 2: relative growth signal and robustness
# ---------------------------------------------------------------------------

growth_main <- read_source("figure3a_lineage_growth_main.tsv") %>%
  mutate(
    lineage = factor(lineage, levels = rev(lineage_order)),
    highlight = lineage == "L1_02.07",
    interval = fmt_interval(median, lower_95, upper_95)
  )
target_growth <- growth_main %>% filter(lineage == "L1_02.07")

pairwise_growth <- read_source("figure3b_l10207_pairwise_growth.tsv") %>%
  select(
    denominator,
    median_main, lower_95_main, upper_95_main,
    median_no_project, lower_95_no_project, upper_95_no_project
  ) %>%
  pivot_longer(
    -denominator,
    names_to = c(".value", "model"),
    names_pattern = "(median|lower_95|upper_95)_(main|no_project)"
  ) %>%
  mutate(
    comparator = factor(denominator, levels = rev(c("L1_01.02", "L1_02.05", "L1_02.06", "Other"))),
    model = recode(model, main = "Project-adjusted", no_project = "No project effects")
  )

growth_robustness <- read_source("figure3c_l10207_growth_robustness.tsv") %>%
  filter(diagnostic_pass) %>%
  mutate(
    country_label = case_when(
      analysis_type == "country_only" ~ recode(country_subset, !!!country_labels),
      analysis_type == "country_omission" ~ recode(omitted_country, !!!country_labels),
      TRUE ~ recode(omitted_project_country, !!!country_labels)
    ),
    model = recode(
      observation_specification,
      project_adjusted = "Project-adjusted",
      no_project = "No project effects"
    )
  )

input_sensitivity <- readr::read_tsv(
  file.path(repo_dir, "results", "model_input_sensitivity_summary", "l10207_input_sensitivity.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    analysis_label = recode(
      analysis_label,
      `Primary: threshold 0.5, lower-bound time` = "Primary",
      `Transition threshold 0.7` = "Threshold 0.7",
      `Transition threshold 0.9` = "Threshold 0.9",
      `Earliest-sample interval midpoint` = "Midpoint time",
      `Earliest-sample interval-uniform` = "Uniform time",
      `Alternative root` = "Alternative root"
    )
  )

scenario_label <- "MT28-associated multiplier set to reference"
scenario_summary <- read_source("figure4d_counterfactual_summary.tsv") %>%
  filter(
    scenario %in% c(
      "no_new_exposure_case_difference_fraction",
      "l10207_growth_scenario_difference_fraction"
    )
  ) %>%
  select(
    country_iso3,
    scenario,
    median_main, lower_95_main, upper_95_main,
    median_no_project, lower_95_no_project, upper_95_no_project
  ) %>%
  pivot_longer(
    cols = c(median_main, lower_95_main, upper_95_main, median_no_project, lower_95_no_project, upper_95_no_project),
    names_to = c(".value", "model"),
    names_pattern = "(median|lower_95|upper_95)_(main|no_project)"
  ) %>%
  mutate(
    country = factor(country_labels[country_iso3], levels = rev(country_labels[country_order])),
    scenario = recode(
      scenario,
      no_new_exposure_case_difference_fraction = "No post-2022 exposure",
      l10207_growth_scenario_difference_fraction = scenario_label
    ),
    model = recode(model, main = "Project-adjusted", no_project = "No project effects")
  )

recovery <- read_source("figure4f_identifiability_recovery.tsv") %>%
  mutate(
    parameter_type = recode(parameter_type, lineage_growth = "Lineage\ngrowth", import_scale = "Import\nscale")
  ) %>%
  pivot_longer(
    c(coverage_95, median_absolute_log_error, correlation_truth_posterior_median),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      coverage_95 = "95% coverage",
      median_absolute_log_error = "Median |log error|",
      correlation_truth_posterior_median = "Truth-estimate r"
    ),
    pass = case_when(
      metric == "coverage_95" ~ value >= 0.8,
      metric == "median_absolute_log_error" & parameter_type == "Lineage\ngrowth" ~ value <= 0.2,
      metric == "median_absolute_log_error" & parameter_type == "Import\nscale" ~ value <= 0.5,
      metric == "correlation_truth_posterior_median" ~ value >= 0.7,
      TRUE ~ FALSE
    ),
    display = case_when(
      metric == "coverage_95" ~ percent(value, accuracy = 0.1),
      TRUE ~ number(value, accuracy = 0.001)
    ),
    metric_label = factor(metric_label, levels = c("95% coverage", "Median |log error|", "Truth-estimate r")),
    parameter_type = factor(parameter_type, levels = c("Lineage\ngrowth", "Import\nscale"))
  )

p2b <- ggplot(growth_main, aes(median, lineage)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.38) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95, colour = highlight),
    orientation = "y",
    width = 0,
    linewidth = 0.68
  ) +
  geom_point(aes(fill = highlight), shape = 21, size = 2.1, stroke = 0.45) +
  scale_colour_manual(values = c(`FALSE` = "#4B4B4B", `TRUE` = "#C75245"), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#C75245"), guide = "none") +
  scale_x_continuous(limits = c(0.88, 1.18), breaks = seq(0.9, 1.15, 0.05)) +
  scale_y_discrete(labels = lineage_display_labels) +
  labs(x = "Relative net-growth multiplier", y = NULL)

p2c <- ggplot(pairwise_growth, aes(median, comparator, colour = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.38) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0,
    position = position_dodge(width = 0.42),
    linewidth = 0.50
  ) +
  geom_point(position = position_dodge(width = 0.42), size = 1.75, stroke = 0.45) +
  scale_colour_manual(values = c(`Project-adjusted` = "#C75245", `No project effects` = "#3B78A8"), name = NULL) +
  scale_shape_manual(values = c(`Project-adjusted` = 16, `No project effects` = 1), name = NULL) +
  scale_x_continuous(limits = c(0.99, 1.27), breaks = c(1, 1.1, 1.2)) +
  scale_y_discrete(labels = lineage_display_labels) +
  labs(x = "MT28-associated / comparator", y = NULL) +
  theme(legend.position = "bottom")

country_only <- growth_robustness %>%
  filter(analysis_type == "country_only") %>%
  mutate(country_label = factor(country_label, levels = rev(country_labels[country_order])))

p2d <- ggplot(country_only, aes(median, country_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.36) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0,
    position = position_dodge(width = 0.40),
    linewidth = 0.48
  ) +
  geom_point(position = position_dodge(width = 0.40), size = 1.75, stroke = 0.45) +
  scale_colour_manual(
    values = c(`Project-adjusted` = "#C75245", `No project effects` = "#3B78A8"),
    name = NULL
  ) +
  scale_shape_manual(values = c(`Project-adjusted` = 16, `No project effects` = 1), name = NULL) +
  scale_x_continuous(limits = c(0.65, 1.55), breaks = c(0.8, 1.0, 1.2, 1.4)) +
  labs(x = "Country-only multiplier", y = NULL) +
  theme(legend.position = "none")

omission_input <- bind_rows(
  growth_robustness %>%
    filter(analysis_type == "country_omission", observation_specification == "project_adjusted") %>%
    transmute(
      row_label = paste("Omit", country_label),
      group = "Omit country",
      median, lower_95, upper_95
    ),
  growth_robustness %>%
    filter(analysis_type == "dominant_project_omission") %>%
    transmute(
      row_label = paste(country_label, "project omitted"),
      group = "Omit dominant project",
      median, lower_95, upper_95
    ),
  input_sensitivity %>%
    transmute(
      row_label = analysis_label,
      group = "Input refit",
      median, lower_95, upper_95
    )
) %>%
  mutate(
    row_label = factor(
      row_label,
      levels = rev(c(
        "Omit Australia", "Omit China", "Omit Japan",
        "Australia project omitted", "China project omitted", "Japan project omitted",
        "Primary", "Threshold 0.7", "Threshold 0.9",
        "Midpoint time", "Uniform time", "Alternative root"
      ))
    )
  )

p2e <- ggplot(omission_input, aes(median, row_label, colour = group)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.36) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0,
    linewidth = 0.45
  ) +
  geom_point(size = 1.45) +
  scale_colour_manual(
    values = c(
      `Omit country` = "#C75245",
      `Omit dominant project` = "#D58A2B",
      `Input refit` = "#6B83A9"
    ),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0.95, 1.38), breaks = c(1.0, 1.1, 1.2, 1.3)) +
  labs(x = "Sensitivity multiplier", y = NULL) +
  theme(axis.text.y = element_text(size = 5.0))

p2f <- ggplot(scenario_summary, aes(median, country, colour = scenario, shape = model)) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.34) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0,
    position = position_dodge(width = 0.50),
    linewidth = 0.48
  ) +
  geom_point(position = position_dodge(width = 0.50), size = 1.55, stroke = 0.42) +
  scale_colour_manual(
    values = c(
      `No post-2022 exposure` = "#D58A2B",
      `MT28-associated multiplier set to reference` = "#C75245"
    ),
    name = NULL
  ) +
  scale_shape_manual(values = c(`Project-adjusted` = 16, `No project effects` = 1), name = NULL) +
  scale_x_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1), labels = label_percent()) +
  labs(x = "Conditional decrease from baseline", y = NULL) +
  theme(legend.position = "bottom")

p2g <- ggplot(recovery, aes(parameter_type, metric_label, fill = pass)) +
  geom_tile(colour = "white", linewidth = 0.65) +
  geom_text(
    aes(label = paste0(display, "\n", ifelse(pass, "PASS", "FAIL"))),
    size = 1.95,
    lineheight = 0.86,
    fontface = "bold",
    colour = ifelse(recovery$pass, "#173D33", "#692C28"),
    family = font_family
  ) +
  scale_fill_manual(values = c(`TRUE` = "#BFDCCE", `FALSE` = "#E8C1BC"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold"),
    panel.border = element_blank()
  )

figure2 <- p2b + p2c + p2d + p2e + p2f + p2g +
  plot_layout(
    design = "
AABB
CCDD
EEFF
",
    widths = c(1, 1, 1, 1),
    heights = c(1.05, 1.15, 1.00),
    guides = "collect",
    tag_level = "new"
  ) +
  publication_annotation() &
  theme(legend.position = "none")

save_letter_figure(figure2, "Figure_2_growth_robustness", 183, 190)

render_manifest <- tibble::tribble(
  ~figure, ~stem, ~width_mm, ~height_mm, ~panels,
  "Figure 1", "Figure_1_observation_structure", 183, 190, "A-E",
  "Figure 2", "Figure_2_growth_robustness", 183, 190, "A-F"
)
readr::write_tsv(render_manifest, file.path(output_dir, "RENDER_MANIFEST.tsv"))

source_manifest <- tibble::tribble(
  ~figure, ~source_file,
  "Figure 1", "figures/source_data/figure2_tree_manifest.tsv",
  "Figure 1", "figures/source_data/figure2a_tree_tip_metadata.tsv",
  "Figure 1", "figures/source_data/figure2b_tip_ancestry_support.tsv",
  "Figure 1", "figures/source_data/figure1a_cases.tsv",
  "Figure 1", "figures/source_data/figure1b_annual_genomes.tsv",
  "Figure 1", "figures/source_data/figure3d_selection_cap_weighted_l10207_shares.tsv",
  "Figure 1", "figures/source_data/figure4e_australia_ct_curve.tsv",
  "Figure 2", "figures/source_data/figure3a_lineage_growth_main.tsv",
  "Figure 2", "figures/source_data/figure3b_l10207_pairwise_growth.tsv",
  "Figure 2", "figures/source_data/figure3c_l10207_growth_robustness.tsv",
  "Figure 2", "figures/source_data/figure4d_counterfactual_summary.tsv",
  "Figure 2", "figures/source_data/figure4f_identifiability_recovery.tsv",
  "Figure 2", "results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv"
)
readr::write_tsv(source_manifest, file.path(output_dir, "LETTER_FIGURE_SOURCE_FILES.tsv"))
