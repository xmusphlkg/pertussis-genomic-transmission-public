#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(scales)
  library(tidyr)
})

theme_set(
  theme_classic(base_size = 7, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      strip.background = element_blank(),
      strip.text = element_text(size = 6.8, face = "bold"),
      plot.title = element_text(size = 7.3, face = "bold"),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.minor = element_blank()
    )
)

out_dir <- "figures/eid"
source_dir <- "figures/source_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

country_levels <- c("CHN", "JPN", "AUS")
country_labels <- c(CHN = "China", JPN = "Japan", AUS = "Australia")
palette <- c(
  possible = "#8FB9D4",
  definite = "#2F6F9F",
  public = "#C75B39",
  comparator = "#707070"
)

cumulative_source <- read_tsv(
  "results/public_availability/eid_cumulative_visibility.tsv",
  show_col_types = FALSE
) %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    event_date = as.Date(event_date),
    resurgence_milestone_date = as.Date(resurgence_milestone_date),
    resurgence_milestone_end = as.Date(resurgence_milestone_end),
    post2022_peak_date = as.Date(post2022_peak_date),
    post2022_peak_end = as.Date(post2022_peak_end),
    country_label = factor(country_labels[country_iso3], levels = country_labels[country_levels])
  )

cumulative <- cumulative_source %>%
  select(
    country_iso3, country_label, event_date, resurgence_milestone_date, resurgence_milestone_end,
    post2022_peak_date, post2022_peak_end,
    n_possibly_collected, n_definitely_collected, n_publicly_available
  ) %>%
  pivot_longer(
    cols = starts_with("n_"),
    names_to = "series",
    values_to = "cumulative_records"
  ) %>%
  mutate(
    series = factor(
      series,
      levels = c("n_possibly_collected", "n_definitely_collected", "n_publicly_available"),
      labels = c("Possibly collected", "Definitely collected", "Public sequence available")
    )
  )

milestones <- cumulative_source %>%
  distinct(country_label, resurgence_milestone_date, resurgence_milestone_end, post2022_peak_date, post2022_peak_end)

case_source <- read_tsv(
  "data/derived/country_month_cases.tsv",
  show_col_types = FALSE
) %>%
  filter(
    country_iso3 %in% country_levels,
    as.Date(model_month) >= as.Date("2023-01-01"),
    as.Date(model_month) <= as.Date("2025-12-01"),
    case_data_available
  ) %>%
  transmute(
    country_iso3,
    country_label = factor(country_labels[country_iso3], levels = country_labels[country_levels]),
    model_month = as.Date(model_month),
    cases = as.numeric(cases),
    source_url,
    source_file,
    data_freeze_date
  )

x_limits <- range(c(cumulative_source$event_date, case_source$model_month), na.rm = TRUE)

shift_source <- read_tsv(
  "results/public_availability/eid_detection_clock_shift.tsv",
  show_col_types = FALSE
) %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    country_label = factor(country_labels[country_iso3], levels = rev(country_labels[country_levels])),
    collection_relative_min = as.numeric(collection_relative_to_resurgence_min_days),
    collection_relative_max = as.numeric(collection_relative_to_resurgence_max_days),
    public_relative_min = as.numeric(public_relative_to_resurgence_min_days),
    public_relative_max = as.numeric(public_relative_to_resurgence_max_days),
    public_relative_mid = (public_relative_min + public_relative_max) / 2,
    collection_label = paste0(collection_relative_min, " to ", collection_relative_max, " d"),
    public_label = paste0(public_relative_min, " to ", public_relative_max, " d")
  )

milestone_row_levels <- c(
  "China Resurgence", "China Peak",
  "Japan Resurgence", "Japan Peak",
  "Australia Resurgence", "Australia Peak"
)

milestone_source <- read_tsv(
  "results/public_availability/eid_milestone_visibility.tsv",
  show_col_types = FALSE
) %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    country_label = factor(country_labels[country_iso3], levels = country_labels[country_levels]),
    milestone_label = recode(
      milestone,
      first_post2022_month_above_2019_max = "Resurgence",
      post2022_peak_month = "Peak"
    ),
    row_label = factor(
      paste(country_labels[country_iso3], milestone_label),
      levels = rev(milestone_row_levels)
    ),
    definitely_collected_by_milestone = as.numeric(definitely_collected_by_milestone),
    possibly_collected_by_milestone = as.numeric(possibly_collected_by_milestone),
    public_by_milestone = as.numeric(public_by_milestone),
    collected_mid = (definitely_collected_by_milestone + possibly_collected_by_milestone) / 2,
    collected_label = if_else(
      definitely_collected_by_milestone == possibly_collected_by_milestone,
      as.character(definitely_collected_by_milestone),
      paste0(definitely_collected_by_milestone, "-", possibly_collected_by_milestone)
    ),
    public_label = as.character(public_by_milestone),
    collected_label_x = if_else(
      definitely_collected_by_milestone == 0 & possibly_collected_by_milestone == 0,
      3,
      collected_mid
    )
  )

comparison_source <- read_tsv(
  "results/public_availability/eid_project_lineage_comparison.tsv",
  show_col_types = FALSE
) %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    country_label = factor(country_labels[country_iso3], levels = country_labels[country_levels]),
    stratum = paste0(country_iso3, " · ", project_id, ", ", collection_year)
  )

comparison <- bind_rows(
  comparison_source %>%
    transmute(
      country_label,
      stratum,
      group = "MT28-associated lineage",
      n = n_target_lineage,
      lag_min = as.numeric(target_median_lag_min_days),
      lag_max = as.numeric(target_median_lag_max_days)
    ),
  comparison_source %>%
    transmute(
      country_label,
      stratum,
      group = "Other lineages",
      n = n_comparator_lineages,
      lag_min = as.numeric(comparator_median_lag_min_days),
      lag_max = as.numeric(comparator_median_lag_max_days)
    )
) %>%
  mutate(
    group = factor(group, levels = c("MT28-associated lineage", "Other lineages")),
    midpoint = (lag_min + lag_max) / 2
  )

write_tsv(cumulative_source, file.path(source_dir, "eid_figure1a_cumulative_visibility.tsv"))
write_tsv(case_source, file.path(source_dir, "eid_figure1a_cases.tsv"))
write_tsv(shift_source, file.path(source_dir, "eid_figure1b_relative_timing.tsv"))
write_tsv(milestone_source, file.path(source_dir, "eid_figure1c_milestone_visibility.tsv"))
write_tsv(comparison_source, file.path(source_dir, "eid_figure1c_project_lineage_comparison.tsv"))

p_a_cases <- ggplot(case_source, aes(model_month, cases, group = country_label)) +
  geom_vline(
    data = milestones,
    aes(xintercept = resurgence_milestone_date),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.3,
    colour = "black"
  ) +
  geom_vline(
    data = milestones,
    aes(xintercept = post2022_peak_date),
    inherit.aes = FALSE,
    linetype = "dotdash",
    linewidth = 0.3,
    colour = "grey45"
  ) +
  geom_area(fill = "grey85", colour = NA) +
  geom_line(linewidth = 0.48, colour = "grey35") +
  facet_wrap(~country_label, ncol = 3, scales = "free_y") +
  scale_x_date(
    limits = x_limits,
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    n.breaks = 3,
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "A",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey90"),
    plot.margin = margin(4, 6, 0, 5)
  )

p_a_cumulative <- ggplot(cumulative, aes(event_date, cumulative_records, colour = series, linetype = series)) +
  geom_vline(
    data = milestones,
    aes(
      xintercept = resurgence_milestone_date,
      colour = "Resurgence milestone",
      linetype = "Resurgence milestone"
    ),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_vline(
    data = milestones,
    aes(
      xintercept = post2022_peak_date,
      colour = "Post-2022 peak",
      linetype = "Post-2022 peak"
    ),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  geom_step(linewidth = 0.62, direction = "hv") +
  facet_wrap(~country_label, ncol = 3, scales = "free_y") +
  scale_colour_manual(
    values = c(
      "Possibly collected" = palette[["possible"]],
      "Definitely collected" = palette[["definite"]],
      "Public sequence available" = palette[["public"]],
      "Resurgence milestone" = "black",
      "Post-2022 peak" = "grey45"
    ),
    breaks = c(
      "Possibly collected", "Definitely collected", "Public sequence available",
      "Resurgence milestone", "Post-2022 peak"
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "Possibly collected" = "22",
      "Definitely collected" = "solid",
      "Public sequence available" = "solid",
      "Resurgence milestone" = "dashed",
      "Post-2022 peak" = "dotdash"
    ),
    breaks = c(
      "Possibly collected", "Definitely collected", "Public sequence available",
      "Resurgence milestone", "Post-2022 peak"
    ),
    name = NULL
  ) +
  scale_x_date(
    limits = x_limits,
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(labels = label_number(accuracy = 1), expand = expansion(mult = c(0, 0.07))) +
  labs(
    x = NULL,
    y = "Cumulative MT28-associated records"
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    axis.text.x = element_text(size = 5.5),
    strip.text = element_blank(),
    plot.margin = margin(0, 6, 1, 5)
  ) +
  guides(
    colour = guide_legend(nrow = 1, byrow = TRUE),
    linetype = guide_legend(nrow = 1, byrow = TRUE)
  )

p_a <- wrap_plots(
  p_a_cases,
  p_a_cumulative,
  ncol = 1,
  heights = c(0.36, 1.0)
)

p_b <- ggplot(shift_source, aes(y = country_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey35") +
  geom_segment(
    aes(
      x = collection_relative_min,
      xend = collection_relative_max,
      yend = country_label,
      colour = "Specimen accumulation"
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = public_relative_min,
      xend = public_relative_max,
      yend = country_label,
      colour = "Public sequence availability"
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  geom_point(aes(x = public_relative_mid, colour = "Public sequence availability"), size = 1.8) +
  geom_text(
    aes(x = (collection_relative_min + collection_relative_max) / 2, label = collection_label),
    nudge_y = 0.17,
    size = 1.75,
    colour = palette[["definite"]]
  ) +
  geom_text(
    aes(x = public_relative_mid, label = public_label),
    nudge_y = -0.17,
    size = 1.75,
    colour = palette[["public"]]
  ) +
  scale_colour_manual(
    values = c(
      "Specimen accumulation" = palette[["definite"]],
      "Public sequence availability" = palette[["public"]]
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = seq(-600, 600, 200),
    labels = label_number(accuracy = 1, style_negative = "minus"),
    expand = expansion(mult = c(0.06, 0.08))
  ) +
  labs(
    title = "B",
    x = "Days relative to resurgence milestone month\nNegative values indicate earlier timing",
    y = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.major.y = element_blank()
  )

p_project_comparison <- ggplot(
  comparison,
  aes(y = reorder(stratum, midpoint), colour = group, shape = group)
) +
  geom_linerange(
    aes(xmin = lag_min, xmax = lag_max),
    orientation = "y",
    position = position_dodge(width = 0.55),
    linewidth = 0.7
  ) +
  geom_point(
    aes(x = midpoint),
    position = position_dodge(width = 0.55),
    size = 1.5,
    fill = "white"
  ) +
  geom_text(
    aes(x = lag_max + 24, label = paste0("n=", n)),
    position = position_dodge(width = 0.55),
    size = 1.55,
    hjust = 0,
    show.legend = FALSE
  ) +
  facet_grid(country_label ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(
    values = c(
      "MT28-associated lineage" = palette[["definite"]],
      "Other lineages" = palette[["comparator"]]
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("MT28-associated lineage" = 21, "Other lineages" = 24),
    name = NULL
  ) +
  scale_x_continuous(
    labels = label_number(suffix = " d"),
    expand = expansion(mult = c(0.06, 0.18))
  ) +
  labs(
    title = "C",
    x = "Median collection-to-public availability interval",
    y = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 5.2),
    strip.text.y = element_blank()
  )

p_c <- ggplot(milestone_source, aes(y = row_label)) +
  geom_segment(
    aes(
      x = definitely_collected_by_milestone,
      xend = possibly_collected_by_milestone,
      yend = row_label
    ),
    linewidth = 1.1,
    lineend = "round",
    colour = palette[["definite"]]
  ) +
  geom_point(
    aes(x = definitely_collected_by_milestone),
    size = 1.55,
    colour = palette[["definite"]]
  ) +
  geom_point(
    aes(x = possibly_collected_by_milestone),
    size = 1.65,
    shape = 21,
    stroke = 0.45,
    fill = "white",
    colour = palette[["definite"]]
  ) +
  geom_point(
    aes(x = public_by_milestone),
    size = 1.65,
    colour = palette[["public"]]
  ) +
  geom_text(
    aes(
      x = collected_label_x,
      label = collected_label
    ),
    nudge_y = 0.17,
    size = 1.55,
    colour = palette[["definite"]]
  ) +
  geom_text(
    aes(
      x = public_by_milestone + 3,
      label = public_label
    ),
    hjust = 0,
    nudge_y = -0.16,
    size = 1.55,
    colour = palette[["public"]]
  ) +
  scale_x_continuous(
    limits = c(-4, 122),
    breaks = c(0, 25, 50, 75, 100),
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(labels = function(x) sub(" ", "\n", x, fixed = TRUE)) +
  labs(
    title = "C",
    x = "MT28-associated records at milestone",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 5.35, lineheight = 0.85),
    axis.text.x = element_text(size = 5.5),
    axis.title.x = element_text(size = 6.2),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey90"),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(4, 5, 4, 4)
  )

bottom <- (p_b | p_c) + plot_layout(widths = c(1.05, 0.95))

combined <- wrap_plots(
  p_a,
  bottom,
  ncol = 1,
  heights = c(1.28, 0.92)
) &
  theme(plot.margin = margin(4, 6, 4, 5))

base <- file.path(out_dir, "Figure_1_release_clock_pertussis_eid")
width_in <- 178 / 25.4
height_in <- 125 / 25.4

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in)
print(combined)
dev.off()

grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_in, height = height_in, family = "sans")
print(combined)
dev.off()

ragg::agg_tiff(
  paste0(base, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600
)
print(combined)
dev.off()

ragg::agg_png(
  paste0(base, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300
)
print(combined)
dev.off()

message("Wrote EID figure exports and source data")
