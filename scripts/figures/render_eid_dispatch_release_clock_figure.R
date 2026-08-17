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
      strip.text = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 7.2, face = "bold"),
      plot.subtitle = element_text(size = 6.1),
      panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey88"),
      panel.grid.minor = element_blank()
    )
)

out_dir <- "figures/eid"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- "figures/source_data"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

country_levels <- c("CHN", "JPN", "AUS")
country_labels <- c(CHN = "China", JPN = "Japan", AUS = "Australia")
clock_cols <- c(collection = "#2F6F9F", public = "#B64A4A")

cases_source <- read_tsv("data/derived/country_month_cases.tsv", show_col_types = FALSE) %>%
  filter(
    country_iso3 %in% country_levels,
    as.character(case_data_available) %in% c("TRUE", "True", "true", "1")
  ) %>%
  mutate(
    model_month = as.Date(model_month),
    cases = as.numeric(cases),
    country_label = country_labels[country_iso3]
  ) %>%
  filter(model_month >= as.Date("2023-01-01"))

cases <- cases_source %>%
  mutate(country_label = factor(country_label, levels = country_labels[country_levels]))

thresholds_source <- read_tsv("results/public_availability/case_thresholds.tsv", show_col_types = FALSE) %>%
  filter(country_iso3 %in% country_levels) %>%
  transmute(
    country_iso3,
    country_label = country_labels[country_iso3],
    cases_2019_max,
    threshold_date = as.Date(first_post2022_month_above_2019_max),
    peak_date = as.Date(post2022_peak_month)
  )

thresholds <- thresholds_source %>%
  mutate(
    country_label = factor(country_label, levels = country_labels[country_levels])
  )

selected_detection_source <- read_tsv(
  "results/public_availability/eid_detection_clock_shift.tsv",
  show_col_types = FALSE
) %>%
  filter(country_iso3 %in% country_levels) %>%
  mutate(
    collection_detection_lower = as.Date(collection_detection_lower),
    collection_detection_upper = as.Date(collection_detection_upper),
    public_detection_date = as.Date(public_detection_date)
  )

selected_detection <- selected_detection_source %>%
  mutate(country_label = factor(country_label, levels = country_labels[country_levels]))

availability_source <- read_tsv("data/derived/public_genome_availability.tsv", show_col_types = FALSE) %>%
  filter(
    country_iso3 %in% country_levels,
    collection_lower >= "2023-01-01",
    !is.na(public_date),
    primary_model_lineage_id == "L1_02.07"
  ) %>%
  mutate(
    country_label = country_labels[country_iso3],
    lag_min_days = as.numeric(lag_min_days),
    lag_max_days = as.numeric(lag_max_days)
  )

availability <- availability_source %>%
  mutate(country_label = factor(country_label, levels = country_labels[country_levels])) %>%
  group_by(country_label) %>%
  arrange(lag_max_days, .by_group = TRUE) %>%
  mutate(
    country_index = as.numeric(country_label),
    x_plot = country_index + seq(-0.18, 0.18, length.out = n())
  ) %>%
  ungroup()

shift_source <- selected_detection_source %>%
  mutate(
    collection_lead_min = as.numeric(collection_lead_to_case_threshold_min_days),
    collection_lead_max = as.numeric(collection_lead_to_case_threshold_max_days),
    public_lead = as.numeric(public_lead_to_case_threshold_days)
  )

shift <- shift_source %>%
  mutate(
    country_label = factor(country_label, levels = rev(country_labels[country_levels])),
    collection_label = paste0(
      if_else(collection_lead_min > 0, "+", ""), collection_lead_min,
      " to ", if_else(collection_lead_max > 0, "+", ""), collection_lead_max, " d"
    ),
    public_label = paste0(if_else(public_lead > 0, "+", ""), public_lead, " d")
  )

write_tsv(cases_source, file.path(source_dir, "eid_figure1a_cases.tsv"))
write_tsv(selected_detection_source, file.path(source_dir, "eid_figure1a_selected_detection.tsv"))
write_tsv(availability_source, file.path(source_dir, "eid_figure1b_release_lags.tsv"))
write_tsv(shift_source, file.path(source_dir, "eid_figure1c_clock_shift.tsv"))

p_a <- ggplot(cases, aes(model_month, cases)) +
  geom_rect(
    data = selected_detection,
    aes(
      xmin = collection_detection_lower,
      xmax = collection_detection_upper,
      ymin = -Inf,
      ymax = Inf,
      fill = "Collection interval"
    ),
    inherit.aes = FALSE,
    alpha = 0.16
  ) +
  geom_line(linewidth = 0.35, colour = "grey30") +
  geom_point(size = 0.6, colour = "grey30") +
  geom_hline(data = thresholds, aes(yintercept = cases_2019_max), linetype = "dotted", linewidth = 0.25, colour = "grey45") +
  geom_vline(data = thresholds, aes(xintercept = threshold_date), linetype = "dashed", linewidth = 0.25, colour = "black") +
  geom_vline(data = thresholds, aes(xintercept = peak_date), linetype = "dotdash", linewidth = 0.25, colour = "grey40") +
  geom_vline(
    data = selected_detection,
    aes(xintercept = public_detection_date, colour = "Public archive"),
    linewidth = 0.5
  ) +
  facet_wrap(~country_label, ncol = 1, scales = "free_y") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y", expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = c("Collection interval" = clock_cols[["collection"]]), name = NULL) +
  scale_colour_manual(values = c("Public archive" = clock_cols[["public"]]), name = NULL) +
  labs(
    title = "A. Cases and genomic timing",
    x = NULL,
    y = "Monthly reported cases"
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left"
  )

p_b_labels <- availability %>%
  group_by(country_label) %>%
  summarise(
    country_index = first(country_index),
    median_min = median(lag_min_days),
    median_max = median(lag_max_days),
    label = paste0("median ", median_min, "–", median_max, " d"),
    y = max(lag_max_days, na.rm = TRUE),
    .groups = "drop"
  )

p_b <- ggplot(availability) +
  geom_segment(
    aes(x = x_plot, xend = x_plot, y = lag_min_days, yend = lag_max_days),
    linewidth = 0.28,
    alpha = 0.38,
    colour = clock_cols[["public"]]
  ) +
  geom_point(aes(x = x_plot, y = lag_min_days), size = 0.35, alpha = 0.45, colour = clock_cols[["public"]]) +
  geom_linerange(
    data = p_b_labels,
    aes(x = country_index, ymin = median_min, ymax = median_max),
    inherit.aes = FALSE,
    linewidth = 1.15,
    colour = "#763232"
  ) +
  geom_text(
    data = p_b_labels,
    aes(x = country_index, y = y, label = label),
    inherit.aes = FALSE,
    vjust = -0.55,
    size = 1.75,
    colour = "grey20"
  ) +
  scale_x_continuous(breaks = seq_along(country_levels), labels = country_labels[country_levels]) +
  scale_y_continuous(labels = label_number(suffix = " d"), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    title = "B. Collection-to-public lag intervals",
    x = NULL,
    y = "Lag"
  ) +
  theme(legend.position = "none")

p_c <- ggplot(shift, aes(y = country_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey35") +
  geom_segment(
    aes(x = collection_lead_min, xend = collection_lead_max, yend = country_label, colour = "Collection interval"),
    linewidth = 1.0,
    lineend = "round"
  ) +
  geom_point(aes(x = public_lead, colour = "Public archive"), size = 1.8) +
  geom_text(
    aes(x = (collection_lead_min + collection_lead_max) / 2, label = collection_label),
    nudge_y = 0.16,
    size = 1.75,
    colour = clock_cols[["collection"]],
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = public_lead, label = public_label),
    nudge_y = -0.16,
    size = 1.75,
    colour = clock_cols[["public"]],
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Collection interval" = clock_cols[["collection"]], "Public archive" = clock_cols[["public"]]),
    name = NULL
  ) +
  scale_x_continuous(labels = label_number(suffix = " d")) +
  labs(
    title = "C. Timing relative to case threshold",
    x = "Days before or after first month above 2019 maximum",
    y = NULL
  ) +
  theme(
    legend.position = "top",
    legend.justification = "left"
  )

combined <- (p_a | (p_b / p_c)) +
  plot_layout(widths = c(1.05, 1.1), heights = c(1, 0.9)) &
  theme(plot.margin = margin(5, 6, 5, 5))

base <- file.path(out_dir, "Figure_1_release_clock_pertussis_eid")
svglite::svglite(paste0(base, ".svg"), width = 178 / 25.4, height = 120 / 25.4)
print(combined)
dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = 178 / 25.4, height = 120 / 25.4, family = "sans")
print(combined)
dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = 178 / 25.4, height = 120 / 25.4, units = "in", res = 600)
print(combined)
dev.off()
ragg::agg_png(paste0(base, ".png"), width = 178 / 25.4, height = 120 / 25.4, units = "in", res = 300)
print(combined)
dev.off()

message("Wrote EID figure exports to ", out_dir)
