#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtree)
  library(grid)
  library(patchwork)
  library(phangorn)
  library(readr)
  library(scales)
  library(systemfonts)
  library(tidyr)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
source_dir <- file.path(repo_dir, "figures", "source_data")
output_dir <- file.path(repo_dir, "figures", "main")
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
  FRA = "France", JPN = "Japan", Background = "Background"
)
country_colors <- c(
  AUS = "#3B78A8", BEL = "#C69232", CHN = "#C75245",
  FRA = "#7C65A8", JPN = "#2F8F83", Background = "#D4D4D4"
)
period_colors <- c(
  prepandemic = "#808080", pre_pandemic = "#808080",
  pandemic = "#9DBDD5", resurgence = "#DF7D6D"
)
lineage_order <- c("L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other")
lineage_colors <- c(
  L1_01.02 = "#3B78A8", L1_02.05 = "#7C65A8",
  L1_02.06 = "#2F8F83", L1_02.07 = "#C75245",
  Other = "#777777", Background = "#D0D0D0"
)
scenario_colors <- c(
  Observed = "#202020", Baseline = "#3B78A8",
  `No new introduction` = "#D58A2B",
  `Equal lineage growth` = "#7C65A8"
)

theme_pub <- function(base_size = 6.6) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#222222"),
      axis.ticks.length = unit(1.2, "mm"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = "#222222"),
      legend.title = element_text(size = base_size - 0.1, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(4.5, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      legend.box.spacing = unit(0.8, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.6, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555"),
      plot.caption = element_text(size = base_size - 0.7, colour = "#555555"),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 5)
    )
}

theme_set(theme_pub())

save_figure <- function(plot, stem, width_mm, height_mm) {
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
    res = 600,
    background = "white"
  )
  print(plot)
  dev.off()
}

read_source <- function(filename) {
  readr::read_tsv(
    file.path(source_dir, filename),
    show_col_types = FALSE,
    progress = FALSE
  )
}

country_factor <- function(x, include_pooled = FALSE) {
  levels <- country_order
  labels <- unname(country_labels[country_order])
  if (include_pooled) {
    levels <- c(levels, "Pooled")
    labels <- c(labels, "Pooled")
  }
  factor(x, levels = levels, labels = labels)
}

publication_annotation <- function() {
  plot_annotation()
}

tag_panel <- function(plot, label) {
  plot +
    labs(tag = label) +
    theme(
      plot.tag = element_text(
        size = 8, face = "bold", family = font_family,
        margin = margin(0, 1.5, 1.5, 0)
      ),
      plot.tag.location = "panel",
      plot.tag.position = c(0.012, 0.982)
    )
}

# ---------------------------------------------------------------------------
# Figure 1: data structure and model
# ---------------------------------------------------------------------------

cases <- read_source("figure1a_cases.tsv") %>%
  mutate(model_month = as.Date(model_month)) %>%
  filter(
    country_iso3 %in% c("AUS", "CHN", "JPN"),
    model_month >= as.Date("2019-01-01"),
    model_month <= as.Date("2025-12-01"),
    case_data_available
  ) %>%
  mutate(country = country_factor(country_iso3))

p1a <- ggplot(cases, aes(model_month, cases, colour = country_iso3)) +
  annotate(
    "rect",
    xmin = as.Date("2020-03-01"), xmax = as.Date("2022-12-31"),
    ymin = -Inf, ymax = Inf, fill = "#DCE9F2", alpha = 0.55
  ) +
  geom_line(linewidth = 0.48, lineend = "round") +
  facet_wrap(~country, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = country_colors, guide = "none") +
  scale_x_date(
    date_breaks = "2 years", date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale()),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(x = NULL, y = "Reported cases per month")

monthly_genomes <- read_source("figure1b_monthly_genomes.tsv") %>%
  mutate(
    model_month = as.Date(model_month),
    country = country_factor(country_iso3),
    epidemic_period = factor(
      epidemic_period,
      levels = c("prepandemic", "pandemic", "resurgence"),
      labels = c("Pre-pandemic", "Pandemic", "Resurgence")
    )
  ) %>%
  filter(model_month >= as.Date("2000-01-01"))

p1b <- ggplot(
  monthly_genomes,
  aes(model_month, n_sampled_genomes, fill = epidemic_period)
) +
  geom_col(width = 26, colour = NA) +
  facet_wrap(~country, ncol = 5, scales = "free_y") +
  scale_fill_manual(
    values = c(
      `Pre-pandemic` = period_colors[["prepandemic"]],
      Pandemic = period_colors[["pandemic"]],
      Resurgence = period_colors[["resurgence"]]
    ),
    name = "Sampling period"
  ) +
  scale_x_date(
    breaks = as.Date(c("2000-01-01", "2010-01-01", "2020-01-01")),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Sampling month", y = "Sampled focal genomes") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

cohort_flow <- read_source("figure1c_cohort_flow.tsv") %>%
  mutate(
    stage = factor(
      stage,
      levels = c("selected", "uniform_qc_pass", "final_alignment_pass"),
      labels = c("Selected", "Uniform QC", "Final tree")
    )
  )

p1c <- ggplot(cohort_flow, aes(stage, n_genomes, fill = stage)) +
  geom_col(width = 0.62, colour = "#333333", linewidth = 0.25) +
  geom_text(aes(label = comma(n_genomes)), vjust = -0.35, size = 2.4) +
  scale_fill_manual(values = c("#C9D6DF", "#7AA6C2", "#3B78A8"), guide = "none") +
  scale_y_continuous(
    limits = c(0, 1320), breaks = c(0, 500, 1000),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Genomes") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

model_nodes <- tibble::tribble(
  ~id, ~x, ~y, ~label, ~family,
  "local", 0.7, 3.0, "Local\ntransmission", "driver",
  "import", 0.7, 1.0, "External ancestry\nexposure", "driver",
  "latent", 2.5, 2.0, "Latent country–month–\nlineage infections", "latent",
  "cases", 4.4, 3.0, "Reported cases", "observed",
  "genomes", 4.4, 1.0, "Project-stratified\ngenome counts", "observed",
  "selection", 2.7, 0.15, "Project selection +\nsequencing success", "filter"
)
model_edges <- tibble::tribble(
  ~from, ~to,
  "local", "latent",
  "import", "latent",
  "latent", "cases",
  "latent", "genomes",
  "selection", "genomes"
) %>%
  left_join(model_nodes %>% select(id, x, y), by = c("from" = "id")) %>%
  rename(x_from = x, y_from = y) %>%
  left_join(model_nodes %>% select(id, x, y), by = c("to" = "id")) %>%
  rename(x_to = x, y_to = y)

p1d <- ggplot() +
  geom_curve(
    data = model_edges,
    aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
    curvature = 0.06, colour = "#6B6B6B", linewidth = 0.42,
    arrow = arrow(type = "closed", length = unit(1.6, "mm"))
  ) +
  geom_label(
    data = model_nodes,
    aes(x, y, label = label, fill = family),
    size = 2.25, linewidth = 0.25, label.padding = unit(1.3, "mm"),
    lineheight = 0.92, family = font_family
  ) +
  scale_fill_manual(
    values = c(
      driver = "#F2D7C1", latent = "#CADDEA",
      observed = "#D8E8DF", filter = "#E3D9EA"
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 5.2), ylim = c(-0.25, 3.5), clip = "off") +
  theme_void(base_family = font_family, base_size = 6.6) +
  theme(
    plot.margin = margin(7, 4, 7, 4)
  )

fig1_design <- "
AAAC
BBBD
"
p1a <- tag_panel(p1a, "A")
p1b <- tag_panel(p1b, "B")
p1c <- tag_panel(p1c, "C")
p1d <- tag_panel(p1d, "D")
figure1 <- p1a + p1b + p1c + p1d +
  plot_layout(
    design = fig1_design,
    widths = c(1, 1, 1, 1.25),
    heights = c(1.6, 1),
    guides = "collect"
  ) +
  publication_annotation() &
  theme(legend.position = "bottom")

save_figure(figure1, "Figure_1_data_and_model", 183, 150)

# ---------------------------------------------------------------------------
# Figure 2: core-SNP phylogeny and resurgence ancestry
# ---------------------------------------------------------------------------

tree_manifest <- read_source("figure2_tree_manifest.tsv")
tree <- phangorn::midpoint(ape::read.tree(tree_manifest$tree_file[[1]]))
tree_meta <- read_source("figure2a_tree_tip_metadata.tsv") %>%
  select(
    tree_sample_id, tree_role, country_iso3, epidemic_period,
    date_lower, date_upper
  )
tip_ancestry <- read_source("figure2b_tip_ancestry_support.tsv") %>%
  select(
    tree_sample_id, primary_model_lineage_id,
    post_reseeding_support, local_persistence_support
  )
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
    period_display = factor(
      epidemic_period,
      levels = c("prepandemic", "pandemic", "resurgence"),
      labels = c("Pre-pandemic", "Pandemic", "Resurgence")
    )
  )

p2a <- ggtree(tree, size = 0.12, colour = "#A7A7A7", ladderize = TRUE) %<+%
  tree_meta
p2a$data <- p2a$data %>%
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

p2a <- p2a +
  geom_point(
    data = p2a$data %>% filter(!isTip, high_support),
    aes(x = x, y = y), inherit.aes = FALSE,
    shape = 16, size = 0.32, colour = "#333333", alpha = 0.75
  ) +
  geom_tippoint(
    aes(
      fill = country_display, colour = lineage_display,
      size = period_display
    ),
    shape = 21, stroke = 0.22, alpha = 0.92
  ) +
  scale_fill_manual(
    values = country_colors,
    breaks = c(country_order, "Background"),
    labels = country_labels[c(country_order, "Background")],
    name = "Country",
    guide = "none"
  ) +
  scale_colour_manual(
    values = lineage_colors,
    breaks = c(lineage_order, "Background"),
    name = "Frozen lineage",
    guide = "none"
  ) +
  scale_size_manual(
    values = c(`Pre-pandemic` = 0.42, Pandemic = 0.62, Resurgence = 0.88),
    name = "Period",
    guide = "none"
  ) +
  geom_treescale(
    x = 0, y = 12, width = 0.001,
    fontsize = 2.1, linesize = 0.35
  ) +
  theme_tree2(base_family = font_family, base_size = 6.4) +
  theme(
    axis.title.x = element_text(size = 6.4),
    axis.text.x = element_text(size = 5.8),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(5, 7, 4, 4)
  ) +
  xlab("Substitutions per site")

ancestry_long <- read_source("figure2b_tip_ancestry_support.tsv") %>%
  filter(epidemic_period == "resurgence", country_iso3 %in% country_order) %>%
  select(
    tree_sample_id, country_iso3,
    `Local persistence` = local_persistence_support,
    Reseeding = post_reseeding_support
  ) %>%
  pivot_longer(
    c(`Local persistence`, Reseeding),
    names_to = "ancestry_type", values_to = "support"
  ) %>%
  mutate(country = country_factor(country_iso3))
ancestry_n <- ancestry_long %>%
  distinct(tree_sample_id, country) %>%
  count(country, name = "n") %>%
  mutate(label = paste0("n=", n))

p2b <- ggplot(
  ancestry_long,
  aes(country, support, fill = ancestry_type)
) +
  geom_violin(
    position = position_dodge(width = 0.72), width = 0.66,
    scale = "width", trim = TRUE, colour = NA, alpha = 0.72
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.72), width = 0.16,
    outlier.shape = NA, linewidth = 0.28, alpha = 0.85
  ) +
  geom_text(
    data = ancestry_n,
    aes(country, 1.075, label = label),
    inherit.aes = FALSE, size = 2.05, colour = "#555555"
  ) +
  scale_fill_manual(
    values = c(`Local persistence` = "#3B78A8", Reseeding = "#D58A2B"),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(0, 1.12), breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Marginal ancestry support") +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ancestry_heat <- read_source("figure2c_country_lineage_ancestry.tsv") %>%
  mutate(
    country = country_factor(country_iso3),
    lineage = factor(primary_model_lineage_id, levels = lineage_order)
  )

p2c <- ggplot(
  ancestry_heat,
  aes(lineage, country, fill = mean_post_import_support)
) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(
    aes(label = paste0(number(mean_post_import_support, accuracy = 0.01), "\n(n=", n_resurgence_tips, ")")),
    size = 1.95, lineheight = 0.88,
    colour = ifelse(ancestry_heat$mean_post_import_support > 0.58, "white", "#222222")
  ) +
  scale_fill_gradient(
    low = "#F7F7F7", high = "#C65B45",
    limits = c(0, 1), labels = label_percent(),
    name = "Mean reseeding\nsupport"
  ) +
  labs(x = "Frozen lineage", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom",
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

p2a <- tag_panel(p2a, "A")
p2b <- tag_panel(p2b, "B")
p2c <- tag_panel(p2c, "C")
figure2 <- p2a | (p2b / p2c) +
  plot_layout(widths = c(1.62, 1), heights = c(1, 1)) +
  publication_annotation()

save_figure(figure2, "Figure_2_phylogeny_and_ancestry", 183, 160)

# ---------------------------------------------------------------------------
# Figure 3: post-reseeding edges and sampled clusters
# ---------------------------------------------------------------------------

edges <- read_source("figure3a_post_reseeding_edges.tsv") %>%
  mutate(
    destination = country_factor(destination_country),
    source = recode(
      top_source_state,
      Europe_other = "Other Europe",
      Unknown_other = "Unknown"
    )
  )
source_levels <- c("AUS", "BEL", "CHN", "FRA", "JPN", "Other Europe", "Unknown")
source_colors <- c(
  country_colors[c("AUS", "BEL", "CHN", "FRA", "JPN")],
  `Other Europe` = "#9DAAB2", Unknown = "#D8D8D8"
)

edge_counts <- edges %>%
  count(destination, source, name = "n") %>%
  mutate(source = factor(source, levels = source_levels))

p3a <- ggplot(edge_counts, aes(n, destination, fill = source)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.25) +
  geom_text(
    data = edges %>% count(destination, name = "total"),
    aes(total + 0.7, destination, label = paste0("n=", total)),
    inherit.aes = FALSE, hjust = 0, size = 2.15
  ) +
  scale_fill_manual(values = source_colors, name = "Top source state") +
  scale_x_continuous(
    limits = c(0, 34), breaks = seq(0, 30, 10),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "Strict post-reseeding edges", y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.position = "bottom")

p3b <- ggplot(
  edges,
  aes(sample_span_days, n_attributed_tips)
) +
  geom_vline(xintercept = 180, linetype = "22", linewidth = 0.35, colour = "#777777") +
  geom_hline(yintercept = 5, linetype = "22", linewidth = 0.35, colour = "#777777") +
  geom_point(
    aes(
      fill = destination_country,
      size = transition_support,
      shape = successful_sampled_cluster
    ),
    colour = "#333333", stroke = 0.32, alpha = 0.86
  ) +
  scale_fill_manual(values = country_colors, guide = "none") +
  scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 24),
    labels = c(`FALSE` = "Other edge", `TRUE` = "Successful sampled cluster"),
    name = NULL
  ) +
  scale_size_continuous(range = c(1.2, 3.0), name = "Transition\nsupport") +
  scale_x_sqrt(
    breaks = c(0, 180, 365, 730),
    labels = c("0", "180", "365", "730")
  ) +
  scale_y_sqrt(breaks = c(1, 5, 10, 20)) +
  labs(
    x = "Sampling span (days; square-root scale)",
    y = "Attributed sampled genomes\n(square-root scale)"
  ) +
  guides(size = "none", shape = "none") +
  theme(legend.position = "bottom")

cluster_probability <- read_source("figure3c_cluster_probability.tsv")
pooled_success <- sum(cluster_probability$n_successful_sampled_clusters)
pooled_events <- sum(cluster_probability$n_high_support_events)
cluster_probability <- bind_rows(
  cluster_probability,
  tibble(
    destination_country = "Pooled",
    n_high_support_events = pooled_events,
    n_successful_sampled_clusters = pooled_success,
    posterior_mean_success_probability = (pooled_success + 1) / (pooled_events + 2),
    lower_95 = qbeta(0.025, pooled_success + 1, pooled_events - pooled_success + 1),
    upper_95 = qbeta(0.975, pooled_success + 1, pooled_events - pooled_success + 1)
  )
) %>%
  mutate(
    country = factor(
      destination_country,
      levels = c("AUS", "BEL", "FRA", "JPN", "Pooled"),
      labels = c("Australia", "Belgium", "France", "Japan", "Pooled")
    ),
    colour_key = if_else(destination_country == "Pooled", "Pooled", destination_country)
  )
cluster_probability_colors <- c(country_colors, Pooled = "#222222")

p3c <- ggplot(
  cluster_probability,
  aes(posterior_mean_success_probability, country, colour = colour_key)
) +
  geom_errorbarh(
    aes(xmin = lower_95, xmax = upper_95),
    height = 0, linewidth = 0.55
  ) +
  geom_point(size = 2.1) +
  geom_text(
    aes(
      x = pmin(upper_95 + 0.03, 0.41),
      label = paste0(n_successful_sampled_clusters, "/", n_high_support_events)
    ),
    hjust = 0, size = 2.0, colour = "#444444"
  ) +
  scale_colour_manual(values = cluster_probability_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, 0.46), breaks = c(0, 0.1, 0.2, 0.3, 0.4),
    labels = label_percent(accuracy = 1)
  ) +
  labs(x = "Sampled-cluster probability (95% CrI)", y = NULL)

thresholds <- read_source("figure3d_threshold_sensitivity.tsv") %>%
  select(
    transition_threshold,
    `Post-reseeding edges` = n_high_support_post_reseeding_edges,
    `Successful clusters` = n_successful_sampled_clusters
  ) %>%
  pivot_longer(
    -transition_threshold, names_to = "quantity", values_to = "count"
  )

p3d <- ggplot(
  thresholds,
  aes(transition_threshold, count, colour = quantity)
) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8) +
  geom_text(aes(label = count), nudge_y = 1.3, size = 2.05, show.legend = FALSE) +
  facet_wrap(~quantity, ncol = 1, scales = "free_y") +
  scale_colour_manual(
    values = c(`Post-reseeding edges` = "#3B78A8", `Successful clusters` = "#C75245"),
    guide = "none"
  ) +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9), limits = c(0.47, 0.93)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = "Transition-support threshold", y = "Count")

fig3_design <- "
AABB
CCDD
"
p3a <- tag_panel(p3a, "A")
p3b <- tag_panel(p3b, "B")
p3c <- tag_panel(p3c, "C")
p3d <- tag_panel(p3d, "D")
figure3 <- p3a + p3b + p3c + p3d +
  plot_layout(design = fig3_design, widths = c(1, 1, 1.12, 0.88)) +
  publication_annotation()

save_figure(figure3, "Figure_3_reseeding_and_clusters", 183, 135)

# ---------------------------------------------------------------------------
# Figure 4: lineage growth and sampling sensitivity
# ---------------------------------------------------------------------------

growth_main <- read_source("figure4a_lineage_growth_main.tsv") %>%
  mutate(
    lineage = factor(lineage, levels = rev(lineage_order)),
    highlight = lineage == "L1_02.07"
  )

p4a <- ggplot(growth_main, aes(median, lineage)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = lower_95, xmax = upper_95, colour = highlight),
    height = 0, linewidth = 0.75
  ) +
  geom_point(aes(fill = highlight), shape = 21, size = 2.5, stroke = 0.45) +
  geom_text(
    aes(
      x = upper_95 + 0.008,
      label = sprintf("%.3f", median),
      colour = highlight
    ),
    hjust = 0, size = 2.05, show.legend = FALSE
  ) +
  scale_colour_manual(values = c(`FALSE` = "#4B4B4B", `TRUE` = "#C75245"), guide = "none") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#C75245"), guide = "none") +
  scale_x_continuous(
    limits = c(0.88, 1.18), breaks = seq(0.9, 1.15, 0.05)
  ) +
  labs(x = "Relative net-growth multiplier (95% CrI)", y = NULL)

growth_sensitivity <- read_source("figure4b_lineage_growth_sensitivity.tsv") %>%
  select(
    lineage,
    `Project-adjusted` = median_main,
    `No project effects` = median_no_project
  ) %>%
  pivot_longer(-lineage, names_to = "model", values_to = "median") %>%
  mutate(lineage = factor(lineage, levels = rev(lineage_order)))

p4b <- ggplot(growth_sensitivity, aes(median, lineage, group = lineage)) +
  geom_vline(xintercept = 1, linetype = "22", colour = "#777777", linewidth = 0.4) +
  geom_line(colour = "#B0B0B0", linewidth = 0.65) +
  geom_point(aes(colour = model, shape = model), size = 2.15, stroke = 0.45) +
  scale_colour_manual(
    values = c(`Project-adjusted` = "#C75245", `No project effects` = "#3B78A8"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(`Project-adjusted` = 16, `No project effects` = 1),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(0.91, 1.13), breaks = seq(0.92, 1.12, 0.04)
  ) +
  labs(x = "Posterior median multiplier", y = NULL) +
  theme(legend.position = "bottom")

latent_shares <- read_source("figure4c_monthly_latent_lineage_shares.tsv") %>%
  mutate(
    model_month = as.Date(model_month),
    country = country_factor(country_iso3),
    lineage = factor(lineage, levels = lineage_order)
  )

p4c <- ggplot(
  latent_shares,
  aes(model_month, median, colour = lineage)
) +
  annotate(
    "rect",
    xmin = as.Date("2020-03-01"), xmax = as.Date("2022-12-31"),
    ymin = -Inf, ymax = Inf, fill = "#E7EDF1", alpha = 0.60
  ) +
  geom_line(linewidth = 0.52, lineend = "round") +
  facet_wrap(~country, ncol = 1) +
  scale_colour_manual(values = lineage_colors, name = "Frozen lineage") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1), labels = label_percent()
  ) +
  labs(x = NULL, y = "Latent lineage share") +
  theme(legend.position = "bottom")

raw_corrected <- read_source("figure4d_raw_vs_corrected_shares.tsv") %>%
  mutate(
    country = country_factor(country_iso3),
    lineage = factor(lineage, levels = lineage_order)
  )

p4d <- ggplot(
  raw_corrected,
  aes(raw_public_tree_share, median, colour = lineage)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#888888", linewidth = 0.4) +
  geom_errorbar(
    aes(ymin = lower_95, ymax = upper_95),
    width = 0, linewidth = 0.42, alpha = 0.8
  ) +
  geom_point(size = 1.8) +
  facet_wrap(~country, ncol = 1) +
  scale_colour_manual(values = lineage_colors, guide = "none") +
  scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
  coord_equal() +
  labs(
    x = "Raw public-tree share",
    y = "Project-adjusted share (95% CrI)"
  )

fig4_design <- "
AABB
CCCD
CCCD
"
p4a <- tag_panel(p4a, "A")
p4b <- tag_panel(p4b, "B")
p4c <- tag_panel(p4c, "C")
p4d <- tag_panel(p4d, "D")
figure4 <- p4a + p4b + p4c + p4d +
  plot_layout(
    design = fig4_design,
    widths = c(1, 1, 1, 1.05),
    heights = c(1, 1.25, 1.25),
    guides = "collect"
  ) +
  publication_annotation() &
  theme(legend.position = "bottom")

save_figure(figure4, "Figure_4_lineage_growth_and_sampling", 183, 165)

# ---------------------------------------------------------------------------
# Figure 5: counterfactuals, Ct calibration, and recovery
# ---------------------------------------------------------------------------

counterfactual_monthly <- read_source("figure5abc_monthly_counterfactuals.tsv") %>%
  mutate(model_month = as.Date(model_month)) %>%
  filter(model_month >= as.Date("2022-01-01"))

make_country_counterfactual <- function(iso3) {
  x <- counterfactual_monthly %>% filter(country_iso3 == iso3)
  lines <- x %>%
    select(
      model_month,
      Observed = observed_cases,
      Baseline = baseline_median,
      `No new introduction` = no_import_median,
      `Equal lineage growth` = equal_lineage_median
    ) %>%
    pivot_longer(-model_month, names_to = "series", values_to = "cases")

  ggplot() +
    geom_ribbon(
      data = x,
      aes(model_month, ymin = fitted_lower_95, ymax = fitted_upper_95),
      fill = "#9CBED5", alpha = 0.28
    ) +
    geom_line(
      data = lines,
      aes(model_month, cases, colour = series, linetype = series),
      linewidth = 0.52, lineend = "round"
    ) +
    scale_colour_manual(values = scenario_colors, name = NULL) +
    scale_linetype_manual(
      values = c(
        Observed = "solid", Baseline = "solid",
        `No new introduction` = "22",
        `Equal lineage growth` = "13"
      ),
      name = NULL
    ) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0.02, 0.10))
    ) +
    labs(x = NULL, y = "Reported cases per month") +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

p5a <- make_country_counterfactual("AUS")
p5b <- make_country_counterfactual("CHN")
p5c <- make_country_counterfactual("JPN")

counterfactual_summary <- read_source("figure5d_counterfactual_summary.tsv") %>%
  filter(
    scenario %in% c(
      "no_new_introduction_case_reduction_fraction",
      "lineage_difference_effect_fraction"
    )
  ) %>%
  mutate(
    country = country_factor(country_iso3),
    scenario_label = recode(
      scenario,
      no_new_introduction_case_reduction_fraction = "No new introduction",
      lineage_difference_effect_fraction = "Equal lineage growth"
    )
  )

p5d <- ggplot(
  counterfactual_summary,
  aes(median_main, country, colour = scenario_label)
) +
  geom_vline(xintercept = 0, colour = "#777777", linewidth = 0.35) +
  geom_errorbarh(
    aes(xmin = lower_95_main, xmax = upper_95_main),
    position = position_dodge(width = 0.42), height = 0, linewidth = 0.55
  ) +
  geom_point(position = position_dodge(width = 0.42), size = 2.0) +
  scale_colour_manual(
    values = c(
      `No new introduction` = scenario_colors[["No new introduction"]],
      `Equal lineage growth` = scenario_colors[["Equal lineage growth"]]
    ),
    name = "Conditional scenario",
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(-0.55, 1.02),
    breaks = c(-0.5, 0, 0.5, 1),
    labels = label_percent()
  ) +
  labs(
    x = "Cumulative post-2022 case reduction (95% CrI)",
    y = NULL
  ) +
  theme(legend.position = "top")

ct_curve <- read_source("figure5e_australia_ct_curve.tsv")

p5e <- ggplot(ct_curve, aes(ct, success_probability)) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper),
    fill = "#9CBED5", alpha = 0.38
  ) +
  geom_line(colour = "#3B78A8", linewidth = 0.8) +
  scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1), labels = label_percent()
  ) +
  labs(
    x = "PCR cycle threshold (Ct)",
    y = "Complete-sequence probability (95% CI)"
  )

recovery <- read_source("figure5f_identifiability_recovery.tsv") %>%
  mutate(
    parameter_type = recode(
      parameter_type,
      lineage_growth = "Lineage growth",
      import_scale = "Import scale"
    )
  ) %>%
  pivot_longer(
    c(coverage_95, median_absolute_log_error, correlation_truth_posterior_median),
    names_to = "metric", values_to = "value"
  ) %>%
  mutate(
    metric_label = recode(
      metric,
      coverage_95 = "95% coverage",
      median_absolute_log_error = "Median |log error|",
      correlation_truth_posterior_median = "Truth–estimate r"
    ),
    pass = case_when(
      metric == "coverage_95" ~ value >= 0.8,
      metric == "median_absolute_log_error" & parameter_type == "Lineage growth" ~ value <= 0.2,
      metric == "median_absolute_log_error" & parameter_type == "Import scale" ~ value <= 0.5,
      metric == "correlation_truth_posterior_median" ~ value >= 0.7,
      TRUE ~ FALSE
    ),
    display = case_when(
      metric == "coverage_95" ~ percent(value, accuracy = 0.1),
      TRUE ~ number(value, accuracy = 0.001)
    ),
    metric_label = factor(
      metric_label,
      levels = c("95% coverage", "Median |log error|", "Truth–estimate r")
    )
  )

p5f <- ggplot(
  recovery,
  aes(parameter_type, metric_label, fill = pass)
) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = paste0(display, "\n", ifelse(pass, "PASS", "FAIL"))),
    size = 2.1, lineheight = 0.88, fontface = "bold",
    colour = ifelse(recovery$pass, "#173D33", "#692C28")
  ) +
  scale_fill_manual(
    values = c(`TRUE` = "#BFDCCE", `FALSE` = "#E8C1BC"),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold"),
    panel.border = element_blank()
  )

fig5_design <- "
ABC
DEF
"
p5a <- tag_panel(p5a, "A")
p5b <- tag_panel(p5b, "B")
p5c <- tag_panel(p5c, "C")
p5d <- tag_panel(p5d, "D")
p5e <- tag_panel(p5e, "E")
p5f <- tag_panel(p5f, "F")
figure5 <- p5a + p5b + p5c + p5d + p5e + p5f +
  plot_layout(
    design = fig5_design,
    heights = c(1.25, 1),
    guides = "collect"
  ) +
  publication_annotation() &
  theme(legend.position = "bottom")

save_figure(figure5, "Figure_5_counterfactuals_and_calibration", 183, 155)

# ---------------------------------------------------------------------------
# Machine-readable render manifest
# ---------------------------------------------------------------------------

render_manifest <- tibble(
  figure = paste0("Figure ", 1:5),
  stem = c(
    "Figure_1_data_and_model",
    "Figure_2_phylogeny_and_ancestry",
    "Figure_3_reseeding_and_clusters",
    "Figure_4_lineage_growth_and_sampling",
    "Figure_5_counterfactuals_and_calibration"
  ),
  width_mm = rep(183, 5),
  height_mm = c(150, 160, 135, 165, 155),
  backend = "R",
  font_family = font_family,
  formats = "PDF-vector;PNG-600dpi",
  rendered_at = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE)
)
write_tsv(render_manifest, file.path(output_dir, "RENDER_MANIFEST.tsv"))

message("Rendered five R multi-panel main figures to: ", output_dir)
