#!/usr/bin/env Rscript
# Figure 1: phylogenetic context and archive observation clocks for China and Japan.
source("scripts/figure_style.R")
suppressPackageStartupMessages({ library(ape); library(scales); library(lubridate) })
root <- normalizePath(".")
outdir <- file.path(root, "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

country_levels <- c("CHN", "JPN")
country_labels <- c(CHN = "China", JPN = "Japan")
country_colors <- c(CHN = "#C75245", JPN = "#2F8F83")

# Panel A: frozen tree with country tips, lineage/year/target strips and node support;
# the 202 resurgence-period target-lineage records are outlined.
tree <- ladderize(read.tree("data/phylogeny/midpoint.tree"))
support_tree <- read.tree("data/phylogeny/support.treefile")
meta <- fread("data/phylogeny/tip_metadata.tsv")
ancestry <- fread("data/phylogeny/lineage_membership.tsv")
stopifnot(length(tree$tip.label) == 989, !anyDuplicated(meta$tree_sample_id), nrow(ancestry) == 774,
          setequal(tree$tip.label, support_tree$tip.label))
ntip <- length(tree$tip.label)
xy <- data.table(node = seq_len(ntip + tree$Nnode), x = node.depth.edgelength(tree), y = NA_real_)
xy[seq_len(ntip), y := seq_len(ntip)]
children <- split(tree$edge[, 2], tree$edge[, 1])
node_y <- function(node) {
  if (is.finite(xy$y[node])) return(xy$y[node])
  ys <- vapply(children[[as.character(node)]], node_y, numeric(1))
  xy$y[node] <<- mean(range(ys)); xy$y[node]
}
invisible(node_y(ntip + 1L))
# Truncate the long midpoint-root branch for display: genetic distances beyond x_cut are shifted left by x_shift.
x_cut <- .0034; x_shift <- .0032
stopifnot(min(xy$x[xy$node != ntip + 1L]) > x_cut)
xy[, x_true := x]; xy[node != ntip + 1L, x := x - x_shift]
# Transfer SH-aLRT/UFBoot support from the unrooted ML tree by matching bipartitions.
split_key <- function(labs, ref) { if (ref %in% labs) labs <- setdiff(tree$tip.label, labs); paste(sort(labs), collapse = "|") }
ref_tip <- sort(tree$tip.label)[1]
sup_parts <- prop.part(support_tree); sup_labs <- attr(sup_parts, "labels")
sup_map <- setNames(support_tree$node.label, vapply(sup_parts, function(i) split_key(sup_labs[i], ref_tip), ""))
tree_parts <- prop.part(tree); tree_labs <- attr(tree_parts, "labels")
xy[, support := NA_character_]
xy[ntip + seq_len(tree$Nnode), support := unname(sup_map[vapply(tree_parts, function(i) split_key(tree_labs[i], ref_tip), "")])]
xy[, high_support := vapply(strsplit(fifelse(is.na(support), "", support), "/", fixed = TRUE), function(z) {
  v <- suppressWarnings(as.numeric(z)); length(v) == 2 && all(is.finite(v)) && all(v >= 95) }, logical(1))]
edge <- data.table(parent = tree$edge[, 1], child = tree$edge[, 2])
edge[, `:=`(x = xy$x[parent], xend = xy$x[child], y = xy$y[child], yend = xy$y[child])]
vertical <- edge[, .(x = xy$x[parent[1]], y = min(y), yend = max(y)), by = parent]
tips <- data.table(tree_sample_id = tree$tip.label, x = xy$x[seq_len(ntip)], y = xy$y[seq_len(ntip)])
tips <- merge(tips, meta[, .(tree_sample_id, country_iso3, tree_role, date_lower, date_upper)], by = "tree_sample_id")
tips <- merge(tips, ancestry[, .(tree_sample_id, epidemic_period, primary_model_lineage_id)], by = "tree_sample_id", all.x = TRUE)
tips[, highlight := country_iso3 %in% country_levels & epidemic_period == "resurgence" & primary_model_lineage_id == "L1_02.07"]
tips[, country_group := fifelse(tree_role == "focal", country_iso3, "Background")]
tips[, lineage_group := fifelse(tree_role == "global_background", "Background", fifelse(is.na(primary_model_lineage_id), "Other", primary_model_lineage_id))]
tips[, sampling_year := decimal_date(as.Date((as.numeric(as.Date(date_lower)) + as.numeric(as.Date(date_upper))) / 2))]
stopifnot(sum(tips$highlight) == 202, nrow(tips) == 989, sum(xy$high_support) > 0)
xmax <- max(xy$x)
all_country_colors <- c(AUS = "#3B78A8", BEL = "#C69232", CHN = "#C75245", FRA = "#7C65A8", JPN = "#2F8F83", Background = "#D4D4D4")
# Lineage strip colours deliberately avoid the country tip hues so the two encodings cannot be confused.
lineage_colors <- c(L1_01.02 = "#E78AC3", L1_02.05 = "#A6D854", L1_02.06 = "#8C6D31", L1_02.07 = "#2B2D42", Other = "#9E9E9E", Background = "#E3E3E3")
# Sampling-year bins are finer after 2010, where most genomes were collected.
year_breaks <- c(-Inf, 2010, 2018, 2022, 2023, 2024, Inf)
year_bin_labels <- c("<2010", "2010\u201317", "2018\u201321", "2022", "2023", "2024+")
tips[, year_bin := cut(sampling_year, year_breaks, labels = year_bin_labels, right = FALSE)]
year_colors <- setNames(c("#3F4A8A", "#4B9FB5", "#A6D96A", "#FDD049", "#F46D43", "#A50026"), year_bin_labels)
strip_w <- xmax * .028
strip_x <- xmax * c(1.04, 1.075, 1.11)
strip_labels <- data.table(x = strip_x, label = c("Lineage", "Year", "Target"))
p_a <- ggplot() +
  geom_segment(data = edge, aes(x, y, xend = xend, yend = yend), colour = "#A7A7A7", linewidth = .13) +
  geom_segment(data = vertical, aes(x, y, xend = x, yend = yend), colour = "#A7A7A7", linewidth = .13) +
  geom_point(data = xy[high_support == TRUE], aes(x, y, shape = "SH-aLRT/UFBoot ≥95%"), size = .3, colour = "#333333", alpha = .75) +
  geom_point(data = tips[highlight == FALSE], aes(x, y, fill = country_group), shape = 21, size = .55, colour = "#555555", stroke = .15, alpha = .9) +
  geom_point(data = tips[highlight == TRUE], aes(x, y, fill = country_group), shape = 21, size = .95, colour = "black", stroke = .3) +
  scale_fill_manual(values = all_country_colors, breaks = names(all_country_colors),
                    labels = c(AUS = "Australia", BEL = "Belgium", CHN = "China", FRA = "France", JPN = "Japan", Background = "Background"),
                    name = "Country", guide = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21, size = 1.6, colour = "#555555", stroke = .2))) +
  scale_shape_manual(values = c(`SH-aLRT/UFBoot ≥95%` = 16), name = "Node support",
                     guide = guide_legend(order = 4, override.aes = list(size = 1.1))) +
  annotate("text", x = .00012, y = range(edge[parent == ntip + 1L]$y), label = "//", size = 2.2, angle = 20, family = ijid_font) +
  annotate("segment", x = .0006, xend = .0011, y = 14, yend = 14, linewidth = .35) +
  annotate("text", x = .00085, y = 34, label = "0.0005", size = 2, family = ijid_font) +
  scale_x_continuous(breaks = c(.004, .005, .006) - x_shift, labels = c("0.004", "0.005", "0.006"), expand = expansion(mult = c(.01, .02))) +
  scale_y_continuous(limits = c(.5, ntip + .5), expand = c(0, 0)) +
  labs(x = "Substitutions per SNP site\n(root branch truncated)", y = NULL) + ijid_theme +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        legend.position = "bottom", plot.margin = margin(5, 1, 5, 5))
strip_theme <- theme_void(base_family = ijid_font, base_size = 6.3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = .5, size = 5.4, colour = "#333333"),
        plot.margin = margin(5, .6, 5, .6), legend.position = "bottom")
strip_plot <- function(fill_scale, label, mapping) {
  ggplot(tips, mapping) + geom_tile(width = .85, height = 1) + fill_scale +
    scale_x_continuous(position = "top", breaks = 1, labels = label, expand = expansion(add = .1)) +
    scale_y_continuous(limits = c(.5, ntip + .5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + strip_theme
}
p_a_lineage <- strip_plot(scale_fill_manual(values = lineage_colors, breaks = names(lineage_colors)[1:5], name = "Lineage",
                                            guide = guide_legend(order = 2, nrow = 1)), "Lineage", aes(1, y, fill = lineage_group))
p_a_year <- strip_plot(scale_fill_manual(values = year_colors, breaks = year_bin_labels, drop = FALSE, na.value = "#B3B3B3", name = "Sampling year",
                                         guide = guide_legend(order = 3, nrow = 1)), "Year", aes(1, y, fill = year_bin))
p_a_target <- strip_plot(scale_fill_manual(values = c(CHN = "#C75245", JPN = "#2F8F83", Other = "#F2F2F2"),
                                           breaks = c("CHN", "JPN"), labels = c(CHN = "China target (n=148)", JPN = "Japan target (n=54)"),
                                           name = "Target records", guide = guide_legend(order = 5, nrow = 1)),
                         "Target", aes(1, y, fill = fifelse(highlight, country_iso3, "Other")))
# Hand-built legend for panel A: one variable per row, keys on a shared column grid so items align across rows.
leg_rows <- list(
  list(title = "Country", type = "point", labels = c("Australia", "Belgium", "China", "France", "Japan", "Global"),
       cols = unname(all_country_colors)),
  list(title = "Lineage", type = "tile", labels = names(lineage_colors)[1:5], cols = unname(lineage_colors[1:5])),
  list(title = "Year", type = "tile", labels = year_bin_labels, cols = unname(year_colors)),
  list(title = "Target", type = "tile", labels = c("China", "Japan"), cols = c("#C75245", "#2F8F83")),
  list(title = "Support", type = "node", labels = "SH-aLRT and UFBoot \u226595%", cols = "#333333"))
leg_x0 <- .15; leg_dx <- .143
leg_keys <- rbindlist(lapply(seq_along(leg_rows), function(i) {
  r <- leg_rows[[i]]
  data.table(row = i, type = r$type, label = r$labels, col = r$cols, x = leg_x0 + (seq_along(r$labels) - 1) * leg_dx)
}))
leg_titles <- data.table(row = seq_along(leg_rows), title = vapply(leg_rows, `[[`, "", "title"))
a_legend <- ggplot() +
  geom_point(data = leg_keys[type == "point"], aes(x, -row, fill = I(col)), shape = 21, size = 1.5, colour = "#555555", stroke = .2) +
  geom_tile(data = leg_keys[type == "tile"], aes(x, -row, fill = I(col)), width = .026, height = .52) +
  geom_point(data = leg_keys[type == "node"], aes(x, -row), size = .8, colour = "#333333") +
  geom_text(data = leg_keys, aes(x + .022, -row, label = label), hjust = 0, size = 5 / .pt, family = ijid_font, colour = "#262626") +
  geom_text(data = leg_titles, aes(.015, -row, label = title), hjust = 0, size = 5.6 / .pt, fontface = "bold", family = ijid_font, colour = "#262626") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-length(leg_rows) - .6, -.4), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(2, 4, 0, 4))
p_a_panel <- wrap_elements(full = p_a + p_a_lineage + p_a_year + p_a_target + plot_layout(widths = c(1, .035, .035, .035)) & theme(legend.position = "none"))

# Shared archive inputs for panels B-D.
cum <- fread("results/reference/cumulative_visibility.tsv")[country_iso3 %in% country_levels]
cum[, event_date := as.Date(event_date)]
cum_long <- melt(cum, id.vars = c("country_iso3", "country_label", "event_date", "resurgence_milestone_date", "post2022_peak_date"),
                 measure.vars = c("n_possibly_collected", "n_definitely_collected", "n_publicly_available"),
                 variable.name = "series", value.name = "n")
cum_long[, series := factor(series, levels = c("n_possibly_collected", "n_definitely_collected", "n_publicly_available"),
                            labels = c("Collection count: upper bound", "Collection count: lower bound", "Publicly available"))]
case <- fread("data/country_month_cases.tsv")[country_iso3 %in% country_levels & case_data_available == TRUE]
case[, model_month := as.Date(model_month)]
case <- case[model_month >= as.Date("2023-01-01") & model_month <= as.Date("2025-12-01")]
case[, country_label := factor(country_labels[country_iso3], levels = country_labels[country_levels])]

# Panel B: reported cases above cumulative collection/public archive clocks.
# Milestone lines carry the facet variable so each country shows only its own threshold and peak.
milestones <- unique(cum[, .(country_iso3, resurgence_milestone_date, post2022_peak_date)])
milestones <- melt(milestones, id.vars = "country_iso3", variable.name = "milestone", value.name = "date")
milestones[, `:=`(date = as.Date(date), country_label = factor(country_labels[country_iso3], levels = country_labels[country_levels]),
                  milestone = factor(milestone, levels = c("resurgence_milestone_date", "post2022_peak_date"),
                                     labels = c("Resurgence threshold", "Post-2022 case peak")))]
stopifnot(nrow(milestones) == 4,
          milestones[country_iso3 == "CHN" & milestone == "Resurgence threshold", date] == as.Date("2023-08-01"),
          milestones[country_iso3 == "CHN" & milestone == "Post-2022 case peak", date] == as.Date("2024-05-01"),
          milestones[country_iso3 == "JPN" & milestone == "Resurgence threshold", date] == as.Date("2025-03-01"),
          milestones[country_iso3 == "JPN" & milestone == "Post-2022 case peak", date] == as.Date("2025-06-01"))
milestone_linetypes <- c("Resurgence threshold" = "dashed", "Post-2022 case peak" = "dotdash")
milestone_colours <- c("Resurgence threshold" = "#222222", "Post-2022 case peak" = "#8A8A8A")
p_b_cases <- ggplot(case, aes(model_month, cases)) +
  geom_area(fill = "#E5E5E5", colour = NA) + geom_line(linewidth = .38, colour = "#555555") +
  geom_vline(data = milestones, aes(xintercept = date, linetype = milestone, colour = milestone), linewidth = .35) +
  facet_wrap(~country_label, nrow = 1, scales = "free_y") +
  scale_linetype_manual(values = milestone_linetypes, name = NULL) + scale_colour_manual(values = milestone_colours, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + scale_y_continuous(labels = label_number(scale_cut = cut_short_scale()), expand = expansion(mult = c(0, .08))) +
  labs(x = NULL, y = "Reported cases") + ijid_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "top", legend.text = element_text(size = 5.4), legend.key.width = grid::unit(6, "mm"), legend.margin = margin(0, 0, 0, 0))
p_b_clock <- ggplot(cum_long, aes(event_date, n, colour = series, linetype = series)) +
  geom_ribbon(data = cum, aes(x = event_date, ymin = n_definitely_collected, ymax = n_possibly_collected), inherit.aes = FALSE, fill = "#A9CFE3", alpha = .55) +
  geom_step(data = cum_long[series == "Collection count: lower bound"], aes(colour = series, linetype = series), linewidth = .48) +
  geom_step(data = cum_long[series == "Collection count: upper bound"], aes(colour = series, linetype = series), linewidth = .48) +
  geom_step(data = cum_long[series == "Publicly available"], aes(colour = series, linetype = series), linewidth = .60) +
  geom_vline(data = milestones[milestone == "Resurgence threshold"], aes(xintercept = date), colour = milestone_colours[[1]], linetype = "dashed", linewidth = .35) +
  geom_vline(data = milestones[milestone == "Post-2022 case peak"], aes(xintercept = date), colour = milestone_colours[[2]], linetype = "dotdash", linewidth = .35) +
  facet_wrap(~country_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("Collection count: upper bound" = "#8FB9D4", "Collection count: lower bound" = "#2F6F9F", "Publicly available" = "#C75B39"), name = NULL) +
  scale_linetype_manual(values = c("Collection count: upper bound" = "22", "Collection count: lower bound" = "solid", "Publicly available" = "solid"), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") + scale_y_continuous(breaks = pretty_breaks(3), expand = expansion(mult = c(0, .08))) +
  labs(x = "Calendar date", y = "Cumulative records") + ijid_theme + theme(legend.position = "bottom", legend.text = element_text(size = 5.4), strip.text = element_blank())
p_b <- wrap_elements(full = p_b_cases / p_b_clock + plot_layout(heights = c(.72, 1)))

# Panel C: project-level concentration of public release (11/11 China/Japan projects shown).
batch <- fread("results/reference/target_project_batch_release.tsv")[country_iso3 %in% country_levels]
batch[, collection := as.Date(collection_lower_min)]; batch[, public := as.Date(modal_public_date)]
batch[, label := sprintf("%s · %s (target n=%d)", country_iso3, project_id, n_target_lineage)]
batch[, label := factor(label, levels = rev(label))]
batch[, all_one_date := modal_public_date_pct == 100]
stopifnot(sum(batch$n_target_lineage) == 202, batch[project_id == "PRJNA1182239", modal_public_date_pct] == 82.9)
p_c <- ggplot(batch, aes(y = label)) +
  geom_segment(aes(x = collection, xend = public, yend = label), colour = "#BFC4C7", linewidth = .7) +
  geom_point(aes(x = collection), colour = "#2F6F9F", size = 1.6) +
  geom_point(aes(x = public, shape = all_one_date), colour = "#C75B39", size = 2) +
  geom_text(data = batch[modal_public_date_pct < 100], aes(x = public + 35, label = sprintf("%.1f%%", modal_public_date_pct)), nudge_y = .18, hjust = 0, size = 2.0, family = ijid_font, colour = "#8F3E2A") +
  scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 0), labels = c(`TRUE` = "All included records released on one date", `FALSE` = "Multiple release dates"), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(.01, .16))) +
  labs(x = "Earliest included collection → modal public-release date", y = NULL) + ijid_theme +
  theme(axis.text.y = element_text(size = 5.8), legend.position = "bottom", legend.text = element_text(size = 5.2))

# Panel D: direct milestone counts, retaining collection-date precision as bounds.
ms <- fread("results/reference/milestone_visibility.tsv")[country_iso3 %in% country_levels]
ms[, milestone_label := fifelse(milestone == "first_post2022_month_above_2019_max", "Resurgence threshold", "Post-2022 peak")]
ms[, row := factor(paste(country_label, milestone_label, sep = " · "), levels = rev(c("China · Resurgence threshold", "China · Post-2022 peak", "Japan · Resurgence threshold", "Japan · Post-2022 peak")))]
ms[, collected_mid := (definitely_collected_by_milestone + possibly_collected_by_milestone) / 2]
ms[, collected_label := ifelse(definitely_collected_by_milestone == possibly_collected_by_milestone, as.character(definitely_collected_by_milestone), paste0(definitely_collected_by_milestone, "–", possibly_collected_by_milestone))]
ms[, public_label := paste0(public_by_milestone, " public")]
ms[, public_type := "Publicly available"]
p_d <- ggplot(ms, aes(y = row)) +
  geom_segment(aes(x = definitely_collected_by_milestone, xend = possibly_collected_by_milestone, yend = row), colour = "#2F6F9F", linewidth = 1.2) +
  geom_point(aes(x = definitely_collected_by_milestone, colour = "Collected"), size = 2) + geom_point(aes(x = possibly_collected_by_milestone, colour = "Collected"), shape = 1, size = 2.5) +
  geom_point(aes(x = public_by_milestone, colour = "Publicly available"), shape = 18, size = 2.1) +
  geom_text(aes(x = pmax(possibly_collected_by_milestone, public_by_milestone) + 5, label = paste0(collected_label, " collected | ", public_by_milestone, " public")), hjust = 0, size = 2.4, family = ijid_font) +
  scale_colour_manual(values = c("Collected" = "#2F6F9F", "Publicly available" = "#C75B39"), name = NULL) +
  scale_x_continuous(limits = c(0, max(ms$possibly_collected_by_milestone) + 75), breaks = pretty_breaks(4)) +
  labs(x = "Number of sampled target-lineage genomes", y = NULL) + ijid_theme +
  theme(axis.text.y = element_text(size = 5.8), legend.position = "bottom", legend.text = element_text(size = 5.2))

# Flat design (nested free() fails in this patchwork version): the A legend sits under the tree,
# freed from alignment with D, and C:D heights are about 2:1.
fig_design <- "
AB
AC
AD
ED"
fig <- p_a_panel + p_b + p_c + p_d + free(a_legend, type = "panel") +
  plot_layout(design = fig_design, widths = c(1.05, 1.25), heights = c(1.6, .74, .07, .30)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", ""))) & ijid_tag_theme
ijid_save(fig, file.path(outdir, "Figure_1_observation_structure"), width_mm = 205, height_mm = 190)

source_out <- rbindlist(list(
  tips[, .(panel = "A", tree_sample_id, country_iso3, x = x + x_shift, y, highlight, primary_model_lineage_id, lineage_group, sampling_year)],
  cum_long[, .(panel = "B", country_iso3, date = as.character(event_date), series = as.character(series), value = n)],
  batch[, .(panel = "C", country_iso3, project_id, collection_date = as.character(collection), public_date = as.character(public), value = modal_public_date_pct)],
  ms[, .(panel = "D", country_iso3, milestone = milestone_label, definite = definitely_collected_by_milestone, possible = possibly_collected_by_milestone, public = public_by_milestone)]
), fill = TRUE)
fwrite(source_out, file.path(outdir, "Figure_1_observation_structure_source.tsv"), sep = "\t")
writeLines(c(
  "Figure 1. Phylogenetic context and archive observation clocks for MT28-associated Bordetella pertussis genomes in China and Japan.",
  "(A) The fixed midpoint-rooted 989-genome core-SNP tree. Tip colours denote country of sampling; the Target annotation identifies the 202 included resurgence-period L1_02.07 records from China (n=148) and Japan (n=54), while the remaining genomes provide phylogenetic context. Strips to the right show the genomic lineage, sampling year (midpoint of the reported date interval, binned more finely after 2010, where most genomes were sampled) and target-record membership of each tip (red, China; green, Japan). Dark internal-node points mark clades with SH-aLRT and UFBoot support both ≥95%. Branch lengths are substitutions per SNP site; the long branch from the midpoint root is truncated (//) for display. Australia, Belgium and France are shown for phylogenetic context only.",
  "(B) Monthly reported cases (upper row) and cumulative target-lineage records (lower row). A light-blue band spans the upper and lower collection-count bounds arising from month/year date precision, with blue boundary lines; orange steps show records with a reproducible public sequence date. The band is a date-precision interval, not a confidence interval. Within each country, black dashed vertical lines mark the resurgence threshold (first post-2022 month above the 2019 maximum; China, August 2023; Japan, March 2025) and grey dot-dash lines mark the post-2022 case peak (China, May 2024; Japan, June 2025).",
  "(C) Earliest included target-lineage specimen collection and modal public-release dates for 11 BioProjects in China and Japan. Labels give the number of included target-lineage records (target n); percentages are shown only when fewer than 100% share the modal release date. Filled orange squares indicate that all included target-lineage records share one modal release date; open squares indicate multiple release dates. Connecting lines represent the interval between the earliest included collection and modal public-release date, not a project-level median delay. The project-level summary in the text additionally includes two Australian projects.",
  "(D) Target-lineage records at the resurgence threshold and subsequent case peak. Filled and open blue points show lower and upper collection bounds; orange diamonds show publicly available records. Labels use the format ‘15–29 collected | 0 public’. The x-axis is the number of sampled target-lineage genomes; bounds reflect collection-date precision rather than statistical confidence intervals.",
  "Source data: Figure_1_observation_structure_source.tsv; frozen tree; results/reference/*.tsv."
), file.path(outdir, "Figure_1_observation_structure_caption.txt"))
