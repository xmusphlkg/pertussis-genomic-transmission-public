# Shared style for the archive accessibility figure.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(systemfonts)
})
ijid_fonts <- unique(systemfonts::system_fonts()$family)
ijid_font <- if ("Arial" %in% ijid_fonts) "Arial" else
  if ("Liberation Sans" %in% ijid_fonts) "Liberation Sans" else "sans"
ijid_palette <- c(ink = "#262626", muted = "#777777", blue = "#3B78A8",
                  teal = "#2F8F83", red = "#C75245", orange = "#C69232")
ijid_theme <- theme_classic(base_size = 7, base_family = ijid_font) +
  theme(axis.line = element_line(linewidth = .32, colour = ijid_palette[["ink"]]),
        axis.ticks = element_line(linewidth = .3), axis.text = element_text(size = 6.3),
        axis.title = element_text(size = 7), plot.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6.3), legend.title = element_text(size = 6.3),
        legend.key.size = grid::unit(3, "mm"), panel.grid = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 6.5, face = "bold"),
        plot.margin = margin(5, 5, 5, 5))
ijid_tag_theme <- theme(plot.tag = element_text(size = 8, face = "bold", family = ijid_font),
                        plot.tag.position = "topleft", plot.tag.location = "margin")
ijid_save <- function(plot, stem, width_mm = 183, height_mm = 155) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h, bg = "white")
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = ijid_font)
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm,
                 units = "mm", res = 600, compression = "lzw", background = "white")
  print(plot); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm,
                units = "mm", res = 300, background = "white")
  print(plot); dev.off()
}
ijid_validate_intervals <- function(dat) {
  stopifnot(all(c("median", "lower_95", "upper_95") %in% names(dat)),
            all(is.finite(dat$median)), all(is.finite(dat$lower_95)),
            all(is.finite(dat$upper_95)), all(dat$lower_95 <= dat$median),
            all(dat$median <= dat$upper_95))
}
