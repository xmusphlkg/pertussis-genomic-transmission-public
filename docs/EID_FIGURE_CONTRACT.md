# EID Figure 1 contract

- **Core conclusion:** cumulative target-lineage specimen accumulation and public sequence availability separated around resurgence milestones in China and Japan, Australia showed the contrasting order, and milestone/peak counts expose the public-archive visibility gap directly.
- **Archetype:** quantitative grid with panel A as the hero evidence.
- **Target/output:** Emerging Infectious Diseases Dispatch; 178 × 125 mm; editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG preview.
- **Backend:** R only (`ggplot2`, `patchwork`, `svglite`, `ragg`).
- **Panel A:** country-specific monthly case curves aligned above cumulative possibly collected, definitely collected, and publicly retrievable target-lineage records, with resurgence and peak markers shared across both layers.
- **Panel B:** specimen accumulation and public sequence dates relative to the resurgence milestone month, using milestone-month intervals for day-offset uncertainty.
- **Panel C:** genomic visibility at epidemiologic milestones, showing definitely-to-possibly collected record intervals and publicly available counts at the resurgence milestone and post-2022 peak.
- **Hero evidence:** separation of accumulation and public-availability curves in China and Japan, with Australia as a contrasting order.
- **Validation evidence:** 10 project- and year-matched lineage strata with group-specific medians and sample sizes, reported outside the main Figure 1 grid unless the editor requests an Appendix figure.
- **Robustness:** collection-date uncertainty remains visible, sequence-count anchors remain interpretive, and sensitivity analyses are reported in the Appendix.
- **Statistics:** descriptive interval medians only; genome-level tests are omitted because records within project release batches are dependent.
- **Source data:** the current `eid_figure1*.tsv` files under `figures/source_data/`, including monthly case series for Panel A, relative timing intervals for Panel B, milestone visibility counts for Panel C, and the Appendix-only matched-lineage comparison source.
- **Image integrity:** no raster source images or selective image adjustment; all panels are generated from tabular data.
- **Reviewer risk:** distinguish external sequence retrievability from local sequencing, lineage assignment, and reporting timestamps.
- **Main-table coordination:** complete; the figure legend and main table align after table columns were split into collection interval, first archive-recorded public date, collection-to-public interval, definitely collected at milestone, possibly collected at milestone, and public at milestone.
