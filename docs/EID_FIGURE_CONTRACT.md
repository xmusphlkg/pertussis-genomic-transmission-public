# EID Figure 1 contract

- **Core conclusion:** target-lineage collection was interval-censored and generally preceded public-archive visibility; collection and public clocks therefore cannot be interpreted as the same surveillance event.
- **Archetype:** quantitative grid with panel A as the hero evidence.
- **Target/output:** Emerging Infectious Diseases Dispatch; 178 × 120 mm; editable SVG/PDF plus 600-dpi TIFF and 300-dpi PNG preview.
- **Backend:** R only (`ggplot2`, `patchwork`, `svglite`, `ragg`).
- **Panel A:** national case curves, collection-detection interval, public-detection date, case threshold, and peak.
- **Panel B:** accession-level minimum-to-maximum collection-to-public lag and country median interval.
- **Panel C:** collection lead-time interval and public-date point relative to the primary case threshold.
- **Hero evidence:** separation between the blue collection interval and red public date in China and Japan.
- **Validation evidence:** accession-level lag intervals.
- **Robustness:** interval width is visible rather than replaced by a midpoint; threshold sensitivity is reported in the Appendix.
- **Statistics:** descriptive order statistics only; `n` is the number of frozen resurgence target-lineage genomes with reproducible public dates.
- **Source data:** the four `eid_figure1*.tsv` files under `figures/source_data/`.
- **Image integrity:** no raster source images or selective image adjustment; all panels are generated from tabular data.
- **Reviewer risk:** retrospective lineage membership must not be presented as evidence that the lineage was analysed or recognized at collection time.
