# Letter figure source data

This directory contains the compact panel-level tables used by the Letter
figures. Authoritative figure-to-file mappings are in:

- `../letter/LETTER_FIGURE_SOURCE_FILES.tsv`
- `../letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv`

Files required by a supplementary panel but stored under `results/` remain in
their analysis-native directory and are referenced from the supplementary
manifest.

The four `eid_figure1*.tsv` files are the panel-level inputs for the EID
public-archive timing figure. Their panel assignments and transformation rules
are documented in [`../../docs/EID_FIGURE_CONTRACT.md`](../../docs/EID_FIGURE_CONTRACT.md).
