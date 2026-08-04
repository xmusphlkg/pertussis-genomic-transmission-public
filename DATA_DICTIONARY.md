# Data Dictionary

TSV files have header rows. JSON files store model inputs, diagnostics, or
validation records. Reader-facing text uses **MT28-associated genomic lineage**;
the frozen machine identifier is `L1_02.07`.

## Main Locations

| Path | Contents |
|---|---|
| `figures/letter/` | Two main Letter figures, three supplementary figures, render manifests, and figure-source maps |
| `figures/source_data/` | Compact panel-level source tables |
| `data/model_inputs/` | Frozen primary Stan input object and its component tables |
| `data/model_sensitivity/` | Frozen inputs for threshold, timing, root, cohort, and initial-state perturbations |
| `data/derived/` | Harmonised public accessions, case series, genome strata, source registry, and lineage display map |
| `results/lineage_nomenclature/` | MT28 label-to-clade crosswalk and validation |
| `results/model_main/` | Primary posterior summaries, pairwise contrasts, scenarios, diagnostics, and recovery |
| `results/model_no_project/` | No-project-effects sensitivity summaries |
| `results/model_growth_robustness/` | Country-only, omission, dominant-project, and cohort-cap analyses |
| `results/model_input_sensitivity_summary/` | Six full input-refit summaries |
| `results/model_anchor_sensitivity/` | Conditional anchored-scenario sensitivity |
| `results/cohort_gate_sensitivity/` | Pre-rescue versus final-cohort audit |
| `results/molecular_characterisation/` | Marker annotation coverage and interpretation limits |
| `results/phylogeny/` | 989-tip tree, core-SNP alignment, metadata, uniform QC, and tree QA |
| `results/cgmlst/`, `results/lineages/`, `results/phylogeography*/` | Extended genomic and sampled-ancestry audit outputs |
| `scripts/figures/` | Letter figure renderers |
| `scripts/model/` | Stan model source |
| `scripts/pipeline/` | Upstream construction and sensitivity scripts |
| `scripts/qa/` | Public-package validation and metadata sanitisation |
| `tests/` | Regression tests for numerical and structural claims |
| `provenance/` | File SHA-256 manifest, excluded-object hashes, R session info, and citation audit |

## Panel Source Files

The 13 released panel tables in `figures/source_data/` are the intended compact
source-data layer. Their figure mapping is recorded in
`figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv` and
`figures/letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv`.

## Naming

Comparator lineages A-C correspond to `L1_01.02`, `L1_02.05`, and `L1_02.06`.
Machine-readable outputs keep internal lineage IDs; display labels are in
`data/derived/model_lineage_display_map.tsv`.

## Redistribution Boundary

Private compute paths are omitted. The Australian direct-specimen source table
is not redistributed. Public BioSample/run accessions used in the tree remain,
but exact Australian Ct values, source-specific aliases, specimen types,
profile-success fields, and article-supplied marker annotations are removed
from per-record released tables. Only the aggregate Ct-dependent recovery curve
is released.
