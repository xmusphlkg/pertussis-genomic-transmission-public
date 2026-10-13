# Pertussis public archive accessibility: data and reproducibility

This repository releases data, computations, analytical figures and numerical results. Submission documents are intentionally excluded.

## Contents

| Directory | Purpose |
| --- | --- |
| `data/` | Accession-level collection intervals, archive release routes and sourced monthly cases |
| `data/archive_sources/` | Cached ENA run reports and project query audit |
| `data/phylogeny/` | Fixed trees, tip metadata and lineage labels for Figure 1A |
| `scripts/` | Numerical reconstruction and standalone figure rendering |
| `results/reference/` | Frozen numerical figure inputs |
| `results/recomputed/` | Recomputed milestones, lag bounds, count displacement, project release patterns and matched comparisons |
| `figures/` | Analytical figures and panel source data |
| `docs/` | Computation process and metadata inventory |
| `provenance/` | Checksums and software information |

## Reproduce

Run from the repository root with Python 3 and pandas:

```bash
python3 scripts/reproduce_results.py
```

Install R packages data.table, ggplot2, patchwork, svglite, ragg, systemfonts, ape, scales and lubridate, then run:

```bash
Rscript scripts/render_figure1.R
```

Numerical checks cover 205 target records (China 148, Japan 54, Australia 3), case milestones, 68/82 (82.9%) modal-date releases, 11/13 single-date projects, 8/9 larger projects, and 10 matched strata comprising 81 target and 44 comparator records. The fixed tree has 989 tips; archive metadata cover 774 focal records. This independently assembled cohort is not a random subset of Zhang et al.'s 8,117 genomes.

The package reproduces archive timing and displays a fixed input tree; it does not rerun genome acquisition or phylogenetic inference. Raw genomes remain at the original archives. See [process documentation](docs/PROCESS.md) for date provenance and interpretation limits. No new repository-wide licence is imposed on third-party records.
