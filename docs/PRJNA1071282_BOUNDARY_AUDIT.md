# PRJNA1071282 boundary audit for the EID Dispatch

Audit updated: 2026-08-17

## Scope and sources

PRJNA1071282 was audited using the frozen-tree metadata, repository-local ENA cache, NAS recovery/QC metadata snapshots, and an accession-level comparison against the complete project. The reproducible snapshot is `data/source_snapshots/eid_nas_prjna1071282_run_audit.tsv`; summary values are generated in `results/public_availability/eid_external_candidate_summary.tsv`.

The complete project contains 734 explicit *Bordetella pertussis* read runs from Shanghai, representing 723 unique sample accessions. ENA collection dates span 2016–2024, but all 734 records are year-level. ENA first-public dates span 2024-09-04 to 2024-09-29.

## Frozen-analysis overlap

The project is not wholly absent from the analysis:

- 16 project records are already represented as tips in the frozen 989-genome tree;
- 6 of those tips were retrospectively assigned to L1_02.07;
- 3 of those target-lineage tips belong to the post-2022 resurgence-period analysis;
- the accession-level public-availability table retains those frozen records according to the main inclusion rules.

Therefore, the manuscript must say: “The frozen cohort included 16 PRJNA1071282 samples; the complete 734-run expansion was not included.” It must not say that PRJNA1071282 was not included.

## Why the complete expansion remains outside the main clock

The EID main analysis compares month-scale collection intervals with monthly national case milestones. Because all 734 project records have year-only collection dates, adding the complete expansion would not provide a month-scale collection clock. Moreover, records outside the frozen tree have not been de-duplicated, quality-controlled, and lineage-placed in the frozen framework. They therefore cannot be counted as L1_02.07 members solely from project membership.

The expansion remains a boundary dataset, not a replacement for the frozen primary cohort. Its public-date range can characterize archive timing, but it cannot support a month-resolved collection threshold or a target-lineage count without additional source dates and phylogenetic placement.

## Re-evaluation requirements

1. obtain month-level or day-level collection dates from a citable source;
2. de-duplicate against the frozen 989 genomes and other recovered records;
3. apply sequence QC and place eligible genomes in the frozen lineage framework;
4. rerun the interval clock as a labelled post-freeze sensitivity analysis, without silently changing the primary cohort.
