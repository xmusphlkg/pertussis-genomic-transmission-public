# EID Dispatch claim–evidence map

## Central argument

This lineage-specific audit separates retrospective archive visibility from contemporaneous external visibility. Among eventually public, archive-derived records, MT28-associated specimens accumulated before national resurgence milestones in China and Japan, while corresponding sequences became publicly retrievable afterward. Australia supplied a descriptive 3-record contrast. Project- and year-matched comparisons connected the observed availability pattern to shared release batches and showed no descriptive evidence of a uniform target-lineage excess interval.

## Major claims

| Claim | Evidence | Primary files | Interpretive boundary |
|---|---|---|---|
| All 205 resurgence-period target-lineage records had a reproducible sequence-public date: China 148, Japan 54, and Australia 3. | Earliest retrievable ENA read or NCBI Assembly date for every included target record. | `data/derived/public_genome_availability.tsv`; `results/public_availability/eid_country_lineage_lag_summary.tsv` | Metadata-only BioSample dates are excluded from the primary endpoint. |
| In China, the 5-record accumulation interval was 2023 Jan 1–Mar 27 and sequence availability reached 5 records on 2024 Oct 4, producing a 557–642-day interval. Relative to the August 2023 milestone month, the public sequence date was 400–430 days later. | Fifth order statistics of collection lower/effective upper bounds and sequence-public dates. | `results/public_availability/eid_detection_clock_shift.tsv`; Figure 1A/1B source data | The count is a descriptive secondary anchor on the cumulative curve. |
| In Japan, the 5-record accumulation interval was 2024 Jan 1–Dec 9 and sequence availability reached 5 records on 2025 May 12, producing a 154–497-day interval. Relative to the March 2025 milestone month, the public sequence date was 42–72 days later. | Same interval estimand as China. | `results/public_availability/eid_detection_clock_shift.tsv`; Figure 1A/1B source data | The reported date precision determines interval width. |
| Australia showed the contrasting order: its 3-record accumulation interval followed the first resurgence milestone, with a 144–174-day accumulation-to-availability interval. Relative to the March 2024 milestone month, the public sequence date was 297–327 days later. | Three target records accumulated during 2024 Aug and were publicly retrievable by 2025 Jan 22. | `results/public_availability/eid_detection_clock_shift.tsv`; Figure 1A/1B source data | Australia is a descriptive 3-record contrast. |
| Cumulative accumulation preceded public visibility at the China and Japan resurgence milestones. | China milestone 15–29 collected/0 public and peak 89–106/0; Japan milestone 38–39/0 and peak 48–49/16; Australia peak 3/0. | `results/public_availability/eid_cumulative_visibility.tsv`; `results/public_availability/eid_milestone_visibility.tsv`; main Table; Figure 1C source data | Definite–possible ranges preserve collection-date precision. |
| Project matching showed no descriptive evidence of a uniform lineage-specific excess interval. | Ten country–BioProject–year strata included 81 target and 44 comparator records; target intervals overlapped in 5 strata, were shorter in 4, and longer in 1. | `results/public_availability/eid_project_lineage_comparison.tsv`; Appendix Table 6 | Descriptive stratum-level comparison accounts for release-batch dependence. |
| Release timing was shared within matched projects. | Target and comparator groups shared the same modal sequence-public date in all 10 matched strata. | `results/public_availability/eid_project_lineage_comparison.tsv` | Shared dates support a project-batched availability pattern. |
| Primary temporal ordering persisted across broad count anchors. | China accumulation preceded the resurgence milestone at counts 1–10; Japan retained the ordering at counts 1–20. | `results/public_availability/eid_threshold_sensitivity.tsv`; Appendix Table 3 | China count 20 spans the resurgence milestone. |
| Median record-level availability remained country-specific. | China 398–446 days; Japan 249.5–261 days; Australia 138–168 days. | `results/public_availability/eid_country_lineage_lag_summary.tsv` | These are intervals derived from reported collection-date resolution. |
| Public release was concentrated by project and route. | Project counts, modal release shares, release spans, and ENA/Assembly route summaries. | `results/public_availability/eid_project_batch_release.tsv`; Appendix Table 5 | The data establish batching patterns and leave causal attribution to stage-specific follow-up data. |
| The month-scale analysis includes 16 frozen PRJNA1071282 genomes, and the full project enters the candidate audit. | 734 explicit runs; 16 frozen tips, 6 target-lineage, 3 resurgence target-lineage; all 734 records have year-only collection dates. | `results/public_availability/eid_external_candidate_summary.tsv`; `data/source_snapshots/eid_nas_prjna1071282_run_audit.tsv`; `docs/PRJNA1071282_BOUNDARY_AUDIT.md` | Month-scale inclusion follows the frozen cohort; the extension remains available for lineage-placement work. |

## Concentrated interpretation boundary

The observed intervals cover collection through external sequence retrievability because sequencing, lineage-assignment, and local-reporting timestamps were unavailable. Retrospective lineage assignment and geographically concentrated sampling confine inference to the included records. Stage-specific timestamps can partition future estimates across specimen processing, sequencing, analysis, reporting, and public release.

## Required wording controls

- Use “specimen accumulation” and “public sequence availability” for the 2 cumulative processes.
- Define sequence-public date as the earliest reproducibly retrievable ENA read or NCBI Assembly date.
- Describe *k* as a descriptive secondary sequence-count anchor on the cumulative curve.
- Present Australia as the contrasting temporal order.
- Report the 10 matched strata, 81 target records, 44 comparator records, 5/4/1 direction pattern, and 10 shared modal dates together.
- Keep 557–642, 154–497, and 144–174 days as collection-to-public intervals; report milestone-relative public timing as 400–430, 42–72, and 297–327 days after the relevant milestone month.
- Place inferential boundaries in the single Conclusions boundary paragraph and the matching appendix interpretation block.
