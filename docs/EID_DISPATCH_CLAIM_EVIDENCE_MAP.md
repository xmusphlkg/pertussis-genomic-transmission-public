# EID Dispatch claim–evidence map

## Central argument

In China and Japan, specimens later assigned to the frozen MT28-associated lineage had been collected before key national case milestones, whereas the public archive reached the corresponding *k*-genome availability thresholds later. Collection, public-archive, and operational-surveillance clocks are distinct; only the first 2 were observed.

## Major claims

| Claim | Evidence | Primary files | Boundary |
|---|---|---|---|
| The China 5-genome collection threshold was 2023 Jan 1–Mar 27 and public detection was 2024 Oct 4, a 557–642-day displacement. | Fifth order statistics of collection lower/effective upper bounds and public dates. | `results/public_availability/eid_detection_clock_shift.tsv`; Figure 1A/1C source data | Retrospective frozen-tree identity; collection is not local analysis. |
| The Japan 5-genome collection threshold was 2024 Jan 1–Dec 9 and public detection was 2025 May 12, a 154–497-day displacement. | Same interval estimand as China. | `results/public_availability/eid_detection_clock_shift.tsv`; Figure 1A/1C source data | National prevalence is not inferred. |
| The Australia 3-genome displacement was 144–174 days. | Only 3 target-lineage genomes were available; the collection interval was 2024 Aug 1–31 and public date 2025 Jan 22. | `results/public_availability/eid_detection_clock_shift.tsv` | Australia is descriptive and uses *k*=3. |
| Collected-but-not-public counts persisted at case milestones. | China: threshold 15–29/0 and peak 89–106/0; Japan: threshold 38–39/0 and peak 48–49/16; Australia peak 3–3/0. | `results/public_availability/eid_milestone_visibility.tsv`; main Table | Counts are definite–possible/public among included genomes, not national sampling denominators. |
| Median target-lineage delays were intervals, not points. | China 398–446 days; Japan 249.5–261; Australia 138–168. | `results/public_availability/eid_country_lineage_lag_summary.tsv`; Figure 1B source data | Bounds derive from reported date resolution. |
| Threshold conclusions are not uniformly robust through *k*=20 in China. | China collection is before the case threshold for *k*=1–10 but spans it at *k*=20; Japan remains before and public remains after for *k*=1–20. | `results/public_availability/eid_threshold_sensitivity.tsv`; Appendix Table 3 | Prevents overclaiming the requested sensitivity direction. |
| Public release was often project-batched and route-dependent. | Project counts, completeness, modal release shares, spans, lag intervals, and ENA/assembly route summaries. | `results/public_availability/eid_project_batch_release.tsv` | Descriptive; no cause of delay is assigned. |
| National case clocks were compared with geographically concentrated samples. | Subnational metadata completeness and composition by country/project/lineage. | `results/public_availability/eid_geography_audit.tsv` | No subnational or national lineage prevalence is estimated. |
| PRJNA1071282 is partly represented in the frozen analysis, but its complete expansion is outside the month-scale analysis. | 734 explicit runs; 16 frozen tips, 6 target-lineage, 3 resurgence target-lineage; all 734 expansion records have year-only collection dates. | `results/public_availability/eid_external_candidate_summary.tsv`; `data/source_snapshots/eid_nas_prjna1071282_run_audit.tsv`; `docs/PRJNA1071282_BOUNDARY_AUDIT.md` | Write “16 frozen samples included; full expansion not included,” never “project not included.” |

## Claims deliberately avoided

- “Lineage information existed” at collection time.
- Public availability proves local real-time identification or warning.
- The target lineage caused national resurgence or had a quantified transmission advantage here.
- Frozen-tree sampling estimates national or subnational lineage prevalence.
- All public-release delays share the same mechanism.
- Post-freeze candidates belong to L1_02.07 without frozen-framework placement.

## Required wording controls

- Use “specimens later assigned to the target lineage had been collected.”
- Use “the public archive reached the corresponding *k*-genome availability threshold.”
- Describe public dates as an optimistic earliest opportunity for external identification.
- Report 557–642, 154–497, and 144–174 days as intervals; do not report their upper bounds as unconditional point estimates.
- State that collection, archive, and operational-surveillance clocks are not interchangeable.
