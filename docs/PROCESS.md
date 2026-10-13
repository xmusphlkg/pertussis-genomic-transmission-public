# Computation and source provenance

1. Select L1_02.07 resurgence records from Australia, China and Japan using the frozen tree assignments. Other resurgence lineages supply project/year-matched comparators.
2. Retain run, BioSample and assembly identifiers, ENA first-public and NCBI assembly release dates, selected public date and route/source labels. Recorded release events are distinct from BioSample-only appearance and historical file downloads.
3. A record is possibly collected when collection_lower is at or before the milestone, and definitely collected when collection_upper_effective is at or before it. These bounds describe date precision. Effective upper dates and temporal-consistency flags are inherited from upstream archive matching.
4. Derive the first post-2022 month above the 2019 maximum or median, the first of two consecutive months above the maximum, and the observed post-2022 peak. Evaluate all counts on the FIRST DAY of the milestone month. Case series retain source URLs and resolution.
5. Calculate paired accession collection-to-release intervals and medians. Separately sort collection bounds and release dates to calculate displacement at cumulative count k; kth collection and kth release need not be the same accession.
6. Calculate target-only modal release dates and percentages by BioProject, breaking ties toward the earliest date. Repeat for projects with at least five target records.
7. Match target and non-target records by country, project and reported collection year, then compare median interval bounds and modal dates. Record counts do not count independent release decisions.
8. Export numerical tables, verify the reported counts and render the fixed phylogenetic/timing figure. Record software versions and checksums.

## Retrieval dates and package freeze

The package was assembled on 13 October 2026. This is a package freeze date, NOT proof that all archive records were queried anew that day. Some ENA project query dates are blank and labelled records_cached_legacy. These gaps are preserved. A uniform archive retrieval date cannot be established from the available query audit. The archived metadata include release events into early 2026; historical download verification is unavailable.

The offline reconstruction checks frozen metadata and numerical results; it does not independently recrawl every accession. A new uniform access date requires dated responses from both archives and a comparison with the frozen inputs. Do not describe package assembly as a complete archive refresh.

## Interpretation

The cohort contains eventually public, geographically/project-selected records and provides no national sequencing denominator. Archive release dates cannot establish sequencing completion, local lineage recognition, public-health communication, deliberate delay or a causal effect on outbreak control.

`metadata_fields.tsv` records accession table fields, inferred types and missingness. The fixed-tree lineage inputs allow audit of cohort membership without redistributing patient-level data or raw genomes.
