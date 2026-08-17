# Appendix. Collection and Public Archive Timing of Pertussis Genomes

## Supplementary Methods

### Study design and frozen analysis boundary

This retrospective analysis compared specimen collection timing, public archive availability, and national pertussis case milestones. It retained the previously frozen 989-genome core-SNP phylogeny, including 774 focal genomes and 215 stratified global-background genomes, and used the existing L1_02.07 target-lineage assignments without re-estimating the tree, lineage definitions, or transmission models. Lineage membership was therefore retrospective. A specimen's collection date does not indicate when its lineage identity was known.

The timing analysis focused on frozen focal genomes from Australia, China, and Japan because these countries had compatible national case series and target-lineage genomes collected during the resurgence period. Projects identified after the phylogenetic freeze were evaluated separately. They were not added to the primary analysis without accession-level de-duplication, date-resolution review, sequence quality control, and placement within the frozen lineage framework.

### Collection intervals and public availability

Reported collection dates were represented by lower and upper bounds. Exact dates had identical bounds; month-level dates spanned the reported calendar month; year-level dates spanned the reported calendar year. The effective collection upper bound was the earlier of the reported upper bound and the first reproducible public date. A public date earlier than the collection lower bound was treated as a temporal conflict and was not silently corrected.

Public availability was defined as the earliest reproducible date from an ENA run or BioSample first-public record or an NCBI Assembly release record. ENA and NCBI dates were retained separately before selecting the earliest route. This endpoint estimates the earliest opportunity for an external archive user to retrieve the record. It does not measure local sequencing completion, bioinformatic analysis, lineage assignment, reporting, or public-health action.

### Detection thresholds and clock displacement

For a cumulative threshold of *k* target-lineage genomes, the collection-detection interval was defined by the *k*th order statistic of collection lower bounds and the *k*th order statistic of effective upper bounds. Public detection was the *k*th order statistic of public dates. Clock displacement ranged from the public date minus the collection-detection upper bound to the public date minus the collection-detection lower bound. The primary thresholds were *k*=5 for China and Japan and *k*=3 for Australia, where only 3 target-lineage genomes were available.

At each epidemiologic milestone, genomes were classified as definitely collected when the effective collection upper bound had passed and possibly collected when the collection lower bound had passed. Publicly available genomes were counted independently from their reproducible public dates. These counts describe the frozen genomic sample and are not estimates of national lineage prevalence.

### Epidemiologic milestones and sensitivity analyses

The primary case milestone was the first post-2022 month in which national reported cases exceeded the country-specific maximum monthly count observed in 2019. The post-2022 peak month was evaluated as a second milestone. Weekly Japanese reports were aggregated to calendar months for the primary comparison. Sensitivity analyses evaluated native reporting resolution, the first 2 consecutive periods above the 2019 maximum, the first month above the 2019 median, and genome thresholds of *k*=1, 3, 5, 10, 20, and 50 when sufficient genomes were available.

Project-level analyses summarized public-date completeness, the concentration of records on the modal release date, the span of release dates, and the first observable public route. Geographic metadata were standardized to subnational units when reported. These data were used to characterize sample composition and metadata completeness, not to calculate subnational or national lineage prevalence.

### Candidate-project and metadata audits

PRJNA1071282 contained 734 runs explicitly annotated as *B. pertussis*. Sixteen were already represented in the frozen tree, including 6 target-lineage genomes and 3 resurgence-period target-lineage genomes. All 734 runs had year-level collection dates; therefore, the full project extension was not eligible for the month-scale primary analysis. The appropriate boundary is that the 16 frozen genomes were included, whereas the complete project extension was not.

The analytic metadata extract retained accession identifiers, BioSample and BioProject identifiers, collection intervals, geography, sequencing technology, separate ENA and NCBI public dates, and record-matching status. It did not contain raw sequence reads, genome assemblies, or identifiable clinical information. Validation included accession uniqueness, cross-source matching, temporal-order checks, monotonicity of detection dates as *k* increased, and regeneration of derived results from the frozen metadata.

## Appendix Table 1. Primary interval-aware detection clock

| Country | Threshold, k | Collection interval | Public date | Clock displacement, d |
| --- | --- | --- | --- | --- |
| Australia | 3 | 2024-08-01 to 2024-08-31 | 2025-01-22 | 144–174 |
| China | 5 | 2023-01-01 to 2023-03-27 | 2024-10-04 | 557–642 |
| Japan | 5 | 2024-01-01 to 2024-12-09 | 2025-05-12 | 154–497 |

## Appendix Table 2. Target-lineage visibility at national case milestones

| Country | Milestone | Date | Definitely collected | Possibly collected | Publicly available |
| --- | --- | --- | --- | --- | --- |
| Australia | First month above 2019 maximum | 2024-03-01 | 0 | 0 | 0 |
| Australia | Post-2022 peak month | 2024-11-01 | 3 | 3 | 0 |
| China | First month above 2019 maximum | 2023-08-01 | 15 | 29 | 0 |
| China | Post-2022 peak month | 2024-05-01 | 89 | 106 | 0 |
| Japan | First month above 2019 maximum | 2025-03-01 | 38 | 39 | 0 |
| Japan | Post-2022 peak month | 2025-06-01 | 48 | 49 | 16 |

## Appendix Table 3. Sensitivity to the cumulative genome threshold

| Country | k | Collection interval | Public date | Displacement, d | Collection timing | Public timing |
| --- | --- | --- | --- | --- | --- | --- |
| Australia | 1 | 2024-08-01 to 2024-08-31 | 2025-01-16 | 138–168 | Collection interval after case threshold | Public on or after case threshold |
| Australia | 3 | 2024-08-01 to 2024-08-31 | 2025-01-22 | 144–174 | Collection interval after case threshold | Public on or after case threshold |
| China | 1 | 2023-01-01 to 2023-02-04 | 2024-09-29 | 603–637 | Collection interval before case threshold | Public on or after case threshold |
| China | 3 | 2023-01-01 to 2023-03-27 | 2024-09-29 | 552–637 | Collection interval before case threshold | Public on or after case threshold |
| China | 5 | 2023-01-01 to 2023-03-27 | 2024-10-04 | 557–642 | Collection interval before case threshold | Public on or after case threshold |
| China | 10 | 2023-01-01 to 2023-06-12 | 2024-10-04 | 480–642 | Collection interval before case threshold | Public on or after case threshold |
| China | 20 | 2023-04-06 to 2023-09-12 | 2024-10-04 | 388–547 | Collection interval spans case threshold | Public on or after case threshold |
| China | 50 | 2023-12-29 to 2023-12-31 | 2025-03-23 | 448–450 | Collection interval after case threshold | Public on or after case threshold |
| Japan | 1 | 2024-01-01 to 2024-09-12 | 2025-05-12 | 242–497 | Collection interval before case threshold | Public on or after case threshold |
| Japan | 3 | 2024-01-01 to 2024-11-13 | 2025-05-12 | 180–497 | Collection interval before case threshold | Public on or after case threshold |
| Japan | 5 | 2024-01-01 to 2024-12-09 | 2025-05-12 | 154–497 | Collection interval before case threshold | Public on or after case threshold |
| Japan | 10 | 2024-11-13 to 2024-12-31 | 2025-05-12 | 132–180 | Collection interval before case threshold | Public on or after case threshold |
| Japan | 20 | 2025-01-21 to 2025-01-22 | 2025-12-10 | 322–323 | Collection interval before case threshold | Public on or after case threshold |
| Japan | 50 | 2025-06-02 to 2025-06-09 | 2026-01-27 | 232–239 | Collection interval after case threshold | Public on or after case threshold |

## Appendix Table 4. Collection-to-public lag intervals for target-lineage genomes

| Country | Target-lineage genomes | Public date available | Median minimum lag, d | Median maximum lag, d | Observed lag range, d |
| --- | --- | --- | --- | --- | --- |
| Australia | 3 | 3 | 138 | 168 | 138–174 |
| China | 148 | 148 | 398.0 | 446.0 | 0–1047 |
| Japan | 54 | 54 | 249.5 | 261.0 | 25–755 |

## Appendix Table 5. Project-level public-release patterns

| Country | BioProject | Resurgence genomes | Target lineage | Public-date span | Batch span, d | Modal-date proportion, % | First observable route |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Australia | PRJNA1178746 | 2 | 2 | 2025-01-16 to 2025-01-22 | 6 | 50.0 | NCBI Assembly only: 2 |
| Australia | PRJNA1199062 | 35 | 1 | 2025-01-16 to 2025-01-16 | 0 | 100.0 | ENA read record only: 35 |
| China | PRJNA1071282 | 4 | 3 | 2024-09-29 to 2024-09-29 | 0 | 100.0 | ENA read record only: 4 |
| China | PRJNA1133929 | 18 | 18 | 2024-11-17 to 2024-11-17 | 0 | 100.0 | ENA read record first: 18 |
| China | PRJNA1143937 | 19 | 19 | 2024-10-04 to 2024-10-04 | 0 | 100.0 | NCBI Assembly only: 19 |
| China | PRJNA1182239 | 83 | 82 | 2025-06-17 to 2025-06-18 | 1 | 83.1 | NCBI Assembly only: 83 |
| China | PRJNA1193776 | 13 | 11 | 2025-03-23 to 2025-03-23 | 0 | 100.0 | ENA read record only: 13 |
| China | PRJNA1295129 | 16 | 15 | 2025-11-13 to 2025-11-13 | 0 | 100.0 | NCBI Assembly only: 16 |
| Japan | PRJDB20292 | 2 | 2 | 2025-05-20 to 2025-05-20 | 0 | 100.0 | NCBI Assembly first: 2 |
| Japan | PRJDB20413 | 15 | 14 | 2025-05-12 to 2025-05-12 | 0 | 100.0 | NCBI Assembly first: 15 |
| Japan | PRJDB34249 | 14 | 11 | 2026-01-25 to 2026-01-25 | 0 | 100.0 | ENA read record only: 14 |
| Japan | PRJDB35593 | 8 | 7 | 2025-12-10 to 2025-12-10 | 0 | 100.0 | NCBI Assembly first: 8 |
| Japan | PRJDB37898 | 21 | 20 | 2026-01-27 to 2026-01-27 | 0 | 100.0 | ENA read record only: 21 |

## Appendix Table 6. Completeness and composition of subnational geographic metadata

| Country | Genomes | With subnational metadata | Completeness, % | Reported locations |
| --- | --- | --- | --- | --- |
| Australia | 37 | 37 | 100.0 | Australian Capital Territory; New South Wales; Queensland; South Australia; Western Australia |
| China | 153 | 153 | 100.0 | Guangzhou; Shanghai; Xinjiang; Zhejiang |
| Japan | 60 | 60 | 100.0 | Okinawa; Osaka; Tottori |

## Appendix Table 7. Sensitivity to the national case-clock definition

| Country | Analysis scale | Source resolution | Milestone definition | Date | Reference case count |
| --- | --- | --- | --- | --- | --- |
| Australia | Harmonized monthly | monthly | First period above 2019 maximum | 2024-03-01 | 1315.0 |
| Australia | Harmonized monthly | monthly | First 2 consecutive periods above 2019 maximum | 2024-03-01 | 1315.0 |
| Australia | Harmonized monthly | monthly | First month above 2019 median | 2024-03-01 | 1002.5 |
| Australia | Harmonized monthly | monthly | Post-2022 peak | 2024-11-01 | 9308.0 |
| China | Harmonized monthly | monthly | First period above 2019 maximum | 2023-08-01 | 4388.0 |
| China | Harmonized monthly | monthly | First 2 consecutive periods above 2019 maximum | 2023-08-01 | 4388.0 |
| China | Harmonized monthly | monthly | First month above 2019 median | 2023-07-01 | 2712.0 |
| China | Harmonized monthly | monthly | Post-2022 peak | 2024-05-01 | 97669.0 |
| Japan | Harmonized monthly | monthly | First period above 2019 maximum | 2025-03-01 | 1389.0 |
| Japan | Harmonized monthly | monthly | First 2 consecutive periods above 2019 maximum | 2025-03-01 | 1389.0 |
| Japan | Harmonized monthly | monthly | First month above 2019 median | 2025-02-01 | 932.0 |
| Japan | Harmonized monthly | monthly | Post-2022 peak | 2025-06-01 | 16156.0 |
| Australia | Native reporting resolution | monthly | First period above 2019 maximum | 2024-03-01 | 1315 |
| Australia | Native reporting resolution | monthly | First 2 consecutive periods above 2019 maximum | 2024-03-01 | 1315 |
| Australia | Native reporting resolution | monthly | Post-2022 peak | 2024-11-01 | 9308 |
| China | Native reporting resolution | monthly | First period above 2019 maximum | 2023-08-01 | 4388 |
| China | Native reporting resolution | monthly | First 2 consecutive periods above 2019 maximum | 2023-08-01 | 4388 |
| China | Native reporting resolution | monthly | Post-2022 peak | 2024-05-01 | 97669 |
| Japan | Native reporting resolution | weekly | First period above 2019 maximum | 2025-02-17 | 353 |
| Japan | Native reporting resolution | weekly | First 2 consecutive periods above 2019 maximum | 2025-03-03 | 353 |
| Japan | Native reporting resolution | weekly | Post-2022 peak | 2025-07-14 | 3908 |

## Appendix Table 8. Projects identified after the phylogenetic freeze

| BioProject | Explicit B. pertussis runs | Country | Collection span | Date resolution | Frozen-tree genomes | Frozen target lineage | Frozen resurgence target | Interpretive status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PRJNA1071282 | 734 | China | 2016-01-01 to 2024-01-01 | Day: 0; month or interval: 0; year only: 734; missing: 0 | 16 | 6 | 3 | Year-only collection dates preclude month-scale inclusion |
| PRJNA1193776 | 62 | China | 2020-04-30 to 2024-05-03 | Day: 60; month or interval: 0; year only: 2; missing: 0 | 21 | 15 | 11 | Partly represented in the frozen cohort; additional records require lineage placement |
| PRJDB39872 | 16 | Japan | 2025-01-01 to 2025-04-01 | Day: 0; month or interval: 16; year only: 0; missing: 0 | 0 | 0 | 0 | Post-freeze project; quality control and lineage placement required |
| PRJNA1455114 | 48 | China | 2017-09-01 to 2024-08-01 | Day: 0; month or interval: 48; year only: 0; missing: 0 | 0 | 0 | 0 | Post-freeze project; quality control and lineage placement required |
| PRJEB88325 | 416 | Belgium | 2014-01-30 to 2023-12-31 | Day: 416; month or interval: 0; year only: 0; missing: 0 | 119 | 7 | 3 | Partly represented in the frozen cohort; additional records require lineage placement |
| PRJNA870170 | 5 | Australia | 2024-03-04 to 2025-03-04 | Day: 5; month or interval: 0; year only: 0; missing: 0 | 0 | 0 | 0 | Post-freeze project; taxonomic review and lineage placement required |

## Appendix Table 9. Metadata components and validation roles

| Data component | Records | Information retained | Analytic role and validation |
| --- | --- | --- | --- |
| China metadata recovery audit | 11 | Project-level metadata recovery and matching summaries | Used to distinguish recovered metadata from unresolved records |
| Three-country focal-genome metadata | 596 | Accessions, collection intervals, geography, sequencing technology, and assembly release dates | Matched to frozen-tree identifiers and checked for accession uniqueness |
| National surveillance observations | 838 | Monthly or weekly reported pertussis cases and source provenance | Checked for temporal completeness and harmonized to monthly scale for the primary analysis |
| PRJNA1071282 run metadata | 734 | Run and sample accessions, organism annotation, collection year, and public date | Restricted to explicit B. pertussis records and reconciled with frozen-tree membership |
| Public-health source registry | 17 | Issuing organization, source type, release information, and country | Used to verify the provenance and reporting context of surveillance series |

## Appendix Table 10. National pertussis surveillance series

| Country | Surveillance source | Native resolution | Observation period | Data freeze |
| --- | --- | --- | --- | --- |
| Australia | Australian National Notifiable Disease Surveillance System | monthly | 2015-01-01 to 2025-12-01 | 2026-04-10 |
| China | Official Chinese national and provincial public-health surveillance reports | monthly | 2015-01-01 to 2025-12-01 | 2026-04-10 |
| Japan | National Institute of Infectious Diseases and Japan Institute for Health Security weekly reports | weekly | 2014-12-29 to 2025-12-22 | 2026-04-10 |

## Definitions and interpretation

The reported collection interval preserves the precision of the source metadata. The effective upper bound enforces temporal coherence with public availability while retaining explicit conflict flags. Minimum and maximum release lags are calculated from the effective upper and lower collection bounds, respectively. Public route identifies whether an ENA read record, an NCBI Assembly record, or both supplied the earliest reproducible date. Subnational metadata describe the locations represented in the genomic sample and should not be interpreted as a population sampling frame.

No newly identified candidate record was assigned to the target lineage unless it was already represented by a frozen-tree identifier. Consequently, all lineage-specific timing estimates preserve the original phylogenetic analysis boundary.
