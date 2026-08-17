# EID public-archive timing data dictionary

Dates use ISO 8601 (`YYYY-MM-DD`). Blank fields denote unavailable or
inapplicable values; they are not zeroes. Country identifiers use ISO3 codes
(`AUS`, `CHN`, `JPN`) with a separate human-readable label where needed.

## Accession-level availability table

`data/derived/public_genome_availability.tsv` contains one row per frozen focal
genome.

| Field group | Fields | Definition |
|---|---|---|
| Identity | `tree_sample_id`, `genome_record_id`, accession fields, `project_id` | Frozen-tree identifier and public archive identifiers |
| Collection interval | `collection_lower`, `collection_upper`, `collection_upper_effective`, `date_resolution` | Earliest and latest dates compatible with reported collection metadata; the effective upper bound is the earlier of the stated upper bound and public date |
| Temporal audit | `temporal_consistency_status` | `ok` when interval ordering is valid; `public_before_collection_lower` flags an unresolved source conflict |
| Public availability | `ena_first_public_date`, `assembly_release_date`, `public_date`, `public_route` | Route-specific dates and their earliest reproducible value |
| Lineage context | `primary_model_lineage_id`, `display_lineage`, `epidemic_period` | Membership from the frozen, retrospective tree and the prespecified period label |
| Geography | `subnational_location`, `location_source`, `location_resolution` | Standardized location and source-resolution audit; blank subnational location means unavailable or country-only metadata |
| Lag interval | `lag_min_days`, `lag_max_days` | Public date minus effective collection upper bound and collection lower bound, respectively |
| Case comparison | threshold/peak date and lead fields | Temporal comparison with national case milestones; positive lead values indicate earlier genome timing |

`public_route` is `ena_run_or_biosample`, `assembly`, `ena_and_assembly_same_date`,
or `not_publicly_dated`. It describes the earliest observed route, not the only
route through which a record could later be found.

## Primary and sensitivity results

| File | Unit and purpose |
|---|---|
| `eid_detection_clock_shift.tsv` | One selected threshold per country; collection interval, public date, and displacement interval |
| `eid_country_lineage_lag_summary.tsv` | Country-lineage summary of accession-level lag intervals |
| `eid_milestone_visibility.tsv` | Definite collection, possible collection, and public availability at case milestones |
| `eid_threshold_sensitivity.tsv` | Results for `k = 1, 3, 5, 10, 20, 50` where sample size permits |
| `eid_case_clock_sensitivity.tsv` | Monthly and native-resolution alternative case milestones |
| `eid_project_batch_release.tsv` | BioProject date completeness, modal release date, batch span, public route, and lag interval |
| `eid_geography_audit.tsv` | Subnational metadata completeness by country, project, and lineage |
| `eid_external_candidate_summary.tsv` | Candidate-project boundary relative to the frozen tree |

For milestone tables, “definitely collected” means
`collection_upper_effective <= milestone`; “possibly collected” means
`collection_lower <= milestone`. Public availability requires
`public_date <= milestone`.

## Released metadata snapshots

The five `data/source_snapshots/eid_nas_*.tsv` files contain only the fields
needed for the EID analysis: public accessions, collection metadata, geography,
sequencing technology, assembly release date, national case observations,
source registry entries, and project-level reconciliation fields. They contain
no reads or assembled sequences. Their record counts and SHA-256 values are in
`provenance/EID_NAS_SNAPSHOT_MANIFEST.tsv`; the public manifest intentionally
omits private storage locations and source-system modification times.

## Candidate-project boundary

`candidate_project_run_metadata.tsv` is a metadata-only ENA audit and does not
assign candidate runs to a lineage. `recommended_tier` describes whether a
project could be considered for a future post-freeze extension; it does not
make a claim about membership in `L1_02.07`.
