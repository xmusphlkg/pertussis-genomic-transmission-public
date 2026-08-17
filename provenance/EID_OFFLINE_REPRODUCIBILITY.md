# EID offline reproducibility verification

Verification date: 2026-08-17

Mode: released metadata snapshots and caches; no private source system or network query used

Pipeline: offline `gtd_40` → `gtd_43` → R figure renderer → public-package QA

The ENA audit covers all 42 BioProjects represented by frozen focal genomes: 23 have cached read-run records and 19 were successfully queried and returned no read-run rows. This distinction is stored in `data/derived/public_availability_ena_project_audit.tsv`, so an empty result is not mistaken for an unqueried project. After this audit was populated, offline execution produced no missing-project warning.

Two consecutive offline executions completed successfully. SHA256 values were compared for 18 generated artifacts: the accession-level availability table, ENA project audit, 3 base clock tables, 8 EID result tables, the academic Appendix, and 4 EID figure source tables. The 2 checksum lists were identical. Word files are not generated. The final public-package QA passed.

SHA256 of the identical 18-file checksum list:

`0bd3d10189e1565e381e0887bfa78beca72952e924a58ca9e9c8b8bee7857606`

Key output checksums:

| File | SHA256 |
|---|---|
| `data/derived/public_availability_ena_project_audit.tsv` | `198317598a2e68e388b103352ff446914c13fe6809c63c420a5ca54382ae89f3` |
| `data/derived/public_genome_availability.tsv` | `6cf27563a70bd79715145b860bd6820fa0ff3a653325e343ae511bc29223d398` |
| `results/public_availability/eid_detection_clock_shift.tsv` | `79874a8340d99183ed270378ca1cd5c2533f93a649bf533705acdf6c33828e42` |
| `results/public_availability/eid_milestone_visibility.tsv` | `6bc195d8badb56ec666904aa0f4ee33a2c414d9c0d58cf148522829aa4c0ef1d` |
| `results/public_availability/eid_threshold_sensitivity.tsv` | `2fb28e14268e448be1839422559ab69ba38d4db7193b529138aa10c508f7d32a` |
| `results/public_availability/eid_project_batch_release.tsv` | `f6a7df2fc05994c7c640f195217b29cc635c9d732b14e43425cfa3a7e7f7c9cf` |
| `results/public_availability/eid_geography_audit.tsv` | `aa5de67bf7cf18cc1e47e42513295b7acbe33678cd4a579199ad8865c2562d10` |
| `results/public_availability/eid_case_clock_sensitivity.tsv` | `db8a9dc8e6202ec835318287d2d5d13252b89f33448eca695919211ba1a64772` |
| `manuscript/eid_dispatch_appendix.md` | `fa9868a68f232310dfc037caf5d70866f6af7f88315423c399d46e18627f5d47` |
