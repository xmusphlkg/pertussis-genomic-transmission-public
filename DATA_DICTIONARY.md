# Data dictionary

This document describes the main public files at directory and table level. Tab-separated files use a header row; JSON files contain validation summaries or model inputs.

## Derived data

- `data/derived/source_registry.tsv`: registry of public source datasets, source URLs, roles, frozen file sizes, and checksums.
- `data/derived/transmission_genome_records.tsv`: harmonised focal-genome records with public accession identifiers, sampling-date resolution, project membership, and analysis eligibility.
- `data/derived/global_tree_background_records.tsv`: harmonised global-background genomes used only for tree topology and ancestry context.
- `data/derived/country_month_cases.tsv`: monthly aggregate pertussis cases, public source URLs, source resolution, and pre-specified reporting eras.
- `data/derived/country_month_case_genome_panel.tsv`: linked country-month case and genome-observation panel.
- `data/derived/country_month_genome_strata.tsv`: project- and lineage-stratified monthly genome counts.
- `data/derived/country_data_gate_summary.tsv`: pre-specified country eligibility gates for the case-genome model.
- `data/derived/preliminary_lineage_screen.tsv`: screening counts and criteria used before formal model-lineage freezing.

Local file-system path fields in released tables are blank because they describe the authors’ private compute environment rather than reusable data. For Australian records linked to public sequence accessions, source-specific sample aliases, exact Ct values, specimen types, profile-success fields, and article-supplied marker annotations are also blank because the source appendix is not redistributed.

## Model inputs

- `data/model_inputs/joint_model_data.json`: complete frozen input object for the Stan model.
- `data/model_inputs/joint_model_data_validation.json`: dimensions and consistency checks for the model input.
- `data/model_inputs/genome_observation_strata.tsv`: country-month-project lineage counts used in the genomic observation likelihood.
- `data/model_inputs/monthly_import_events.tsv`: lineage-weighted phylogeographic event inputs.
- `data/model_inputs/project_index.tsv`: stable project identifiers used by the model.
- `data/model_inputs/transmission_spline_basis.tsv`: monthly spline basis used for country-specific transmission functions.

## Results

- `results/phylogeny/`: final 989-tip tree, 11,550-position core-SNP alignment, tip metadata, uniform-QC results, and tree-QA reports.
- `results/lineages/`: primary and alternative rhierBAPS assignments, frozen model-lineage assignments, sensitivity comparisons, and validation summaries.
- `results/phylogeography/` and `results/phylogeography_alternative_root/`: marginal geographic-state reconstructions, country-entry edges, persistence and reseeding support, threshold sensitivities, and alternative-root analyses.
- `results/cgmlst/`: cgMLST tree and concordance with the core-SNP analysis.
- `results/model_main/`: project-adjusted posterior summaries, diagnostics, counterfactuals, sampling corrections, and recovery results.
- `results/model_no_project/`: corresponding summaries from the no-project observation-model sensitivity analysis.

## Figures and provenance

- `figures/source_data/`: one or more tabular source files for every main figure panel and a provenance index.
- `figures/main/`: main figures in PNG and PDF formats.
- `figures/supplementary/`: supplementary figures in PNG format.
- `provenance/FILE_MANIFEST_SHA256.tsv`: relative path, byte size, and SHA-256 digest for every released file except the manifest itself.
- `provenance/EXTERNAL_LARGE_OBJECTS.tsv`: checksums and sizes of large fitted objects that were not deposited.
- `provenance/R_SESSION_INFO.txt`: R version, platform, and package-session information.
- `provenance/CITATION_AUDIT.tsv`: bibliographic verification record for manuscript references.
