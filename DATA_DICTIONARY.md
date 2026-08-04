# Data dictionary

This dictionary is organised around the Letter rather than the former
long-article figure sequence. Tab-separated files contain a header row; JSON
files contain model inputs, diagnostics, or machine-readable validation
records.

Reader-facing documents use **MT28-associated genomic lineage**. Its frozen
machine identifier is `L1_02.07`. Comparator lineages A–C correspond to
`L1_01.02`, `L1_02.05`, and `L1_02.06`; machine-readable outputs retain those
identifiers.

## Letter figures

- `figures/letter/Figure_1_observation_structure.{png,pdf}`: phylogenetic
  definition plus case, genome-timing, cohort-cap, and Ct-calibration context.
- `figures/letter/Figure_2_growth_robustness.{png,pdf}`: shared and pairwise
  growth, country and omission analyses, input refits, conditional scenarios,
  and recovery boundaries.
- `figures/letter/Supplementary_Figure_S1_genomic_validation.png`: cgMLST,
  nearest-neighbour, and lineage-partition validation.
- `figures/letter/Supplementary_Figure_S2_posterior_predictive_checks.png`:
  fitted trajectories and posterior-predictive metrics.
- `figures/letter/Supplementary_Figure_S3_ancestry_input_sensitivity.png`:
  sampled-ancestry rank and model-input sensitivity.
- `figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv` and
  `LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv`: figure-to-source-file maps.

## Panel-level source data

The compact Letter source tables are under `figures/source_data/`:

- `figure2_tree_manifest.tsv`, `figure2a_tree_tip_metadata.tsv`, and
  `figure2b_tip_ancestry_support.tsv`: Figure 1 phylogeny and annotations.
- `figure1a_cases.tsv` and `figure1b_annual_genomes.tsv`: reported cases and
  public-genome timing.
- `figure3d_selection_cap_weighted_l10207_shares.tsv`: unweighted and
  deterministic selection-cap-weighted sampled-tree shares.
- `figure4e_australia_ct_curve.tsv`: aggregate Ct-dependent complete-profile
  recovery curve.
- `figure3a_lineage_growth_main.tsv`: shared lineage multipliers.
- `figure3b_l10207_pairwise_growth.tsv`: direct MT28-associated-to-comparator
  posterior ratios.
- `figure3c_l10207_growth_robustness.tsv`: country, country-omission, and
  dominant-project-omission estimates.
- `figure4d_counterfactual_summary.tsv`: conditional anchored scenario
  summaries.
- `figure4f_identifiability_recovery.tsv`: lineage-growth and import-scale
  recovery metrics.
- `figure4abc_monthly_counterfactuals.tsv`: monthly trajectories used in
  Supplementary Figure S2.

## Derived public-data tables

- `data/derived/source_registry.tsv`: public source datasets, URLs, roles, and
  frozen checksums.
- `data/derived/transmission_genome_records.tsv`: harmonised focal-genome
  records, accessions, sampling-date resolution, project, and eligibility.
- `data/derived/global_tree_background_records.tsv`: stratified background
  genomes used only for topology and sampled ancestry.
- `data/derived/country_month_cases.tsv`: monthly case series for Australia,
  China, and Japan.
- `data/derived/country_month_case_genome_panel.tsv`: linked case and genome
  observation panel.
- `data/derived/country_month_genome_strata.tsv`: project- and
  lineage-stratified genome counts.
- `data/derived/country_data_gate_summary.tsv`: country eligibility gates.
- `data/derived/model_lineage_display_map.tsv`: internal-to-reader-facing
  lineage mapping.
- `data/derived/preliminary_lineage_screen.tsv`: pre-freeze screening counts;
  these labels did not define model lineages.

Local private-compute paths are not reusable data and are omitted from released
tables. The Australian direct-specimen source table is not redistributed.
Public BioSample/run accessions used in the tree remain, but exact Australian
Ct values, source-specific aliases, specimen types, profile-success fields, and
article-supplied marker annotations are removed from per-record released
tables. Only the aggregate Ct-dependent complete-profile recovery curve is
released.

## Frozen model inputs

- `data/model_inputs/joint_model_data.json`: complete input for the primary
  Stan model.
- `data/model_inputs/joint_model_data_validation.json`: dimensions and
  consistency checks.
- `data/model_inputs/genome_observation_strata.tsv`: country-month-project
  lineage counts used in the Dirichlet-multinomial observation model.
- `data/model_inputs/initial_lineage_prior.tsv`: strictly pre-model initial
  lineage-state prior.
- `data/model_inputs/monthly_import_events.tsv`: fixed sampled-ancestry
  covariates. The legacy filename is retained for reproducibility; these rows
  are not observed imported-case counts.
- `data/model_inputs/project_index.tsv`: stable project identifiers.
- `data/model_inputs/transmission_spline_basis.tsv`: monthly spline basis.

The corresponding frozen perturbation inputs are under
`data/model_sensitivity/`: transition thresholds 0.7 and 0.9, midpoint and
interval-uniform event timing, alternative rooting, historical and symmetric
initial states, and the pre-rescue cohort.

## Core Letter results

- `results/lineage_nomenclature/`: published MT28-label counts, tip-level
  crosswalk, MRCA match, and validation of the display name.
- `results/model_main/`: primary project-adjusted posterior summaries,
  diagnostics, pairwise contrasts, conditional scenarios, and recovery.
- `results/model_no_project/`: corresponding no-project summaries.
- `results/model_growth_robustness/`: country-only, leave-one-country-out,
  dominant-project-omission, and selection-cap checks. Files beginning
  `legacy_omit_country_JPN_` document the failed initial run and do not define
  the reported estimate.
- `results/model_input_sensitivity_summary/`: six full input-refit summaries.
- `results/model_anchor_sensitivity/`: conditional scenario sensitivity across
  four anchor dates.
- `results/cohort_gate_sensitivity/`: pre-rescue versus final-cohort results and
  input checksums.
- `results/model_initial_state_prior_sensitivity/`: historical versus symmetric
  initial-state checks.
- `results/model_temporal_holdout_2025/`: leakage-audited temporal stress test;
  this is an extended validation, not a Letter forecast claim.
- `results/molecular_characterisation/`: marker composition, annotation
  coverage, project confounding, and interpretation limits.

## Extended audit results

- `results/phylogeny/`: 989-tip tree, 11,550-position core-SNP alignment,
  metadata, uniform QC, and tree QA.
- `results/lineages/`: primary and alternative rhierBAPS assignments and
  partition sensitivity.
- `results/cgmlst/`: cgMLST tree and concordance with the core-SNP analysis.
- `results/phylogeography/` and `results/phylogeography_alternative_root/`:
  sampled geographic-state reconstruction, threshold checks, and root
  sensitivity.
- `results/phylogeography_sensitivity_summary/`: compact cross-topology and
  sampled-reference sensitivity summaries.
- `data/model_sensitivity/cohort_gate_pre_rescue/`: compact input and output
  package for the frozen pre-rescue comparison.

## Code and provenance

- `scripts/model/gtd_joint_transmission_sampling.stan`: joint model source.
- `scripts/figures/`: Letter figure renderers.
- `scripts/lib/`: shared anchored-scenario functions.
- `scripts/pipeline/`: upstream construction and sensitivity programs.
- `scripts/qa/sanitise_public_metadata.py`: removes non-redistributable
  specimen-level fields from released TSVs.
- `scripts/qa/validate_public_package.py`: Letter-package validation.
- `tests/`: numerical and structural regression tests.
- `provenance/FILE_MANIFEST_SHA256.tsv`: released-file size and SHA-256 index.
- `provenance/EXTERNAL_LARGE_OBJECTS.tsv`: excluded fitted objects and hashes;
  private workspace paths are not released.
- `provenance/R_SESSION_INFO.txt`: principal R environment.
- `provenance/CITATION_AUDIT.tsv`: bibliographic verification record.
