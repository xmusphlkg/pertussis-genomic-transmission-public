# Letter evidence map

This map links each central Letter claim to the smallest public evidence set.
It is a navigation aid, not a substitute for the manuscript's methods and
limitations.

## Claims and evidence

| Letter claim | Primary public evidence | Supporting audit |
|---|---|---|
| The frozen tree contains 989 genomes: 774 focal and 215 background | `figures/source_data/figure2_tree_manifest.tsv`; `results/phylogeny/tree_tip_metadata.tsv` | `results/phylogeny/TREE_QA_REPORT.md` |
| All 99 source-published MT28-labelled tips map to one 288-tip monophyletic population | `results/lineage_nomenclature/mt28_lineage_nomenclature_validation.json` | `mt28_lineage_tip_crosswalk.tsv`; `published_mt28_label_counts.tsv` in the same directory |
| Case dynamics, genomic composition, and sampled ancestry have distinct inferential roles | `data/model_inputs/joint_model_data.json`; `scripts/model/gtd_joint_transmission_sampling.stan` | `data/model_inputs/genome_observation_strata.tsv`; `monthly_import_events.tsv` |
| Deterministic cohort caps do not generate the China/Japan sampled-tree share shifts | `figures/source_data/figure3d_selection_cap_weighted_l10207_shares.tsv` | `results/model_growth_robustness/selection_cap_weighted_l1_02_07_shares.tsv` |
| Ct is associated with complete-profile recovery in the Australian calibration | `figures/source_data/figure4e_australia_ct_curve.tsv` | Aggregate output only; the third-party specimen table is not redistributed |
| The shared MT28-associated multiplier is above the lineage reference | `results/model_main/lineage_relative_transmission.tsv` | `figures/source_data/figure3a_lineage_growth_main.tsv` |
| The MT28-associated lineage exceeds every comparator directly | `results/model_main/l10207_pairwise_growth.tsv`; `results/model_no_project/l10207_pairwise_growth.tsv` | `figures/source_data/figure3b_l10207_pairwise_growth.tsv` |
| China and Japan are above reference; Australia is uninformative | `results/model_growth_robustness/l1_02_07_growth_robustness.tsv` | `figures/source_data/figure3c_l10207_growth_robustness.tsv` |
| Country and dominant-project omissions retain the signal | `results/model_growth_robustness/l1_02_07_growth_robustness.tsv` | `fit_diagnostics.tsv` and `run_configuration.tsv` in the same directory |
| Six full input refits retain complete intervals above one | `results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv` | `model_input_sensitivity_diagnostics.tsv` in the same directory |
| Conditional anchored scenarios are decompositions, not causal effects | `figures/source_data/figure4d_counterfactual_summary.tsv` | `results/model_anchor_sensitivity/anchor_scenario_sensitivity.tsv` |
| Lineage growth passes the recovery gate; import scale does not | `results/model_main/recovery_summary.tsv`; `identifiability_domain_interpretation.json` | `figures/source_data/figure4f_identifiability_recovery.tsv` |
| Available molecular annotations do not identify an expansion mechanism | `results/molecular_characterisation/interpretation_limits.tsv` | `annotation_coverage_by_project.tsv`; `molecular_characterisation_validation.json` |

## Figure 1 panel map

| Panel | Evidence |
|---|---|
| A — core-SNP tree and MT28-associated clade | `figure2_tree_manifest.tsv`, `figure2a_tree_tip_metadata.tsv`, `figure2b_tip_ancestry_support.tsv` |
| B — monthly reported cases | `figure1a_cases.tsv` |
| C — annual public genomes | `figure1b_annual_genomes.tsv` |
| D — unweighted and cohort-cap-weighted shares | `figure3d_selection_cap_weighted_l10207_shares.tsv` |
| E — Ct-dependent complete-profile recovery | `figure4e_australia_ct_curve.tsv` |

All panel files in this section are under `figures/source_data/`.

## Figure 2 panel map

| Panel | Evidence |
|---|---|
| A — shared lineage multipliers | `figures/source_data/figure3a_lineage_growth_main.tsv` |
| B — direct pairwise ratios | `figures/source_data/figure3b_l10207_pairwise_growth.tsv` |
| C — country-only fits | `figures/source_data/figure3c_l10207_growth_robustness.tsv` |
| D — country/project/input robustness | preceding robustness table plus `results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv` |
| E — conditional anchored scenarios | `figures/source_data/figure4d_counterfactual_summary.tsv` |
| F — recovery boundary | `figures/source_data/figure4f_identifiability_recovery.tsv` |

The authoritative machine-readable maps are
`figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv` and
`figures/letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv`.
