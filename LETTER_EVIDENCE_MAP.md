# Letter Evidence Map

Smallest public file set for each central Letter claim.

| Claim | Primary files |
|---|---|
| 989-genome tree, with 774 focal and 215 background genomes | `figures/source_data/figure2_tree_manifest.tsv`; `results/phylogeny/tree_tip_metadata.tsv` |
| All 99 source-published MT28-labelled tips fall in one 288-tip population | `results/lineage_nomenclature/mt28_lineage_nomenclature_validation.json`; `results/lineage_nomenclature/mt28_lineage_tip_crosswalk.tsv` |
| Case data, genome composition, and sampled ancestry are separate model inputs | `data/model_inputs/joint_model_data.json`; `data/model_inputs/genome_observation_strata.tsv`; `data/model_inputs/monthly_import_events.tsv` |
| Cohort-cap weighting does not explain the China/Japan tree-share shifts | `figures/source_data/figure3d_selection_cap_weighted_l10207_shares.tsv`; `results/model_growth_robustness/selection_cap_weighted_l1_02_07_shares.tsv` |
| Australian Ct calibration is aggregate only | `figures/source_data/figure4e_australia_ct_curve.tsv` |
| MT28-associated relative growth is above the lineage reference | `results/model_main/lineage_relative_transmission.tsv`; `figures/source_data/figure3a_lineage_growth_main.tsv` |
| MT28-associated lineage exceeds each comparator | `results/model_main/l10207_pairwise_growth.tsv`; `results/model_no_project/l10207_pairwise_growth.tsv` |
| China and Japan support the signal; Australia is sparse | `results/model_growth_robustness/l1_02_07_growth_robustness.tsv` |
| Country, project, and input perturbations retain the signal | `results/model_growth_robustness/l1_02_07_growth_robustness.tsv`; `results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv` |
| Anchored scenarios are conditional decompositions | `figures/source_data/figure4d_counterfactual_summary.tsv`; `results/model_anchor_sensitivity/anchor_scenario_sensitivity.tsv` |
| Lineage-growth recovery passes; import scale is not identifiable | `results/model_main/recovery_summary.tsv`; `results/model_main/identifiability_domain_interpretation.json` |
| Available molecular annotations do not identify a mechanism | `results/molecular_characterisation/interpretation_limits.tsv`; `results/molecular_characterisation/molecular_characterisation_validation.json` |

Authoritative figure-source maps are in `figures/letter/`.
