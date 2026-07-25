#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
python_bin="${PYTHON_BIN:-python3}"
analysis_root="${repo_root}/analysis/genomic_transmission_dynamics"
cases="${analysis_root}/derived/country_month_cases.tsv"
lineages="${analysis_root}/lineages_gate_rescue4_depth/primary_finalized/model_lineage_assignments.tsv"
phylo="${analysis_root}/phylogeography_gate_rescue4_depth/events_thr0_5"
tip_support="${phylo}/tip_persistence_reseeding_support.tsv"
events="${phylo}/introduction_events.tsv"
model_root="${analysis_root}/model_gate_rescue4_depth"
stan="${analysis_root}/model/gtd_joint_transmission_sampling.stan"

for required in "${cases}" "${lineages}" "${tip_support}" "${events}" "${stan}"; do
  if [[ ! -s "${required}" ]]; then
    echo "Required final-analysis input is missing: ${required}" >&2
    exit 2
  fi
done

"${python_bin}" \
  "${analysis_root}/bin/gtd_12_build_joint_model_data.py" \
  "${cases}" "${lineages}" "${tip_support}" "${model_root}/data"

"${python_bin}" \
  "${analysis_root}/bin/gtd_12_build_joint_model_data.py" \
  "${cases}" "${lineages}" "${tip_support}" "${model_root}/data_no_project" \
  --disable-project-effects

Rscript "${analysis_root}/bin/gtd_13_fit_joint_model.R" \
  "${stan}" "${model_root}/data/joint_model_data.json" \
  "${model_root}/fit_main" 4 2000 20260725

Rscript "${analysis_root}/bin/gtd_13_fit_joint_model.R" \
  "${stan}" "${model_root}/data_no_project/joint_model_data.json" \
  "${model_root}/fit_no_project" 4 2000 20260731

Rscript "${analysis_root}/bin/gtd_15_summarise_joint_model.R" \
  "${model_root}/fit_main/posterior_outputs.rds" \
  "${model_root}/data/joint_model_data.json" \
  "${model_root}/data/genome_observation_strata.tsv" \
  "${events}" \
  "${analysis_root}/results_gate_rescue4_depth/main"

Rscript "${analysis_root}/bin/gtd_15_summarise_joint_model.R" \
  "${model_root}/fit_no_project/posterior_outputs.rds" \
  "${model_root}/data_no_project/joint_model_data.json" \
  "${model_root}/data_no_project/genome_observation_strata.tsv" \
  "${events}" \
  "${analysis_root}/results_gate_rescue4_depth/no_project"

Rscript "${analysis_root}/bin/gtd_16_identifiability_recovery.R" \
  "${model_root}/fit_main/compiled_model.rds" \
  "${model_root}/fit_main/joint_model_fit.rds" \
  "${model_root}/data/joint_model_data.json" \
  "${model_root}/identifiability_recovery" \
  6 2 1200 20260801

echo "Final joint-model fits, counterfactuals, and recovery audit completed."
