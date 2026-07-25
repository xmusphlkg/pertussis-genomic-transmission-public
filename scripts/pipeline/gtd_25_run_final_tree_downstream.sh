#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
bio_python="${PYTHON_BIN:-python3}"
pastml="${PASTML_BIN:-pastml}"

tree_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/results/gate_rescue4_depth"
tree="${tree_dir}/gtd_gate_rescue4_depth_core_snp.treefile"
iqtree_report="${tree_dir}/gtd_gate_rescue4_depth_core_snp.iqtree"
metadata="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/qc_gate_rescue4/uniform_sequence_qc.tsv"
lineages="${repo_root}/analysis/genomic_transmission_dynamics/lineages_gate_rescue4_depth/primary_finalized/model_lineage_assignments.tsv"
phylogeography="${repo_root}/analysis/genomic_transmission_dynamics/phylogeography_gate_rescue4_depth"
pastml_dir="${phylogeography}/pastml_primary"

if [[ ! -s "${tree}" || ! -s "${iqtree_report}" ]]; then
  echo "Final IQ-TREE outputs are incomplete; refusing downstream inference." >&2
  exit 2
fi

"${bio_python}" \
  "${repo_root}/analysis/genomic_transmission_dynamics/bin/gtd_08_tree_qa.py" \
  --tree "${tree}" \
  --qc "${metadata}" \
  --output-dir "${tree_dir}" \
  --permutations 1000 \
  --seed 20260725

Rscript \
  "${repo_root}/analysis/genomic_transmission_dynamics/bin/gtd_10_prepare_phylogeography.R" \
  "${tree}" \
  "${tree_dir}/tree_tip_metadata.tsv" \
  "${phylogeography}"

mkdir -p "${pastml_dir}"
"${pastml}" \
  --tree "${phylogeography}/gtd_primary_final_midpoint_rooted.tree" \
  --data "${phylogeography}/tip_geography.tsv" \
  --columns geo_state \
  --prediction_method MPPA \
  --model F81 \
  --out_data "${pastml_dir}/ancestral_states.tsv" \
  --work_dir "${pastml_dir}" \
  --threads 2

marginals="${pastml_dir}/marginal_probabilities.character_geo_state.model_F81.tab"
if [[ ! -s "${marginals}" ]]; then
  echo "PastML marginal-probability output is missing." >&2
  exit 2
fi

for threshold in 0.5 0.7 0.9; do
  label="${threshold/./_}"
  "${bio_python}" \
    "${repo_root}/analysis/genomic_transmission_dynamics/bin/gtd_11_infer_phylogeographic_events.py" \
    "${phylogeography}/gtd_primary_final_midpoint_rooted.tree" \
    "${marginals}" \
    "${tree_dir}/tree_tip_metadata.tsv" \
    "${lineages}" \
    "${phylogeography}/events_thr${label}" \
    --transition-threshold "${threshold}"
done

echo "Final tree QA and midpoint-rooted phylogeography completed."
