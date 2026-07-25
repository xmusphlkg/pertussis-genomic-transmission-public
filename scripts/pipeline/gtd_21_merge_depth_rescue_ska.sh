#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
ska_bin="${SKA_BIN:-ska}"
python_bin="${PYTHON_BIN:-python3}"
source_skf="${repo_root}/pertussis_data/pertussis_gene/genomic_transmission_dynamics/phylogeny/gate_rescue4/gate_rescue4.skf"
reference="${repo_root}/pertussis_data/pertussis_gene/workflow/phylo/core.ref.fa"
rescue_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/depth_rescue"
work_dir="${repo_root}/pertussis_data/pertussis_gene/genomic_transmission_dynamics/phylogeny/gate_rescue4_depth"
result_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/results/gate_rescue4_depth"

mkdir -p "${work_dir}" "${result_dir}"

if [[ ! -s "${work_dir}/depth_replacements.skf" ]]; then
  "${ska_bin}" build \
    -o "${work_dir}/depth_replacements" \
    -k 31 --min-count 3 --min-qual 20 --qual-filter middle --threads 28 \
    -f "${rescue_dir}/completed_depth_rescue_ska_input.tsv" \
    >"${result_dir}/ska_depth_build.log" 2>&1
fi

if [[ ! -s "${work_dir}/base_without_depth_replacements.skf" ]]; then
  mapfile -t rescue_sample_ids \
    < "${rescue_dir}/completed_depth_rescue_sample_ids.txt"
  "${ska_bin}" delete \
    --skf-file "${source_skf}" \
    -o "${work_dir}/base_without_depth_replacements" \
    "${rescue_sample_ids[@]}" \
    >"${result_dir}/ska_delete_old_depth.log" 2>&1
fi

"${ska_bin}" merge \
  -o "${work_dir}/gate_rescue4_depth" \
  "${work_dir}/base_without_depth_replacements.skf" \
  "${work_dir}/depth_replacements.skf" \
  >"${result_dir}/ska_merge.log" 2>&1

"${ska_bin}" map \
  --repeat-mask --ambig-mask \
  -o "${work_dir}/core.full.aln" \
  "${reference}" "${work_dir}/gate_rescue4_depth.skf" \
  >"${result_dir}/ska_map.log" 2>&1

"${python_bin}" \
  "${repo_root}/analysis/genomic_transmission_dynamics/bin/gtd_06_alignment_metrics.py" \
  "${work_dir}/core.full.aln" \
  "${result_dir}/alignment_missingness_initial.tsv" \
  --max-missing 0.20 \
  --exclude-list "${work_dir}/exclude_high_missing.txt"
