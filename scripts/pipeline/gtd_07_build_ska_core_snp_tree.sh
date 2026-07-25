#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
ska_bin="${SKA_BIN:-ska}"
python_bin="${PYTHON_BIN:-python3}"
snp_sites_bin="${SNP_SITES_BIN:-snp-sites}"
iqtree_bin="${IQTREE_BIN:-iqtree2}"
qc_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/qc"
run_label="${1:-primary}"
quality_filter="${3:-middle}"
work_dir="${repo_root}/pertussis_data/pertussis_gene/genomic_transmission_dynamics/phylogeny/${run_label}"
result_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/results/${run_label}"
reference="${repo_root}/pertussis_data/pertussis_gene/workflow/phylo/core.ref.fa"
input_list="${2:-${qc_dir}/ska_primary_input.tsv}"

mkdir -p "${work_dir}" "${result_dir}"

if [[ ! -s "${work_dir}/${run_label}.skf" ]]; then
  "${ska_bin}" build \
    -o "${work_dir}/${run_label}" \
    -k 31 \
    --min-count 3 \
    --min-qual 20 \
    --qual-filter "${quality_filter}" \
    --threads 32 \
    -f "${input_list}" \
    2>&1 | tee "${result_dir}/ska_build.log"
fi

if [[ "${SKIP_NK:-0}" != "1" ]]; then
  "${ska_bin}" nk "${work_dir}/${run_label}.skf" \
    > "${result_dir}/ska_kmer_summary.txt"
fi

if [[ ! -s "${work_dir}/core.full.aln" ]]; then
  "${ska_bin}" map \
    --repeat-mask \
    --ambig-mask \
    --threads 32 \
    -o "${work_dir}/core.full.aln" \
    "${reference}" \
    "${work_dir}/${run_label}.skf" \
    2>&1 | tee "${result_dir}/ska_map.log"
fi

"${python_bin}" "${script_dir}/gtd_06_alignment_metrics.py" \
  "${work_dir}/core.full.aln" \
  "${result_dir}/alignment_missingness_initial.tsv" \
  --max-missing 0.20 \
  --exclude-list "${work_dir}/exclude_high_missing.txt"

if [[ -s "${work_dir}/exclude_high_missing.txt" ]]; then
  mapfile -t excluded_samples < "${work_dir}/exclude_high_missing.txt"
  "${ska_bin}" delete \
    --skf-file "${work_dir}/${run_label}.skf" \
    -o "${work_dir}/${run_label}.filtered" \
    "${excluded_samples[@]}"
  "${ska_bin}" map \
    --repeat-mask \
    --ambig-mask \
    --threads 32 \
    -o "${work_dir}/core.filtered.full.aln" \
    "${reference}" \
    "${work_dir}/${run_label}.filtered.skf"
else
  cp "${work_dir}/core.full.aln" "${work_dir}/core.filtered.full.aln"
  cp "${work_dir}/${run_label}.skf" "${work_dir}/${run_label}.filtered.skf"
fi

"${python_bin}" "${script_dir}/gtd_06_alignment_metrics.py" \
  "${work_dir}/core.filtered.full.aln" \
  "${result_dir}/alignment_missingness_final.tsv" \
  --max-missing 0.20

"${snp_sites_bin}" \
  -o "${work_dir}/core.filtered.snps.aln" \
  "${work_dir}/core.filtered.full.aln"

"${iqtree_bin}" \
  -s "${work_dir}/core.filtered.snps.aln" \
  -m MFP+ASC \
  -B 1000 \
  --alrt 1000 \
  -T 16 \
  --prefix "${result_dir}/gtd_${run_label}_core_snp" \
  2>&1 | tee "${result_dir}/iqtree.log"

cp "${work_dir}/core.filtered.snps.aln" "${result_dir}/core.filtered.snps.aln"

"${ska_bin}" distance \
  --min-freq 0.90 \
  --threads 32 \
  -o "${result_dir}/ska_pairwise_distances.tsv" \
  "${work_dir}/${run_label}.filtered.skf"
