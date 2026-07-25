#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
python_bin="${PYTHON_BIN:-python3}"
snp_sites_bin="${SNP_SITES_BIN:-snp-sites}"
iqtree_bin="${IQTREE_BIN:-iqtree2}"
work_dir="${repo_root}/pertussis_data/pertussis_gene/genomic_transmission_dynamics/phylogeny/gate_rescue4_depth"
result_dir="${repo_root}/analysis/genomic_transmission_dynamics/phylogeny/results/gate_rescue4_depth"

if [[ ! -s "${work_dir}/core.full.aln" ]]; then
  echo "Missing completed full alignment: ${work_dir}/core.full.aln" >&2
  exit 2
fi
if [[ ! -s "${work_dir}/exclude_high_missing.txt" ]]; then
  echo "Missing exclusion list: ${work_dir}/exclude_high_missing.txt" >&2
  exit 2
fi

seqkit grep \
  --invert-match \
  --pattern-file "${work_dir}/exclude_high_missing.txt" \
  "${work_dir}/core.full.aln" \
  > "${work_dir}/core.filtered.full.aln"

"${python_bin}" \
  "${repo_root}/analysis/genomic_transmission_dynamics/bin/gtd_06_alignment_metrics.py" \
  "${work_dir}/core.filtered.full.aln" \
  "${result_dir}/alignment_missingness_final.tsv" \
  --max-missing 0.20

"${snp_sites_bin}" \
  -o "${work_dir}/core.filtered.snps.aln" \
  "${work_dir}/core.filtered.full.aln"

"${iqtree_bin}" \
  -s "${work_dir}/core.filtered.snps.aln" \
  -m GTR+F+R2+ASC \
  -B 1000 \
  --alrt 1000 \
  -T 16 \
  --prefix "${result_dir}/gtd_gate_rescue4_depth_core_snp" \
  2>&1 | tee "${result_dir}/iqtree.log"

cp "${work_dir}/core.filtered.snps.aln" \
  "${result_dir}/core.filtered.snps.aln"
