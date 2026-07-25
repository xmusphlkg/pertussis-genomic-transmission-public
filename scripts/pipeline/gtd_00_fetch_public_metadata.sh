#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analysis_dir="$(cd "${script_dir}/.." && pwd)"
raw_dir="${analysis_dir}/inputs/raw"
mkdir -p "${raw_dir}"

curl -L --fail --retry 3 \
  --output "${raw_dir}/australia_fong_2026_appendix2.xlsx" \
  "https://ars.els-cdn.com/content/image/1-s2.0-S2666524725002149-mmc2.xlsx"

fetch_ena_project() {
  local accession="$1"
  local output_name="$2"
  curl -L --fail --retry 3 --get \
    "https://www.ebi.ac.uk/ena/portal/api/filereport" \
    --data-urlencode "accession=${accession}" \
    --data-urlencode "result=read_run" \
    --data-urlencode \
    "fields=run_accession,sample_accession,study_accession,collection_date,country,first_public,fastq_ftp" \
    --data-urlencode "format=tsv" \
    --output "${raw_dir}/${output_name}"
}

fetch_ena_project "PRJNA1199062" "australia_prjna1199062_ena_runs.tsv"
fetch_ena_project "PRJEB18624" "japan_prjeb18624_ena_runs.tsv"
fetch_ena_project "PRJDB20292" "japan_prjdb20292_ena_runs.tsv"
fetch_ena_project "PRJDB20413" "japan_prjdb20413_ena_runs.tsv"
fetch_ena_project "PRJDB34249" "japan_prjdb34249_ena_runs.tsv"
fetch_ena_project "PRJDB35593" "japan_prjdb35593_ena_runs.tsv"
fetch_ena_project "PRJDB37898" "japan_prjdb37898_ena_runs.tsv"
