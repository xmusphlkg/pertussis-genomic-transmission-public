#!/usr/bin/env python3
"""Scout external ENA studies for release-date-aware pertussis analyses.

This is a metadata-only discovery step. It queries ENA Portal for
*Bordetella pertussis* read runs that became public from 2023 onward, groups
them by study accession, and labels whether each study is already represented
in the current release-date analysis or might be worth follow-up.

The script does not download reads, assemblies, or large sequence files.
"""

from __future__ import annotations

import csv
import re
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


REPO = Path.cwd()
OUTPUT = REPO / "results" / "public_availability" / "external_ena_study_scout.tsv"
PROJECT_INDEX = REPO / "data" / "model_inputs" / "project_index.tsv"

FOCUS_COUNTRIES = {
    "Australia",
    "Belgium",
    "China",
    "France",
    "Japan",
    "Netherlands",
    "United Kingdom",
    "USA",
    "United States",
}

LAB_TITLE_PATTERN = re.compile(
    r"RNA-seq|RNAseq|transcript|gene expression|CHIPseq|knockout|mutant|"
    r"passage|Tohama|D420|BPSM|cultivated|autoagglutination|RACE",
    flags=re.IGNORECASE,
)


def read_known_projects() -> set[str]:
    with PROJECT_INDEX.open(newline="") as handle:
        return {row["project_id"] for row in csv.DictReader(handle, delimiter="\t")}


def fetch_ena_rows() -> list[dict[str, str]]:
    fields = (
        "run_accession,study_accession,study_title,collection_date,country,"
        "first_public"
    )
    params = {
        "result": "read_run",
        "query": 'tax_eq(520) AND first_public>="2023-01-01"',
        "fields": fields,
        "format": "tsv",
        "limit": "0",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(params)
    request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def country_root(value: str) -> str:
    return (value or "").split(":")[0].strip()


def has_post2022_collection(value: str) -> bool:
    match = re.search(r"(20[2-9][3-9])", value or "")
    return bool(match)


def collection_resolution(value: str) -> str:
    value = (value or "").strip().lower()
    if not value or value in {"missing", "not provided", "not applicable", "not collected"}:
        return "missing_or_not_applicable"
    if re.match(r"^\d{4}-\d{2}-\d{2}$", value):
        return "day"
    if re.match(r"^\d{4}-\d{2}$", value):
        return "month"
    if re.match(r"^\d{4}$", value):
        return "year"
    if re.search(r"\d{4}-\d{2}/\d{4}-\d{2}", value):
        return "month_interval"
    return "other"


def classify_study(status: str, countries: set[str], title: str, rows: list[dict[str, str]]) -> str:
    if status == "known_current_cohort":
        return status
    if LAB_TITLE_PATTERN.search(title):
        return "exclude_or_discussion_lab_study"
    if not (countries & FOCUS_COUNTRIES):
        return "discussion_only_nonfocus_geography"

    resolutions = [collection_resolution(row.get("collection_date", "")) for row in rows]
    n_month_or_day = sum(resolution in {"day", "month", "month_interval"} for resolution in resolutions)
    n_post2022 = sum(has_post2022_collection(row.get("collection_date", "")) for row in rows)
    n_with_first_public = sum(bool(row.get("first_public")) for row in rows)

    if n_post2022 == 0 or n_with_first_public == 0:
        return "discussion_only_missing_resurgence_or_release_clock"
    if len(rows) >= 5 and n_month_or_day > 0:
        return "main_analysis_candidate"
    return "sensitivity_candidate"


def main() -> None:
    known = read_known_projects()
    rows = fetch_ena_rows()

    by_study: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_study[row["study_accession"]].append(row)

    output_rows: list[dict[str, object]] = []
    for study, study_rows in sorted(by_study.items()):
        status = "known_current_cohort" if study in known else "new_to_current_project_index"
        countries = {country_root(row.get("country", "")) for row in study_rows if row.get("country")}
        countries.discard("")
        collection_dates = [row.get("collection_date", "") for row in study_rows if row.get("collection_date")]
        first_public = [row.get("first_public", "") for row in study_rows if row.get("first_public")]
        titles = sorted({row.get("study_title", "") for row in study_rows if row.get("study_title")})
        title = titles[0] if titles else ""
        resolutions = [collection_resolution(row.get("collection_date", "")) for row in study_rows]
        tier = classify_study(status, countries, title, study_rows)

        output_rows.append(
            {
                "study_accession": study,
                "project_status": status,
                "candidate_tier": tier,
                "n_runs": len(study_rows),
                "countries": ";".join(sorted(countries)),
                "study_title": title,
                "min_collection_date": min(collection_dates) if collection_dates else "",
                "max_collection_date": max(collection_dates) if collection_dates else "",
                "first_public_min": min(first_public) if first_public else "",
                "first_public_max": max(first_public) if first_public else "",
                "n_collection_day": sum(resolution == "day" for resolution in resolutions),
                "n_collection_month_or_interval": sum(
                    resolution in {"month", "month_interval"} for resolution in resolutions
                ),
                "n_collection_year_only": sum(resolution == "year" for resolution in resolutions),
                "n_collection_missing_or_not_applicable": sum(
                    resolution == "missing_or_not_applicable" for resolution in resolutions
                ),
                "metadata_probe_url": "https://www.ebi.ac.uk/ena/portal/api/filereport?"
                + urllib.parse.urlencode(
                    {
                        "accession": study,
                        "result": "read_run",
                        "fields": "run_accession,sample_accession,study_accession,collection_date,country,first_public",
                        "format": "tsv",
                    }
                ),
            }
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "study_accession",
        "project_status",
        "candidate_tier",
        "n_runs",
        "countries",
        "study_title",
        "min_collection_date",
        "max_collection_date",
        "first_public_min",
        "first_public_max",
        "n_collection_day",
        "n_collection_month_or_interval",
        "n_collection_year_only",
        "n_collection_missing_or_not_applicable",
        "metadata_probe_url",
    ]
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"scouted {len(rows)} ENA runs across {len(by_study)} studies")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
