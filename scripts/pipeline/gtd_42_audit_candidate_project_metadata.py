#!/usr/bin/env python3
"""Audit shortlisted public pertussis projects for possible data extension.

This is a metadata-only follow-up to ``gtd_41_scout_external_ena_studies.py``.
It fetches ENA read-run metadata for shortlisted BioProjects/studies, compares
run/sample/study identifiers with the current 989-genome release-clock cohort,
and writes a compact decision table.

The script deliberately does not download FASTQ, FASTA, or assemblies. It also
does not perform lineage placement; any detection-clock implication is therefore
reported as "possible if lineage-compatible" rather than as a lineage claim.
"""

from __future__ import annotations

import csv
import datetime as dt
import os
import re
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CURRENT_COHORT = REPO / "data" / "derived" / "public_genome_availability.tsv"
CASE_THRESHOLDS = REPO / "results" / "public_availability" / "case_thresholds.tsv"
RUN_OUTPUT = REPO / "results" / "public_availability" / "candidate_project_run_metadata.tsv"
SUMMARY_OUTPUT = REPO / "results" / "public_availability" / "candidate_project_metadata_audit.tsv"

SHORTLIST = [
    "PRJNA1071282",
    "PRJNA1193776",
    "PRJDB39872",
    "PRJNA1455114",
    "PRJEB88325",
    "PRJNA870170",
]

FOCUS_COUNTRY_TO_ISO3 = {
    "Australia": "AUS",
    "Belgium": "BEL",
    "China": "CHN",
    "Japan": "JPN",
}


def fetch_ena_read_run(project_id: str) -> list[dict[str, str]]:
    fields = (
        "run_accession,sample_accession,study_accession,study_title,"
        "tax_id,scientific_name,collection_date,country,first_public"
    )
    params = {
        "accession": project_id,
        "result": "read_run",
        "fields": fields,
        "format": "tsv",
    }
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + urllib.parse.urlencode(params)
    text = ""
    last_error: Exception | None = None
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                text = response.read().decode("utf-8")
            break
        except Exception as exc:  # pragma: no cover - network branch
            last_error = exc
            time.sleep(1.5 * (attempt + 1))
    if not text.strip() and last_error is not None:
        print(f"warning: ENA fetch failed for {project_id}: {last_error}")
    if not text.strip():
        return []
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    for row in rows:
        row["query_project_id"] = project_id
        row["metadata_probe_url"] = url
    return rows


def is_explicit_pertussis(row: dict[str, str]) -> bool:
    tax_id = (row.get("tax_id") or "").strip()
    scientific_name = (row.get("scientific_name") or "").strip().lower()
    return tax_id == "520" or scientific_name == "bordetella pertussis"


def read_current_ids() -> dict[str, set[str]]:
    ids = {
        "project_id": set(),
        "run_accession": set(),
        "biosample_accession": set(),
        "assembly_accession": set(),
    }
    with CURRENT_COHORT.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for field in ids:
                value = (row.get(field) or "").strip()
                if value:
                    ids[field].add(value)
    return ids


def read_case_thresholds() -> dict[str, dict[str, str]]:
    if not CASE_THRESHOLDS.exists():
        return {}
    with CASE_THRESHOLDS.open(newline="") as handle:
        return {row["country_iso3"]: row for row in csv.DictReader(handle, delimiter="\t")}


def country_root(value: str) -> str:
    return (value or "").split(":")[0].strip()


def collection_resolution(value: str) -> str:
    value = (value or "").strip().lower()
    if not value or value in {"missing", "not provided", "not applicable", "not collected"}:
        return "missing_or_not_applicable"
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return "day"
    if re.fullmatch(r"\d{4}-\d{2}", value):
        return "month"
    if re.fullmatch(r"\d{4}", value):
        return "year"
    if re.search(r"\d{4}-\d{2}/\d{4}-\d{2}", value):
        return "month_interval"
    return "other"


def parse_collection_lower(value: str) -> dt.date | None:
    value = (value or "").strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return dt.date.fromisoformat(value)
    if re.fullmatch(r"\d{4}-\d{2}", value):
        return dt.date.fromisoformat(value + "-01")
    if re.fullmatch(r"\d{4}", value):
        return dt.date(int(value), 1, 1)
    match = re.match(r"^(\d{4}-\d{2})/(\d{4}-\d{2})$", value)
    if match:
        return dt.date.fromisoformat(match.group(1) + "-01")
    return None


def parse_public_date(value: str) -> dt.date | None:
    value = (value or "").strip()
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return dt.date.fromisoformat(value)
    return None


def min_date(values: list[dt.date]) -> str:
    return min(values).isoformat() if values else ""


def max_date(values: list[dt.date]) -> str:
    return max(values).isoformat() if values else ""


def lead_days(public_date: dt.date | None, threshold: str) -> str:
    if public_date is None or not threshold:
        return ""
    return str((dt.date.fromisoformat(threshold) - public_date).days)


def classify_summary(project_id: str, rows: list[dict[str, str]], overlap_any: int) -> tuple[str, str]:
    countries = {country_root(row.get("country", "")) for row in rows if row.get("country")}
    resolutions = [collection_resolution(row.get("collection_date", "")) for row in rows]
    n_day_or_month = sum(resolution in {"day", "month", "month_interval"} for resolution in resolutions)
    n_post2022 = sum(
        (parse_collection_lower(row.get("collection_date", "")) or dt.date(1900, 1, 1))
        >= dt.date(2023, 1, 1)
        for row in rows
    )

    if not rows:
        return "exclude", "No explicit B. pertussis ENA read_run rows returned."
    if project_id == "PRJNA1071282":
        return (
            "boundary_dataset_pending_date_resolution",
            "Sixteen frozen tips are already represented; the complete project expansion has year-level collection dates and is unsuitable for the month-level detection-clock analysis.",
        )
    if overlap_any == len(rows):
        return "partial_current_project_extension", "All explicit B. pertussis rows overlap current cohort identifiers."
    if overlap_any:
        return "partial_current_project_extension", "Some explicit B. pertussis rows overlap current cohort identifiers."
    if project_id == "PRJEB88325":
        return (
            "supplement_sensitivity",
            "Large Belgium dataset; useful for release-lag generalisability, but outside current three-country case-clock model.",
        )
    if project_id == "PRJNA870170":
        return (
            "discussion_or_sensitivity_after_taxon_check",
            "Small Australia candidate; study title is broad, so inspect taxon/sample metadata before any claim.",
        )
    if n_day_or_month and n_post2022 and countries & {"China", "Japan", "Australia"}:
        return (
            "supplement_sensitivity_then_lineage_place",
            "Metadata are good enough for a release-lag sensitivity table, but detection-date claims require lineage placement.",
        )
    return (
        "discussion_only_pending_metadata",
        "Current metadata are insufficient for primary detection-clock analysis.",
    )


def main() -> None:
    current_ids = read_current_ids()
    thresholds = read_case_thresholds()

    all_rows: list[dict[str, str]] = []
    if os.environ.get("GTD42_OFFLINE_CACHE") == "1":
        if not RUN_OUTPUT.exists():
            raise SystemExit(f"Missing offline candidate cache: {RUN_OUTPUT}")
        with RUN_OUTPUT.open(newline="") as handle:
            cached_rows = list(csv.DictReader(handle, delimiter="\t"))
        all_rows = [
            row
            for row in cached_rows
            if row.get("query_project_id") in SHORTLIST and is_explicit_pertussis(row)
        ]
    else:
        for project_id in SHORTLIST:
            all_rows.extend(row for row in fetch_ena_read_run(project_id) if is_explicit_pertussis(row))

    for row in all_rows:
        country = country_root(row.get("country", ""))
        row["country_root"] = country
        row["country_iso3_guess"] = FOCUS_COUNTRY_TO_ISO3.get(country, "")
        row["collection_resolution"] = collection_resolution(row.get("collection_date", ""))
        row["overlap_run_accession"] = str(row.get("run_accession", "") in current_ids["run_accession"])
        row["overlap_sample_accession"] = str(row.get("sample_accession", "") in current_ids["biosample_accession"])
        row["overlap_study_project"] = str(row.get("study_accession", "") in current_ids["project_id"])

    by_project: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        by_project[row["query_project_id"]].append(row)

    summary_rows: list[dict[str, str]] = []
    for project_id in SHORTLIST:
        rows = by_project.get(project_id, [])
        countries = sorted({row["country_root"] for row in rows if row["country_root"]})
        titles = sorted({row.get("study_title", "") for row in rows if row.get("study_title")})
        resolutions = Counter(row["collection_resolution"] for row in rows)
        collection_lowers = [
            parsed
            for parsed in (parse_collection_lower(row.get("collection_date", "")) for row in rows)
            if parsed is not None
        ]
        public_dates = [
            parsed
            for parsed in (parse_public_date(row.get("first_public", "")) for row in rows)
            if parsed is not None
        ]
        overlap_run = sum(row["overlap_run_accession"] == "True" for row in rows)
        overlap_sample = sum(row["overlap_sample_accession"] == "True" for row in rows)
        overlap_project = sum(row["overlap_study_project"] == "True" for row in rows)
        overlap_any = sum(
            row["overlap_run_accession"] == "True"
            or row["overlap_sample_accession"] == "True"
            or row["overlap_study_project"] == "True"
            for row in rows
        )
        tier, decision_note = classify_summary(project_id, rows, overlap_any)

        iso3s = sorted({row["country_iso3_guess"] for row in rows if row["country_iso3_guess"]})
        threshold_values = [
            thresholds[iso3]["first_post2022_month_above_2019_max"]
            for iso3 in iso3s
            if iso3 in thresholds and thresholds[iso3].get("first_post2022_month_above_2019_max")
        ]
        peak_values = [
            thresholds[iso3]["post2022_peak_month"]
            for iso3 in iso3s
            if iso3 in thresholds and thresholds[iso3].get("post2022_peak_month")
        ]
        earliest_public = min(public_dates) if public_dates else None

        summary_rows.append(
            {
                "project_id": project_id,
                "n_runs": str(len(rows)),
                "countries": ";".join(countries),
                "study_title": titles[0] if titles else "",
                "collection_lower_min": min_date(collection_lowers),
                "collection_lower_max": max_date(collection_lowers),
                "first_public_min": min_date(public_dates),
                "first_public_max": max_date(public_dates),
                "n_collection_day": str(resolutions["day"]),
                "n_collection_month_or_interval": str(resolutions["month"] + resolutions["month_interval"]),
                "n_collection_year_only": str(resolutions["year"]),
                "n_collection_missing_or_not_applicable": str(resolutions["missing_or_not_applicable"]),
                "overlap_run_count": str(overlap_run),
                "overlap_sample_count": str(overlap_sample),
                "overlap_project_count": str(overlap_project),
                "overlap_any_count": str(overlap_any),
                "earliest_public_lead_to_focus_threshold_days": ";".join(
                    lead_days(earliest_public, threshold) for threshold in threshold_values
                ),
                "earliest_public_lead_to_focus_peak_days": ";".join(
                    lead_days(earliest_public, peak) for peak in peak_values
                ),
                "recommended_tier": tier,
                "decision_note": decision_note,
                "metadata_probe_url": rows[0]["metadata_probe_url"] if rows else "",
            }
        )

    RUN_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    run_fields = [
        "query_project_id",
        "run_accession",
        "sample_accession",
        "study_accession",
        "study_title",
        "tax_id",
        "scientific_name",
        "collection_date",
        "collection_resolution",
        "country",
        "country_root",
        "country_iso3_guess",
        "first_public",
        "overlap_run_accession",
        "overlap_sample_accession",
        "overlap_study_project",
        "metadata_probe_url",
    ]
    with RUN_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=run_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in run_fields} for row in all_rows)

    summary_fields = [
        "project_id",
        "n_runs",
        "countries",
        "study_title",
        "collection_lower_min",
        "collection_lower_max",
        "first_public_min",
        "first_public_max",
        "n_collection_day",
        "n_collection_month_or_interval",
        "n_collection_year_only",
        "n_collection_missing_or_not_applicable",
        "overlap_run_count",
        "overlap_sample_count",
        "overlap_project_count",
        "overlap_any_count",
        "earliest_public_lead_to_focus_threshold_days",
        "earliest_public_lead_to_focus_peak_days",
        "recommended_tier",
        "decision_note",
        "metadata_probe_url",
    ]
    with SUMMARY_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=summary_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"wrote {RUN_OUTPUT}")
    print(f"wrote {SUMMARY_OUTPUT}")


if __name__ == "__main__":
    main()
