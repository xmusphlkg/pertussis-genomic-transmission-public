#!/usr/bin/env python3
"""Build interval-aware genome availability tables for the EID Dispatch.

This script is intentionally downstream of the frozen JOI analysis. It does not
change lineage definitions, trees, model inputs, or fitted quantities. It adds
an external-visibility layer: the earliest date on which an ENA read record or
an NCBI assembly could be reproducibly retrieved for each included genome.
"""

from __future__ import annotations

import csv
import datetime as dt
import os
import statistics
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
DATA_DERIVED = REPO / "data" / "derived"
RESULTS = REPO / "results" / "public_availability"
ENA_CACHE = DATA_DERIVED / "public_availability_ena_runs.tsv"
ENA_PROJECT_AUDIT = DATA_DERIVED / "public_availability_ena_project_audit.tsv"
PUBLIC_TABLE = DATA_DERIVED / "public_genome_availability.tsv"
SUMMARY_TABLE = RESULTS / "public_availability_summary.tsv"
CUMULATIVE_TABLE = RESULTS / "cumulative_detection_dates.tsv"
CASE_THRESHOLDS_TABLE = RESULTS / "case_thresholds.tsv"
NAS_FOCAL_SNAPSHOT = REPO / "data" / "source_snapshots" / "eid_nas_focal_genome_metadata.tsv"

MODEL_COUNTRIES = {"AUS", "CHN", "JPN"}
EID_FOCUS_LINEAGE = "L1_02.07"
DISPLAY_LINEAGE = {
    "L1_01.02": "Comparator genomic lineage A",
    "L1_02.05": "Comparator genomic lineage B",
    "L1_02.06": "Comparator genomic lineage C",
    "L1_02.07": "MT28-associated genomic lineage",
    "Other": "Other genomic lineages",
}


def parse_date(value: str | None) -> dt.date | None:
    value = (value or "").strip()
    if not value or value in {"missing", "not applicable", "NA", "NaN"}:
        return None
    try:
        if len(value) == 4:
            return dt.date(int(value), 1, 1)
        if len(value) == 7:
            year, month = map(int, value.split("-"))
            return dt.date(year, month, 1)
        return dt.date.fromisoformat(value[:10])
    except ValueError:
        return None


def date_text(value: dt.date | None) -> str:
    return value.isoformat() if value else ""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def fetch_ena_runs(projects: list[str]) -> list[dict[str, str]]:
    projects = sorted(set(project for project in projects if project))
    fieldnames = [
        "query_project_id",
        "run_accession",
        "sample_accession",
        "study_accession",
        "collection_date",
        "country",
        "first_public",
    ]

    cached_rows = read_tsv(ENA_CACHE) if ENA_CACHE.exists() else []
    cached_by_project: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cached_rows:
        cached_by_project[row.get("query_project_id", "")].append(row)

    audit_fields = ["query_project_id", "query_status", "n_rows", "last_successful_query_date"]
    audit_rows = read_tsv(ENA_PROJECT_AUDIT) if ENA_PROJECT_AUDIT.exists() else []
    audit_by_project = {row.get("query_project_id", ""): row for row in audit_rows}
    for project, rows_for_project in cached_by_project.items():
        if project and project not in audit_by_project:
            audit_by_project[project] = {
                "query_project_id": project,
                "query_status": "records_cached_legacy",
                "n_rows": len(rows_for_project),
                "last_successful_query_date": "",
            }
    completed_projects = set(cached_by_project) | {
        project
        for project, row in audit_by_project.items()
        if row.get("query_status") == "queried_no_read_run_rows"
    }

    if os.environ.get("GTD40_OFFLINE_CACHE") == "1":
        missing = sorted(set(projects) - completed_projects)
        if missing:
            print(
                "warning: using offline ENA cache; missing projects: " + ",".join(missing),
                file=sys.stderr,
            )
        return cached_rows

    if os.environ.get("GTD40_FORCE_REFRESH") == "1":
        fetch_projects = projects
        rows: list[dict[str, str]] = []
        audit_by_project = {}
    else:
        fetch_projects = [project for project in projects if project not in completed_projects]
        if not fetch_projects:
            write_tsv(ENA_PROJECT_AUDIT, list(sorted(audit_by_project.values(), key=lambda row: row["query_project_id"])), audit_fields)
            return cached_rows
        fetch_project_set = set(fetch_projects)
        rows = [
            row
            for row in cached_rows
            if row.get("query_project_id", "") in set(projects) - fetch_project_set
        ]

    base = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    fields = "run_accession,sample_accession,study_accession,collection_date,country,first_public"
    for project in projects:
        if project not in fetch_projects:
            continue
        query = urllib.parse.urlencode(
            {
                "accession": project,
                "result": "read_run",
                "fields": fields,
                "format": "tsv",
            }
        )
        url = f"{base}?{query}"
        text = ""
        request_succeeded = False
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                with urllib.request.urlopen(request, timeout=45) as response:
                    text = response.read().decode("utf-8")
                request_succeeded = True
                break
            except Exception as exc:  # pragma: no cover - network branch
                if attempt == 2:
                    print(f"warning: ENA fetch failed for {project}: {exc}", file=sys.stderr)
                    rows.extend(cached_by_project.get(project, []))
                time.sleep(1.5 * (attempt + 1))
        if not request_succeeded:
            continue
        rows = [row for row in rows if row.get("query_project_id") != project]
        response_rows = list(csv.DictReader(text.splitlines(), delimiter="\t")) if text.strip() else []
        for row in response_rows:
            row["query_project_id"] = project
            rows.append(row)
        audit_by_project[project] = {
            "query_project_id": project,
            "query_status": "records_cached" if response_rows else "queried_no_read_run_rows",
            "n_rows": len(response_rows),
            "last_successful_query_date": dt.date.today().isoformat(),
        }
        write_tsv(ENA_CACHE, rows, fieldnames)
        write_tsv(
            ENA_PROJECT_AUDIT,
            list(sorted(audit_by_project.values(), key=lambda row: row["query_project_id"])),
            audit_fields,
        )

    write_tsv(ENA_CACHE, rows, fieldnames)
    write_tsv(
        ENA_PROJECT_AUDIT,
        list(sorted(audit_by_project.values(), key=lambda row: row["query_project_id"])),
        audit_fields,
    )
    return rows


def focal_project_ids() -> set[str]:
    projects: set[str] = set()
    for row in read_tsv(REPO / "results" / "phylogeny" / "tree_tip_metadata.tsv"):
        if row.get("tree_role") != "focal":
            continue
        project = (row.get("project_id") or "").strip()
        if project:
            projects.add(project)
    return projects


def load_nas_snapshot() -> tuple[dict[str, dict[str, str]], dict[str, dt.date]]:
    if not NAS_FOCAL_SNAPSHOT.exists():
        raise SystemExit(
            f"Missing {NAS_FOCAL_SNAPSHOT}. Generate it with "
            "python3 scripts/pipeline/gtd_39_snapshot_eid_nas_inputs.py while the NAS is mounted."
        )
    rows = read_tsv(NAS_FOCAL_SNAPSHOT)
    by_tree_id = {row["tree_sample_id"]: row for row in rows}
    release_by_accession: dict[str, dt.date] = {}
    for row in rows:
        release_date = parse_date(row.get("assembly_release_date"))
        if not release_date:
            continue
        for column in ("assembly_accession",):
            accession = (row.get(column) or "").strip()
            if not accession:
                continue
            release_by_accession[accession] = release_date
            release_by_accession[accession.split(".")[0]] = release_date
    return by_tree_id, release_by_accession


def normalize_location(country_iso3: str, raw_value: str) -> tuple[str, str]:
    value = " ".join((raw_value or "").strip().split())
    if not value or value.lower() in {"missing", "na", "not applicable"}:
        return "", "missing"
    if ":" not in value:
        return "", "country"
    _, suffix = value.split(":", 1)
    suffix = " ".join(suffix.strip().split())
    if not suffix:
        return "", "country"
    return suffix, "subnational"


def public_route(ena_date: dt.date | None, assembly_date: dt.date | None) -> str:
    if ena_date and assembly_date:
        if ena_date == assembly_date:
            return "ena_read_and_ncbi_assembly_same_day"
        return "ena_read_first" if ena_date < assembly_date else "ncbi_assembly_first"
    if ena_date:
        return "ena_read_only"
    if assembly_date:
        return "ncbi_assembly_only"
    return "not_publicly_dated"


def lineage_map() -> dict[str, str]:
    rows = read_tsv(REPO / "results" / "lineages" / "primary_finalized" / "model_lineage_assignments.tsv")
    mapping: dict[str, str] = {}
    for row in rows:
        lineage = row.get("primary_model_lineage_id") or row.get("model_lineage_id") or ""
        mapping[row["tree_sample_id"]] = lineage
    return mapping


def case_thresholds() -> dict[str, dict[str, object]]:
    rows = read_tsv(DATA_DERIVED / "country_month_cases.tsv")
    by_country: dict[str, list[tuple[dt.date, float]]] = defaultdict(list)
    for row in rows:
        if row.get("case_data_available") != "True":
            continue
        country = row.get("country_iso3", "")
        month = parse_date(row.get("model_month"))
        if not country or not month:
            continue
        by_country[country].append((month, float(row["cases"])))

    thresholds: dict[str, dict[str, object]] = {}
    output_rows: list[dict[str, object]] = []
    for country, values in sorted(by_country.items()):
        values = sorted(values)
        year_2019 = [case for month, case in values if month.year == 2019]
        if not year_2019:
            continue
        max_2019 = max(year_2019)
        median_2019 = statistics.median(year_2019)
        post = [(month, case) for month, case in values if month >= dt.date(2023, 1, 1)]
        first_above_max = next(((month, case) for month, case in post if case > max_2019), (None, None))
        first_above_median = next(((month, case) for month, case in post if case > median_2019), (None, None))
        first_two_above_max: tuple[dt.date | None, float | None] = (None, None)
        for index in range(len(post) - 1):
            month, case = post[index]
            next_month, next_case = post[index + 1]
            expected_next = dt.date(month.year + (month.month == 12), 1 if month.month == 12 else month.month + 1, 1)
            if next_month == expected_next and case > max_2019 and next_case > max_2019:
                first_two_above_max = (month, case)
                break
        peak_month, peak_cases = max(post, key=lambda item: item[1]) if post else (None, None)
        record = {
            "country_iso3": country,
            "cases_2019_max": max_2019,
            "cases_2019_median": median_2019,
            "first_post2022_month_above_2019_max": first_above_max[0],
            "first_post2022_cases_above_2019_max": first_above_max[1],
            "first_post2022_month_above_2019_median": first_above_median[0],
            "first_post2022_cases_above_2019_median": first_above_median[1],
            "first_two_consecutive_months_above_2019_max": first_two_above_max[0],
            "first_two_consecutive_months_first_cases": first_two_above_max[1],
            "post2022_peak_month": peak_month,
            "post2022_peak_cases": peak_cases,
        }
        thresholds[country] = record
        output_rows.append({key: date_text(value) if isinstance(value, dt.date) else value for key, value in record.items()})

    write_tsv(
        CASE_THRESHOLDS_TABLE,
        output_rows,
        [
            "country_iso3",
            "cases_2019_max",
            "cases_2019_median",
            "first_post2022_month_above_2019_max",
            "first_post2022_cases_above_2019_max",
            "first_post2022_month_above_2019_median",
            "first_post2022_cases_above_2019_median",
            "first_two_consecutive_months_above_2019_max",
            "first_two_consecutive_months_first_cases",
            "post2022_peak_month",
            "post2022_peak_cases",
        ],
    )
    return thresholds


def public_availability_rows() -> list[dict[str, object]]:
    project_rows = read_tsv(REPO / "data" / "model_inputs" / "project_index.tsv")
    projects = sorted({row["project_id"] for row in project_rows if row.get("project_id")} | focal_project_ids())
    ena_rows = fetch_ena_runs(projects)
    nas_by_tree_id, assembly_release = load_nas_snapshot()
    lineages = lineage_map()
    thresholds = case_thresholds()

    first_by_run: dict[str, dt.date] = {}
    first_by_sample: dict[str, dt.date] = {}
    ena_row_by_run: dict[str, dict[str, str]] = {}
    ena_row_by_sample: dict[str, dict[str, str]] = {}
    for row in ena_rows:
        first_public = parse_date(row.get("first_public"))
        if not first_public:
            continue
        run = (row.get("run_accession") or "").strip()
        sample = (row.get("sample_accession") or "").strip()
        if run:
            first_by_run[run] = min(first_by_run.get(run, first_public), first_public)
            if run not in ena_row_by_run or first_public <= parse_date(ena_row_by_run[run].get("first_public")):
                ena_row_by_run[run] = row
        if sample:
            first_by_sample[sample] = min(first_by_sample.get(sample, first_public), first_public)
            if sample not in ena_row_by_sample or first_public <= parse_date(ena_row_by_sample[sample].get("first_public")):
                ena_row_by_sample[sample] = row

    public_rows: list[dict[str, object]] = []
    for row in read_tsv(REPO / "results" / "phylogeny" / "tree_tip_metadata.tsv"):
        if row.get("tree_role") != "focal":
            continue

        collection_lower = parse_date(row.get("date_lower"))
        collection_upper = parse_date(row.get("date_upper"))
        nas_metadata = nas_by_tree_id.get(row["tree_sample_id"], {})
        lineage = lineages.get(row["tree_sample_id"], "")
        candidates: list[tuple[str, dt.date]] = []

        run = (row.get("run_accession") or "").strip()
        ena_metadata = ena_row_by_run.get(run, {})
        if run and run in first_by_run:
            candidates.append(("ena_read_first_public_run", first_by_run[run]))

        for sample_column in ("sample_accession", "biosample_accession"):
            sample = (row.get(sample_column) or "").strip()
            if sample and sample in first_by_sample:
                candidates.append((f"ena_read_first_public_via_{sample_column}", first_by_sample[sample]))
                if not ena_metadata:
                    ena_metadata = ena_row_by_sample.get(sample, {})

        assembly = (row.get("assembly_accession") or "").strip()
        for accession in (assembly, assembly.split(".")[0] if assembly else ""):
            if accession and accession in assembly_release:
                candidates.append(("ncbi_assembly_release", assembly_release[accession]))

        ena_first_public = min(
            (date for source, date in candidates if source.startswith("ena_read_")),
            default=None,
        )
        assembly_release_date = min(
            (date for source, date in candidates if source == "ncbi_assembly_release"),
            default=None,
        )
        public_date = min((date for _, date in candidates), default=None)
        public_sources = sorted(source for source, date in candidates if date == public_date)
        temporal_status = "ok"
        collection_upper_effective = collection_upper
        if public_date and collection_lower and public_date < collection_lower:
            temporal_status = "public_before_collection_lower"
            collection_upper_effective = None
        elif public_date and collection_upper:
            collection_upper_effective = min(collection_upper, public_date)

        ena_geo = (ena_metadata.get("country") or "").strip()
        nas_geo = (nas_metadata.get("geo_raw") or "").strip()
        ena_location, ena_resolution = normalize_location(row.get("country_iso3", ""), ena_geo)
        nas_location, nas_resolution = normalize_location(row.get("country_iso3", ""), nas_geo)
        if ena_resolution == "subnational":
            subnational_location, location_resolution, location_source = ena_location, ena_resolution, "ena_country"
        elif nas_resolution == "subnational":
            subnational_location, location_resolution, location_source = nas_location, nas_resolution, "nas_focal_snapshot"
        elif ena_geo:
            subnational_location, location_resolution, location_source = "", ena_resolution, "ena_country"
        elif nas_geo:
            subnational_location, location_resolution, location_source = "", nas_resolution, "nas_focal_snapshot"
        else:
            subnational_location, location_resolution, location_source = "", "missing", ""
        threshold = thresholds.get(row.get("country_iso3", ""), {})
        threshold_max_date = threshold.get("first_post2022_month_above_2019_max")
        peak_date = threshold.get("post2022_peak_month")

        public_rows.append(
            {
                "tree_sample_id": row["tree_sample_id"],
                "genome_record_id": row.get("genome_record_id", ""),
                "country_iso3": row.get("country_iso3", ""),
                "country_name": row.get("country_name", ""),
                "collection_lower": date_text(collection_lower),
                "collection_upper": date_text(collection_upper),
                "collection_upper_effective": date_text(collection_upper_effective),
                "temporal_consistency_status": temporal_status,
                "date_resolution": row.get("date_resolution", ""),
                "year": row.get("year", ""),
                "month": row.get("month", ""),
                "project_id": row.get("project_id", ""),
                "run_accession": run,
                "biosample_accession": row.get("biosample_accession", ""),
                "assembly_accession": assembly,
                "sequence_acquisition": row.get("sequence_acquisition", ""),
                "primary_model_lineage_id": lineage,
                "display_lineage": DISPLAY_LINEAGE.get(lineage, lineage),
                "epidemic_period": row.get("epidemic_period", ""),
                "public_date": date_text(public_date),
                "public_date_source": ";".join(public_sources),
                "sequence_public_date": date_text(public_date),
                "sequence_public_date_source": ";".join(public_sources),
                "ena_first_public_date": date_text(ena_first_public),
                "assembly_release_date": date_text(assembly_release_date),
                "public_route": public_route(ena_first_public, assembly_release_date),
                "subnational_location": subnational_location,
                "location_source": location_source,
                "location_resolution": location_resolution,
                "nas_snapshot_record_id": nas_metadata.get("nas_snapshot_record_id", ""),
                "lag_from_collection_lower_days": (public_date - collection_lower).days if public_date and collection_lower else "",
                "lag_from_collection_upper_days": (public_date - collection_upper).days if public_date and collection_upper else "",
                "lag_min_days": (public_date - collection_upper_effective).days
                if public_date and collection_upper_effective
                else "",
                "lag_max_days": (public_date - collection_lower).days if public_date and collection_lower else "",
                "first_post2022_month_above_2019_max": date_text(threshold_max_date if isinstance(threshold_max_date, dt.date) else None),
                "post2022_peak_month": date_text(peak_date if isinstance(peak_date, dt.date) else None),
                "public_lead_to_2019max_threshold_days": (threshold_max_date - public_date).days
                if public_date and isinstance(threshold_max_date, dt.date)
                else "",
                "public_lead_to_post2022_peak_days": (peak_date - public_date).days
                if public_date and isinstance(peak_date, dt.date)
                else "",
                "collection_lead_to_2019max_threshold_days": (threshold_max_date - collection_lower).days
                if collection_lower and isinstance(threshold_max_date, dt.date)
                else "",
                "collection_lead_to_post2022_peak_days": (peak_date - collection_lower).days
                if collection_lower and isinstance(peak_date, dt.date)
                else "",
            }
        )
    return public_rows


def summarise(public_rows: list[dict[str, object]]) -> None:
    focus = [
        row
        for row in public_rows
        if row["country_iso3"] in MODEL_COUNTRIES
        and row["collection_lower"] >= "2023-01-01"
        and row["primary_model_lineage_id"]
    ]

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in focus:
        grouped[(str(row["country_iso3"]), str(row["primary_model_lineage_id"]))].append(row)

    summary_rows: list[dict[str, object]] = []
    for (country, lineage), rows in sorted(grouped.items()):
        public_subset = [row for row in rows if row["public_date"]]
        lag_mins = sorted(int(row["lag_min_days"]) for row in public_subset if row["lag_min_days"] != "")
        lag_maxs = sorted(int(row["lag_max_days"]) for row in public_subset if row["lag_max_days"] != "")
        first_public = min((str(row["public_date"]) for row in public_subset), default="")
        first_collection = min(str(row["collection_lower"]) for row in rows)
        first_public_rows = [row for row in public_subset if row["public_date"] == first_public]
        lead_to_threshold = first_public_rows[0]["public_lead_to_2019max_threshold_days"] if first_public_rows else ""
        lead_to_peak = first_public_rows[0]["public_lead_to_post2022_peak_days"] if first_public_rows else ""
        summary_rows.append(
            {
                "country_iso3": country,
                "primary_model_lineage_id": lineage,
                "display_lineage": DISPLAY_LINEAGE.get(lineage, lineage),
                "n_resurgence_focal_tips": len(rows),
                "n_with_public_date": len(public_subset),
                "first_collection_lower": first_collection,
                "first_public_date": first_public,
                "median_lag_min_days": statistics.median(lag_mins) if lag_mins else "",
                "median_lag_max_days": statistics.median(lag_maxs) if lag_maxs else "",
                "min_lag_min_days": min(lag_mins) if lag_mins else "",
                "max_lag_max_days": max(lag_maxs) if lag_maxs else "",
                "median_lag_from_collection_lower_days": statistics.median(lag_maxs) if lag_maxs else "",
                "max_lag_from_collection_lower_days": max(lag_maxs) if lag_maxs else "",
                "first_public_lead_to_2019max_threshold_days": lead_to_threshold,
                "first_public_lead_to_post2022_peak_days": lead_to_peak,
            }
        )

    write_tsv(
        SUMMARY_TABLE,
        summary_rows,
        [
            "country_iso3",
            "primary_model_lineage_id",
            "display_lineage",
            "n_resurgence_focal_tips",
            "n_with_public_date",
            "first_collection_lower",
            "first_public_date",
            "median_lag_min_days",
            "median_lag_max_days",
            "min_lag_min_days",
            "max_lag_max_days",
            "median_lag_from_collection_lower_days",
            "max_lag_from_collection_lower_days",
            "first_public_lead_to_2019max_threshold_days",
            "first_public_lead_to_post2022_peak_days",
        ],
    )

    cumulative_rows: list[dict[str, object]] = []
    thresholds = case_thresholds()
    for (country, lineage), rows in sorted(grouped.items()):
        for clock, date_column in (
            ("collection_lower", "collection_lower"),
            ("collection_upper_effective", "collection_upper_effective"),
            ("public_date", "public_date"),
        ):
            dated = [row for row in rows if row.get(date_column)]
            dated.sort(key=lambda item: str(item[date_column]))
            for minimum_tips in (1, 3, 5, 10, 20, 50):
                if len(dated) < minimum_tips:
                    continue
                detection_date = parse_date(str(dated[minimum_tips - 1][date_column]))
                threshold_date = thresholds.get(country, {}).get("first_post2022_month_above_2019_max")
                peak_date = thresholds.get(country, {}).get("post2022_peak_month")
                cumulative_rows.append(
                    {
                        "country_iso3": country,
                        "primary_model_lineage_id": lineage,
                        "display_lineage": DISPLAY_LINEAGE.get(lineage, lineage),
                        "clock": clock,
                        "minimum_cumulative_tips": minimum_tips,
                        "detection_date": date_text(detection_date),
                        "lead_to_2019max_threshold_days": (threshold_date - detection_date).days
                        if isinstance(threshold_date, dt.date) and detection_date
                        else "",
                        "lead_to_post2022_peak_days": (peak_date - detection_date).days
                        if isinstance(peak_date, dt.date) and detection_date
                        else "",
                    }
                )

    write_tsv(
        CUMULATIVE_TABLE,
        cumulative_rows,
        [
            "country_iso3",
            "primary_model_lineage_id",
            "display_lineage",
            "clock",
            "minimum_cumulative_tips",
            "detection_date",
            "lead_to_2019max_threshold_days",
            "lead_to_post2022_peak_days",
        ],
    )


def main() -> None:
    public_rows = public_availability_rows()
    write_tsv(
        PUBLIC_TABLE,
        public_rows,
        [
            "tree_sample_id",
            "genome_record_id",
            "country_iso3",
            "country_name",
            "collection_lower",
            "collection_upper",
            "collection_upper_effective",
            "temporal_consistency_status",
            "date_resolution",
            "year",
            "month",
            "project_id",
            "run_accession",
            "biosample_accession",
            "assembly_accession",
            "sequence_acquisition",
            "primary_model_lineage_id",
            "display_lineage",
            "epidemic_period",
            "public_date",
            "public_date_source",
            "sequence_public_date",
            "sequence_public_date_source",
            "ena_first_public_date",
            "assembly_release_date",
            "public_route",
            "subnational_location",
            "location_source",
            "location_resolution",
            "nas_snapshot_record_id",
            "lag_from_collection_lower_days",
            "lag_from_collection_upper_days",
            "lag_min_days",
            "lag_max_days",
            "first_post2022_month_above_2019_max",
            "post2022_peak_month",
            "public_lead_to_2019max_threshold_days",
            "public_lead_to_post2022_peak_days",
            "collection_lead_to_2019max_threshold_days",
            "collection_lead_to_post2022_peak_days",
        ],
    )
    summarise(public_rows)
    print(f"wrote {PUBLIC_TABLE}")
    print(f"wrote {SUMMARY_TABLE}")
    print(f"wrote {CUMULATIVE_TABLE}")
    print(f"wrote {CASE_THRESHOLDS_TABLE}")


if __name__ == "__main__":
    main()
