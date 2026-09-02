#!/usr/bin/env python3
"""Build EID Dispatch tables for the lineage-specific visibility audit.

The script is downstream of the frozen 989-genome tree and lineage assignment.
It reconstructs cumulative specimen accumulation and public sequence
availability, uses selected k values as interpretive anchors, and compares the
target lineage with other lineages within country, BioProject, and collection
year. Collection dates remain interval-censored.
"""

from __future__ import annotations

import csv
import datetime as dt
import calendar
import statistics
from collections import Counter, defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results" / "public_availability"
PUBLIC_TABLE = REPO / "data" / "derived" / "public_genome_availability.tsv"
CUMULATIVE_TABLE = RESULTS / "cumulative_detection_dates.tsv"
SUMMARY_TABLE = RESULTS / "public_availability_summary.tsv"
CANDIDATE_TABLE = RESULTS / "candidate_project_metadata_audit.tsv"
CASE_THRESHOLDS = RESULTS / "case_thresholds.tsv"
NATIVE_CASES = REPO / "data" / "derived" / "figure1a_native_resolution_surveillance.tsv"
SNAPSHOT_MANIFEST = REPO / "provenance" / "EID_NAS_SNAPSHOT_MANIFEST.tsv"
HIGHRES_CASE_SNAPSHOT = REPO / "data" / "source_snapshots" / "eid_nas_highres_cases.tsv"

SHIFT_OUTPUT = RESULTS / "eid_detection_clock_shift.tsv"
MILESTONE_OUTPUT = RESULTS / "eid_milestone_visibility.tsv"
LAG_SUMMARY_OUTPUT = RESULTS / "eid_country_lineage_lag_summary.tsv"
CANDIDATE_OUTPUT = RESULTS / "eid_external_candidate_summary.tsv"
THRESHOLD_SENSITIVITY_OUTPUT = RESULTS / "eid_threshold_sensitivity.tsv"
PROJECT_BATCH_OUTPUT = RESULTS / "eid_project_batch_release.tsv"
GEOGRAPHY_OUTPUT = RESULTS / "eid_geography_audit.tsv"
CASE_SENSITIVITY_OUTPUT = RESULTS / "eid_case_clock_sensitivity.tsv"
CUMULATIVE_VISIBILITY_OUTPUT = RESULTS / "eid_cumulative_visibility.tsv"
PROJECT_LINEAGE_COMPARISON_OUTPUT = RESULTS / "eid_project_lineage_comparison.tsv"
APPENDIX_MD = REPO / "manuscript" / "eid_dispatch_appendix.md"

TARGET_COUNTRIES = ("AUS", "CHN", "JPN")
COUNTRY_LABELS = {"AUS": "Australia", "CHN": "China", "JPN": "Japan"}
TARGET_LINEAGE = "L1_02.07"
SELECTED_K = {"AUS": 3, "CHN": 5, "JPN": 5}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def parse_date(value: str | None) -> dt.date | None:
    text = (value or "").strip()
    if not text:
        return None
    return dt.date.fromisoformat(text[:10])


def date_text(value: dt.date | None) -> str:
    return value.isoformat() if value else ""


def month_end(value: dt.date | None) -> dt.date | None:
    if not value:
        return None
    return dt.date(value.year, value.month, calendar.monthrange(value.year, value.month)[1])


def date_range_text(start: dt.date | None, end: dt.date | None) -> str:
    if not start:
        return ""
    if not end or start == end:
        return date_text(start)
    return f"{date_text(start)} to {date_text(end)}"


def number_text(value: int | float | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def pct_text(numerator: int, denominator: int) -> str:
    return "NA" if denominator == 0 else f"{100 * numerator / denominator:.1f}"


def counter_text(values: list[str]) -> str:
    counts = Counter(value for value in values if value)
    return ";".join(f"{key}={counts[key]}" for key in sorted(counts))


def relative_range_text(min_value: object, max_value: object) -> str:
    if min_value in {None, ""} or max_value in {None, ""}:
        return ""
    if min_value == max_value:
        return str(min_value)
    return f"{min_value}–{max_value}"


def target_resurgence_rows() -> list[dict[str, str]]:
    return [
        row
        for row in read_tsv(PUBLIC_TABLE)
        if row.get("country_iso3") in TARGET_COUNTRIES
        and row.get("primary_model_lineage_id") == TARGET_LINEAGE
        and row.get("epidemic_period") == "resurgence"
    ]


def sequence_public_date(row: dict[str, str]) -> dt.date | None:
    return parse_date(row.get("sequence_public_date") or row.get("public_date"))


def build_threshold_sensitivity() -> list[dict[str, object]]:
    rows = [
        row
        for row in read_tsv(CUMULATIVE_TABLE)
        if row.get("country_iso3") in TARGET_COUNTRIES
        and row.get("primary_model_lineage_id") == TARGET_LINEAGE
    ]
    milestone_starts = {
        row["country_iso3"]: parse_date(row.get("first_post2022_month_above_2019_max"))
        for row in read_tsv(CASE_THRESHOLDS)
        if row.get("country_iso3") in TARGET_COUNTRIES
    }
    lookup = {
        (row["country_iso3"], row["clock"], int(row["minimum_cumulative_tips"])): row
        for row in rows
    }
    output_rows: list[dict[str, object]] = []
    for country in TARGET_COUNTRIES:
        for threshold in (1, 3, 5, 10, 20, 50):
            lower = lookup.get((country, "collection_lower", threshold))
            upper = lookup.get((country, "collection_upper_effective", threshold))
            public = lookup.get((country, "public_date", threshold))
            if not lower or not upper or not public:
                continue
            lower_date = parse_date(lower["detection_date"])
            upper_date = parse_date(upper["detection_date"])
            public_date = parse_date(public["detection_date"])
            milestone_start = milestone_starts.get(country)
            milestone_end = month_end(milestone_start)
            shift_min = (public_date - upper_date).days if public_date and upper_date else None
            shift_max = (public_date - lower_date).days if public_date and lower_date else None
            if milestone_start and milestone_end and lower_date and upper_date:
                collection_relative_min = (lower_date - milestone_end).days
                collection_relative_max = (upper_date - milestone_start).days
                lead_min = -collection_relative_max
                lead_max = -collection_relative_min
            else:
                lead_min = int(upper["lead_to_2019max_threshold_days"])
                lead_max = int(lower["lead_to_2019max_threshold_days"])
                collection_relative_min = -lead_max
                collection_relative_max = -lead_min
            if milestone_start and milestone_end and public_date:
                public_relative_min = (public_date - milestone_end).days
                public_relative_max = (public_date - milestone_start).days
                public_lead_min = -public_relative_max
                public_lead_max = -public_relative_min
            else:
                public_lead = int(public["lead_to_2019max_threshold_days"])
                public_relative_min = -public_lead
                public_relative_max = -public_lead
                public_lead_min = public_lead
                public_lead_max = public_lead
            if lead_min > 0:
                collection_class = "interval_before_case_threshold"
            elif lead_max < 0:
                collection_class = "interval_after_case_threshold"
            else:
                collection_class = "interval_spans_case_threshold"
            public_class = (
                "before_case_threshold"
                if public_relative_max < 0
                else "on_or_after_case_threshold"
            )
            output_rows.append(
                {
                    "country_iso3": country,
                    "country_label": COUNTRY_LABELS[country],
                    "target_lineage": TARGET_LINEAGE,
                    "cumulative_genome_threshold": threshold,
                    "collection_detection_lower": date_text(lower_date),
                    "collection_detection_upper": date_text(upper_date),
                    "public_detection_date": date_text(public_date),
                    "resurgence_milestone_start": date_text(milestone_start),
                    "resurgence_milestone_end": date_text(milestone_end),
                    "clock_shift_min_days": number_text(shift_min),
                    "clock_shift_max_days": number_text(shift_max),
                    "collection_lead_to_case_threshold_min_days": lead_min,
                    "collection_lead_to_case_threshold_max_days": lead_max,
                    "public_lead_to_case_threshold_min_days": public_lead_min,
                    "public_lead_to_case_threshold_max_days": public_lead_max,
                    "public_lead_to_case_threshold_days": public_lead_max,
                    "collection_relative_to_resurgence_min_days": collection_relative_min,
                    "collection_relative_to_resurgence_max_days": collection_relative_max,
                    "public_relative_to_resurgence_min_days": public_relative_min,
                    "public_relative_to_resurgence_max_days": public_relative_max,
                    "public_relative_to_resurgence_days": public_relative_max,
                    "collection_timing_class": collection_class,
                    "public_timing_class": public_class,
                }
            )
    return output_rows


def build_detection_shift(threshold_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in threshold_rows
        if int(row["cumulative_genome_threshold"]) == SELECTED_K[str(row["country_iso3"])]
    ]


def build_milestone_visibility() -> list[dict[str, object]]:
    rows = target_resurgence_rows()
    output_rows: list[dict[str, object]] = []
    for country in TARGET_COUNTRIES:
        country_rows = [row for row in rows if row["country_iso3"] == country]
        if not country_rows:
            continue
        milestones = {
            "first_post2022_month_above_2019_max": parse_date(
                country_rows[0]["first_post2022_month_above_2019_max"]
            ),
            "post2022_peak_month": parse_date(country_rows[0]["post2022_peak_month"]),
        }
        for milestone, milestone_date in milestones.items():
            if not milestone_date:
                continue
            milestone_period_end = month_end(milestone_date)
            possible = sum(
                (parse_date(row.get("collection_lower")) or dt.date.max) <= milestone_date
                for row in country_rows
            )
            definite = sum(
                (parse_date(row.get("collection_upper_effective")) or dt.date.max) <= milestone_date
                for row in country_rows
            )
            public = sum(
                bool(row.get("public_date"))
                and (parse_date(row.get("public_date")) or dt.date.max) <= milestone_date
                for row in country_rows
            )
            output_rows.append(
                {
                    "country_iso3": country,
                    "country_label": COUNTRY_LABELS[country],
                    "target_lineage": TARGET_LINEAGE,
                    "milestone": milestone,
                    "milestone_date": date_text(milestone_date),
                    "milestone_period_start": date_text(milestone_date),
                    "milestone_period_end": date_text(milestone_period_end),
                    "milestone_period": date_range_text(milestone_date, milestone_period_end),
                    "total_resurgence_target_genomes": len(country_rows),
                    "definitely_collected_by_milestone": definite,
                    "possibly_collected_by_milestone": possible,
                    "public_by_milestone": public,
                    "definitely_collected_not_public": max(definite - public, 0),
                    "possibly_collected_not_public": max(possible - public, 0),
                    "public_among_possible_collected_pct": pct_text(public, possible),
                }
            )
    return output_rows


def build_cumulative_visibility() -> list[dict[str, object]]:
    rows = target_resurgence_rows()
    output_rows: list[dict[str, object]] = []
    for country in TARGET_COUNTRIES:
        country_rows = [row for row in rows if row["country_iso3"] == country]
        if not country_rows:
            continue
        milestone_date = parse_date(country_rows[0].get("first_post2022_month_above_2019_max"))
        peak_date = parse_date(country_rows[0].get("post2022_peak_month"))
        milestone_end = month_end(milestone_date)
        peak_end = month_end(peak_date)
        event_dates: set[dt.date] = {dt.date(2023, 1, 1)}
        for row in country_rows:
            for value in (
                parse_date(row.get("collection_lower")),
                parse_date(row.get("collection_upper_effective")),
                sequence_public_date(row),
            ):
                if value:
                    event_dates.add(value)
        for value in (milestone_date, peak_date):
            if value:
                event_dates.add(value)
        for event_date in sorted(event_dates):
            output_rows.append(
                {
                    "country_iso3": country,
                    "country_label": COUNTRY_LABELS[country],
                    "target_lineage": TARGET_LINEAGE,
                    "event_date": date_text(event_date),
                    "n_possibly_collected": sum(
                        (parse_date(row.get("collection_lower")) or dt.date.max) <= event_date
                        for row in country_rows
                    ),
                    "n_definitely_collected": sum(
                        (parse_date(row.get("collection_upper_effective")) or dt.date.max) <= event_date
                        for row in country_rows
                    ),
                    "n_publicly_available": sum(
                        bool(sequence_public_date(row))
                        and (sequence_public_date(row) or dt.date.max) <= event_date
                        for row in country_rows
                    ),
                    "total_target_lineage_records": len(country_rows),
                    "resurgence_milestone_date": date_text(milestone_date),
                    "resurgence_milestone_end": date_text(milestone_end),
                    "post2022_peak_date": date_text(peak_date),
                    "post2022_peak_end": date_text(peak_end),
                }
            )
    return output_rows


def build_project_lineage_comparison() -> list[dict[str, object]]:
    rows = [
        row
        for row in read_tsv(PUBLIC_TABLE)
        if row.get("country_iso3") in TARGET_COUNTRIES
        and row.get("epidemic_period") == "resurgence"
        and row.get("project_id")
        and row.get("year")
    ]
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["country_iso3"], row["project_id"], row["year"])].append(row)

    def lag_values(group: list[dict[str, str]], field: str) -> list[int]:
        return [int(row[field]) for row in group if row.get(field) not in {None, ""}]

    def modal_public(group: list[dict[str, str]]) -> tuple[str, int]:
        counts = Counter(
            date_text(sequence_public_date(row))
            for row in group
            if sequence_public_date(row)
        )
        if not counts:
            return "", 0
        modal_n = max(counts.values())
        return min(date for date, count in counts.items() if count == modal_n), modal_n

    output_rows: list[dict[str, object]] = []
    for (country, project, year), group in sorted(grouped.items()):
        target = [row for row in group if row.get("primary_model_lineage_id") == TARGET_LINEAGE]
        comparator = [row for row in group if row.get("primary_model_lineage_id") != TARGET_LINEAGE]
        if not target or not comparator:
            continue
        target_min = lag_values(target, "lag_min_days")
        target_max = lag_values(target, "lag_max_days")
        comparator_min = lag_values(comparator, "lag_min_days")
        comparator_max = lag_values(comparator, "lag_max_days")
        if not all((target_min, target_max, comparator_min, comparator_max)):
            continue
        target_median_min = statistics.median(target_min)
        target_median_max = statistics.median(target_max)
        comparator_median_min = statistics.median(comparator_min)
        comparator_median_max = statistics.median(comparator_max)
        if target_median_max < comparator_median_min:
            direction = "target_shorter"
        elif target_median_min > comparator_median_max:
            direction = "target_longer"
        else:
            direction = "intervals_overlap"
        project_modal_date, project_modal_n = modal_public(group)
        target_modal_date, _ = modal_public(target)
        comparator_modal_date, _ = modal_public(comparator)
        output_rows.append(
            {
                "country_iso3": country,
                "country_label": COUNTRY_LABELS[country],
                "project_id": project,
                "collection_year": year,
                "n_target_lineage": len(target),
                "n_comparator_lineages": len(comparator),
                "target_median_lag_min_days": number_text(target_median_min),
                "target_median_lag_max_days": number_text(target_median_max),
                "comparator_median_lag_min_days": number_text(comparator_median_min),
                "comparator_median_lag_max_days": number_text(comparator_median_max),
                "target_minus_comparator_min_days": number_text(
                    target_median_min - comparator_median_max
                ),
                "target_minus_comparator_max_days": number_text(
                    target_median_max - comparator_median_min
                ),
                "comparison_direction": direction,
                "project_modal_public_date": project_modal_date,
                "project_modal_public_date_pct": pct_text(project_modal_n, len(group)),
                "target_modal_public_date": target_modal_date,
                "comparator_modal_public_date": comparator_modal_date,
                "groups_share_modal_public_date": target_modal_date == comparator_modal_date,
            }
        )
    return output_rows


def build_lag_summary() -> list[dict[str, object]]:
    rows = [row for row in read_tsv(SUMMARY_TABLE) if row.get("country_iso3") in TARGET_COUNTRIES]
    return [
        {
            "country_iso3": row["country_iso3"],
            "country_label": COUNTRY_LABELS[row["country_iso3"]],
            "lineage": row["primary_model_lineage_id"],
            "display_lineage": row["display_lineage"],
            "n_resurgence_focal_tips": row["n_resurgence_focal_tips"],
            "n_with_public_date": row["n_with_public_date"],
            "first_collection_lower": row["first_collection_lower"],
            "first_public_date": row["first_public_date"],
            "median_lag_min_days": row["median_lag_min_days"],
            "median_lag_max_days": row["median_lag_max_days"],
            "min_lag_min_days": row["min_lag_min_days"],
            "max_lag_max_days": row["max_lag_max_days"],
        }
        for row in rows
    ]


def build_project_batch_release() -> list[dict[str, object]]:
    rows = [
        row
        for row in read_tsv(PUBLIC_TABLE)
        if row.get("country_iso3") in TARGET_COUNTRIES and row.get("epidemic_period") == "resurgence"
    ]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["country_iso3"], row.get("project_id") or "missing_project")].append(row)

    output_rows: list[dict[str, object]] = []
    for (country, project), group in sorted(grouped.items()):
        public_dates = [row["public_date"] for row in group if row.get("public_date")]
        public_date_values = [parse_date(value) for value in public_dates]
        public_date_values = [value for value in public_date_values if value is not None]
        date_counts = Counter(public_dates)
        modal_date = ""
        modal_n = 0
        if date_counts:
            modal_n = max(date_counts.values())
            modal_date = min(date for date, count in date_counts.items() if count == modal_n)
        lag_mins = [int(row["lag_min_days"]) for row in group if row.get("lag_min_days") != ""]
        lag_maxs = [int(row["lag_max_days"]) for row in group if row.get("lag_max_days") != ""]
        output_rows.append(
            {
                "country_iso3": country,
                "country_label": COUNTRY_LABELS[country],
                "project_id": project,
                "n_focal_resurgence_genomes": len(group),
                "n_target_lineage": sum(row.get("primary_model_lineage_id") == TARGET_LINEAGE for row in group),
                "n_with_public_date": len(public_dates),
                "public_date_completeness_pct": pct_text(len(public_dates), len(group)),
                "collection_lower_min": min(row["collection_lower"] for row in group if row.get("collection_lower")),
                "collection_upper_effective_max": max(
                    row["collection_upper_effective"] for row in group if row.get("collection_upper_effective")
                ),
                "public_date_min": min(public_dates, default=""),
                "public_date_max": max(public_dates, default=""),
                "public_batch_span_days": (
                    (max(public_date_values) - min(public_date_values)).days if public_date_values else ""
                ),
                "modal_public_date": modal_date,
                "modal_public_date_n": modal_n,
                "modal_public_date_pct": pct_text(modal_n, len(public_dates)),
                "median_lag_min_days": number_text(statistics.median(lag_mins) if lag_mins else None),
                "median_lag_max_days": number_text(statistics.median(lag_maxs) if lag_maxs else None),
                "date_resolution_counts": counter_text([row.get("date_resolution", "") for row in group]),
                "public_route_counts": counter_text([row.get("public_route", "") for row in group]),
                "subnational_locations": ";".join(
                    sorted({row["subnational_location"] for row in group if row.get("subnational_location")})
                ),
            }
        )
    return output_rows


def build_geography_audit() -> list[dict[str, object]]:
    rows = [
        row
        for row in read_tsv(PUBLIC_TABLE)
        if row.get("country_iso3") in TARGET_COUNTRIES and row.get("epidemic_period") == "resurgence"
    ]
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["country_iso3"], row.get("project_id", ""), row.get("primary_model_lineage_id", ""))].append(row)

    output_rows: list[dict[str, object]] = []
    for (country, project, lineage), group in sorted(grouped.items()):
        subnational = [row for row in group if row.get("location_resolution") == "subnational"]
        output_rows.append(
            {
                "audit_level": "country_project_lineage",
                "country_iso3": country,
                "country_label": COUNTRY_LABELS[country],
                "project_id": project,
                "lineage": lineage,
                "n_genomes": len(group),
                "n_with_subnational_location": len(subnational),
                "subnational_location_completeness_pct": pct_text(len(subnational), len(group)),
                "subnational_locations": ";".join(
                    sorted({row["subnational_location"] for row in subnational if row.get("subnational_location")})
                ),
                "location_source_counts": counter_text([row.get("location_source", "") for row in group]),
            }
        )
    for country in TARGET_COUNTRIES:
        group = [row for row in rows if row["country_iso3"] == country]
        subnational = [row for row in group if row.get("location_resolution") == "subnational"]
        output_rows.append(
            {
                "audit_level": "country_total",
                "country_iso3": country,
                "country_label": COUNTRY_LABELS[country],
                "project_id": "ALL",
                "lineage": "ALL",
                "n_genomes": len(group),
                "n_with_subnational_location": len(subnational),
                "subnational_location_completeness_pct": pct_text(len(subnational), len(group)),
                "subnational_locations": ";".join(
                    sorted({row["subnational_location"] for row in subnational if row.get("subnational_location")})
                ),
                "location_source_counts": counter_text([row.get("location_source", "") for row in group]),
            }
        )
    return sorted(
        output_rows,
        key=lambda row: (str(row["country_iso3"]), str(row["audit_level"]), str(row["project_id"]), str(row["lineage"])),
    )


def first_two_consecutive(values: list[tuple[dt.date, float]], threshold: float) -> tuple[dt.date | None, float | None]:
    for index in range(len(values) - 1):
        date, value = values[index]
        next_date, next_value = values[index + 1]
        if value > threshold and next_value > threshold:
            return date, value
    return None, None


def build_case_clock_sensitivity() -> list[dict[str, object]]:
    threshold_rows = {row["country_iso3"]: row for row in read_tsv(CASE_THRESHOLDS)}
    output_rows: list[dict[str, object]] = []
    monthly_definitions = (
        ("first_above_2019_max", "first_post2022_month_above_2019_max", "cases_2019_max"),
        ("first_two_consecutive_above_2019_max", "first_two_consecutive_months_above_2019_max", "cases_2019_max"),
        ("first_above_2019_median", "first_post2022_month_above_2019_median", "cases_2019_median"),
        ("post2022_peak", "post2022_peak_month", "post2022_peak_cases"),
    )
    for country in TARGET_COUNTRIES:
        row = threshold_rows[country]
        for definition, date_column, value_column in monthly_definitions:
            output_rows.append(
                {
                    "country_iso3": country,
                    "country_label": COUNTRY_LABELS[country],
                    "analysis_scale": "harmonized_monthly",
                    "source_resolution": "monthly",
                    "milestone_definition": definition,
                    "milestone_date": row.get(date_column, ""),
                    "reference_value": row.get(value_column, ""),
                }
            )

    native_by_country: dict[str, list[tuple[dt.date, float, str]]] = defaultdict(list)
    for row in read_tsv(NATIVE_CASES):
        country = row.get("country_iso3", "")
        if country not in TARGET_COUNTRIES or row.get("value_status") not in {"reported", "reported_provisional"}:
            continue
        date = parse_date(row.get("observation_date"))
        if date:
            native_by_country[country].append((date, float(row["value"]), row.get("time_resolution", "")))
    for country in TARGET_COUNTRIES:
        values = sorted(native_by_country[country])
        values_2019 = [(date, value) for date, value, _ in values if date.year == 2019]
        post = [(date, value) for date, value, _ in values if date >= dt.date(2023, 1, 1)]
        if not values_2019 or not post:
            continue
        threshold = max(value for _, value in values_2019)
        first = next(((date, value) for date, value in post if value > threshold), (None, None))
        consecutive = first_two_consecutive(post, threshold)
        peak = max(post, key=lambda item: item[1])
        resolution = next((item[2] for item in values if item[2]), "native")
        for definition, result in (
            ("first_above_2019_max", first),
            ("first_two_consecutive_above_2019_max", consecutive),
            ("post2022_peak", peak),
        ):
            output_rows.append(
                {
                    "country_iso3": country,
                    "country_label": COUNTRY_LABELS[country],
                    "analysis_scale": "native_resolution",
                    "source_resolution": resolution,
                    "milestone_definition": definition,
                    "milestone_date": date_text(result[0]),
                    "reference_value": number_text(threshold if definition != "post2022_peak" else result[1]),
                }
            )
    return output_rows


def build_candidate_summary() -> list[dict[str, object]]:
    if not CANDIDATE_TABLE.exists():
        return []
    frozen_rows = read_tsv(PUBLIC_TABLE)
    frozen_by_project: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in frozen_rows:
        frozen_by_project[row.get("project_id", "")].append(row)
    output_rows: list[dict[str, object]] = []
    for row in read_tsv(CANDIDATE_TABLE):
        project = row["project_id"]
        overlap = frozen_by_project.get(project, [])
        output_rows.append(
            {
                "project_id": project,
                "n_explicit_pertussis_runs": row["n_runs"],
                "countries": row["countries"],
                "collection_span": f"{row['collection_lower_min']} to {row['collection_lower_max']}",
                "first_public_span": f"{row['first_public_min']} to {row['first_public_max']}",
                "date_resolution_summary": (
                    f"day={row['n_collection_day']};month_or_interval={row['n_collection_month_or_interval']};"
                    f"year_only={row['n_collection_year_only']};missing={row['n_collection_missing_or_not_applicable']}"
                ),
                "n_frozen_tree_tips": len(overlap),
                "n_frozen_target_lineage": sum(item.get("primary_model_lineage_id") == TARGET_LINEAGE for item in overlap),
                "n_frozen_resurgence_target": sum(
                    item.get("primary_model_lineage_id") == TARGET_LINEAGE
                    and item.get("epidemic_period") == "resurgence"
                    for item in overlap
                ),
                "recommended_tier": row["recommended_tier"],
            }
        )
    return output_rows


def markdown_table(rows: list[dict[str, object]], fieldnames: list[str]) -> str:
    lines = ["| " + " | ".join(fieldnames) + " |", "| " + " | ".join(["---"] * len(fieldnames)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(field, "")) for field in fieldnames) + " |")
    return "\n".join(lines)


def write_appendix(
    shift_rows: list[dict[str, object]],
    milestone_rows: list[dict[str, object]],
    lag_rows: list[dict[str, object]],
    candidate_rows: list[dict[str, object]],
    threshold_rows: list[dict[str, object]],
    project_rows: list[dict[str, object]],
    project_lineage_rows: list[dict[str, object]],
    geography_rows: list[dict[str, object]],
    case_rows: list[dict[str, object]],
) -> None:
    target_lags = [row for row in lag_rows if row["lineage"] == TARGET_LINEAGE]
    country_geography = [row for row in geography_rows if row["audit_level"] == "country_total"]
    snapshot_rows = [
        row for row in read_tsv(SNAPSHOT_MANIFEST) if row.get("record_type") == "repository_snapshot"
    ]
    country_name = COUNTRY_LABELS
    milestone_labels = {
        "first_post2022_month_above_2019_max": "First month above 2019 maximum",
        "post2022_peak_month": "Post-2022 peak month",
    }
    case_definition_labels = {
        "first_above_2019_max": "First period above 2019 maximum",
        "first_two_consecutive_above_2019_max": "First 2 consecutive periods above 2019 maximum",
        "first_above_2019_median": "First month above 2019 median",
        "post2022_peak": "Post-2022 peak",
    }
    timing_labels = {
        "interval_before_case_threshold": "Collection interval before case threshold",
        "interval_spans_case_threshold": "Collection interval spans case threshold",
        "interval_after_case_threshold": "Collection interval after case threshold",
        "before_case_threshold": "Public before case threshold",
        "on_or_after_case_threshold": "Public on or after case threshold",
    }
    route_labels = {
        "ena_read_only": "ENA read record only",
        "ena_read_first": "ENA read record first",
        "ncbi_assembly_only": "NCBI Assembly only",
        "ncbi_assembly_first": "NCBI Assembly first",
        "same_date": "ENA and NCBI on same date",
        "not_publicly_dated": "No reproducible public date",
    }
    candidate_labels = {
        "boundary_dataset_pending_date_resolution": "Year-only collection dates preclude month-scale inclusion",
        "partial_current_project_extension": "Partly represented in the frozen cohort; additional records require lineage placement",
        "supplement_sensitivity_then_lineage_place": "Post-freeze project; quality control and lineage placement required",
        "discussion_or_sensitivity_after_taxon_check": "Post-freeze project; taxonomic review and lineage placement required",
    }

    def display_counts(value: object) -> str:
        parts: list[str] = []
        for item in str(value or "").split(";"):
            if not item:
                continue
            label, count = item.rsplit("=", 1)
            parts.append(f"{route_labels.get(label, label)}: {count}")
        return "; ".join(parts)

    def display_date_resolution(value: object) -> str:
        parsed = {
            item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in str(value or "").split(";")
            if "=" in item
        }
        return (
            f"Day: {parsed.get('day', '0')}; month or interval: {parsed.get('month_or_interval', '0')}; "
            f"year only: {parsed.get('year_only', '0')}; missing: {parsed.get('missing', '0')}"
        )

    primary_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "Threshold, k": row["cumulative_genome_threshold"],
            "Collection interval": f"{row['collection_detection_lower']} to {row['collection_detection_upper']}",
            "Sequence-public date": row["public_detection_date"],
            "Clock displacement, d": f"{row['clock_shift_min_days']}–{row['clock_shift_max_days']}",
            "Public timing after milestone month, d": relative_range_text(
                row.get("public_relative_to_resurgence_min_days"),
                row.get("public_relative_to_resurgence_max_days"),
            ),
        }
        for row in shift_rows
    ]
    milestone_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "Milestone": milestone_labels[str(row["milestone"])],
            "Milestone month": row["milestone_period"],
            "Definitely collected": row["definitely_collected_by_milestone"],
            "Possibly collected": row["possibly_collected_by_milestone"],
            "Public sequence available": row["public_by_milestone"],
        }
        for row in milestone_rows
    ]
    threshold_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "k": row["cumulative_genome_threshold"],
            "Collection interval": f"{row['collection_detection_lower']} to {row['collection_detection_upper']}",
            "Sequence-public date": row["public_detection_date"],
            "Displacement, d": f"{row['clock_shift_min_days']}–{row['clock_shift_max_days']}",
            "Collection timing relative to milestone month, d": relative_range_text(
                row.get("collection_relative_to_resurgence_min_days"),
                row.get("collection_relative_to_resurgence_max_days"),
            ),
            "Public timing after milestone month, d": relative_range_text(
                row.get("public_relative_to_resurgence_min_days"),
                row.get("public_relative_to_resurgence_max_days"),
            ),
            "Collection timing": timing_labels[str(row["collection_timing_class"])],
            "Public timing": timing_labels[str(row["public_timing_class"])],
        }
        for row in threshold_rows
    ]
    lag_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "Target-lineage genomes": row["n_resurgence_focal_tips"],
            "Sequence-public date available": row["n_with_public_date"],
            "Median minimum lag, d": row["median_lag_min_days"],
            "Median maximum lag, d": row["median_lag_max_days"],
            "Observed lag range, d": f"{row['min_lag_min_days']}–{row['max_lag_max_days']}",
        }
        for row in target_lags
    ]
    project_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "BioProject": row["project_id"],
            "Resurgence genomes": row["n_focal_resurgence_genomes"],
            "Target lineage": row["n_target_lineage"],
            "Sequence-public date span": f"{row['public_date_min']} to {row['public_date_max']}",
            "Batch span, d": row["public_batch_span_days"],
            "Modal sequence-public date proportion, %": row["modal_public_date_pct"],
            "First observable sequence-data route": display_counts(row["public_route_counts"]),
        }
        for row in project_rows
    ]
    project_lineage_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "BioProject": row["project_id"],
            "Collection year": row["collection_year"],
            "Target n": row["n_target_lineage"],
            "Comparator n": row["n_comparator_lineages"],
            "Target median interval, d": (
                f"{row['target_median_lag_min_days']}–{row['target_median_lag_max_days']}"
            ),
            "Comparator median interval, d": (
                f"{row['comparator_median_lag_min_days']}–{row['comparator_median_lag_max_days']}"
            ),
            "Direction": str(row["comparison_direction"]).replace("_", " "),
            "Shared modal public date": row["groups_share_modal_public_date"],
        }
        for row in project_lineage_rows
    ]
    geography_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "Genomes": row["n_genomes"],
            "With subnational metadata": row["n_with_subnational_location"],
            "Completeness, %": row["subnational_location_completeness_pct"],
            "Reported locations": str(row["subnational_locations"]).replace(";", "; "),
        }
        for row in country_geography
    ]
    case_display = [
        {
            "Country": country_name[str(row["country_iso3"])],
            "Analysis scale": "Harmonized monthly" if row["analysis_scale"] == "harmonized_monthly" else "Native reporting resolution",
            "Source resolution": row["source_resolution"],
            "Milestone definition": case_definition_labels[str(row["milestone_definition"])],
            "Date": row["milestone_date"],
            "Reference case count": row["reference_value"],
        }
        for row in case_rows
    ]
    candidate_display = [
        {
            "BioProject": row["project_id"],
            "Explicit B. pertussis runs": row["n_explicit_pertussis_runs"],
            "Country": row["countries"],
            "Collection span": row["collection_span"],
            "Date resolution": display_date_resolution(row["date_resolution_summary"]),
            "Frozen-tree genomes": row["n_frozen_tree_tips"],
            "Frozen target lineage": row["n_frozen_target_lineage"],
            "Frozen resurgence target": row["n_frozen_resurgence_target"],
            "Interpretive status": candidate_labels.get(str(row["recommended_tier"]), str(row["recommended_tier"])),
        }
        for row in candidate_rows
    ]

    metadata_descriptions = {
        "eid_nas_focal_genome_metadata.tsv": (
            "Three-country focal-genome metadata",
            "Accessions, collection intervals, geography, sequencing technology, and assembly release dates",
            "Matched to frozen-tree identifiers and checked for accession uniqueness",
        ),
        "eid_nas_highres_cases.tsv": (
            "National surveillance observations",
            "Monthly or weekly reported pertussis cases and source provenance",
            "Checked for temporal completeness and harmonized to monthly scale for the primary analysis",
        ),
        "eid_nas_public_health_sources.tsv": (
            "Public-health source registry",
            "Issuing organization, source type, release information, and country",
            "Used to verify the provenance and reporting context of surveillance series",
        ),
        "eid_nas_prjna1071282_run_audit.tsv": (
            "PRJNA1071282 run metadata",
            "Run and sample accessions, organism annotation, collection year, and public date",
            "Restricted to explicit B. pertussis records and reconciled with frozen-tree membership",
        ),
        "eid_nas_china_recovery_summary.tsv": (
            "China metadata recovery audit",
            "Project-level metadata recovery and matching summaries",
            "Used to distinguish recovered metadata from unresolved records",
        ),
    }
    metadata_display: list[dict[str, object]] = []
    for row in snapshot_rows:
        key = Path(str(row["path"])).name
        if key not in metadata_descriptions:
            continue
        domain, content, validation = metadata_descriptions[key]
        metadata_display.append(
            {
                "Data component": domain,
                "Records": row["n_records"],
                "Information retained": content,
                "Analytic role and validation": validation,
            }
        )

    case_data = [row for row in read_tsv(HIGHRES_CASE_SNAPSHOT) if row.get("country_iso3") in TARGET_COUNTRIES]
    case_source_names = {
        "AUS": "Australian National Notifiable Disease Surveillance System",
        "CHN": "Official Chinese national and provincial public-health surveillance reports",
        "JPN": "National Institute of Infectious Diseases and Japan Institute for Health Security weekly reports",
    }
    case_source_display: list[dict[str, object]] = []
    for country in TARGET_COUNTRIES:
        rows = [row for row in case_data if row["country_iso3"] == country]
        dates = sorted(row["date"] for row in rows if row.get("date"))
        case_source_display.append(
            {
                "Country": country_name[country],
                "Surveillance source": case_source_names[country],
                "Native resolution": ", ".join(sorted({row["time_resolution"] for row in rows})),
                "Observation period": f"{dates[0]} to {dates[-1]}",
                "Data freeze": max(row["data_freeze_date"] for row in rows),
            }
        )

    appendix = f"""# Appendix. Public-Archive Timing for MT28-Associated Pertussis Genomes

## Supplementary Methods

### Study design and frozen analysis boundary

This retrospective audit compared specimen accumulation, public sequence availability, and national pertussis resurgence milestones. It retained the previously frozen 989-genome core-SNP phylogeny, including 774 focal genomes and 215 stratified global-background genomes. Existing L1_02.07 assignments supplied retrospective lineage labels; the phylogeny was not rebuilt, lineage definitions were not re-estimated, and fitted transmission models were outside this audit.

The timing analysis focused on frozen focal genomes from Australia, China, and Japan because these countries had compatible national case series and target-lineage genomes collected during the resurgence period. Projects identified after the phylogenetic freeze entered a separate accession, date-resolution, sequence-quality, and lineage-placement audit.

### Collection intervals and public sequence availability

Reported collection dates were represented by lower and upper bounds. Exact dates had identical bounds; month-level dates spanned the reported calendar month; year-level dates spanned the reported calendar year. The effective collection upper bound was the earlier of the reported upper bound and the first reproducible public date. A public date earlier than the collection lower bound retained an explicit temporal-conflict flag.

Public sequence availability was the earliest archive-recorded first-public or release date for an ENA read record or NCBI Assembly record. ENA and NCBI dates were retained separately before selecting the earliest sequence-data route. Metadata-only BioSample first-public dates were excluded from the primary endpoint and interpreted, when reviewed, only as metadata visibility rather than sequence-data availability. This endpoint measures archive-recorded external sequence visibility, not a historical replay of database search results.

### Cumulative visibility curves and interpretive anchors

Cumulative curves reported possibly collected, definitely collected, and publicly retrievable target-lineage records at every observed event date. Milestone-month counts are reported as counts at the first day of the milestone month, whereas relative day offsets use the full milestone-month interval to avoid assigning day-level precision to monthly surveillance data. The *k*=5 sequence count summarized the China and Japan curves as a descriptive secondary anchor. Australia used its 3 available target-lineage genomes as a descriptive contrast. The collection-to-public availability interval ranged from the public sequence date minus the accumulation upper bound to the public sequence date minus the accumulation lower bound.

At each epidemiologic milestone, genomes were classified as definitely collected when the effective collection upper bound had passed and possibly collected when the collection lower bound had passed. Publicly available genomes were counted independently from their reproducible public dates. These counts characterize accumulation and external visibility within the frozen genomic sample.

### Epidemiologic milestones and sensitivity analyses

The primary case milestone was the first post-2022 calendar month in which national reported cases exceeded the country-specific maximum monthly count observed in 2019. The highest observed post-2022 monthly count through the case-data cutoff was evaluated as a second marker. Monthly milestones were treated as calendar-month intervals for day-offset calculations. Weekly Japanese reports were aggregated to calendar months for the primary comparison and retained at native weekly resolution in sensitivity analyses. National case milestones supplied contextual epidemiologic clocks only and should not be interpreted as evidence that the included lineage caused, represented, or predicted national resurgence. Sensitivity analyses evaluated native reporting resolution, the first 2 consecutive periods above the 2019 maximum, the first month above the 2019 median, and genome thresholds of *k*=1, 3, 5, 10, 20, and 50 when sufficient genomes were available.

Project-level analyses summarized sequence-public date completeness, modal release concentration, release-date span, and first observable sequence-data route. Target and comparator lineages were additionally compared within country, BioProject, and collection year. We reported interval medians and their direction without genome-level significance tests because records within release batches were dependent. Geographic metadata were standardized to subnational units and characterized sample composition.

### Candidate-project and metadata audits

PRJNA1071282 contained 734 runs explicitly annotated as *B. pertussis*. Sixteen were represented in the frozen tree, including 6 target-lineage genomes and 3 resurgence-period target-lineage genomes. All 734 runs had year-level collection dates. The month-scale analysis used the 16 frozen genomes, and the complete extension entered the candidate-project audit.

The analytic metadata extract retained accession identifiers, BioSample and BioProject identifiers, collection intervals, geography, sequencing technology, separate ENA read and NCBI Assembly public dates, and record-matching status. Its analytic scope was accession-level metadata for records already represented in the frozen tree. Validation included accession uniqueness, cross-source matching, temporal-order checks, monotonicity of availability dates as *k* increased, and regeneration of derived results from the frozen metadata.

## Appendix Table 1. Primary interval-aware sequence-availability clock

{markdown_table(primary_display, list(primary_display[0]))}

## Appendix Table 2. Target-lineage visibility at national case milestones

{markdown_table(milestone_display, list(milestone_display[0]))}

## Appendix Table 3. Sensitivity to the cumulative genome threshold

{markdown_table(threshold_display, list(threshold_display[0]))}

## Appendix Table 4. Collection-to-public lag intervals for target-lineage genomes

{markdown_table(lag_display, list(lag_display[0]))}

## Appendix Table 5. Project-level public-release patterns

{markdown_table(project_display, list(project_display[0]))}

## Appendix Table 6. Project- and year-matched target versus comparator lineages

{markdown_table(project_lineage_display, list(project_lineage_display[0]))}

## Appendix Table 7. Completeness and composition of subnational geographic metadata

{markdown_table(geography_display, list(geography_display[0]))}

## Appendix Table 8. Sensitivity to the national resurgence-milestone definition

{markdown_table(case_display, list(case_display[0]))}

## Appendix Table 9. Projects identified after the phylogenetic freeze

{markdown_table(candidate_display, list(candidate_display[0]))}

## Appendix Table 10. Metadata components and validation roles

{markdown_table(metadata_display, list(metadata_display[0]))}

## Appendix Table 11. National pertussis surveillance series

{markdown_table(case_source_display, list(case_source_display[0]))}

## Definitions and interpretation

The reported collection interval preserves source precision. The effective upper bound enforces temporal coherence with public sequence availability while retaining explicit conflict flags. Minimum and maximum availability lags are calculated from the effective upper and lower collection bounds. Public route identifies whether an ENA read record, an NCBI Assembly record, or both supplied the earliest reproducible date. Subnational metadata describe locations represented in the genomic sample.

Lineage-specific analyses assigned target status to frozen-tree identifiers and thereby preserved the original phylogenetic analysis boundary.
"""
    APPENDIX_MD.write_text(appendix)


def main() -> None:
    threshold_rows = build_threshold_sensitivity()
    shift_rows = build_detection_shift(threshold_rows)
    milestone_rows = build_milestone_visibility()
    cumulative_visibility_rows = build_cumulative_visibility()
    project_lineage_rows = build_project_lineage_comparison()
    lag_rows = build_lag_summary()
    project_rows = build_project_batch_release()
    geography_rows = build_geography_audit()
    case_rows = build_case_clock_sensitivity()
    candidate_rows = build_candidate_summary()

    outputs = [
        (SHIFT_OUTPUT, shift_rows, [
            "country_iso3", "country_label", "target_lineage", "cumulative_genome_threshold",
            "collection_detection_lower", "collection_detection_upper", "public_detection_date",
            "resurgence_milestone_start", "resurgence_milestone_end",
            "clock_shift_min_days", "clock_shift_max_days",
            "collection_lead_to_case_threshold_min_days", "collection_lead_to_case_threshold_max_days",
            "public_lead_to_case_threshold_min_days", "public_lead_to_case_threshold_max_days",
            "public_lead_to_case_threshold_days",
            "collection_relative_to_resurgence_min_days", "collection_relative_to_resurgence_max_days",
            "public_relative_to_resurgence_min_days", "public_relative_to_resurgence_max_days",
            "public_relative_to_resurgence_days", "collection_timing_class", "public_timing_class",
        ]),
        (MILESTONE_OUTPUT, milestone_rows, [
            "country_iso3", "country_label", "target_lineage", "milestone", "milestone_date",
            "milestone_period_start", "milestone_period_end", "milestone_period",
            "total_resurgence_target_genomes", "definitely_collected_by_milestone",
            "possibly_collected_by_milestone", "public_by_milestone",
            "definitely_collected_not_public", "possibly_collected_not_public",
            "public_among_possible_collected_pct",
        ]),
        (LAG_SUMMARY_OUTPUT, lag_rows, [
            "country_iso3", "country_label", "lineage", "display_lineage",
            "n_resurgence_focal_tips", "n_with_public_date", "first_collection_lower", "first_public_date",
            "median_lag_min_days", "median_lag_max_days", "min_lag_min_days", "max_lag_max_days",
        ]),
        (CANDIDATE_OUTPUT, candidate_rows, [
            "project_id", "n_explicit_pertussis_runs", "countries", "collection_span", "first_public_span",
            "date_resolution_summary", "n_frozen_tree_tips", "n_frozen_target_lineage",
            "n_frozen_resurgence_target", "recommended_tier",
        ]),
        (THRESHOLD_SENSITIVITY_OUTPUT, threshold_rows, [
            "country_iso3", "country_label", "target_lineage", "cumulative_genome_threshold",
            "collection_detection_lower", "collection_detection_upper", "public_detection_date",
            "resurgence_milestone_start", "resurgence_milestone_end",
            "clock_shift_min_days", "clock_shift_max_days",
            "collection_lead_to_case_threshold_min_days", "collection_lead_to_case_threshold_max_days",
            "public_lead_to_case_threshold_min_days", "public_lead_to_case_threshold_max_days",
            "public_lead_to_case_threshold_days",
            "collection_relative_to_resurgence_min_days", "collection_relative_to_resurgence_max_days",
            "public_relative_to_resurgence_min_days", "public_relative_to_resurgence_max_days",
            "public_relative_to_resurgence_days", "collection_timing_class", "public_timing_class",
        ]),
        (PROJECT_BATCH_OUTPUT, project_rows, [
            "country_iso3", "country_label", "project_id", "n_focal_resurgence_genomes", "n_target_lineage",
            "n_with_public_date", "public_date_completeness_pct", "collection_lower_min",
            "collection_upper_effective_max", "public_date_min", "public_date_max", "public_batch_span_days", "modal_public_date",
            "modal_public_date_n", "modal_public_date_pct", "median_lag_min_days", "median_lag_max_days",
            "date_resolution_counts", "public_route_counts", "subnational_locations",
        ]),
        (GEOGRAPHY_OUTPUT, geography_rows, [
            "audit_level", "country_iso3", "country_label", "project_id", "lineage", "n_genomes",
            "n_with_subnational_location", "subnational_location_completeness_pct", "subnational_locations",
            "location_source_counts",
        ]),
        (CASE_SENSITIVITY_OUTPUT, case_rows, [
            "country_iso3", "country_label", "analysis_scale", "source_resolution",
            "milestone_definition", "milestone_date", "reference_value",
        ]),
        (CUMULATIVE_VISIBILITY_OUTPUT, cumulative_visibility_rows, [
            "country_iso3", "country_label", "target_lineage", "event_date",
            "n_possibly_collected", "n_definitely_collected", "n_publicly_available",
            "total_target_lineage_records", "resurgence_milestone_date", "resurgence_milestone_end",
            "post2022_peak_date", "post2022_peak_end",
        ]),
        (PROJECT_LINEAGE_COMPARISON_OUTPUT, project_lineage_rows, [
            "country_iso3", "country_label", "project_id", "collection_year",
            "n_target_lineage", "n_comparator_lineages",
            "target_median_lag_min_days", "target_median_lag_max_days",
            "comparator_median_lag_min_days", "comparator_median_lag_max_days",
            "target_minus_comparator_min_days", "target_minus_comparator_max_days",
            "comparison_direction", "project_modal_public_date", "project_modal_public_date_pct",
            "target_modal_public_date", "comparator_modal_public_date",
            "groups_share_modal_public_date",
        ]),
    ]
    for path, rows, fields in outputs:
        write_tsv(path, rows, fields)
    write_appendix(
        shift_rows,
        milestone_rows,
        lag_rows,
        candidate_rows,
        threshold_rows,
        project_rows,
        project_lineage_rows,
        geography_rows,
        case_rows,
    )
    for path, _, _ in outputs:
        print(f"wrote {path.relative_to(REPO)}")
    print(f"wrote {APPENDIX_MD.relative_to(REPO)}")


if __name__ == "__main__":
    main()
