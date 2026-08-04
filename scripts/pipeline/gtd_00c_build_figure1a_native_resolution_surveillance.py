#!/usr/bin/env python3
"""Build native-resolution surveillance series used in Figure 1A.

Figure 1A is descriptive rather than a harmonised incidence comparison.  This
builder therefore preserves the native time resolution and surveillance
measure of each source:

* Australia: monthly notifications;
* China: monthly reported cases;
* Japan: weekly reported cases;
* Belgium: monthly NRC-confirmed cases digitised from the vector bars in the
  official 2024 surveillance report; and
* France: monthly PCR positivity from 2019 to 2024 and monthly PCR-positive
  test counts for 2023-2024 from official 3Labos figures.

The Belgium and France series are not used in the case likelihood.  Approximate
digitised values are marked as such, and PCR-positive tests are not treated as
deduplicated patients or compatible national case denominators.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_JAPAN_INPUT = ROOT / "data/raw/pertussis_incidence_timeseries.csv"
START_YEAR = 2019
END_YEAR = 2025

BELGIUM_SOURCE_URL = (
    "https://www.sciensano.be/sites/default/files/"
    "kinkhoest-epidemiologie-jaarrapport-2024-nl.pdf"
)
FRANCE_SOURCE_URL = "https://www.santepubliquefrance.fr/coqueluche/donnees"
FRANCE_REPORT_URL = (
    "https://www.santepubliquefrance.fr/sites/default/files/rdd/document/"
    "bullnat_coqueluche_20251013.pdf"
)

# Values read from the vector bar geometry in Figure 2 of the Sciensano report.
# The 2024 bars sum to 3,398, whereas the later report text states 3,404; the
# bar-derived monthly values are retained without forcing them to the annual
# textual total.
BELGIUM_MONTHLY_NRC_CASES = {
    2019: [65, 56, 52, 46, 58, 46, 61, 82, 94, 40, 35, 50],
    2023: [13, 17, 35, 29, 38, 71, 92, 161, 163, 138, 137, 154],
    2024: [199, 230, 335, 437, 443, 436, 354, 351, 278, 146, 126, 63],
}

# Values for 2019-2022 were digitised to the nearest 0.1 percentage point from
# the monthly positivity line on the official data page.  Values for 2023-2024
# were transcribed from labels in Figure 13 of the final 2024 national report.
FRANCE_MONTHLY_PCR_POSITIVITY = {
    2019: [6.1, 5.7, 9.8, 10.9, 13.5, 16.6, 19.0, 20.5, 12.5, 7.3, 6.6, 4.3],
    2020: [6.8, 6.4, 10.0, 16.0, 9.3, 7.0, 4.6, 2.1, 0.4, 0.4, 0.1, 0.4],
    2021: [0.9, 0.4, 0.5, 1.1, 0.6, 0.2, 0.4, 0.1, 0.1, 0.1, 0.2, 0.4],
    2022: [0.6, 0.9, 1.4, 2.5, 3.7, 7.9, 7.5, 5.1, 1.2, 1.2, 1.2, 0.6],
    2023: [0.5, 1.0, 1.3, 1.0, 1.4, 2.4, 4.7, 5.8, 7.8, 4.7, 5.2, 5.9],
    2024: [7.6, 11.4, 22.2, 24.3, 27.9, 21.4, 23.9, 28.8, 14.4, 7.6, 8.1, 5.4],
}

# Exact labels in Figure 13.  Annual sums reproduce the official Figure 12
# totals exactly: 518 positive tests in 2023 and 38,847 in 2024.
FRANCE_MONTHLY_PCR_POSITIVE_TESTS = {
    2023: [7, 11, 14, 9, 10, 21, 49, 58, 86, 65, 82, 106],
    2024: [146, 237, 710, 1665, 5626, 15469, 5915, 5073, 2318, 883, 484, 321],
}

FIELDNAMES = [
    "country_iso3",
    "country_name",
    "observation_date",
    "period_start",
    "period_end",
    "year",
    "month",
    "week",
    "value",
    "time_resolution",
    "metric",
    "y_axis_title",
    "surveillance_scope",
    "extraction_method",
    "value_status",
    "source_url",
    "source_detail",
    "data_freeze_date",
]


def next_month_start(year: int, month: int) -> date:
    if month == 12:
        return date(year + 1, 1, 1)
    return date(year, month + 1, 1)


def read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def monthly_rows(
    source_path: Path,
    data_freeze_date: str,
) -> tuple[list[dict[str, object]], dict[tuple[str, int, int], int]]:
    source = read_delimited(source_path, "\t")
    output: list[dict[str, object]] = []
    monthly_lookup: dict[tuple[str, int, int], int] = {}
    metadata = {
        "AUS": {
            "country_name": "Australia",
            "metric": "notified_cases",
            "y_axis_title": "Notifications per month",
            "surveillance_scope": "National notifiable-disease notifications",
        },
        "CHN": {
            "country_name": "China",
            "metric": "reported_cases",
            "y_axis_title": "Reported cases per month",
            "surveillance_scope": "National reported pertussis cases",
        },
    }

    for row in source:
        iso3 = row["country_iso3"]
        if iso3 not in metadata:
            continue
        year = int(row["year"])
        month = int(row["month"])
        if not START_YEAR <= year <= END_YEAR:
            continue
        if row["case_data_available"].lower() != "true":
            raise ValueError(f"Missing monthly value for {iso3} {year}-{month:02d}")
        value_float = float(row["cases"])
        if value_float < 0 or not value_float.is_integer():
            raise ValueError(
                f"Invalid monthly case count for {iso3} {year}-{month:02d}: "
                f"{row['cases']!r}"
            )
        value = int(value_float)
        period_start = date(year, month, 1)
        country_meta = metadata[iso3]
        output.append(
            {
                "country_iso3": iso3,
                "country_name": country_meta["country_name"],
                "observation_date": period_start.isoformat(),
                "period_start": period_start.isoformat(),
                "period_end": next_month_start(year, month).isoformat(),
                "year": year,
                "month": month,
                "week": "",
                "value": value,
                "time_resolution": "monthly",
                "metric": country_meta["metric"],
                "y_axis_title": country_meta["y_axis_title"],
                "surveillance_scope": country_meta["surveillance_scope"],
                "extraction_method": "direct_tabular",
                "value_status": "reported",
                "source_url": row["source_url"],
                "source_detail": row["source_file"],
                "data_freeze_date": data_freeze_date,
            }
        )
        monthly_lookup[(iso3, year, month)] = value

    for iso3 in metadata:
        n_rows = sum(row["country_iso3"] == iso3 for row in output)
        if n_rows != 84:
            raise ValueError(f"Expected 84 monthly rows for {iso3}, found {n_rows}")
    return output, monthly_lookup


def japan_weekly_rows(
    source_path: Path,
    monthly_lookup: dict[tuple[str, int, int], int],
    data_freeze_date: str,
) -> list[dict[str, object]]:
    source = read_delimited(source_path, ",")
    output: list[dict[str, object]] = []
    monthly_sums: defaultdict[tuple[str, int, int], int] = defaultdict(int)

    for row in source:
        if row["iso3"] != "JPN":
            continue
        year = int(row["Year"])
        if not START_YEAR <= year <= END_YEAR:
            continue
        if row["reporting_frequency"] != "weekly":
            raise ValueError(f"Unexpected Japan frequency: {row['reporting_frequency']!r}")
        month = int(row["Month"])
        value_float = float(row["Cases"])
        week_float = float(row["Week"])
        if (
            value_float < 0
            or not value_float.is_integer()
            or not week_float.is_integer()
        ):
            raise ValueError(f"Invalid Japan weekly row: {row!r}")
        value = int(value_float)
        week = int(week_float)
        source_period_start = datetime.strptime(
            row["period_start"], "%Y-%m-%d"
        ).date()
        reporting_year = source_period_start.isocalendar().year
        # Reconstruct the Monday boundary from the published epidemiological
        # week.  This corrects a known date-offset artefact in several late
        # 2025 rows of the harmonised input while leaving reported counts
        # unchanged.
        period_start = date.fromisocalendar(reporting_year, week, 1)
        period_end = date.fromisocalendar(reporting_year, week, 7)
        period_end = date.fromordinal(period_end.toordinal() + 1)
        official_index_url = (
            "https://id-info.jihs.go.jp/en/surveillance/idwr/rapid/"
            f"{reporting_year}/index.html"
        )
        output.append(
            {
                "country_iso3": "JPN",
                "country_name": "Japan",
                "observation_date": period_start.isoformat(),
                "period_start": period_start.isoformat(),
                "period_end": period_end.isoformat(),
                "year": period_start.year,
                "month": period_start.month,
                "week": week,
                "value": value,
                "time_resolution": "weekly",
                "metric": "reported_cases",
                "y_axis_title": "Reported cases per week",
                "surveillance_scope": (
                    "National notifiable-disease surveillance; weekly "
                    "provisional reports"
                ),
                "extraction_method": "direct_tabular",
                "value_status": "reported_provisional",
                "source_url": official_index_url,
                "source_detail": (
                    "IDWR weekly surveillance table; weekly rows aggregate "
                    "exactly to the frozen monthly model series"
                ),
                "data_freeze_date": data_freeze_date,
            }
        )
        monthly_sums[("JPN", year, month)] += value

    if len(output) != 364:
        raise ValueError(f"Expected 364 Japan weekly rows, found {len(output)}")
    expected_months = {
        key: value
        for key, value in monthly_lookup.items()
        if key[0] == "JPN"
    }
    if monthly_sums != expected_months:
        mismatches = {
            key: (expected_months.get(key), monthly_sums.get(key))
            for key in set(expected_months) | set(monthly_sums)
            if expected_months.get(key) != monthly_sums.get(key)
        }
        raise ValueError(
            "Japan weekly values do not reproduce the frozen monthly series: "
            f"{mismatches}"
        )
    return output


def digitised_monthly_rows(data_freeze_date: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for year, values in BELGIUM_MONTHLY_NRC_CASES.items():
        if len(values) != 12:
            raise ValueError(f"Belgium {year} does not contain 12 monthly values")
        for month, value in enumerate(values, start=1):
            period_start = date(year, month, 1)
            output.append(
                {
                    "country_iso3": "BEL",
                    "country_name": "Belgium",
                    "observation_date": period_start.isoformat(),
                    "period_start": period_start.isoformat(),
                    "period_end": next_month_start(year, month).isoformat(),
                    "year": year,
                    "month": month,
                    "week": "",
                    "value": value,
                    "time_resolution": "monthly",
                    "metric": "nrc_confirmed_cases",
                    "y_axis_title": "NRC-confirmed cases per month",
                    "surveillance_scope": (
                        "National Reference Centre laboratory-confirmed "
                        "surveillance; non-exhaustive trend series"
                    ),
                    "extraction_method": (
                        "digitised_from_official_vector_figure"
                    ),
                    "value_status": "approximate",
                    "source_url": BELGIUM_SOURCE_URL,
                    "source_detail": (
                        "Figure 2, monthly NRC reports; 2020-2022 monthly "
                        "bars are not displayed and remain missing"
                    ),
                    "data_freeze_date": data_freeze_date,
                }
            )

    for year, values in FRANCE_MONTHLY_PCR_POSITIVITY.items():
        if len(values) != 12:
            raise ValueError(f"France {year} does not contain 12 monthly values")
        labelled_in_report = year >= 2023
        for month, value in enumerate(values, start=1):
            period_start = date(year, month, 1)
            output.append(
                {
                    "country_iso3": "FRA",
                    "country_name": "France",
                    "observation_date": period_start.isoformat(),
                    "period_start": period_start.isoformat(),
                    "period_end": next_month_start(year, month).isoformat(),
                    "year": year,
                    "month": month,
                    "week": "",
                    "value": f"{value:.1f}",
                    "time_resolution": "monthly",
                    "metric": "pcr_positivity_percent",
                    "y_axis_title": "PCR positivity per month (%)",
                    "surveillance_scope": (
                        "3Labos (Cerba and Eurofins-Biomnis) PCR testing "
                        "activity; positivity is not a national case count"
                    ),
                    "extraction_method": (
                        "transcribed_from_official_figure_labels"
                        if labelled_in_report
                        else "digitised_from_official_raster_figure"
                    ),
                    "value_status": (
                        "reported" if labelled_in_report else "approximate"
                    ),
                    "source_url": (
                        FRANCE_REPORT_URL
                        if labelled_in_report
                        else FRANCE_SOURCE_URL
                    ),
                    "source_detail": (
                        "Figure 13, monthly positivity labels"
                        if labelled_in_report
                        else (
                            "Monthly positivity line in the official 3Labos "
                            "figure; digitised to 0.1 percentage point"
                        )
                    ),
                    "data_freeze_date": data_freeze_date,
                }
            )

    for year, values in FRANCE_MONTHLY_PCR_POSITIVE_TESTS.items():
        if len(values) != 12:
            raise ValueError(
                f"France {year} does not contain 12 monthly positive-test values"
            )
        for month, value in enumerate(values, start=1):
            period_start = date(year, month, 1)
            output.append(
                {
                    "country_iso3": "FRA",
                    "country_name": "France",
                    "observation_date": period_start.isoformat(),
                    "period_start": period_start.isoformat(),
                    "period_end": next_month_start(year, month).isoformat(),
                    "year": year,
                    "month": month,
                    "week": "",
                    "value": value,
                    "time_resolution": "monthly",
                    "metric": "pcr_positive_tests",
                    "y_axis_title": "PCR-positive tests per month",
                    "surveillance_scope": (
                        "3Labos (Cerba and Eurofins-Biomnis) PCR testing "
                        "activity; positive tests are not deduplicated "
                        "patients or nationally notified cases"
                    ),
                    "extraction_method": (
                        "transcribed_from_official_figure_labels"
                    ),
                    "value_status": "reported",
                    "source_url": FRANCE_REPORT_URL,
                    "source_detail": (
                        "Figure 13, monthly PCR-positive test labels; annual "
                        "sums match Figure 12 (2023: 518; 2024: 38,847)"
                    ),
                    "data_freeze_date": data_freeze_date,
                }
            )

    france_annual_sums = {
        year: sum(values)
        for year, values in FRANCE_MONTHLY_PCR_POSITIVE_TESTS.items()
    }
    if france_annual_sums != {2023: 518, 2024: 38847}:
        raise ValueError(
            "France monthly PCR-positive tests do not reproduce official "
            f"annual totals: {france_annual_sums}"
        )
    return output


def write_tsv(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDNAMES,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--monthly-input",
        type=Path,
        default=ROOT / "data/derived/country_month_cases.tsv",
    )
    parser.add_argument(
        "--japan-input",
        type=Path,
        default=DEFAULT_JAPAN_INPUT,
        help=(
            "Harmonised source table retaining Japan's weekly rows; the "
            "builder verifies that they reproduce the frozen monthly series"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "data/derived/figure1a_native_resolution_surveillance.tsv"
        ),
    )
    parser.add_argument(
        "--data-freeze-date",
        default=date.today().isoformat(),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.japan_input.is_file():
        raise FileNotFoundError(
            "Japan weekly source table was not found. Supply --japan-input "
            "with the harmonised weekly incidence table."
        )
    monthly, monthly_lookup = monthly_rows(
        args.monthly_input,
        args.data_freeze_date,
    )

    # Add the frozen monthly Japan values only to the validation lookup; they
    # are not written because Figure 1A displays the weekly source resolution.
    for row in read_delimited(args.monthly_input, "\t"):
        if (
            row["country_iso3"] == "JPN"
            and START_YEAR <= int(row["year"]) <= END_YEAR
            and row["case_data_available"].lower() == "true"
        ):
            monthly_lookup[
                ("JPN", int(row["year"]), int(row["month"]))
            ] = int(float(row["cases"]))

    rows = (
        monthly
        + digitised_monthly_rows(args.data_freeze_date)
        + japan_weekly_rows(
            args.japan_input,
            monthly_lookup,
            args.data_freeze_date,
        )
    )
    country_order = {"AUS": 0, "BEL": 1, "CHN": 2, "FRA": 3, "JPN": 4}
    rows.sort(
        key=lambda row: (
            country_order[str(row["country_iso3"])],
            str(row["observation_date"]),
        )
    )
    if len(rows) != 664:
        raise ValueError(f"Expected 664 native-resolution rows, found {len(rows)}")
    write_tsv(rows, args.output)
    print(f"wrote {len(rows)} native-resolution observations to {args.output}")


if __name__ == "__main__":
    main()
