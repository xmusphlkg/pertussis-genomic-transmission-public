#!/usr/bin/env python3
"""Fetch Belgium and France annual pertussis counts for source audit.

Figure 1A displays a separate native-resolution descriptive table.  The WHO
annual observations remain useful as a source audit, but must not be
disaggregated into artificial monthly counts or substituted for the
high-frequency surveillance measures shown in the figure.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import date
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
INDICATOR_CODE = "WHS3_43"
API_ENDPOINT = f"https://ghoapi.azureedge.net/api/{INDICATOR_CODE}"
INDICATOR_URL = (
    "https://www.who.int/data/gho/data/indicators/indicator-details/GHO/"
    "pertussis-number-of-reported-cases"
)
COUNTRIES = {"BEL": "Belgium", "FRA": "France"}
START_YEAR = 2015
END_YEAR = 2025
SURVEILLANCE_SCOPE = {
    "BEL": (
        "Country-reported annual count; Belgium uses voluntary "
        "sentinel-laboratory surveillance"
    ),
    "FRA": (
        "Country-reported annual count; France uses hospital-based "
        "sentinel surveillance"
    ),
}
COMPARABILITY_NOTE = (
    "Descriptive surveillance series only; not an independent national "
    "case-model replication and not directly comparable as population-wide "
    "incidence. Annual totals must not be disaggregated into monthly counts."
)


def build_api_url() -> str:
    query = {
        "$filter": (
            "(SpatialDim eq 'BEL' or SpatialDim eq 'FRA') "
            f"and TimeDim ge {START_YEAR} and TimeDim le {END_YEAR}"
        ),
        "$select": "SpatialDim,TimeDim,NumericValue,Value,Date",
        "$orderby": "SpatialDim,TimeDim",
        "$format": "json",
    }
    return f"{API_ENDPOINT}?{urlencode(query)}"


def fetch_payload(api_url: str) -> dict[str, object]:
    request = Request(
        api_url,
        headers={
            "Accept": "application/json",
            "User-Agent": "pertussis-genomic-transmission-figure-source/1.0",
        },
    )
    with urlopen(request, timeout=90) as response:
        return json.load(response)


def normalise_rows(
    payload: dict[str, object],
    api_url: str,
    data_freeze_date: str,
) -> list[dict[str, object]]:
    records = payload.get("value")
    if not isinstance(records, list):
        raise ValueError("WHO GHO response does not contain a list-valued 'value' field")

    by_country_year: dict[tuple[str, int], dict[str, object]] = {}
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("WHO GHO response contains a non-object record")
        country = str(record.get("SpatialDim", ""))
        year = int(record.get("TimeDim", 0))
        if country not in COUNTRIES or not START_YEAR <= year <= END_YEAR:
            raise ValueError(f"Unexpected WHO GHO record: country={country!r}, year={year}")
        key = (country, year)
        if key in by_country_year:
            raise ValueError(f"Duplicate WHO GHO record for {country} {year}")
        by_country_year[key] = record

    rows: list[dict[str, object]] = []
    for country, country_name in COUNTRIES.items():
        for year in range(START_YEAR, END_YEAR + 1):
            record = by_country_year.get((country, year), {})
            numeric_value = record.get("NumericValue")
            available = numeric_value is not None
            if available:
                numeric_float = float(numeric_value)
                if numeric_float < 0 or not numeric_float.is_integer():
                    raise ValueError(
                        f"Invalid reported-case count for {country} {year}: "
                        f"{numeric_value!r}"
                    )
                reported_cases: int | str = int(numeric_float)
            else:
                reported_cases = ""

            rows.append(
                {
                    "country_iso3": country,
                    "country_name": country_name,
                    "year": year,
                    "reported_cases": reported_cases,
                    "case_data_available": str(available),
                    "time_resolution": "annual",
                    "surveillance_scope": SURVEILLANCE_SCOPE[country],
                    "comparability_note": COMPARABILITY_NOTE,
                    "indicator_code": INDICATOR_CODE,
                    "source": (
                        "WHO/UNICEF Joint Reporting Form via WHO Global "
                        "Health Observatory"
                    ),
                    "source_url": INDICATOR_URL,
                    "source_api_url": api_url,
                    "api_last_updated": record.get("Date", ""),
                    "data_freeze_date": data_freeze_date,
                }
            )
    return rows


def write_tsv(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "data/derived/figure1a_belgium_france_annual_cases.tsv"
        ),
    )
    parser.add_argument(
        "--data-freeze-date",
        default=date.today().isoformat(),
        help="YYYY-MM-DD date recorded in the frozen output",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    api_url = build_api_url()
    payload = fetch_payload(api_url)
    rows = normalise_rows(payload, api_url, args.data_freeze_date)
    write_tsv(rows, args.output)

    summary = {
        country: {
            "available_years": sum(
                row["country_iso3"] == country
                and row["case_data_available"] == "True"
                for row in rows
            ),
            "first_year": START_YEAR,
            "last_year": END_YEAR,
        }
        for country in COUNTRIES
    }
    print(json.dumps({"output": str(args.output), "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
