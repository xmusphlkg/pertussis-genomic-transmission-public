#!/usr/bin/env python3
"""Build sampling-aware genomic transmission model tables.

This stage deliberately stops before phylogeny construction or transmission
model fitting. It freezes:

1. one row per focal or background genome candidate;
2. one row per Australian specimen in the Ct/sequencing calibration cohort;
3. one row per country-month of case surveillance;
4. one row per country-month-project-preliminary-lineage genome stratum; and
5. explicit country-level GO/NO-GO metadata gates.

The 8,117-genome compilation is tree background only. Its rows never enter the
country genome observation denominator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from calendar import monthrange
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd


FOCAL_COUNTRIES = {
    "CHN": "China",
    "AUS": "Australia",
    "JPN": "Japan",
    "BEL": "Belgium",
    "FRA": "France",
}
DEVELOPMENT_COUNTRIES = ("CHN", "AUS", "JPN")
DATA_FREEZE_DATE = "2026-07-24"
PANDEMIC_START = pd.Timestamp("2020-01-01")
POST_PANDEMIC_START = pd.Timestamp("2023-01-01")
FONG_SOURCE_DATASET = "australia_fong_2026_direct_specimens"
FONG_REDACTED_NOTE = "Third-party specimen metadata not redistributed"
FONG_RESTRICTED_PUBLIC_COLUMNS = {
    "sample_id",
    "genome_qc_status",
    "sampling_process_observed",
    "sequencing_success",
    "ct_value",
    "specimen_type",
    "preliminary_lineage_id",
    "lineage_definition_status",
    "published_branch",
    "published_lineage",
    "published_sublineage",
    "ptxP_label",
    "fim3_label",
    "marker_23s_status",
}

COUNTRY_OVERRIDES = {
    "Australia": "AUS",
    "Belgium": "BEL",
    "China": "CHN",
    "France": "FRA",
    "Japan": "JPN",
    "United States": "USA",
    "United States of America": "USA",
    "United Kingdom": "GBR",
    "New Zealand": "NZL",
    "South Korea": "KOR",
    "Republic of Korea": "KOR",
    "Taiwan": "TWN",
    "Russia": "RUS",
    "Vietnam": "VNM",
    "Iran": "IRN",
    "Czech Republic": "CZE",
    "The Netherlands": "NLD",
}

IDENTIFIER_COLUMNS = (
    "biosample_accession",
    "assembly_accession",
    "run_accession",
    "sample_id",
)

GENOME_COLUMNS = [
    "genome_record_id",
    "sample_id",
    "biosample_accession",
    "assembly_accession",
    "run_accession",
    "country_iso3",
    "country_name",
    "collection_date_raw",
    "date_lower",
    "date_upper",
    "year",
    "month",
    "quarter",
    "date_resolution",
    "date_month_or_quarter",
    "project_id",
    "project_id_source",
    "source_dataset",
    "source_role",
    "model_denominator_role",
    "genome_candidate",
    "genome_qc_status",
    "genome_observation_eligible",
    "sampling_process_observed",
    "sequencing_success",
    "ct_value",
    "specimen_type",
    "preliminary_lineage_id",
    "lineage_definition_status",
    "published_branch",
    "published_lineage",
    "published_sublineage",
    "ptxP_label",
    "fim3_label",
    "marker_23s_status",
    "duplicate_of_focal",
    "duplicate_group_size",
    "data_freeze_date",
    "provenance_note",
]


@dataclass(frozen=True)
class DateInterval:
    raw: str
    lower: str
    upper: str
    year: int | None
    month: int | None
    quarter: int | None
    resolution: str
    month_or_quarter: bool


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def clean_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    text = str(value).strip()
    if text.casefold() in {"nan", "nat", "none", "missing"}:
        return ""
    if re.fullmatch(r"\d{4}\.0", text):
        return text[:-2]
    return text


def clean_accession(value: object) -> str:
    return clean_text(value).split(".")[0] if clean_text(value) else ""


def parse_date_interval(value: object, fallback_year: object = "") -> DateInterval:
    raw = clean_text(value)
    year_fallback = clean_text(fallback_year)
    if re.fullmatch(r"\d{4}", year_fallback):
        fallback = int(year_fallback)
    else:
        fallback = None

    day_match = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})(?:[ T].*)?", raw)
    if day_match:
        year, month, day = map(int, day_match.groups())
        try:
            date = pd.Timestamp(year=year, month=month, day=day)
        except ValueError:
            date = None
        if date is not None:
            iso = date.date().isoformat()
            return DateInterval(raw, iso, iso, year, month, ((month - 1) // 3) + 1, "day", True)

    month_match = re.fullmatch(r"(\d{4})-(\d{2})", raw)
    if month_match:
        year, month = map(int, month_match.groups())
        if 1 <= month <= 12:
            lower = f"{year:04d}-{month:02d}-01"
            upper = f"{year:04d}-{month:02d}-{monthrange(year, month)[1]:02d}"
            return DateInterval(raw, lower, upper, year, month, ((month - 1) // 3) + 1, "month", True)

    quarter_match = re.fullmatch(r"(\d{4})[- ]?[Qq]([1-4])", raw)
    if quarter_match:
        year, quarter = map(int, quarter_match.groups())
        first_month = (quarter - 1) * 3 + 1
        last_month = first_month + 2
        lower = f"{year:04d}-{first_month:02d}-01"
        upper = f"{year:04d}-{last_month:02d}-{monthrange(year, last_month)[1]:02d}"
        return DateInterval(raw, lower, upper, year, None, quarter, "quarter", True)

    year_match = re.fullmatch(r"(\d{4})", raw)
    if year_match:
        year = int(year_match.group(1))
        return DateInterval(raw, f"{year:04d}-01-01", f"{year:04d}-12-31", year, None, None, "year", False)

    slash_years = [int(token) for token in re.findall(r"(?<!\d)((?:19|20)\d{2})(?!\d)", raw)]
    if slash_years:
        year = max(slash_years)
        return DateInterval(raw, f"{year:04d}-01-01", f"{year:04d}-12-31", year, None, None, "year", False)

    if fallback is not None:
        return DateInterval(raw, f"{fallback:04d}-01-01", f"{fallback:04d}-12-31", fallback, None, None, "year", False)

    return DateInterval(raw, "", "", None, None, None, "unknown", False)


def normalize_marker(value: object, prefix: str) -> str:
    text = clean_text(value)
    if not text:
        return ""
    match = re.search(rf"{re.escape(prefix)}[_ (]?([0-9]+)", text, flags=re.IGNORECASE)
    if not match:
        return text
    separator = "-" if prefix.casefold() == "fim3" else ""
    return f"{prefix}{separator}{match.group(1)}"


def marker_lineage(ptxp: object, fim3: object) -> str:
    ptxp_value = normalize_marker(ptxp, "ptxP")
    fim3_value = normalize_marker(fim3, "fim3")
    if ptxp_value and fim3_value:
        return f"{ptxp_value}/{fim3_value}"
    return ptxp_value or fim3_value


def country_to_iso3(country: object) -> str:
    text = clean_text(country)
    if not text:
        return ""
    text = text.split(":", 1)[0].strip()
    if text in COUNTRY_OVERRIDES:
        return COUNTRY_OVERRIDES[text]
    try:
        import pycountry

        return pycountry.countries.lookup(text).alpha_3
    except (ImportError, LookupError):
        return ""


def empty_genome_row() -> dict[str, object]:
    return {column: "" for column in GENOME_COLUMNS if column != "genome_record_id"}


def apply_date_fields(row: dict[str, object], value: object, fallback_year: object = "") -> None:
    parsed = parse_date_interval(value, fallback_year)
    row.update(
        {
            "collection_date_raw": parsed.raw,
            "date_lower": parsed.lower,
            "date_upper": parsed.upper,
            "year": parsed.year if parsed.year is not None else "",
            "month": parsed.month if parsed.month is not None else "",
            "quarter": parsed.quarter if parsed.quarter is not None else "",
            "date_resolution": parsed.resolution,
            "date_month_or_quarter": parsed.month_or_quarter,
        }
    )


def source_paths(root: Path) -> dict[str, Path]:
    base = root / "analysis/genomic_transmission_dynamics"
    raw = base / "inputs/raw"
    return {
        "archive": root / "analysis/prn_vaccine_evolution/derived/frozen_archive_isolates.tsv",
        "extended": root / "pertussis_data/pertussis_gene/step1_ingest/bp_extended_metadata.tsv",
        "cases": root / "pertussis_data/pertussis_gene/public_health/outputs/ph_highres_cases.tsv",
        "bel_fra": root / "analysis/post_covid_prn_turnover/derived/public_wgs_run_inventory.tsv",
        "global_8117": root / "archive/1-s2.0-S0163445326000435-mmc2.xlsx",
        "aus_workbook": raw / "australia_fong_2026_appendix2.xlsx",
        "aus_ena": raw / "australia_prjna1199062_ena_runs.tsv",
        "jpn_historical": raw / "japan_prjeb18624_ena_runs.tsv",
        "jpn_20292": raw / "japan_prjdb20292_ena_runs.tsv",
        "jpn_20413": raw / "japan_prjdb20413_ena_runs.tsv",
        "jpn_34249": raw / "japan_prjdb34249_ena_runs.tsv",
        "jpn_35593": raw / "japan_prjdb35593_ena_runs.tsv",
        "jpn_37898": raw / "japan_prjdb37898_ena_runs.tsv",
    }


def ensure_inputs(paths: dict[str, Path]) -> None:
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing required inputs:\n" + "\n".join(missing))


def load_global_8117(path: Path) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name="Global_Dataset", dtype=str).fillna("")
    if len(frame) != 8117:
        raise ValueError(f"expected 8,117 global rows, observed {len(frame):,}")
    return frame


def global_annotation_maps(global_frame: pd.DataFrame) -> dict[tuple[str, str], dict[str, str]]:
    maps: dict[tuple[str, str], dict[str, str]] = {}
    for _, source in global_frame.iterrows():
        annotation = {
            "published_branch": clean_text(source.get("Branch(I_II)", "")),
            "published_lineage": clean_text(source.get("lineage(ptxP/fim3)", "")),
            "published_sublineage": clean_text(source.get("Sublineages", "")),
            "ptxP_label": normalize_marker(source.get("ptxP", ""), "ptxP"),
            "fim3_label": normalize_marker(source.get("fim3", ""), "fim3"),
            "marker_23s_status": clean_text(source.get("23S_rRNA", "")),
        }
        identifiers = {
            "biosample_accession": clean_accession(source.get("Biosample", "")),
            "assembly_accession": clean_accession(source.get("Accession_number", "")),
            "run_accession": clean_accession(source.get("SRA-runinfo", "")),
            "sample_id": clean_text(source.get("Label", "")),
        }
        for field, value in identifiers.items():
            if value:
                maps.setdefault((field, value), annotation)
    return maps


def enrich_with_global_annotation(row: dict[str, object], maps: dict[tuple[str, str], dict[str, str]]) -> None:
    annotation: dict[str, str] | None = None
    for field in IDENTIFIER_COLUMNS:
        value = clean_accession(row.get(field, "")) if field != "sample_id" else clean_text(row.get(field, ""))
        if value and (field, value) in maps:
            annotation = maps[(field, value)]
            break
    if annotation is None:
        return
    for field, value in annotation.items():
        if not clean_text(row.get(field, "")) and value:
            row[field] = value
    if not clean_text(row.get("preliminary_lineage_id", "")):
        row["preliminary_lineage_id"] = (
            clean_text(row.get("published_sublineage", ""))
            or clean_text(row.get("published_lineage", ""))
            or marker_lineage(row.get("ptxP_label", ""), row.get("fim3_label", ""))
        )


def build_archive_rows(paths: dict[str, Path]) -> list[dict[str, object]]:
    archive = pd.read_csv(paths["archive"], sep="\t", dtype=str).fillna("")
    archive = archive[archive["country_iso3"].isin(FOCAL_COUNTRIES)].copy()
    extended = pd.read_csv(paths["extended"], sep="\t", dtype=str).fillna("")
    extended = extended[
        [
            "Assembly Accession",
            "Assembly BioProject Accession",
            "Assembly BioSample Collection date",
        ]
    ].drop_duplicates("Assembly Accession")
    archive = archive.merge(
        extended,
        how="left",
        left_on="assembly_accession",
        right_on="Assembly Accession",
    ).fillna("")

    rows: list[dict[str, object]] = []
    for _, source in archive.iterrows():
        row = empty_genome_row()
        year = clean_text(source.get("year", ""))
        date_raw = (
            clean_text(source.get("Assembly BioSample Collection date", ""))
            or clean_text(source.get("collection_date_raw", ""))
        )
        project = clean_text(source.get("Assembly BioProject Accession", ""))
        project_source = "assembly_metadata"
        if not project:
            project = clean_text(source.get("study_accession", ""))
            project_source = "study_accession"
        if not project:
            base_block = clean_text(source.get("base_block_id", ""))
            if base_block and not base_block.startswith("singleton:"):
                project = base_block
                project_source = "study_block"
        if not project:
            project = clean_text(source.get("bioproject_accession", ""))
            project_source = "archive_bioproject"

        row.update(
            {
                "sample_id": clean_text(source.get("sample_id_canonical", "")),
                "biosample_accession": clean_accession(source.get("biosample_accession", "")),
                "assembly_accession": clean_accession(source.get("assembly_accession", "")),
                "run_accession": "",
                "country_iso3": clean_text(source.get("country_iso3", "")),
                "country_name": FOCAL_COUNTRIES.get(clean_text(source.get("country_iso3", "")), ""),
                "project_id": project,
                "project_id_source": project_source if project else "unresolved",
                "source_dataset": "current_2406_archive",
                "source_role": "historical_tree_anchor" if year and float(year) <= 2019 else "focal_model_observation",
                "model_denominator_role": "focal_genome_observation",
                "genome_candidate": True,
                "genome_qc_status": "archive_analysis_genome",
                "genome_observation_eligible": True,
                "sampling_process_observed": False,
                "sequencing_success": True,
                "ct_value": clean_text(source.get("ct_or_dna_input", "")),
                "specimen_type": clean_text(source.get("specimen_type", "")),
                "published_lineage": clean_text(source.get("published_lineage_label", "")),
                "published_sublineage": clean_text(source.get("published_sublineage_label", "")),
                "ptxP_label": normalize_marker(source.get("ptxP_label", ""), "ptxP"),
                "fim3_label": normalize_marker(source.get("fim3_label", ""), "fim3"),
                "marker_23s_status": clean_text(source.get("marker_23s_status", "")),
                "lineage_definition_status": "provisional_marker_or_published_annotation",
                "duplicate_of_focal": False,
                "data_freeze_date": DATA_FREEZE_DATE,
                "provenance_note": "Legacy PRN endpoints ignored; record retained only for genome/date/project/marker metadata.",
            }
        )
        row["preliminary_lineage_id"] = (
            clean_text(row["published_sublineage"])
            or clean_text(row["published_lineage"])
            or marker_lineage(row["ptxP_label"], row["fim3_label"])
        )
        apply_date_fields(row, date_raw, year)
        rows.append(row)
    return rows


def build_public_inventory_rows(path: Path) -> list[dict[str, object]]:
    inventory = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    rows: list[dict[str, object]] = []
    for _, source in inventory.iterrows():
        country_iso3 = clean_text(source.get("country_iso3", ""))
        row = empty_genome_row()
        row.update(
            {
                "sample_id": clean_text(source.get("sample_accession", "")),
                "biosample_accession": clean_accession(source.get("sample_accession", "")),
                "assembly_accession": "",
                "run_accession": clean_accession(source.get("run_accession", "")),
                "country_iso3": country_iso3,
                "country_name": FOCAL_COUNTRIES.get(country_iso3, clean_text(source.get("country", ""))),
                "project_id": clean_text(source.get("study_accession", "")),
                "project_id_source": "ena_study_accession",
                "source_dataset": "targeted_belgium_france_public_inventory",
                "source_role": "historical_tree_anchor"
                if parse_date_interval(source.get("collection_date", "")).year
                and parse_date_interval(source.get("collection_date", "")).year <= 2019
                else "focal_model_observation",
                "model_denominator_role": "focal_genome_observation",
                "genome_candidate": True,
                "genome_qc_status": "public_raw_reads_pending_uniform_qc",
                "genome_observation_eligible": True,
                "sampling_process_observed": country_iso3 == "BEL",
                "sequencing_success": True,
                "lineage_definition_status": "pending_new_phylogeny",
                "duplicate_of_focal": False,
                "data_freeze_date": DATA_FREEZE_DATE,
                "provenance_note": clean_text(source.get("phenotype_linkage_status", "")),
            }
        )
        apply_date_fields(row, source.get("collection_date", ""))
        rows.append(row)
    return rows


def load_ena_runs(paths: Iterable[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if frame.empty:
            continue
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).drop_duplicates("run_accession")


def build_japan_rows(paths: dict[str, Path]) -> list[dict[str, object]]:
    recent_paths = [
        paths["jpn_20292"],
        paths["jpn_20413"],
        paths["jpn_34249"],
        paths["jpn_35593"],
        paths["jpn_37898"],
    ]
    historical = load_ena_runs([paths["jpn_historical"]])
    historical["source_role"] = "historical_tree_anchor"
    recent = load_ena_runs(recent_paths)
    recent["source_role"] = "focal_model_observation"
    combined = pd.concat([historical, recent], ignore_index=True).drop_duplicates("run_accession")

    rows: list[dict[str, object]] = []
    for _, source in combined.iterrows():
        row = empty_genome_row()
        row.update(
            {
                "sample_id": clean_text(source.get("sample_accession", "")),
                "biosample_accession": clean_accession(source.get("sample_accession", "")),
                "assembly_accession": "",
                "run_accession": clean_accession(source.get("run_accession", "")),
                "country_iso3": "JPN",
                "country_name": "Japan",
                "project_id": clean_text(source.get("study_accession", "")),
                "project_id_source": "ena_study_accession",
                "source_dataset": "japan_targeted_ena_projects",
                "source_role": clean_text(source.get("source_role", "")),
                "model_denominator_role": "focal_genome_observation",
                "genome_candidate": True,
                "genome_qc_status": "public_raw_reads_pending_uniform_qc",
                "genome_observation_eligible": True,
                "sampling_process_observed": False,
                "sequencing_success": True,
                "lineage_definition_status": "pending_new_phylogeny",
                "duplicate_of_focal": False,
                "data_freeze_date": DATA_FREEZE_DATE,
                "provenance_note": "Targeted Japanese historical or recent public project.",
            }
        )
        apply_date_fields(row, source.get("collection_date", ""))
        rows.append(row)
    return rows


def build_australia_tables(paths: dict[str, Path]) -> tuple[list[dict[str, object]], pd.DataFrame]:
    specimens = pd.read_excel(
        paths["aus_workbook"],
        sheet_name="Full Collection (n=255)",
        header=1,
    ).fillna("")
    if len(specimens) != 255:
        raise ValueError(f"expected 255 Australian specimens, observed {len(specimens)}")
    ena = load_ena_runs([paths["aus_ena"]])
    ena = ena.rename(
        columns={
            "run_accession": "SRA Accession",
            "sample_accession": "ena_sample_accession",
            "study_accession": "ena_study_accession",
            "collection_date": "ena_collection_date",
            "country": "ena_country",
            "first_public": "ena_first_public",
        }
    )
    ena_columns = [
        "SRA Accession",
        "ena_sample_accession",
        "ena_study_accession",
        "ena_collection_date",
        "ena_country",
        "ena_first_public",
    ]
    specimens = specimens.merge(ena[ena_columns], on="SRA Accession", how="left").fillna("")

    sampling_rows: list[dict[str, object]] = []
    genome_rows: list[dict[str, object]] = []
    for _, source in specimens.iterrows():
        run = clean_accession(source.get("SRA Accession", ""))
        profile_status = clean_text(source.get("Profile Status", ""))
        date = parse_date_interval(source.get("ena_collection_date", ""), source.get("Year", ""))
        ptxp = normalize_marker(source.get("ptxP", ""), "ptxP")
        fim3 = normalize_marker(source.get("fim3", ""), "fim3")
        raw_sequence_public = bool(run)
        profile_complete = profile_status.casefold() == "complete"
        sampling_rows.append(
            {
                "sampling_record_id": f"AUS_FONG_2026::{clean_text(source.get('Patient Number', ''))}",
                "patient_number": clean_text(source.get("Patient Number", "")),
                "probe_wgs_id": clean_text(source.get("Probe WGS ID", "")),
                "culture_wgs_id": clean_text(source.get("Culture WGS ID", "")),
                "biosample_accession": clean_accession(
                    source.get("BioSample", "") or source.get("ena_sample_accession", "")
                ),
                "run_accession": run,
                "country_iso3": "AUS",
                "state": clean_text(source.get("State", "")),
                "specimen_type": clean_text(source.get("Sample type", "")),
                "collection_date_raw": date.raw,
                "date_lower": date.lower,
                "date_upper": date.upper,
                "year": date.year if date.year is not None else "",
                "month": date.month if date.month is not None else "",
                "date_resolution": date.resolution,
                "ct_is481": clean_text(source.get("IS481", "")),
                "pcr_species_result": clean_text(source.get("PCR Result", "")),
                "mapped_coverage": clean_text(source.get("Bor. sp. Mapped coverage", "")),
                "mapped_depth": clean_text(source.get("Bor. sp. Mapped depth", "")),
                "profile_status": profile_status,
                "raw_sequence_public": raw_sequence_public,
                "profile_complete": profile_complete,
                "ptxP_label": ptxp,
                "fim3_label": fim3,
                "preliminary_lineage_id": marker_lineage(ptxp, fim3),
                "resistance": clean_text(source.get("Resistance", "")),
                "project_id": clean_text(source.get("ena_study_accession", "")) or "PRJNA1199062",
                "data_freeze_date": DATA_FREEZE_DATE,
            }
        )
        if not raw_sequence_public:
            continue
        row = empty_genome_row()
        row.update(
            {
                "sample_id": clean_text(source.get("Probe WGS ID", "")),
                "biosample_accession": clean_accession(
                    source.get("BioSample", "") or source.get("ena_sample_accession", "")
                ),
                "assembly_accession": "",
                "run_accession": run,
                "country_iso3": "AUS",
                "country_name": "Australia",
                "project_id": clean_text(source.get("ena_study_accession", "")) or "PRJNA1199062",
                "project_id_source": "ena_study_accession",
                "source_dataset": "australia_fong_2026_direct_specimens",
                "source_role": "historical_tree_anchor" if date.year and date.year <= 2019 else "focal_model_observation",
                "model_denominator_role": "focal_genome_observation",
                "genome_candidate": True,
                "genome_qc_status": f"fong_profile_{profile_status.casefold() or 'not_assigned'}",
                "genome_observation_eligible": profile_complete,
                "sampling_process_observed": True,
                "sequencing_success": profile_complete,
                "ct_value": clean_text(source.get("IS481", "")),
                "specimen_type": clean_text(source.get("Sample type", "")),
                "preliminary_lineage_id": marker_lineage(ptxp, fim3),
                "lineage_definition_status": "provisional_marker_lineage",
                "ptxP_label": ptxp,
                "fim3_label": fim3,
                "marker_23s_status": clean_text(source.get("Resistance", "")),
                "duplicate_of_focal": False,
                "data_freeze_date": DATA_FREEZE_DATE,
                "provenance_note": "Direct-specimen sequencing cohort; Ct and failures retained in sampling-process table.",
            }
        )
        apply_date_fields(row, source.get("ena_collection_date", ""), source.get("Year", ""))
        genome_rows.append(row)
    return genome_rows, pd.DataFrame(sampling_rows)


def identifier_keys(row: dict[str, object] | pd.Series) -> list[tuple[str, str]]:
    keys: list[tuple[str, str]] = []
    for field in IDENTIFIER_COLUMNS:
        raw = row.get(field, "")
        value = clean_accession(raw) if field != "sample_id" else clean_text(raw)
        if value:
            keys.append((field, value))
    return keys


def merge_focal_duplicates(rows: list[dict[str, object]]) -> pd.DataFrame:
    parent = list(range(len(rows)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    seen: dict[tuple[str, str], int] = {}
    for index, row in enumerate(rows):
        for key in identifier_keys(row):
            if key in seen:
                union(index, seen[key])
            else:
                seen[key] = index

    groups: dict[int, list[int]] = {}
    for index in range(len(rows)):
        groups.setdefault(find(index), []).append(index)

    priority = {
        "australia_fong_2026_direct_specimens": 100,
        "targeted_belgium_france_public_inventory": 90,
        "japan_targeted_ena_projects": 80,
        "current_2406_archive": 70,
    }
    merged_rows: list[dict[str, object]] = []
    for indices in groups.values():
        ordered = sorted(indices, key=lambda item: priority.get(clean_text(rows[item].get("source_dataset", "")), 0), reverse=True)
        primary = dict(rows[ordered[0]])
        sources = []
        notes = []
        for index in ordered:
            candidate = rows[index]
            source_name = clean_text(candidate.get("source_dataset", ""))
            if source_name and source_name not in sources:
                sources.append(source_name)
            note = clean_text(candidate.get("provenance_note", ""))
            if note and note not in notes:
                notes.append(note)
            for column in GENOME_COLUMNS:
                if column in {"source_dataset", "provenance_note", "duplicate_group_size", "genome_record_id"}:
                    continue
                if not clean_text(primary.get(column, "")) and clean_text(candidate.get(column, "")):
                    primary[column] = candidate[column]
            for boolean_field in (
                "genome_candidate",
                "genome_observation_eligible",
                "sampling_process_observed",
                "sequencing_success",
                "date_month_or_quarter",
            ):
                primary[boolean_field] = any(bool(rows[item].get(boolean_field, False)) for item in ordered)
        primary["source_dataset"] = ";".join(sources)
        primary["provenance_note"] = " | ".join(notes)
        primary["duplicate_group_size"] = len(indices)
        merged_rows.append(primary)

    for row in merged_rows:
        preferred = (
            clean_text(row.get("biosample_accession", ""))
            or clean_text(row.get("run_accession", ""))
            or clean_text(row.get("assembly_accession", ""))
            or clean_text(row.get("sample_id", ""))
        )
        row["genome_record_id"] = f"{clean_text(row.get('country_iso3', 'UNK'))}::{preferred}"
    frame = pd.DataFrame(merged_rows)
    if frame["genome_record_id"].duplicated().any():
        duplicates = frame.loc[frame["genome_record_id"].duplicated(False), "genome_record_id"].tolist()
        raise ValueError(f"duplicate focal genome_record_id values: {duplicates[:5]}")
    return frame[GENOME_COLUMNS]


def sanitise_fong_public_records(frame: pd.DataFrame) -> pd.DataFrame:
    released = frame.copy()
    if "source_dataset" not in released.columns:
        return released
    mask = released["source_dataset"].eq(FONG_SOURCE_DATASET)
    for column in FONG_RESTRICTED_PUBLIC_COLUMNS.intersection(released.columns):
        released.loc[mask, column] = ""
    if "provenance_note" in released.columns:
        released.loc[mask, "provenance_note"] = FONG_REDACTED_NOTE
    return released


def build_global_background_rows(
    global_frame: pd.DataFrame,
    focal: pd.DataFrame,
) -> pd.DataFrame:
    focal_keys = {key for _, row in focal.iterrows() for key in identifier_keys(row)}
    rows: list[dict[str, object]] = []
    for _, source in global_frame.iterrows():
        row = empty_genome_row()
        identifiers = {
            "sample_id": clean_text(source.get("Label", "")),
            "biosample_accession": clean_accession(source.get("Biosample", "")),
            "assembly_accession": clean_accession(source.get("Accession_number", "")),
            "run_accession": clean_accession(source.get("SRA-runinfo", "")),
        }
        duplicate = any(key in focal_keys for key in identifier_keys(identifiers))
        country_name = clean_text(source.get("Country", ""))
        country_iso3 = country_to_iso3(country_name)
        ptxp = normalize_marker(source.get("ptxP", ""), "ptxP")
        fim3 = normalize_marker(source.get("fim3", ""), "fim3")
        row.update(
            {
                **identifiers,
                "country_iso3": country_iso3,
                "country_name": country_name,
                "project_id": "",
                "project_id_source": "not_available_in_global_compilation",
                "source_dataset": "global_8117_compilation",
                "source_role": "global_tree_background",
                "model_denominator_role": "global_tree_background_only",
                "genome_candidate": True,
                "genome_qc_status": "published_global_compilation_qc",
                "genome_observation_eligible": False,
                "sampling_process_observed": False,
                "sequencing_success": True,
                "preliminary_lineage_id": clean_text(source.get("Sublineages", ""))
                or clean_text(source.get("lineage(ptxP/fim3)", ""))
                or marker_lineage(ptxp, fim3),
                "lineage_definition_status": "published_background_annotation_not_model_lineage",
                "published_branch": clean_text(source.get("Branch(I_II)", "")),
                "published_lineage": clean_text(source.get("lineage(ptxP/fim3)", "")),
                "published_sublineage": clean_text(source.get("Sublineages", "")),
                "ptxP_label": ptxp,
                "fim3_label": fim3,
                "marker_23s_status": clean_text(source.get("23S_rRNA", "")),
                "duplicate_of_focal": duplicate,
                "duplicate_group_size": 1,
                "data_freeze_date": DATA_FREEZE_DATE,
                "provenance_note": "Excluded from all national case/genome observation denominators; eligible only for stratified tree background.",
            }
        )
        apply_date_fields(row, source.get("Collection_date", ""))
        preferred = (
            clean_text(row["biosample_accession"])
            or clean_text(row["run_accession"])
            or clean_text(row["assembly_accession"])
            or clean_text(row["sample_id"])
        )
        row["genome_record_id"] = f"BG::{preferred}"
        rows.append(row)
    background = pd.DataFrame(rows)[GENOME_COLUMNS]
    return background


def build_case_months(path: Path) -> pd.DataFrame:
    cases = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    cases = cases[cases["country_iso3"].isin(DEVELOPMENT_COUNTRIES)].copy()
    cases["date"] = pd.to_datetime(cases["date"])
    cases["cases"] = pd.to_numeric(cases["cases"], errors="coerce")
    cases["model_month"] = cases["date"].dt.to_period("M").dt.to_timestamp()

    metadata_columns = [
        "country_name",
        "source_url",
        "source_file",
        "reporting_era_record_iso3",
        "reporting_era_scope_type",
        "reporting_era_match_type",
        "reporting_era_confidence",
        "pcr_lab_guideline_year",
        "reporting_case_definition_change_year",
        "surveillance_platform_change_year",
    ]
    aggregations: dict[str, object] = {
        "cases": "sum",
        "date": "size",
        "time_resolution": lambda values: ";".join(sorted(set(values))),
    }
    aggregations.update({column: "first" for column in metadata_columns})
    monthly = (
        cases.groupby(["country_iso3", "model_month"], as_index=False)
        .agg(aggregations)
        .rename(columns={"date": "n_source_intervals", "time_resolution": "source_time_resolution"})
    )

    skeleton_rows = []
    for country_iso3, country_name in FOCAL_COUNTRIES.items():
        for month in pd.date_range("2015-01-01", "2025-12-01", freq="MS"):
            skeleton_rows.append(
                {
                    "country_iso3": country_iso3,
                    "country_name_skeleton": country_name,
                    "model_month": month,
                }
            )
    skeleton = pd.DataFrame(skeleton_rows)
    monthly = skeleton.merge(monthly, on=["country_iso3", "model_month"], how="left")
    monthly["country_name"] = monthly["country_name"].fillna(monthly["country_name_skeleton"])
    monthly = monthly.drop(columns="country_name_skeleton")
    monthly["case_data_available"] = monthly["cases"].notna()
    monthly["year"] = monthly["model_month"].dt.year
    monthly["month"] = monthly["model_month"].dt.month
    monthly["period"] = "post_pandemic_resurgence"
    monthly.loc[monthly["model_month"] < POST_PANDEMIC_START, "period"] = "pandemic_interruption"
    monthly.loc[monthly["model_month"] < PANDEMIC_START, "period"] = "pre_pandemic"
    monthly["model_month"] = monthly["model_month"].dt.date.astype(str)
    monthly["data_freeze_date"] = DATA_FREEZE_DATE
    ordered = [
        "country_iso3",
        "country_name",
        "model_month",
        "year",
        "month",
        "period",
        "cases",
        "case_data_available",
        "n_source_intervals",
        "source_time_resolution",
        *metadata_columns[1:],
        "data_freeze_date",
    ]
    return monthly[ordered].sort_values(["country_iso3", "model_month"])


def build_genome_strata(focal: pd.DataFrame) -> pd.DataFrame:
    eligible = focal[
        focal["model_denominator_role"].eq("focal_genome_observation")
        & focal["source_role"].eq("focal_model_observation")
        & focal["date_month_or_quarter"].astype(bool)
    ].copy()
    eligible["model_month"] = pd.to_datetime(eligible["date_lower"]).dt.to_period("M").dt.to_timestamp()
    eligible["project_id"] = eligible["project_id"].replace("", "UNRESOLVED_PROJECT")
    eligible["preliminary_lineage_id"] = eligible["preliminary_lineage_id"].replace("", "PENDING_LINEAGE_TYPING")
    grouped = (
        eligible.groupby(
            ["country_iso3", "model_month", "project_id", "preliminary_lineage_id"],
            as_index=False,
        )
        .agg(
            n_genome_candidates=("genome_record_id", "size"),
            n_genome_observation_eligible=("genome_observation_eligible", "sum"),
            n_sampling_process_observed=("sampling_process_observed", "sum"),
        )
    )
    grouped["model_month"] = grouped["model_month"].dt.date.astype(str)
    grouped["lineage_definition_status"] = "preliminary_only_pending_new_phylogeny"
    grouped["data_freeze_date"] = DATA_FREEZE_DATE
    return grouped.sort_values(
        ["country_iso3", "model_month", "project_id", "preliminary_lineage_id"]
    )


def build_integrated_panel(case_months: pd.DataFrame, strata: pd.DataFrame) -> pd.DataFrame:
    panel = case_months.merge(strata, on=["country_iso3", "model_month"], how="left", suffixes=("", "_genome"))
    no_genome = panel["project_id"].isna()
    panel.loc[no_genome, "project_id"] = "NO_GENOME_OBSERVED"
    panel.loc[no_genome, "preliminary_lineage_id"] = "NO_GENOME_OBSERVED"
    panel.loc[no_genome, "lineage_definition_status"] = "no_genome_observed"
    for column in (
        "n_genome_candidates",
        "n_genome_observation_eligible",
        "n_sampling_process_observed",
    ):
        panel[column] = panel[column].fillna(0).astype(int)
    panel["n_genome_strata_in_country_month"] = panel.groupby(
        ["country_iso3", "model_month"]
    )["project_id"].transform("size")
    panel["case_value_repeated_across_genome_strata"] = (
        panel["case_data_available"].astype(bool)
        & panel["n_genome_strata_in_country_month"].gt(1)
    )
    panel["model_lineage_id"] = ""
    panel["model_lineage_status"] = "LOCKED_PENDING_NEW_PHYLOGENY"
    panel["data_freeze_date"] = DATA_FREEZE_DATE
    return panel.sort_values(
        ["country_iso3", "model_month", "project_id", "preliminary_lineage_id"]
    )


def build_country_gates(focal: pd.DataFrame, case_months: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for country_iso3 in DEVELOPMENT_COUNTRIES:
        genomes = focal[
            focal["country_iso3"].eq(country_iso3)
            & focal["model_denominator_role"].eq("focal_genome_observation")
        ].copy()
        genomes["year_numeric"] = pd.to_numeric(genomes["year"], errors="coerce")
        recent = genomes[genomes["year_numeric"].ge(2023)].copy()
        historical = genomes[genomes["year_numeric"].le(2019)].copy()
        cases = case_months[
            case_months["country_iso3"].eq(country_iso3)
            & case_months["case_data_available"].astype(bool)
        ].copy()
        pre_case_months = int(sum(pd.to_datetime(cases["model_month"]) < PANDEMIC_START))
        post_case_months = int(sum(pd.to_datetime(cases["model_month"]) >= POST_PANDEMIC_START))
        precise_recent = int(recent["date_month_or_quarter"].astype(bool).sum())
        recent_n = len(recent)
        precise_fraction = precise_recent / recent_n if recent_n else 0.0
        resolved_project = recent["project_id"].fillna("").ne("")
        project_fraction = float(resolved_project.mean()) if recent_n else 0.0
        project_counts = recent.loc[resolved_project, "project_id"].value_counts()
        dominant_project = clean_text(project_counts.index[0]) if len(project_counts) else ""
        dominant_share = float(project_counts.iloc[0] / recent_n) if recent_n else 0.0
        sample_gate = recent_n >= 30 and len(historical) >= 50
        date_gate = precise_fraction >= 0.70
        case_gate = pre_case_months >= 24 and post_case_months >= 18
        project_gate = project_fraction >= 0.90
        metadata_gate = sample_gate and date_gate and case_gate and project_gate
        sampling_status = "pass"
        if not project_gate:
            sampling_status = "fail_project_resolution"
        elif dominant_share >= 0.70:
            sampling_status = "pass_with_project_dominance_sensitivity_required"
        rows.append(
            {
                "country_iso3": country_iso3,
                "country_name": FOCAL_COUNTRIES[country_iso3],
                "n_focal_genome_candidates": len(genomes),
                "n_historical_genomes_2019_or_earlier": len(historical),
                "n_recent_genomes_2023_or_later": recent_n,
                "n_recent_genomes_observation_eligible": int(
                    recent["genome_observation_eligible"].astype(bool).sum()
                ),
                "n_recent_month_or_quarter_dates": precise_recent,
                "recent_month_or_quarter_date_fraction": round(precise_fraction, 4),
                "recent_project_id_fraction": round(project_fraction, 4),
                "dominant_recent_project_id": dominant_project,
                "dominant_recent_project_share": round(dominant_share, 4),
                "n_pre_pandemic_case_months": pre_case_months,
                "n_post_pandemic_case_months": post_case_months,
                "sample_size_gate_pass": sample_gate,
                "date_resolution_gate_pass": date_gate,
                "case_time_coverage_gate_pass": case_gate,
                "project_encoding_gate_pass": project_gate,
                "sampling_gate_status": sampling_status,
                "temporal_metadata_gate_pass": metadata_gate,
                "phylogenetic_clock_gate_status": "PENDING_NEW_TRANSMISSION_TREE_AND_DATE_RANDOMISATION",
                "formal_lineage_gate_status": "PENDING_PREDEFINED_PHYLOGENETIC_LINEAGES",
                "development_gate_status": (
                    "PROVISIONAL_GO_METADATA_ONLY"
                    if metadata_gate
                    else "NO_GO_UNTIL_METADATA_GAPS_RESOLVED"
                ),
                "data_freeze_date": DATA_FREEZE_DATE,
            }
        )
    return pd.DataFrame(rows)


def build_lineage_screen(focal: pd.DataFrame) -> pd.DataFrame:
    frame = focal[
        focal["model_denominator_role"].eq("focal_genome_observation")
        & focal["preliminary_lineage_id"].fillna("").ne("")
    ].copy()
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "preliminary_lineage_id",
                "n_genomes",
                "n_countries",
                "countries",
                "n_periods",
                "periods",
                "screening_minimum_20_pass",
                "cross_period_or_country_pass",
                "formal_model_lineage_eligible",
                "reason",
            ]
        )
    frame["year_numeric"] = pd.to_numeric(frame["year"], errors="coerce")
    frame["period"] = "post_pandemic_resurgence"
    frame.loc[frame["year_numeric"] < 2023, "period"] = "pandemic_interruption"
    frame.loc[frame["year_numeric"] < 2020, "period"] = "pre_pandemic"
    rows = []
    for lineage, group in frame.groupby("preliminary_lineage_id"):
        countries = sorted(set(group["country_iso3"]) - {""})
        periods = sorted(set(group["period"]) - {""})
        minimum = len(group) >= 20
        cross = len(countries) >= 2 or len(periods) >= 2
        rows.append(
            {
                "preliminary_lineage_id": lineage,
                "n_genomes": len(group),
                "n_countries": len(countries),
                "countries": ";".join(countries),
                "n_periods": len(periods),
                "periods": ";".join(periods),
                "screening_minimum_20_pass": minimum,
                "cross_period_or_country_pass": cross,
                "formal_model_lineage_eligible": False,
                "reason": "Preliminary published/marker label only; model lineage remains locked until the new phylogeny is defined without outcome peeking.",
            }
        )
    return pd.DataFrame(rows).sort_values(["n_genomes", "preliminary_lineage_id"], ascending=[False, True])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_source_registry(paths: dict[str, Path], root: Path) -> pd.DataFrame:
    urls = {
        "aus_workbook": "https://ars.els-cdn.com/content/image/1-s2.0-S2666524725002149-mmc2.xlsx",
        "aus_ena": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1199062",
        "jpn_historical": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB18624",
        "jpn_20292": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB20292",
        "jpn_20413": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB20413",
        "jpn_34249": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB34249",
        "jpn_35593": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB35593",
        "jpn_37898": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB37898",
    }
    roles = {
        "global_8117": "global_tree_background_only",
        "cases": "national_case_observation",
        "archive": "focal_genome_seed",
        "extended": "project_and_date_enrichment",
        "bel_fra": "focal_genome_seed",
        "aus_workbook": "sampling_process_calibration_and_focal_genomes",
        "aus_ena": "date_and_project_enrichment",
        "jpn_historical": "historical_tree_anchor",
        "jpn_20292": "recent_focal_genomes",
        "jpn_20413": "recent_focal_genomes",
        "jpn_34249": "recent_focal_genomes",
        "jpn_35593": "recent_focal_genomes",
        "jpn_37898": "recent_focal_genomes",
    }
    rows = []
    for source_id, path in paths.items():
        try:
            display_path = str(path.relative_to(root))
        except ValueError:
            display_path = str(path)
        rows.append(
            {
                "source_id": source_id,
                "path": display_path,
                "source_url": urls.get(source_id, ""),
                "role": roles.get(source_id, ""),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
                "data_freeze_date": DATA_FREEZE_DATE,
            }
        )
    return pd.DataFrame(rows)


def validate_outputs(
    focal: pd.DataFrame,
    background: pd.DataFrame,
    sampling: pd.DataFrame,
    cases: pd.DataFrame,
    gates: pd.DataFrame,
) -> dict[str, object]:
    checks = {
        "focal_genome_record_id_unique": bool(not focal["genome_record_id"].duplicated().any()),
        "global_8117_rows_preserved": bool(len(background) == 8117),
        "global_background_never_in_model_denominator": bool(
            background["genome_observation_eligible"].eq(False).all()
        ),
        "australia_sampling_rows_255": bool(len(sampling) == 255),
        "australia_public_sequence_rows_154": bool(
            int(sampling["raw_sequence_public"].astype(bool).sum()) == 154
        ),
        "case_country_month_key_unique": bool(
            not cases[["country_iso3", "model_month"]].duplicated().any()
        ),
        "three_development_countries_audited": bool(
            set(gates["country_iso3"]) == set(DEVELOPMENT_COUNTRIES)
        ),
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ValueError("output validation failed: " + ", ".join(failed))
    return {
        "status": "PASS",
        "checks": checks,
        "n_focal_genomes": len(focal),
        "n_global_background_rows": len(background),
        "n_global_background_duplicates_of_focal": int(background["duplicate_of_focal"].astype(bool).sum()),
        "n_australia_sampling_records": len(sampling),
        "n_case_month_rows": len(cases),
        "development_gate_status": dict(zip(gates["country_iso3"], gates["development_gate_status"])),
        "hard_stops_remaining": [
            "No formal transmission lineage is released until the new tree is built and lineages are locked.",
            "No molecular-clock inference is released until root-to-tip and date-randomisation checks pass.",
            "No joint import-versus-growth model is fit until identifiability simulations recover both quantities.",
        ],
        "data_freeze_date": DATA_FREEZE_DATE,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root() / "analysis/genomic_transmission_dynamics/derived",
    )
    parser.add_argument(
        "--audit-dir",
        type=Path,
        default=repo_root() / "analysis/genomic_transmission_dynamics/audit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = repo_root()
    paths = source_paths(root)
    ensure_inputs(paths)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.audit_dir.mkdir(parents=True, exist_ok=True)

    global_frame = load_global_8117(paths["global_8117"])
    annotation_maps = global_annotation_maps(global_frame)

    focal_rows = build_archive_rows(paths)
    focal_rows.extend(build_public_inventory_rows(paths["bel_fra"]))
    focal_rows.extend(build_japan_rows(paths))
    australia_rows, sampling = build_australia_tables(paths)
    focal_rows.extend(australia_rows)
    for row in focal_rows:
        enrich_with_global_annotation(row, annotation_maps)

    focal = merge_focal_duplicates(focal_rows)
    background = build_global_background_rows(global_frame, focal)
    cases = build_case_months(paths["cases"])
    strata = build_genome_strata(focal)
    panel = build_integrated_panel(cases, strata)
    gates = build_country_gates(focal, cases)
    lineage_screen = build_lineage_screen(focal)
    registry = build_source_registry(paths, root)

    released_focal = sanitise_fong_public_records(focal)
    released_focal.to_csv(args.output_dir / "transmission_genome_records.tsv", sep="\t", index=False)
    background.to_csv(args.output_dir / "global_tree_background_records.tsv", sep="\t", index=False)
    cases.to_csv(args.output_dir / "country_month_cases.tsv", sep="\t", index=False)
    strata.to_csv(args.output_dir / "country_month_genome_strata.tsv", sep="\t", index=False)
    panel.to_csv(args.output_dir / "country_month_case_genome_panel.tsv", sep="\t", index=False)
    gates.to_csv(args.output_dir / "country_data_gate_summary.tsv", sep="\t", index=False)
    lineage_screen.to_csv(args.output_dir / "preliminary_lineage_screen.tsv", sep="\t", index=False)
    registry.to_csv(args.output_dir / "source_registry.tsv", sep="\t", index=False)

    report = validate_outputs(focal, background, sampling, cases, gates)
    with (args.audit_dir / "build_validation.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
