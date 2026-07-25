#!/usr/bin/env python3
"""Freeze focal and stratified background cohorts for the transmission tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


FOCAL_COUNTRIES = {"CHN", "AUS", "JPN", "BEL", "FRA"}
CONTINENT_BY_FOCAL = {
    "CHN": "Asia",
    "AUS": "Oceania",
    "JPN": "Asia",
    "BEL": "Europe",
    "FRA": "Europe",
}
DATA_FREEZE_DATE = "2026-07-24"
BACKGROUND_REPLICATES = ("primary", "replicate_b", "replicate_c")


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    value_text = str(value).strip()
    return "" if value_text.casefold() in {"nan", "none", "nat"} else value_text


def accession(value: object) -> str:
    value_text = text(value)
    return value_text.split(".")[0] if value_text else ""


def stable_hash(seed: str, value: str) -> str:
    return hashlib.sha256(f"{seed}|{value}".encode()).hexdigest()


def period_bin(year: object) -> str:
    numeric = pd.to_numeric(pd.Series([year]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "unknown"
    value = int(numeric)
    if value < 2000:
        return "pre2000"
    if value <= 2009:
        return "2000_2009"
    if value <= 2014:
        return "2010_2014"
    if value <= 2019:
        return "2015_2019"
    if value <= 2022:
        return "2020_2022"
    if value <= 2025:
        return "2023_2025"
    return "post2025"


def five_year_bin(year: object) -> str:
    numeric = pd.to_numeric(pd.Series([year]), errors="coerce").iloc[0]
    if pd.isna(numeric):
        return "unknown"
    start = (int(numeric) // 5) * 5
    return f"{start}_{start + 4}"


def modern_time_bin(row: pd.Series) -> str:
    year = pd.to_numeric(pd.Series([row.get("year", "")]), errors="coerce").iloc[0]
    if pd.isna(year):
        return "unknown"
    if text(row.get("date_resolution", "")) in {"day", "month", "quarter"}:
        lower = pd.to_datetime(row.get("date_lower", ""), errors="coerce")
        if not pd.isna(lower):
            return lower.to_period("M").strftime("%Y-%m")
    return str(int(year))


def safe_tree_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")
    return clean[:80] or "unnamed"


def build_catalog_maps(catalog: pd.DataFrame) -> dict[str, dict[str, dict[str, str]]]:
    maps: dict[str, dict[str, dict[str, str]]] = {
        "biosample_accession": {},
        "assembly_accession": {},
        "sample_id": {},
    }
    for _, row in catalog.iterrows():
        payload = {
            "catalog_assembly_accession": text(row.get("assembly_accession", "")),
            "catalog_run_accession": text(row.get("ena_run_accession", ""))
            or text(row.get("sra_run_accession", "")),
            "local_fasta_path": text(row.get("primary_fasta_path", "")),
            "total_sequence_length": text(row.get("total_sequence_length", "")),
            "gc_percent": text(row.get("gc_percent", "")),
            "n_contigs": text(row.get("n_contigs", "")),
            "contig_n50": text(row.get("contig_n50", "")),
            "catalog_qc_status": text(row.get("primary_fasta_status", "")),
        }
        values = {
            "biosample_accession": accession(row.get("biosample_accession", "")),
            "assembly_accession": accession(row.get("assembly_accession", "")),
            "sample_id": text(row.get("sample_id_canonical", "")),
        }
        for field, value in values.items():
            if value:
                maps[field].setdefault(value, payload)
    return maps


def annotate_local_sequence(frame: pd.DataFrame, maps: dict[str, dict[str, dict[str, str]]]) -> pd.DataFrame:
    payload_rows = []
    for _, row in frame.iterrows():
        payload: dict[str, str] = {}
        candidates = [
            ("biosample_accession", accession(row.get("biosample_accession", ""))),
            ("assembly_accession", accession(row.get("assembly_accession", ""))),
            ("sample_id", text(row.get("sample_id", ""))),
        ]
        for field, value in candidates:
            if value and value in maps[field]:
                payload = maps[field][value]
                break
        payload_rows.append(payload)
    payload_frame = pd.DataFrame(payload_rows)
    for target_column, source_column in (
        ("assembly_accession", "catalog_assembly_accession"),
        ("run_accession", "catalog_run_accession"),
    ):
        catalog_values = payload_frame.get(
            source_column, pd.Series([""] * len(frame))
        ).fillna("")
        missing = frame[target_column].map(text).eq("")
        frame.loc[missing, target_column] = catalog_values.loc[missing].values
    for column in (
        "local_fasta_path",
        "total_sequence_length",
        "gc_percent",
        "n_contigs",
        "contig_n50",
        "catalog_qc_status",
    ):
        frame[column] = payload_frame.get(column, pd.Series([""] * len(frame))).fillna("").values
    return frame


def global_continent_map(path: Path) -> dict[tuple[str, str], str]:
    global_frame = pd.read_excel(path, sheet_name="Global_Dataset", dtype=str).fillna("")
    mapping: dict[tuple[str, str], str] = {}
    source_columns = {
        "biosample_accession": "Biosample",
        "assembly_accession": "Accession_number",
        "run_accession": "SRA-runinfo",
        "sample_id": "Label",
    }
    for _, row in global_frame.iterrows():
        continent = text(row.get("Continent", "")) or "unknown"
        for target, source in source_columns.items():
            value = text(row.get(source, ""))
            value = accession(value) if target != "sample_id" else value
            if value:
                mapping.setdefault((target, value), continent)
    return mapping


def annotate_continent(frame: pd.DataFrame, mapping: dict[tuple[str, str], str]) -> pd.DataFrame:
    continents = []
    for _, row in frame.iterrows():
        country = text(row.get("country_iso3", ""))
        if country in CONTINENT_BY_FOCAL:
            continents.append(CONTINENT_BY_FOCAL[country])
            continue
        continent = ""
        for field in ("biosample_accession", "assembly_accession", "run_accession", "sample_id"):
            value = text(row.get(field, ""))
            value = accession(value) if field != "sample_id" else value
            if value and (field, value) in mapping:
                continent = mapping[(field, value)]
                break
        continents.append(continent or "unknown")
    frame["continent"] = continents
    return frame


def add_acquisition_route(frame: pd.DataFrame) -> pd.DataFrame:
    frame["sequence_acquisition"] = "unresolved"
    frame.loc[frame["local_fasta_path"].ne(""), "sequence_acquisition"] = "local_fasta"
    frame.loc[
        frame["local_fasta_path"].eq("") & frame["assembly_accession"].ne(""),
        "sequence_acquisition",
    ] = "ncbi_assembly_download"
    frame.loc[
        frame["local_fasta_path"].eq("")
        & frame["assembly_accession"].eq("")
        & frame["run_accession"].ne(""),
        "sequence_acquisition",
    ] = "ena_fastq_download"
    return frame


def select_focal(
    focal: pd.DataFrame,
    seed: str = "primary",
    historical_cell_cap: int = 2,
    historical_country_cap: int = 60,
    historical_all_countries: set[str] | None = None,
) -> pd.DataFrame:
    focal = focal.copy()
    focal["year_numeric"] = pd.to_numeric(focal["year"], errors="coerce")
    focal["lineage_stratum"] = focal["preliminary_lineage_id"].replace("", "UNKNOWN_LINEAGE")
    focal["time_stratum"] = focal.apply(modern_time_bin, axis=1)
    historical_all = focal[focal["year_numeric"].le(2019)].copy()
    rescue_countries = historical_all_countries or set()
    historical_rescue = historical_all[
        historical_all["country_iso3"].isin(rescue_countries)
    ].copy()
    historical = historical_all[
        ~historical_all["country_iso3"].isin(rescue_countries)
    ].copy()
    historical["time_stratum"] = historical["year"].map(five_year_bin)
    historical["rank_hash"] = historical["genome_record_id"].map(lambda value: stable_hash(seed, value))
    historical = historical.sort_values("rank_hash")
    historical = historical.groupby(
        ["country_iso3", "time_stratum", "project_id", "lineage_stratum"],
        dropna=False,
        group_keys=False,
    ).head(historical_cell_cap)
    historical = (
        historical.sort_values("rank_hash")
        .groupby("country_iso3", group_keys=False)
        .head(historical_country_cap)
    )
    historical["selection_reason"] = (
        f"historical_stratified_anchor_cap{historical_cell_cap}_cell_"
        f"cap{historical_country_cap}_country"
    )
    if not historical_rescue.empty:
        historical_rescue["rank_hash"] = historical_rescue["genome_record_id"].map(
            lambda value: stable_hash(seed, value)
        )
        historical_rescue["selection_reason"] = (
            "historical_all_available_gate_rescue"
        )
        historical = pd.concat([historical, historical_rescue], ignore_index=True)

    modern = focal[focal["year_numeric"].ge(2020)].copy()
    modern["rank_hash"] = modern["genome_record_id"].map(lambda value: stable_hash(seed, value))
    modern = modern.sort_values("rank_hash")
    modern = modern.groupby(
        ["country_iso3", "time_stratum", "project_id", "lineage_stratum"],
        dropna=False,
        group_keys=False,
    ).head(8)
    modern["selection_reason"] = "modern_country_time_project_lineage_cap8"
    selected = pd.concat([historical, modern], ignore_index=True)
    return selected.drop_duplicates("genome_record_id")


def select_background(background: pd.DataFrame, replicate: str) -> pd.DataFrame:
    eligible = background[
        ~background["duplicate_of_focal"].astype(str).str.casefold().eq("true")
        & ~background["country_iso3"].isin(FOCAL_COUNTRIES)
        & background["country_iso3"].ne("")
        & background["sequence_acquisition"].ne("unresolved")
    ].copy()
    eligible["period_stratum"] = eligible["year"].map(period_bin)
    eligible["lineage_stratum"] = eligible["preliminary_lineage_id"].replace("", "UNKNOWN_LINEAGE")
    eligible["rank_hash"] = eligible["genome_record_id"].map(
        lambda value: stable_hash(replicate, value)
    )
    route_priority = {
        "local_fasta": 0,
        "ncbi_assembly_download": 1,
        "ena_fastq_download": 2,
    }
    eligible["route_priority"] = eligible["sequence_acquisition"].map(route_priority).fillna(9)
    eligible["date_priority"] = ~eligible["date_resolution"].isin({"day", "month", "quarter"})
    eligible = eligible.sort_values(["route_priority", "date_priority", "rank_hash"])
    selected = eligible.groupby(
        ["country_iso3", "period_stratum", "lineage_stratum"],
        dropna=False,
        group_keys=False,
    ).head(2).copy()
    selected["selection_reason"] = (
        "global_background_country_period_lineage_cap2_prefer_local_and_precise_date"
    )
    selected["background_replicate"] = replicate
    return selected


def finalise_manifest(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["sample_key"] = frame.apply(
        lambda row: text(row.get("biosample_accession", ""))
        or text(row.get("run_accession", ""))
        or text(row.get("assembly_accession", ""))
        or text(row.get("sample_id", ""))
        or text(row.get("genome_record_id", "")),
        axis=1,
    )
    frame["tree_sample_id"] = frame.apply(
        lambda row: safe_tree_name(
            f"{text(row.get('country_iso3', 'UNK'))}__{text(row.get('sample_key', ''))}"
        ),
        axis=1,
    )
    duplicates = frame["tree_sample_id"].duplicated(keep=False)
    if duplicates.any():
        frame.loc[duplicates, "tree_sample_id"] = frame.loc[duplicates].apply(
            lambda row: safe_tree_name(
                f"{row['tree_sample_id']}__{stable_hash('tree-id', row['genome_record_id'])[:8]}"
            ),
            axis=1,
        )
    frame["tree_include"] = True
    frame["uniform_qc_status"] = "PENDING"
    frame["fastq_r1"] = ""
    frame["fastq_r2"] = ""
    frame["sequence_input_path"] = frame["local_fasta_path"]
    frame["data_freeze_date"] = DATA_FREEZE_DATE
    columns = [
        "tree_sample_id",
        "genome_record_id",
        "tree_role",
        "selection_reason",
        "background_replicate",
        "country_iso3",
        "country_name",
        "continent",
        "date_lower",
        "date_upper",
        "date_resolution",
        "year",
        "month",
        "project_id",
        "preliminary_lineage_id",
        "lineage_definition_status",
        "sample_id",
        "biosample_accession",
        "assembly_accession",
        "run_accession",
        "sequence_acquisition",
        "local_fasta_path",
        "fastq_r1",
        "fastq_r2",
        "sequence_input_path",
        "total_sequence_length",
        "gc_percent",
        "n_contigs",
        "contig_n50",
        "catalog_qc_status",
        "uniform_qc_status",
        "tree_include",
        "data_freeze_date",
    ]
    for column in columns:
        if column not in frame:
            frame[column] = ""
    return frame[columns].sort_values(["tree_role", "country_iso3", "year", "tree_sample_id"])


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "analysis/genomic_transmission_dynamics/phylogeny/cohort",
    )
    parser.add_argument("--historical-cell-cap", type=int, default=2)
    parser.add_argument("--historical-country-cap", type=int, default=60)
    parser.add_argument(
        "--historical-all-countries",
        default="",
        help="Comma-separated countries for which all pre-2020 candidates are retained",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = repo_root()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    focal_path = root / "analysis/genomic_transmission_dynamics/derived/transmission_genome_records.tsv"
    background_path = root / "analysis/genomic_transmission_dynamics/derived/global_tree_background_records.tsv"
    catalog_path = root / "state/manifest/genome_catalog.tsv"
    global_path = root / "archive/1-s2.0-S0163445326000435-mmc2.xlsx"

    focal = pd.read_csv(focal_path, sep="\t", dtype=str).fillna("")
    background = pd.read_csv(background_path, sep="\t", dtype=str).fillna("")
    catalog = pd.read_csv(catalog_path, sep="\t", dtype=str).fillna("")
    maps = build_catalog_maps(catalog)
    continent_map = global_continent_map(global_path)
    focal = add_acquisition_route(annotate_continent(annotate_local_sequence(focal, maps), continent_map))
    background = add_acquisition_route(
        annotate_continent(annotate_local_sequence(background, maps), continent_map)
    )

    focal_selected = select_focal(
        focal,
        historical_cell_cap=args.historical_cell_cap,
        historical_country_cap=args.historical_country_cap,
        historical_all_countries={
            value.strip()
            for value in args.historical_all_countries.split(",")
            if value.strip()
        },
    )
    focal_selected["tree_role"] = "focal"
    focal_selected["background_replicate"] = ""
    focal_selected.to_csv(
        args.output_dir / "focal_phylogeny_selection.tsv",
        sep="\t",
        index=False,
    )

    replicate_frames = []
    for replicate in BACKGROUND_REPLICATES:
        selected = select_background(background, replicate)
        replicate_frames.append(selected)
        selected.to_csv(
            args.output_dir / f"background_selection_{replicate}.tsv",
            sep="\t",
            index=False,
        )
    all_background = pd.concat(replicate_frames, ignore_index=True)
    all_background.to_csv(
        args.output_dir / "background_selection_all_replicates.tsv",
        sep="\t",
        index=False,
    )

    primary_background = all_background[all_background["background_replicate"].eq("primary")].copy()
    primary_background["tree_role"] = "global_background"
    manifest = finalise_manifest(pd.concat([focal_selected, primary_background], ignore_index=True))
    manifest.to_csv(args.output_dir / "primary_phylogeny_manifest.tsv", sep="\t", index=False)

    overlap = (
        all_background.assign(selected=True)
        .pivot_table(
            index="genome_record_id",
            columns="background_replicate",
            values="selected",
            aggfunc="max",
            fill_value=False,
        )
        .reset_index()
    )
    overlap.to_csv(args.output_dir / "background_replicate_membership.tsv", sep="\t", index=False)

    counts = {
        "n_focal_selected": int(len(focal_selected)),
        "n_background_primary": int(len(primary_background)),
        "n_primary_tree_total": int(len(manifest)),
        "primary_by_role_country": {
            f"{role}:{country}": int(value)
            for (role, country), value in manifest.groupby(["tree_role", "country_iso3"]).size().items()
        },
        "primary_by_acquisition": {
            str(key): int(value) for key, value in manifest["sequence_acquisition"].value_counts().items()
        },
        "background_replicate_sizes": {
            replicate: int(sum(all_background["background_replicate"].eq(replicate)))
            for replicate in BACKGROUND_REPLICATES
        },
        "selection_rules": {
            "historical_focal": (
                f"cap {args.historical_cell_cap} per "
                "country-five-year-project-preliminary-lineage cell, then cap "
                f"{args.historical_country_cap} per country"
            ),
            "historical_all_countries": [
                value.strip()
                for value in args.historical_all_countries.split(",")
                if value.strip()
            ],
            "modern_focal": "cap 8 per country-month-or-year-project-preliminary-lineage cell",
            "global_background": "cap 2 per non-focal-country-period-preliminary-lineage cell; prefer local FASTA and precise dates",
        },
        "data_freeze_date": DATA_FREEZE_DATE,
    }
    with (args.output_dir / "cohort_freeze_report.json").open("w", encoding="utf-8") as handle:
        json.dump(counts, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps(counts, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
