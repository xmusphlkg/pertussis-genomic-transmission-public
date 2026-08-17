#!/usr/bin/env python3
"""Create minimal, deterministic NAS snapshots for the EID release-clock audit.

The generated snapshots contain only metadata needed by the frozen 989-genome
analysis.  Raw reads and genome sequences are never copied.  Downstream EID
scripts read these repository-local snapshots and therefore do not require the
NAS mount or a neighbouring backup Git repository.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
NAS_ROOT_VALUE = os.environ.get("EID_NAS_PERTUSSIS_ROOT", "").strip()
NAS_ROOT = Path(NAS_ROOT_VALUE) if NAS_ROOT_VALUE else REPO / "__EID_NAS_ROOT_NOT_CONFIGURED__"
SNAPSHOT_DIR = REPO / "data" / "source_snapshots"
PROVENANCE_DIR = REPO / "provenance"

TREE_METADATA = REPO / "results" / "phylogeny" / "tree_tip_metadata.tsv"
CANDIDATE_RUNS = REPO / "results" / "public_availability" / "candidate_project_run_metadata.tsv"

EXTENDED_METADATA = NAS_ROOT / "step1_ingest" / "bp_extended_metadata.tsv"
PUBLIC_MANIFEST = NAS_ROOT / "step1_ingest" / "outputs" / "bp_public_genome_manifest.tsv"
HIGHRES_CASES = NAS_ROOT / "public_health" / "outputs" / "ph_highres_cases.tsv"
PUBLIC_HEALTH_SOURCES = NAS_ROOT / "public_health" / "outputs" / "ph_source_registry.tsv"
RECOVERY_SUMMARY = NAS_ROOT / "step6_epi_transmission" / "outputs" / "bp_focal_country_recovery_summary.tsv"
RAW_QC = NAS_ROOT / "step1_ingest" / "outputs" / "bp_raw_read_assembly_qc_pass_combined.tsv"
PAPER_COMPARISON = (
    NAS_ROOT
    / "snapshots"
    / "repo_root_outputs_legacy_20260423"
    / "paper_dataset_compare_20260330"
    / "paper_included_comparison.tsv"
)

FOCUS_COUNTRIES = {"AUS", "CHN", "JPN"}
PRJNA1071282 = "PRJNA1071282"


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def accession_root(value: str) -> str:
    return (value or "").strip().split(".")[0]


def first_nonmissing(*values: str) -> str:
    for value in values:
        text = (value or "").strip()
        if text and text.lower() not in {"missing", "na", "not applicable"}:
            return text
    return ""


def index_rows(
    rows: list[dict[str, str]],
    assembly_columns: tuple[str, ...],
    biosample_column: str,
    run_column: str = "",
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    by_assembly: dict[str, dict[str, str]] = {}
    by_biosample: dict[str, dict[str, str]] = {}
    by_run: dict[str, dict[str, str]] = {}
    for row in rows:
        for column in assembly_columns:
            accession = accession_root(row.get(column, ""))
            if accession:
                by_assembly.setdefault(accession, row)
        biosample = (row.get(biosample_column) or "").strip()
        if biosample:
            by_biosample.setdefault(biosample, row)
        if run_column:
            run = (row.get(run_column) or "").strip()
            if run:
                by_run.setdefault(run, row)
    return by_assembly, by_biosample, by_run


def lookup_metadata(
    row: dict[str, str],
    by_assembly: dict[str, dict[str, str]],
    by_biosample: dict[str, dict[str, str]],
    by_run: dict[str, dict[str, str]],
) -> dict[str, str] | None:
    assembly = accession_root(row.get("assembly_accession", ""))
    biosample = (row.get("biosample_accession") or "").strip()
    run = (row.get("run_accession") or "").strip()
    return by_assembly.get(assembly) or by_biosample.get(biosample) or by_run.get(run)


def build_focal_metadata_snapshot() -> tuple[Path, int]:
    tree_rows = [
        row
        for row in read_tsv(TREE_METADATA)
        if row.get("tree_role") == "focal" and row.get("country_iso3") in FOCUS_COUNTRIES
    ]
    extended_rows = read_tsv(EXTENDED_METADATA)
    manifest_rows = read_tsv(PUBLIC_MANIFEST)
    ext_assembly, ext_biosample, ext_run = index_rows(
        extended_rows,
        ("Assembly Accession", "Current Accession"),
        "Assembly BioSample Accession",
    )
    pub_assembly, pub_biosample, pub_run = index_rows(
        manifest_rows,
        ("assembly_accession", "current_accession"),
        "biosample_accession",
        "sra_run_accession",
    )

    output_rows: list[dict[str, object]] = []
    for tree in sorted(tree_rows, key=lambda item: item["tree_sample_id"]):
        ext = lookup_metadata(tree, ext_assembly, ext_biosample, ext_run)
        public = lookup_metadata(tree, pub_assembly, pub_biosample, pub_run)
        match_sources = [name for name, match in (("extended_metadata", ext), ("public_manifest", public)) if match]
        output_rows.append(
            {
                "nas_snapshot_record_id": tree["tree_sample_id"],
                "tree_sample_id": tree["tree_sample_id"],
                "country_iso3": tree.get("country_iso3", ""),
                "project_id": tree.get("project_id", ""),
                "run_accession": tree.get("run_accession", ""),
                "biosample_accession": tree.get("biosample_accession", ""),
                "assembly_accession": tree.get("assembly_accession", ""),
                "assembly_release_date": first_nonmissing(
                    (ext or {}).get("Assembly Release Date", ""),
                    (public or {}).get("assembly_release_date", ""),
                ),
                "collection_date_raw": first_nonmissing(
                    (ext or {}).get("Assembly BioSample Collection date", ""),
                    (public or {}).get("collection_date_raw", ""),
                    tree.get("collection_date", ""),
                ),
                "geo_raw": first_nonmissing(
                    (ext or {}).get("Assembly BioSample Geographic location", ""),
                    (public or {}).get("geo_raw", ""),
                    tree.get("country", ""),
                ),
                "sequencing_tech": first_nonmissing(
                    (ext or {}).get("Assembly Sequencing Tech", ""),
                    (public or {}).get("sequencing_tech", ""),
                ),
                "source_database": first_nonmissing(
                    (ext or {}).get("Source Database", ""),
                    (public or {}).get("source_database", ""),
                    tree.get("sequence_acquisition", ""),
                ),
                "match_status": ";".join(match_sources) if match_sources else "unmatched_frozen_record",
            }
        )

    output = SNAPSHOT_DIR / "eid_nas_focal_genome_metadata.tsv"
    fields = [
        "nas_snapshot_record_id",
        "tree_sample_id",
        "country_iso3",
        "project_id",
        "run_accession",
        "biosample_accession",
        "assembly_accession",
        "assembly_release_date",
        "collection_date_raw",
        "geo_raw",
        "sequencing_tech",
        "source_database",
        "match_status",
    ]
    write_tsv(output, output_rows, fields)
    return output, len(output_rows)


def build_public_health_snapshots() -> list[tuple[Path, int]]:
    cases = [row for row in read_tsv(HIGHRES_CASES) if row.get("country_iso3") in FOCUS_COUNTRIES]
    cases.sort(key=lambda row: (row.get("country_iso3", ""), row.get("date", "")))
    case_output = SNAPSHOT_DIR / "eid_nas_highres_cases.tsv"
    case_fields = [
        "country_iso3",
        "country_name",
        "time_resolution",
        "date",
        "year",
        "month",
        "week",
        "cases",
        "source_url",
        "source_file",
        "data_freeze_date",
        "notes",
    ]
    write_tsv(case_output, cases, case_fields)

    sources = [row for row in read_tsv(PUBLIC_HEALTH_SOURCES) if row.get("country_iso3") in FOCUS_COUNTRIES]
    sources.sort(key=lambda row: (row.get("country_iso3", ""), row.get("source_id", "")))
    source_output = SNAPSHOT_DIR / "eid_nas_public_health_sources.tsv"
    source_fields = [
        "source_id",
        "source_name",
        "source_url",
        "source_kind",
        "source_domain",
        "country_iso3",
        "source_release_date",
        "source_access_date",
        "freeze_policy",
        "notes",
    ]
    write_tsv(source_output, sources, source_fields)
    return [(case_output, len(cases)), (source_output, len(sources))]


def build_prjna1071282_snapshots() -> list[tuple[Path, int]]:
    candidate_rows = [
        row for row in read_tsv(CANDIDATE_RUNS) if row.get("query_project_id") == PRJNA1071282
    ]
    run_ids = {row.get("run_accession", "") for row in candidate_rows}

    comparison_by_run = {
        row.get("paper_sra_runinfo", ""): row
        for row in read_tsv(PAPER_COMPARISON)
        if row.get("paper_sra_runinfo", "") in run_ids
    }
    qc_by_run = {
        row.get("run_accession", ""): row
        for row in read_tsv(RAW_QC)
        if row.get("run_accession", "") in run_ids
    }

    output_rows: list[dict[str, object]] = []
    for candidate in sorted(candidate_rows, key=lambda row: row.get("run_accession", "")):
        run = candidate.get("run_accession", "")
        comparison = comparison_by_run.get(run, {})
        qc = qc_by_run.get(run, {})
        matches = [name for name, value in (("paper_comparison", comparison), ("raw_read_qc", qc)) if value]
        output_rows.append(
            {
                "project_id": PRJNA1071282,
                "run_accession": run,
                "sample_accession": candidate.get("sample_accession", ""),
                "collection_date": candidate.get("collection_date", ""),
                "country": candidate.get("country", ""),
                "first_public": candidate.get("first_public", ""),
                "comparison_status": comparison.get("comparison_status", ""),
                "gapfill_success_n": comparison.get("gapfill_success_n", ""),
                "gapfill_success_runs": comparison.get("gapfill_success_runs", ""),
                "qc_decision": qc.get("qc_decision", ""),
                "qc_reason": qc.get("qc_reason", ""),
                "checkm_completeness": qc.get("checkm_completeness", ""),
                "checkm_contamination": qc.get("checkm_contamination", ""),
                "nas_match_status": ";".join(matches) if matches else "no_nas_recovery_record",
            }
        )

    project_output = SNAPSHOT_DIR / "eid_nas_prjna1071282_run_audit.tsv"
    project_fields = [
        "project_id",
        "run_accession",
        "sample_accession",
        "collection_date",
        "country",
        "first_public",
        "comparison_status",
        "gapfill_success_n",
        "gapfill_success_runs",
        "qc_decision",
        "qc_reason",
        "checkm_completeness",
        "checkm_contamination",
        "nas_match_status",
    ]
    write_tsv(project_output, output_rows, project_fields)

    recovery = [row for row in read_tsv(RECOVERY_SUMMARY) if row.get("country_iso3") == "CHN"]
    recovery.sort(key=lambda row: int(row.get("year") or 0))
    recovery_output = SNAPSHOT_DIR / "eid_nas_china_recovery_summary.tsv"
    recovery_fields = [
        "country_iso3",
        "country_name",
        "year",
        "current_interpretable",
        "current_total_genomes",
        "success_runs_reconciled",
        "success_runs_missing_manifest",
        "planned_only_remaining",
        "reconciled_total_genomes_est",
        "recovery_status",
        "manifest_backfill_required",
        "notes",
    ]
    write_tsv(recovery_output, recovery, recovery_fields)
    return [(project_output, len(output_rows)), (recovery_output, len(recovery))]


def write_manifest(outputs: list[tuple[Path, int]]) -> Path:
    source_specs = [
        (EXTENDED_METADATA, "identifier match to frozen AUS/CHN/JPN focal records"),
        (PUBLIC_MANIFEST, "identifier match to frozen AUS/CHN/JPN focal records"),
        (HIGHRES_CASES, "country_iso3 in AUS, CHN, JPN"),
        (PUBLIC_HEALTH_SOURCES, "country_iso3 in AUS, CHN, JPN"),
        (PAPER_COMPARISON, "paper_sra_runinfo in PRJNA1071282 ENA run list"),
        (RAW_QC, "run_accession in PRJNA1071282 ENA run list"),
        (RECOVERY_SUMMARY, "country_iso3 equals CHN"),
    ]
    manifest_rows: list[dict[str, object]] = []
    for path, filter_rule in source_specs:
        stat = path.stat()
        manifest_rows.append(
            {
                "record_type": "nas_source",
                "path": str(path),
                "size_bytes": stat.st_size,
                "mtime_epoch_seconds": int(stat.st_mtime),
                "sha256": sha256(path),
                "filter_rule": filter_rule,
                "n_records": "",
            }
        )
    for path, n_records in outputs:
        stat = path.stat()
        manifest_rows.append(
            {
                "record_type": "repository_snapshot",
                "path": str(path.relative_to(REPO)),
                "size_bytes": stat.st_size,
                "mtime_epoch_seconds": "",
                "sha256": sha256(path),
                "filter_rule": "deterministic output of gtd_39_snapshot_eid_nas_inputs.py",
                "n_records": n_records,
            }
        )
    manifest_rows.sort(key=lambda row: (str(row["record_type"]), str(row["path"])))
    output = PROVENANCE_DIR / "EID_NAS_SNAPSHOT_MANIFEST.tsv"
    write_tsv(
        output,
        manifest_rows,
        ["record_type", "path", "size_bytes", "mtime_epoch_seconds", "sha256", "filter_rule", "n_records"],
    )
    return output


def main() -> None:
    if not NAS_ROOT_VALUE:
        raise SystemExit(
            "Set EID_NAS_PERTUSSIS_ROOT to the mounted source-data root before regenerating snapshots."
        )
    required = [
        TREE_METADATA,
        CANDIDATE_RUNS,
        EXTENDED_METADATA,
        PUBLIC_MANIFEST,
        HIGHRES_CASES,
        PUBLIC_HEALTH_SOURCES,
        RECOVERY_SUMMARY,
        RAW_QC,
        PAPER_COMPARISON,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing required snapshot inputs:\n" + "\n".join(missing))

    outputs = [build_focal_metadata_snapshot()]
    outputs.extend(build_public_health_snapshots())
    outputs.extend(build_prjna1071282_snapshots())
    manifest = write_manifest(outputs)
    for path, n_records in outputs:
        print(f"wrote {path.relative_to(REPO)} ({n_records} records)")
    print(f"wrote {manifest.relative_to(REPO)}")


if __name__ == "__main__":
    main()
