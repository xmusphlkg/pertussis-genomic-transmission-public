#!/usr/bin/env python3
"""Remove non-redistributable specimen-level fields from public TSV outputs."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FONG_SOURCE_DATASET = "australia_fong_2026_direct_specimens"
REDACTED_NOTE = "Third-party specimen metadata not redistributed"

TARGET_TABLES = (
    "data/derived/transmission_genome_records.tsv",
    "results/phylogeny/focal_phylogeny_selection.tsv",
    "results/phylogeny/primary_phylogeny_manifest.tsv",
    "results/phylogeny/tree_tip_metadata.tsv",
    "results/phylogeny/uniform_sequence_qc.tsv",
    "figures/source_data/figure2a_tree_tip_metadata.tsv",
)

RESTRICTED_COLUMNS = {
    "sample_id",
    "genome_qc_status",
    "sampling_process_observed",
    "sequencing_success",
    "ct_value",
    "specimen_type",
    "preliminary_lineage_id",
    "lineage_definition_status",
    "lineage_stratum",
    "published_branch",
    "published_lineage",
    "published_sublineage",
    "ptxP_label",
    "fim3_label",
    "marker_23s_status",
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return reader.fieldnames or [], rows


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with tempfile.NamedTemporaryFile(
        "w",
        delete=False,
        dir=path.parent,
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        temp_path = Path(handle.name)
    temp_path.replace(path)


def fong_biosamples() -> set[str]:
    _, rows = read_tsv(ROOT / "data/derived/transmission_genome_records.tsv")
    return {
        row["biosample_accession"]
        for row in rows
        if row.get("source_dataset") == FONG_SOURCE_DATASET
        and row.get("biosample_accession")
    }


def is_fong_row(row: dict[str, str], biosamples: set[str]) -> bool:
    return (
        row.get("source_dataset") == FONG_SOURCE_DATASET
        or row.get("biosample_accession") in biosamples
        or row.get("sample_accession") in biosamples
        or row.get("sample_id", "").startswith("24-BPE-")
    )


def sanitise_table(relative_path: str, biosamples: set[str]) -> int:
    path = ROOT / relative_path
    if not path.is_file():
        return 0

    fieldnames, rows = read_tsv(path)
    changed_rows = 0
    for row in rows:
        if not is_fong_row(row, biosamples):
            continue
        changed_rows += 1
        for column in RESTRICTED_COLUMNS.intersection(fieldnames):
            row[column] = ""
        if "provenance_note" in fieldnames:
            row["provenance_note"] = REDACTED_NOTE

    write_tsv(path, fieldnames, rows)
    return changed_rows


def main() -> None:
    biosamples = fong_biosamples()
    total = 0
    for relative_path in TARGET_TABLES:
        total += sanitise_table(relative_path, biosamples)
    print(
        f"Sanitised {total} Fong-linked rows across {len(TARGET_TABLES)} public TSVs"
    )


if __name__ == "__main__":
    main()
