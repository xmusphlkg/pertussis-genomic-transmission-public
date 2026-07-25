#!/usr/bin/env python3
"""Resolve reproducible NCBI assembly and ENA FASTQ inputs for the tree."""

from __future__ import annotations

import argparse
import io
import json
import os
import time
from pathlib import Path

import pandas as pd
import requests


DATA_FREEZE_DATE = "2026-07-24"
ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def query_ena_runs(run_accessions: list[str], cache_path: Path | None = None) -> pd.DataFrame:
    fields = [
        "run_accession",
        "fastq_ftp",
        "fastq_bytes",
        "library_layout",
        "read_count",
        "base_count",
        "study_accession",
        "sample_accession",
        "collection_date",
        "country",
    ]
    frames = []
    cached_runs: set[str] = set()
    if cache_path and cache_path.is_file():
        cached = pd.read_csv(cache_path, sep="\t", dtype=str).fillna("")
        frames.append(cached[cached["run_accession"].isin(run_accessions)])
        cached_runs = set(cached["run_accession"])
    missing_runs = sorted(set(run_accessions) - cached_runs)
    for batch in chunks(missing_runs, 40):
        query = " OR ".join(f'run_accession="{run}"' for run in batch)
        response = None
        for attempt in range(5):
            try:
                response = requests.get(
                    ENA_SEARCH_URL,
                    params={
                        "result": "read_run",
                        "query": query,
                        "fields": ",".join(fields),
                        "format": "tsv",
                        "limit": "0",
                    },
                    timeout=120,
                )
                response.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 4:
                    raise
                time.sleep(2**attempt)
        assert response is not None
        if response.text.strip():
            frames.append(pd.read_csv(io.StringIO(response.text), sep="\t", dtype=str).fillna(""))
    if not frames:
        return pd.DataFrame(columns=fields)
    return pd.concat(frames, ignore_index=True).drop_duplicates("run_accession")


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/cohort/primary_phylogeny_manifest.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "analysis/genomic_transmission_dynamics/phylogeny/acquisition",
    )
    parser.add_argument(
        "--sequence-data-root",
        type=Path,
        default=root
        / "pertussis_data/pertussis_gene/genomic_transmission_dynamics",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    reads_root = args.sequence_data_root / "reads/raw"
    assemblies_root = args.sequence_data_root / "assemblies"
    reads_root.mkdir(parents=True, exist_ok=True)
    assemblies_root.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    # Canonical assembly files are not uniformly readable on the shared NFS.
    # Re-materialise every public assembly accession into a project-owned path.
    has_assembly = manifest["assembly_accession"].str.match(
        r"^GC[AF]_\d+(?:\.\d+)?$"
    )
    has_run = manifest["run_accession"].ne("")
    readable_local = manifest["local_fasta_path"].map(
        lambda value: bool(value) and Path(value).is_file()
        and Path(value).stat().st_mode & 0o444 != 0
    )
    manifest.loc[has_assembly, "sequence_acquisition"] = "ncbi_assembly_download"
    manifest.loc[~has_assembly & readable_local, "sequence_acquisition"] = "local_fasta"
    manifest.loc[
        ~has_assembly & ~readable_local & has_run, "sequence_acquisition"
    ] = "ena_fastq_download"
    manifest.loc[
        ~has_assembly & ~has_run & ~readable_local, "sequence_acquisition"
    ] = "unresolved"
    legacy_root = repo_root() / "pertussis_data/bp_step4/outputs/assemblies"
    for index in manifest.index[
        manifest["sequence_acquisition"].eq("ena_fastq_download")
    ]:
        run = manifest.at[index, "run_accession"]
        candidates = [
            legacy_root / server / run / "contigs.fa"
            for server in ("server1", "server2", "server3")
        ]
        for candidate in candidates:
            try:
                usable = (
                    candidate.is_file()
                    and os.access(candidate, os.R_OK)
                    and candidate.stat().st_size > 1_000_000
                )
            except OSError:
                usable = False
            if usable:
                manifest.at[index, "sequence_acquisition"] = "legacy_read_assembly"
                manifest.at[index, "local_fasta_path"] = str(candidate.resolve())
                break

    required_runs = manifest.loc[
        manifest["sequence_acquisition"].eq("ena_fastq_download"), "run_accession"
    ].tolist()
    ena_cache = args.output_dir / "ena_run_resolution.tsv"
    ena = query_ena_runs(required_runs, ena_cache)
    ena.to_csv(ena_cache, sep="\t", index=False)
    manifest = manifest.merge(ena, on="run_accession", how="left", suffixes=("", "_ena")).fillna("")

    manifest["fastq_ftp_list"] = manifest["fastq_ftp"].astype(str)
    manifest["fastq_bytes_list"] = manifest["fastq_bytes"].astype(str)
    manifest["fastq_total_bytes"] = 0
    manifest["source_resolution_status"] = "unresolved_sequence_source"
    manifest.loc[
        manifest["sequence_acquisition"].eq("local_fasta"),
        "source_resolution_status",
    ] = "resolved_local_fasta"
    manifest.loc[
        manifest["sequence_acquisition"].eq("ncbi_assembly_download"),
        "source_resolution_status",
    ] = "pending_ncbi_assembly"
    manifest.loc[
        manifest["sequence_acquisition"].eq("legacy_read_assembly"),
        "source_resolution_status",
    ] = "resolved_legacy_read_assembly"
    manifest.loc[
        manifest["sequence_acquisition"].eq("ena_fastq_download"),
        "source_resolution_status",
    ] = "unresolved_ena_run"

    download_rows = []
    assembly_rows = []
    for index, row in manifest.iterrows():
        if row["sequence_acquisition"] == "ncbi_assembly_download":
            accession = row["assembly_accession"]
            local_path = assemblies_root / f"{accession}.fasta"
            manifest.at[index, "local_fasta_path"] = str(local_path)
            manifest.at[index, "sequence_input_path"] = str(local_path)
            if local_path.is_file() and local_path.stat().st_size > 0:
                manifest.at[index, "source_resolution_status"] = "resolved_ncbi_assembly"
            assembly_rows.append(
                {
                    "tree_sample_id": row["tree_sample_id"],
                    "assembly_accession": accession,
                    "local_path": str(local_path),
                    "already_present": local_path.is_file()
                    and local_path.stat().st_size > 0,
                }
            )
            continue
        if row["sequence_acquisition"] != "ena_fastq_download":
            manifest.at[index, "sequence_input_path"] = row["local_fasta_path"]
            continue
        all_urls = [value for value in str(row["fastq_ftp"]).split(";") if value]
        paired_urls = [
            value
            for value in all_urls
            if value.endswith("_1.fastq.gz") or value.endswith("_2.fastq.gz")
        ]
        urls = paired_urls if len(paired_urls) == 2 else all_urls
        byte_values = []
        for value in str(row["fastq_bytes"]).split(";"):
            try:
                byte_values.append(int(value))
            except ValueError:
                pass
        manifest.at[index, "fastq_total_bytes"] = sum(byte_values)
        if row["library_layout"] != "PAIRED" or len(urls) != 2:
            manifest.at[index, "source_resolution_status"] = "blocked_not_two_file_paired_fastq"
            continue
        run = row["run_accession"]
        local_paths = []
        for read_number, url in enumerate(urls, start=1):
            filename = f"{run}_{read_number}.fastq.gz"
            local_path = reads_root / filename
            local_paths.append(str(local_path))
            expected_bytes = byte_values[read_number - 1] if read_number <= len(byte_values) else ""
            download_rows.append(
                {
                    "tree_sample_id": row["tree_sample_id"],
                    "run_accession": run,
                    "read_number": read_number,
                    "url": f"https://{url}",
                    "local_path": str(local_path),
                    "expected_bytes": expected_bytes,
                    "already_present": local_path.exists(),
                    "present_bytes": local_path.stat().st_size if local_path.exists() else 0,
                }
            )
        manifest.at[index, "fastq_r1"] = local_paths[0]
        manifest.at[index, "fastq_r2"] = local_paths[1]
        manifest.at[index, "sequence_input_path"] = ";".join(local_paths)
        manifest.at[index, "source_resolution_status"] = "resolved_ena_paired_fastq"

    download_plan = pd.DataFrame(download_rows)
    download_plan.to_csv(args.output_dir / "fastq_download_plan.tsv", sep="\t", index=False)
    assembly_plan = pd.DataFrame(assembly_rows)
    assembly_plan.to_csv(args.output_dir / "assembly_download_plan.tsv", sep="\t", index=False)
    (args.output_dir / "assembly_accessions.txt").write_text(
        "\n".join(
            sorted(
                assembly_plan.loc[
                    ~assembly_plan["already_present"], "assembly_accession"
                ].drop_duplicates()
            )
        )
        + ("\n" if (~assembly_plan["already_present"]).any() else ""),
        encoding="utf-8",
    )
    manifest.to_csv(
        args.output_dir / "primary_phylogeny_manifest_resolved.tsv",
        sep="\t",
        index=False,
    )

    aria2_lines = []
    for _, row in download_plan.iterrows():
        if bool(row["already_present"]) and int(row["present_bytes"]) == int(row["expected_bytes"]):
            continue
        local_path = Path(row["local_path"])
        aria2_lines.extend(
            [
                str(row["url"]),
                f"  dir={local_path.parent}",
                f"  out={local_path.name}",
            ]
        )
    (args.output_dir / "aria2_fastq_downloads.txt").write_text(
        "\n".join(aria2_lines) + ("\n" if aria2_lines else ""),
        encoding="utf-8",
    )

    report = {
        "n_tree_samples": int(len(manifest)),
        "n_local_fasta": int(manifest["sequence_acquisition"].eq("local_fasta").sum()),
        "n_ncbi_assembly": int(
            manifest["sequence_acquisition"].eq("ncbi_assembly_download").sum()
        ),
        "n_legacy_read_assembly": int(
            manifest["sequence_acquisition"].eq("legacy_read_assembly").sum()
        ),
        "n_resolved_ncbi_assembly": int(
            manifest["source_resolution_status"].eq("resolved_ncbi_assembly").sum()
        ),
        "n_ena_fastq": int(manifest["sequence_acquisition"].eq("ena_fastq_download").sum()),
        "n_resolved_ena_fastq": int(
            manifest["source_resolution_status"].eq("resolved_ena_paired_fastq").sum()
        ),
        "n_blocked_sequence_sources": int(
            manifest["source_resolution_status"].str.startswith("blocked").sum()
            + manifest["source_resolution_status"].str.startswith("unresolved").sum()
        ),
        "expected_fastq_bytes": int(pd.to_numeric(download_plan["expected_bytes"], errors="coerce").sum()),
        "expected_fastq_gib": round(
            float(pd.to_numeric(download_plan["expected_bytes"], errors="coerce").sum())
            / (1024**3),
            2,
        ),
        "download_files": int(len(download_plan)),
        "data_freeze_date": DATA_FREEZE_DATE,
    }
    with (args.output_dir / "sequence_source_resolution_report.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
