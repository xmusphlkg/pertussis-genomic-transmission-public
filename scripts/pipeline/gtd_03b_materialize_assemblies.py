#!/usr/bin/env python3
"""Download selected NCBI assemblies in batches and materialise one FASTA per accession."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plan",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/acquisition/assembly_download_plan.tsv",
    )
    parser.add_argument("--batch-size", type=int, default=100)
    return parser.parse_args()


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def materialise_batch(batch: list[str], targets: dict[str, Path]) -> dict[str, str]:
    statuses: dict[str, str] = {}
    with tempfile.TemporaryDirectory(prefix="gtd_ncbi_") as temp_name:
        temp = Path(temp_name)
        accessions = temp / "accessions.txt"
        accessions.write_text("\n".join(batch) + "\n", encoding="utf-8")
        package = temp / "ncbi_dataset.zip"
        datasets_binary = Path(sys.executable).with_name("datasets")
        command = [
            str(datasets_binary),
            "download",
            "genome",
            "accession",
            "--inputfile",
            str(accessions),
            "--include",
            "genome",
            "--filename",
            str(package),
            "--no-progressbar",
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if result.returncode != 0 or not package.is_file():
            message = (result.stderr or result.stdout).strip().replace("\t", " ")[:500]
            return {accession: f"FAIL_DATASETS:{message}" for accession in batch}
        with zipfile.ZipFile(package) as archive:
            candidates: dict[str, str] = {}
            for name in archive.namelist():
                if not name.endswith(("_genomic.fna", "_genomic.fna.gz")):
                    continue
                parts = Path(name).parts
                for accession in batch:
                    if accession in parts or accession in name:
                        candidates.setdefault(accession, name)
                        break
            for accession in batch:
                member = candidates.get(accession)
                if not member:
                    statuses[accession] = "FAIL_NO_GENOMIC_FASTA_IN_PACKAGE"
                    continue
                target = targets[accession]
                target.parent.mkdir(parents=True, exist_ok=True)
                temporary_target = target.with_suffix(".fasta.partial")
                with archive.open(member) as source, temporary_target.open("wb") as destination:
                    shutil.copyfileobj(source, destination)
                if temporary_target.stat().st_size == 0:
                    temporary_target.unlink(missing_ok=True)
                    statuses[accession] = "FAIL_EMPTY_FASTA"
                    continue
                temporary_target.replace(target)
                statuses[accession] = "PASS"
    return statuses


def main() -> None:
    args = parse_args()
    plan = pd.read_csv(args.plan, sep="\t", dtype=str).fillna("")
    plan = plan.drop_duplicates("assembly_accession").copy()
    targets = {
        row["assembly_accession"]: Path(row["local_path"])
        for _, row in plan.iterrows()
    }
    pending = [
        accession
        for accession, target in targets.items()
        if not target.is_file() or target.stat().st_size == 0
    ]
    statuses = {
        accession: "PASS_EXISTING"
        for accession, target in targets.items()
        if target.is_file() and target.stat().st_size > 0
    }
    batches = chunks(pending, args.batch_size)
    for batch_number, batch in enumerate(batches, start=1):
        batch_status = materialise_batch(batch, targets)
        statuses.update(batch_status)
        passed = sum(value == "PASS" for value in batch_status.values())
        print(f"{batch_number}/{len(batches)}: {passed}/{len(batch)} assemblies materialised")

    plan["materialisation_status"] = plan["assembly_accession"].map(statuses).fillna(
        "FAIL_NOT_ATTEMPTED"
    )
    report_path = args.plan.parent / "assembly_materialisation_status.tsv"
    plan.to_csv(report_path, sep="\t", index=False)
    report = {
        "n_assembly_accessions": int(len(targets)),
        "n_pass": int(plan.drop_duplicates("assembly_accession")[
            "materialisation_status"
        ].str.startswith("PASS").sum()),
        "n_fail": int(plan.drop_duplicates("assembly_accession")[
            "materialisation_status"
        ].str.startswith("FAIL").sum()),
        "status_file": str(report_path),
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
