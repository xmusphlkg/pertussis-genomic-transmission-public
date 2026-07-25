#!/usr/bin/env python3
"""Prepare focal, uniformly QC-passed assemblies for cgMLST validation."""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import subprocess
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniform-qc", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--skesa", default="skesa")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--cores-per-job", type=int, default=2)
    parser.add_argument(
        "--assembly-inputs-only",
        action="store_true",
        help="Restrict validation to uniformly QC-passed pre-existing assemblies.",
    )
    return parser.parse_args()


def fasta_metrics(path: Path) -> tuple[int, int, int]:
    lengths: list[int] = []
    current = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    if not lengths:
        return 0, 0, 0
    ordered = sorted(lengths, reverse=True)
    total = sum(ordered)
    cumulative = 0
    n50 = 0
    for length in ordered:
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return total, len(lengths), n50


def assemble(row: dict[str, str], args: argparse.Namespace) -> dict[str, object]:
    sample = row["tree_sample_id"]
    if row["qc_input_type"] == "assembly":
        target = Path(row["sequence_input_path"])
        status = "PASS_EXISTING" if target.is_file() and target.stat().st_size else "FAIL_MISSING"
    else:
        target = args.output_root / f"{sample}.fasta"
        if target.is_file() and target.stat().st_size:
            status = "PASS_EXISTING_CGMLST_ASSEMBLY"
        else:
            temporary = target.with_suffix(".fasta.partial")
            temporary.unlink(missing_ok=True)
            command = [
                args.skesa,
                "--fastq",
                f"{row['r1_path']},{row['r2_path']}",
                "--contigs_out",
                str(temporary),
                "--cores",
                str(args.cores_per_job),
            ]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            if result.returncode == 0 and temporary.is_file() and temporary.stat().st_size:
                temporary.replace(target)
                status = "PASS_ASSEMBLED"
            else:
                temporary.unlink(missing_ok=True)
                return {
                    "tree_sample_id": sample,
                    "cgmlst_assembly_path": str(target),
                    "assembly_status": "FAIL_SKESA",
                    "assembly_error": (result.stderr or result.stdout)[-1000:].replace("\t", " "),
                    "assembly_length": 0,
                    "assembly_contigs": 0,
                    "assembly_n50": 0,
                }
    length, contigs, n50 = fasta_metrics(target) if status.startswith("PASS") else (0, 0, 0)
    if status.startswith("PASS") and not (3_200_000 <= length <= 4_500_000):
        status = "FAIL_CGMLST_ASSEMBLY_LENGTH"
    return {
        "tree_sample_id": sample,
        "cgmlst_assembly_path": str(target),
        "assembly_status": status,
        "assembly_error": "",
        "assembly_length": length,
        "assembly_contigs": contigs,
        "assembly_n50": n50,
    }


def main() -> None:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    qc = pd.read_csv(args.uniform_qc, sep="\t", dtype=str).fillna("")
    focal = qc[
        qc["tree_role"].eq("focal") & qc["uniform_qc_status"].eq("PASS")
    ].copy()
    if args.assembly_inputs_only:
        focal = focal[focal["qc_input_type"].eq("assembly")].copy()
    rows = focal.to_dict(orient="records")
    results: list[dict[str, object]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(assemble, row, args): row for row in rows}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            results.append(result)
            print(
                f"{index}/{len(rows)} {result['tree_sample_id']} "
                f"{result['assembly_status']}",
                flush=True,
            )
    status = pd.DataFrame(results)
    manifest = focal.merge(status, on="tree_sample_id", how="left", validate="one_to_one")
    manifest.to_csv(args.report_dir / "cgmlst_assembly_manifest.tsv", sep="\t", index=False)
    report = {
        "n_focal_uniform_qc_pass": len(focal),
        "assembly_inputs_only": args.assembly_inputs_only,
        "n_cgmlst_assembly_pass": int(
            manifest["assembly_status"].str.startswith("PASS").sum()
        ),
        "n_cgmlst_assembly_fail": int(
            (~manifest["assembly_status"].str.startswith("PASS")).sum()
        ),
        "status_counts": manifest["assembly_status"].value_counts().to_dict(),
    }
    (args.report_dir / "cgmlst_assembly_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
