#!/usr/bin/env python3
"""Apply source-independent assembly/read QC and write SKA2 inputs."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def fasta_metrics(path: Path) -> dict[str, object]:
    lengths = []
    total = 0
    gc = 0
    n_bases = 0
    current = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            sequence = line.strip().upper()
            current += len(sequence)
            total += len(sequence)
            gc += sequence.count("G") + sequence.count("C")
            n_bases += sequence.count("N")
    if current:
        lengths.append(current)
    sorted_lengths = sorted(lengths, reverse=True)
    halfway = total / 2
    cumulative = 0
    n50 = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= halfway:
            n50 = length
            break
    return {
        "observed_length": total,
        "observed_gc_percent": 100 * gc / total if total else 0,
        "observed_n_fraction": n_bases / total if total else 1,
        "observed_contigs": len(lengths),
        "observed_n50": n50,
    }


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/primary_phylogeny_manifest_with_fastq.tsv",
    )
    parser.add_argument(
        "--fastq-qc",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/fastq_subset_500000_qc.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "analysis/genomic_transmission_dynamics/phylogeny/qc",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    fastq_qc = (
        pd.read_csv(args.fastq_qc, sep="\t", dtype=str).fillna("")
        if args.fastq_qc.exists()
        else pd.DataFrame(columns=["run_accession"])
    )
    fastq_lookup = fastq_qc.set_index("run_accession").to_dict("index")
    rows = []
    for _, source in manifest.iterrows():
        route = source["sequence_acquisition"]
        result = source.to_dict()
        result.update(
            {
                "qc_input_type": "",
                "observed_length": "",
                "observed_gc_percent": "",
                "observed_n_fraction": "",
                "observed_contigs": "",
                "observed_n50": "",
                "paired_records": "",
                "read_q30_min": "",
                "uniform_qc_status": "FAIL_UNRESOLVED_INPUT",
                "uniform_qc_reason": "",
            }
        )
        if route in {
            "local_fasta",
            "ncbi_assembly_download",
            "legacy_read_assembly",
        }:
            path = Path(source["local_fasta_path"])
            result["qc_input_type"] = "assembly"
            if not path.exists():
                result["uniform_qc_reason"] = "local_fasta_missing"
            else:
                try:
                    metrics = fasta_metrics(path)
                except (OSError, PermissionError) as error:
                    result["uniform_qc_reason"] = f"assembly_unreadable:{type(error).__name__}"
                    rows.append(result)
                    continue
                result.update(metrics)
                pass_flags = {
                    "length": 3_500_000 <= int(metrics["observed_length"]) <= 4_500_000,
                    "gc": 65.0 <= float(metrics["observed_gc_percent"]) <= 70.0,
                    "n_fraction": float(metrics["observed_n_fraction"]) <= 0.05,
                    "contigs": int(metrics["observed_contigs"]) <= 500,
                    "n50": int(metrics["observed_n50"]) >= 5_000,
                }
                failed = [name for name, passed in pass_flags.items() if not passed]
                result["uniform_qc_status"] = "PASS" if not failed else "FAIL_ASSEMBLY_QC"
                result["uniform_qc_reason"] = ";".join(failed)
                result["sequence_input_path"] = str(path)
        elif route == "ena_fastq_download":
            result["qc_input_type"] = "paired_reads"
            metrics = fastq_lookup.get(source["run_accession"], {})
            paired_records = int(float(metrics.get("paired_records", 0) or 0))
            q30_values = [
                float(metrics.get("r1_q30_fraction", 0) or 0),
                float(metrics.get("r2_q30_fraction", 0) or 0),
            ]
            result["paired_records"] = paired_records
            result["read_q30_min"] = min(q30_values)
            result["fastq_r1"] = metrics.get("r1_path", source.get("fastq_r1", ""))
            result["fastq_r2"] = metrics.get("r2_path", source.get("fastq_r2", ""))
            pass_flags = {
                "subset_status": metrics.get("status", "") == "PASS",
                "paired_records": paired_records >= 50_000,
                # Do not use modern-platform Q30 yield as a country/project
                # selection filter. Catastrophic-quality runs fail here;
                # strict Q20 split-kmer filtering and alignment missingness
                # remain the sequence-level hard filters.
                "q30": min(q30_values) >= 0.45,
                "files": bool(result["fastq_r1"])
                and bool(result["fastq_r2"])
                and Path(result["fastq_r1"]).exists()
                and Path(result["fastq_r2"]).exists(),
            }
            failed = [name for name, passed in pass_flags.items() if not passed]
            result["uniform_qc_status"] = "PASS" if not failed else "FAIL_READ_QC"
            result["uniform_qc_reason"] = ";".join(failed)
            result["sequence_input_path"] = (
                f"{result['fastq_r1']};{result['fastq_r2']}" if not failed else ""
            )
        rows.append(result)

    qc = pd.DataFrame(rows)
    qc["tree_include_after_uniform_qc"] = qc["uniform_qc_status"].eq("PASS")
    qc.to_csv(args.output_dir / "uniform_sequence_qc.tsv", sep="\t", index=False)

    passed = qc[qc["tree_include_after_uniform_qc"]].copy()
    ska_lines = []
    for _, row in passed.iterrows():
        if row["qc_input_type"] == "assembly":
            ska_lines.append(f"{row['tree_sample_id']}\t{row['local_fasta_path']}")
        else:
            ska_lines.append(
                f"{row['tree_sample_id']}\t{row['fastq_r1']}\t{row['fastq_r2']}"
            )
    (args.output_dir / "ska_primary_input.tsv").write_text(
        "\n".join(ska_lines) + ("\n" if ska_lines else ""),
        encoding="utf-8",
    )
    report = {
        "n_manifest": int(len(qc)),
        "n_uniform_qc_pass": int(qc["uniform_qc_status"].eq("PASS").sum()),
        "n_uniform_qc_fail": int((~qc["uniform_qc_status"].eq("PASS")).sum()),
        "status_counts": {
            str(key): int(value) for key, value in qc["uniform_qc_status"].value_counts().items()
        },
        "pass_by_role_country": {
            f"{role}:{country}": int(value)
            for (role, country), value in passed.groupby(["tree_role", "country_iso3"]).size().items()
        },
    }
    report["n_uniform_qc_fail"] = int(
        (~qc["uniform_qc_status"].eq("PASS")).sum()
    )
    (args.output_dir / "uniform_sequence_qc_report.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
