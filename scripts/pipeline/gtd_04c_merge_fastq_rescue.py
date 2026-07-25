#!/usr/bin/env python3
"""Merge primary and deeper-rescue FASTQ QC, preferring successful rescue rows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primary",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/fastq_subset_100000_qc.tsv",
    )
    parser.add_argument(
        "--rescue",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/rescue250k/fastq_subset_250000_qc.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/fastq_combined_qc.tsv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primary = pd.read_csv(args.primary, sep="\t", dtype=str).fillna("")
    rescue = pd.read_csv(args.rescue, sep="\t", dtype=str).fillna("")
    successful_rescue = rescue[rescue["status"].eq("PASS")].copy()
    combined = pd.concat(
        [
            primary[~primary["run_accession"].isin(successful_rescue["run_accession"])],
            successful_rescue,
        ],
        ignore_index=True,
    ).sort_values("run_accession")
    combined.to_csv(args.output, sep="\t", index=False)
    report = {
        "n_primary_runs": int(len(primary)),
        "n_rescue_attempted": int(len(rescue)),
        "n_rescue_successful": int(len(successful_rescue)),
        "n_combined_runs": int(len(combined)),
        "combined_status_counts": {
            str(key): int(value)
            for key, value in combined["status"].value_counts().items()
        },
    }
    args.output.with_suffix(".report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
