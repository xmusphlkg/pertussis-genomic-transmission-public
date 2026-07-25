#!/usr/bin/env python3
"""Select read inputs that fail reference missingness for deeper streaming rescue."""

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
        "--manifest",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/uniform_sequence_qc.tsv",
    )
    parser.add_argument(
        "--alignment-metrics",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/results/primary/alignment_missingness_initial.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/rescue250k_manifest.tsv",
    )
    parser.add_argument("--max-missing", type=float, default=0.20)
    parser.add_argument("--target-pairs", type=int, default=250000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    metrics = pd.read_csv(args.alignment_metrics, sep="\t", dtype=str).fillna("")
    metrics["missing_fraction"] = pd.to_numeric(
        metrics["missing_fraction"], errors="coerce"
    )
    rescue_ids = set(
        metrics.loc[
            metrics["missing_fraction"].gt(args.max_missing), "tree_sample_id"
        ]
    )
    rescue = manifest[
        manifest["tree_sample_id"].isin(rescue_ids)
        & manifest["sequence_acquisition"].eq("ena_fastq_download")
        & manifest["uniform_qc_status"].eq("PASS")
    ].copy()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rescue.to_csv(args.output, sep="\t", index=False)
    report = {
        "n_alignment_fail": int(len(rescue_ids)),
        "n_read_rescue": int(len(rescue)),
        "by_country": {
            str(key): int(value)
            for key, value in rescue["country_iso3"].value_counts().items()
        },
        "rescue_target_pairs": args.target_pairs,
        "unchanged_hard_missingness_threshold": args.max_missing,
    }
    report_path = args.output.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
