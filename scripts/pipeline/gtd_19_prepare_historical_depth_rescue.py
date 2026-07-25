#!/usr/bin/env python3
"""Freeze historical paired-read samples failing the 100k-read alignment gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--uniform-qc", type=Path, required=True)
    parser.add_argument("--alignment-missingness", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    qc = pd.read_csv(args.uniform_qc, sep="\t", dtype=str).fillna("")
    missing = pd.read_csv(args.alignment_missingness, sep="\t", dtype=str).fillna("")
    merged = qc.merge(
        missing[["tree_sample_id", "missing_fraction", "alignment_qc_pass"]],
        on="tree_sample_id",
        how="inner",
        validate="one_to_one",
    )
    years = pd.to_numeric(merged["year"], errors="coerce")
    rescue = merged[
        merged["tree_role"].eq("focal")
        & merged["country_iso3"].isin(["AUS", "CHN", "JPN"])
        & years.lt(2020)
        & merged["uniform_qc_status"].eq("PASS")
        & merged["qc_input_type"].eq("paired_reads")
        & merged["alignment_qc_pass"].eq("False")
    ].copy()
    rescue["depth_rescue_reason"] = "pre2020_alignment_coverage_below_80pct_at_100k"
    rescue.to_csv(
        args.output_dir / "historical_depth_rescue_manifest.tsv",
        sep="\t",
        index=False,
    )
    report = {
        "n_depth_rescue_samples": len(rescue),
        "by_country": rescue["country_iso3"].value_counts().sort_index().to_dict(),
        "source_read_depth": 100000,
        "target_read_depth": 500000,
    }
    (args.output_dir / "historical_depth_rescue_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
