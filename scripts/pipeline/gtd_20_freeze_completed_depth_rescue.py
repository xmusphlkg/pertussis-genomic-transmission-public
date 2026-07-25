#!/usr/bin/env python3
"""Freeze successfully completed high-depth replacements for the final SKA tree."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rescue-manifest", type=Path, required=True)
    parser.add_argument("--read-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.rescue_manifest, sep="\t", dtype=str).fillna("")
    rows: list[dict[str, str]] = []
    for _, row in manifest.iterrows():
        done = args.read_root / f"{row['run_accession']}.subset500000.done.json"
        if not done.is_file():
            continue
        payload = json.loads(done.read_text(encoding="utf-8"))
        if payload.get("status") != "PASS":
            continue
        record = row.to_dict()
        record["depth_rescue_r1"] = payload["r1_path"]
        record["depth_rescue_r2"] = payload["r2_path"]
        record["depth_rescue_records"] = str(payload["paired_records"])
        rows.append(record)
    frozen = pd.DataFrame(rows)
    frozen.to_csv(
        args.output_dir / "completed_depth_rescue_manifest.tsv",
        sep="\t",
        index=False,
    )
    with (args.output_dir / "completed_depth_rescue_ska_input.tsv").open(
        "w", encoding="utf-8"
    ) as handle:
        for _, row in frozen.iterrows():
            handle.write(
                f"{row['tree_sample_id']}\t{row['depth_rescue_r1']}\t"
                f"{row['depth_rescue_r2']}\n"
            )
    (args.output_dir / "completed_depth_rescue_sample_ids.txt").write_text(
        "\n".join(frozen["tree_sample_id"]) + "\n", encoding="utf-8"
    )
    report = {
        "n_completed_depth_rescue": len(frozen),
        "by_country": frozen["country_iso3"].value_counts().sort_index().to_dict(),
        "read_pairs_per_file": 500000,
    }
    (args.output_dir / "completed_depth_rescue_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
