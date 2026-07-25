#!/usr/bin/env python3
"""Audit whether rebuilding the transmission tree changes frozen lineages."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("previous_assignments", type=Path)
    parser.add_argument("rebuilt_assignments", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    columns = [
        "tree_sample_id",
        "model_sublineage_id",
        "primary_model_lineage_id",
    ]
    old = pd.read_csv(args.previous_assignments, sep="\t", usecols=columns)
    new = pd.read_csv(args.rebuilt_assignments, sep="\t", usecols=columns)
    shared = old.merge(
        new,
        on="tree_sample_id",
        suffixes=("_previous", "_rebuilt"),
        validate="one_to_one",
    )
    if shared.empty:
        raise RuntimeError("The previous and rebuilt cohorts have no shared samples")

    contingency = pd.crosstab(
        shared["primary_model_lineage_id_previous"],
        shared["primary_model_lineage_id_rebuilt"],
        dropna=False,
    )
    contingency.to_csv(
        args.output_dir / "primary_lineage_rebuild_contingency.tsv", sep="\t"
    )
    shared.to_csv(
        args.output_dir / "shared_sample_lineage_comparison.tsv",
        sep="\t",
        index=False,
    )

    exact = shared["primary_model_lineage_id_previous"].eq(
        shared["primary_model_lineage_id_rebuilt"]
    )
    result = {
        "n_previous": int(len(old)),
        "n_rebuilt": int(len(new)),
        "n_shared": int(len(shared)),
        "level2_adjusted_rand_index": float(
            adjusted_rand_score(
                shared["model_sublineage_id_previous"],
                shared["model_sublineage_id_rebuilt"],
            )
        ),
        "primary_lineage_exact_agreement": float(exact.mean()),
        "n_primary_lineage_disagreements": int((~exact).sum()),
        "interpretation": (
            "A value of 1 denotes identical formal model-lineage labels among "
            "shared samples; expanded-cohort samples are audited separately."
        ),
    }
    (args.output_dir / "lineage_rebuild_validation.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
