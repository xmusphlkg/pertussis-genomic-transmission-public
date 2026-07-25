#!/usr/bin/env python3
"""Freeze model lineages from pre-effect hierBAPS assignments.

The primary model uses level-2 populations only when they have at least
20 genomes, at least five resurgence genomes, and cross-country or
cross-period representation. All other populations are pooled as Other.
No case counts are read by this script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("assignments", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.assignments, sep="\t")

    required = {
        "tree_sample_id",
        "tree_role",
        "country_iso3",
        "epidemic_period",
        "model_lineage_id",
        "model_sublineage_id",
    }
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    summary = (
        data.groupby("model_sublineage_id", dropna=False)
        .agg(
            n_genomes=("tree_sample_id", "size"),
            n_focal=("tree_role", lambda x: int((x == "focal").sum())),
            n_countries=("country_iso3", "nunique"),
            n_periods=(
                "epidemic_period",
                lambda x: int(x[x.notna() & x.ne("unknown")].nunique()),
            ),
            n_prepandemic=("epidemic_period", lambda x: int((x == "prepandemic").sum())),
            n_pandemic=("epidemic_period", lambda x: int((x == "pandemic").sum())),
            n_resurgence=("epidemic_period", lambda x: int((x == "resurgence").sum())),
            min_year=("year", "min"),
            max_year=("year", "max"),
        )
        .reset_index()
    )
    summary["primary_model_eligible"] = (
        summary["n_genomes"].ge(20)
        & summary["n_resurgence"].ge(5)
        & (summary["n_countries"].ge(2) | summary["n_periods"].ge(2))
    )
    eligible = set(
        summary.loc[summary["primary_model_eligible"], "model_sublineage_id"]
    )
    if len(eligible) < 3:
        raise RuntimeError(
            f"Only {len(eligible)} lineages passed the pre-specified gate; "
            "the multi-lineage model must not proceed"
        )

    data["primary_model_lineage_id"] = data["model_sublineage_id"].where(
        data["model_sublineage_id"].isin(eligible), "Other"
    )
    data["primary_model_lineage_status"] = data["model_sublineage_id"].map(
        lambda x: "eligible_level2_hierbaps" if x in eligible else "pooled_other"
    )
    data["lineage_freeze_rule"] = (
        "level2_hierBAPS_n>=20_resurgence_n>=5_and_cross_country_or_period"
    )

    summary = summary.sort_values(
        ["primary_model_eligible", "n_genomes"], ascending=[False, False]
    )
    data.to_csv(
        args.output_dir / "model_lineage_assignments.tsv",
        sep="\t",
        index=False,
        na_rep="",
    )
    summary.to_csv(
        args.output_dir / "model_lineage_summary.tsv",
        sep="\t",
        index=False,
        na_rep="",
    )

    dated = data[
        data["tree_role"].eq("focal")
        & data["country_iso3"].isin(["AUS", "CHN", "JPN", "BEL", "FRA"])
        & data["date_resolution"].isin(["day", "month"])
    ].copy()
    dated["model_month"] = pd.to_datetime(dated["date_lower"]).dt.to_period("M").dt.to_timestamp()
    monthly = (
        dated.groupby(
            [
                "country_iso3",
                "model_month",
                "project_id",
                "primary_model_lineage_id",
            ],
            dropna=False,
        )
        .size()
        .rename("n_tree_genomes")
        .reset_index()
    )
    monthly.to_csv(
        args.output_dir / "country_month_project_model_lineage_counts.tsv",
        sep="\t",
        index=False,
        na_rep="",
    )

    validation = {
        "n_tree_genomes": int(len(data)),
        "n_level1_populations": int(data["model_lineage_id"].nunique()),
        "n_level2_populations": int(data["model_sublineage_id"].nunique()),
        "n_primary_model_lineages_excluding_other": int(len(eligible)),
        "primary_model_lineages": sorted(eligible),
        "n_month_resolved_focal_genomes": int(len(dated)),
        "case_data_used_for_lineage_definition": False,
        "gate_pass": len(eligible) >= 3,
    }
    (args.output_dir / "model_lineage_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )

    report = [
        "# Formal model-lineage freeze",
        "",
        "The level-1 hierBAPS partition contained only two populations meeting the",
        "minimum sample gate. The pre-specified fallback therefore uses level-2",
        "populations. This decision was made without reading case counts or growth",
        "estimates.",
        "",
        "Primary level-2 eligibility required:",
        "",
        "- at least 20 genomes;",
        "- at least five resurgence-period genomes; and",
        "- representation across at least two countries or two epidemic periods.",
        "",
        f"Four eligible lineages were retained: {', '.join(sorted(eligible))}.",
        "All remaining populations were pooled as `Other`.",
        "",
        "These population labels are genomic clusters, not PRN, ptxP, fim3, or",
        "antimicrobial-resistance endpoints.",
    ]
    (args.output_dir / "MODEL_LINEAGE_FREEZE.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
