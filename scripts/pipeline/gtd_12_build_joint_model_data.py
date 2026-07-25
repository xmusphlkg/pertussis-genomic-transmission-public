#!/usr/bin/env python3
"""Build frozen Stan inputs for the modular transmission model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import dmatrix

COUNTRIES = ["AUS", "CHN", "JPN"]
LINEAGES = ["L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other"]
START = pd.Timestamp("2019-01-01")
END = pd.Timestamp("2025-12-01")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("cases", type=Path)
    p.add_argument("lineage_assignments", type=Path)
    p.add_argument("tip_attribution", type=Path)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--disable-project-effects", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    months = pd.date_range(START, END, freq="MS")
    country_index = {x: i + 1 for i, x in enumerate(COUNTRIES)}
    lineage_index = {x: i + 1 for i, x in enumerate(LINEAGES)}

    case = pd.read_csv(args.cases, sep="\t")
    case["model_month"] = pd.to_datetime(case["model_month"])
    case = case[
        case["country_iso3"].isin(COUNTRIES)
        & case["model_month"].between(START, END)
    ].copy()
    case_grid = (
        pd.MultiIndex.from_product(
            [COUNTRIES, months], names=["country_iso3", "model_month"]
        )
        .to_frame(index=False)
        .merge(
            case[["country_iso3", "model_month", "cases"]],
            on=["country_iso3", "model_month"],
            how="left",
            validate="one_to_one",
        )
    )
    if case_grid["cases"].isna().any():
        raise ValueError("The frozen 2019-2025 case grid contains missing values")
    case_matrix = (
        case_grid.pivot(index="country_iso3", columns="model_month", values="cases")
        .reindex(index=COUNTRIES, columns=months)
        .round()
        .astype(int)
        .to_numpy()
    )

    lineages = pd.read_csv(args.lineage_assignments, sep="\t")
    lineages["primary_model_lineage_id"] = lineages[
        "primary_model_lineage_id"
    ].fillna("Other")
    lineages["model_month"] = (
        pd.to_datetime(lineages["date_lower"], errors="coerce")
        .dt.to_period("M")
        .dt.to_timestamp()
    )

    # Initial lineage proportions use all pre-2020 focal tree tips, including
    # interval-censored year-only samples. These are priors, not monthly counts.
    pre = lineages[
        lineages["country_iso3"].isin(COUNTRIES)
        & lineages["year"].le(2019)
    ]
    initial_alpha = np.full((len(COUNTRIES), len(LINEAGES)), 0.5)
    for (country, lineage), n in pre.groupby(
        ["country_iso3", "primary_model_lineage_id"]
    ).size().items():
        if lineage not in lineage_index:
            lineage = "Other"
        initial_alpha[country_index[country] - 1, lineage_index[lineage] - 1] += int(n)

    # Monthly genomic observations are restricted to exact day/month dates and
    # are grouped by country, month and public project.
    observed = lineages[
        lineages["country_iso3"].isin(COUNTRIES)
        & lineages["date_resolution"].isin(["day", "month"])
        & lineages["model_month"].between(START, END)
    ].copy()
    observed["project_id"] = observed["project_id"].fillna("UNSPECIFIED")
    observed.loc[
        ~observed["primary_model_lineage_id"].isin(LINEAGES),
        "primary_model_lineage_id",
    ] = "Other"
    projects = sorted(observed["project_id"].unique())
    project_index = {x: i + 1 for i, x in enumerate(projects)}

    obs_rows = []
    for (country, month, project), group in observed.groupby(
        ["country_iso3", "model_month", "project_id"], sort=True
    ):
        counts = (
            group["primary_model_lineage_id"]
            .value_counts()
            .reindex(LINEAGES, fill_value=0)
            .astype(int)
            .tolist()
        )
        obs_rows.append(
            {
                "country_iso3": country,
                "model_month": month,
                "project_id": project,
                "country_id": country_index[country],
                "month_id": months.get_loc(month) + 1,
                "project_id_numeric": project_index[project],
                **{f"n_{lineage}": counts[i] for i, lineage in enumerate(LINEAGES)},
            }
        )
    obs = pd.DataFrame(obs_rows)
    y_genome = obs[[f"n_{x}" for x in LINEAGES]].astype(int).to_numpy()

    attribution = pd.read_csv(args.tip_attribution, sep="\t")
    tip_fields = lineages[
        [
            "tree_sample_id",
            "date_lower",
            "country_iso3",
            "primary_model_lineage_id",
            "epidemic_period",
        ]
    ]
    attribution = attribution.drop(
        columns=[
            x
            for x in [
                "country_iso3",
                "primary_model_lineage_id",
                "epidemic_period",
            ]
            if x in attribution.columns
        ]
    ).merge(tip_fields, on="tree_sample_id", how="left", validate="one_to_one")
    attribution["date_lower"] = pd.to_datetime(attribution["date_lower"], errors="coerce")
    event_tips = attribution[
        attribution["country_iso3"].isin(COUNTRIES)
        & attribution["epidemic_period"].eq("resurgence")
        & attribution["strongest_post_event_id"].fillna("").ne("")
        & attribution["strongest_post_transition_support"].ge(0.5)
    ].copy()

    resurgence_support = attribution[
        attribution["country_iso3"].isin(COUNTRIES)
        & attribution["epidemic_period"].eq("resurgence")
    ].copy()
    country_fallback = (
        resurgence_support.groupby("country_iso3")["local_persistence_support"].mean()
    )
    persistence_support = np.zeros((len(COUNTRIES), len(LINEAGES)))
    for country in COUNTRIES:
        for lineage in LINEAGES:
            values = resurgence_support.loc[
                resurgence_support["country_iso3"].eq(country)
                & resurgence_support["primary_model_lineage_id"].fillna("Other").eq(lineage),
                "local_persistence_support",
            ]
            value = (
                float(values.mean())
                if len(values)
                else float(country_fallback.get(country, 0.5))
            )
            persistence_support[
                country_index[country] - 1, lineage_index[lineage] - 1
            ] = np.clip(value, 0.01, 0.99)

    # One event contributes its transition support once, split across model
    # lineages according to attributed descendants. The first observed month
    # supplies an interval-censored introduction proxy, not an exact event date.
    event_rows = []
    for (event_id, country), group in event_tips.groupby(
        ["strongest_post_event_id", "country_iso3"]
    ):
        month = group["date_lower"].min().to_period("M").to_timestamp()
        if month < START or month > END:
            continue
        support = float(group["strongest_post_transition_support"].max())
        counts = group["primary_model_lineage_id"].fillna("Other").value_counts()
        total = int(counts.sum())
        for lineage, n in counts.items():
            lineage = lineage if lineage in lineage_index else "Other"
            event_rows.append(
                {
                    "event_id": event_id,
                    "country_iso3": country,
                    "model_month": month,
                    "primary_model_lineage_id": lineage,
                    "transition_support": support,
                    "lineage_event_weight": support * int(n) / total,
                    "n_attributed_resurgence_tips": int(n),
                }
            )
    event_table = pd.DataFrame(event_rows)
    exposure = np.zeros((len(COUNTRIES), len(months), len(LINEAGES)))
    for row in event_table.itertuples(index=False):
        exposure[
            country_index[row.country_iso3] - 1,
            months.get_loc(row.model_month),
            lineage_index[row.primary_model_lineage_id] - 1,
        ] += row.lineage_event_weight

    # Cubic B-spline basis for the shared country transmission intensity.
    time = np.linspace(0.0, 1.0, len(months))
    basis = np.asarray(
        dmatrix(
            "bs(x, df=14, degree=3, include_intercept=True) - 1",
            {"x": time},
            return_type="dataframe",
        )
    )

    # Known reporting-platform indicator. Only China changes inside this
    # analysis window (2024); AUS/JPN changes predate it.
    reporting_change = np.zeros((len(COUNTRIES), len(months)))
    reporting_change[country_index["CHN"] - 1, months >= pd.Timestamp("2024-01-01")] = 1

    stan_data = {
        "C": len(COUNTRIES),
        "T": len(months),
        "L": len(LINEAGES),
        "P": len(projects),
        "J": len(obs),
        "K": int(basis.shape[1]),
        "use_project_effects": 0 if args.disable_project_effects else 1,
        "cases": case_matrix.tolist(),
        "B": basis.tolist(),
        "reporting_change": reporting_change.tolist(),
        "initial_alpha": initial_alpha.tolist(),
        "persistence_support": persistence_support.tolist(),
        "import_exposure": exposure.tolist(),
        "obs_country": obs["country_id"].astype(int).tolist(),
        "obs_month": obs["month_id"].astype(int).tolist(),
        "obs_project": obs["project_id_numeric"].astype(int).tolist(),
        "y_genome": y_genome.tolist(),
    }
    (args.output_dir / "joint_model_data.json").write_text(
        json.dumps(stan_data, separators=(",", ":"))
    )
    obs.to_csv(args.output_dir / "genome_observation_strata.tsv", sep="\t", index=False)
    event_table.to_csv(args.output_dir / "monthly_import_events.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "model_month": months,
            **{f"basis_{i + 1}": basis[:, i] for i in range(basis.shape[1])},
        }
    ).to_csv(args.output_dir / "transmission_spline_basis.tsv", sep="\t", index=False)
    pd.DataFrame(
        {
            "project_id_numeric": range(1, len(projects) + 1),
            "project_id": projects,
        }
    ).to_csv(args.output_dir / "project_index.tsv", sep="\t", index=False)

    validation = {
        "countries": COUNTRIES,
        "lineages": LINEAGES,
        "start_month": str(START.date()),
        "end_month": str(END.date()),
        "n_case_months_per_country": len(months),
        "n_genome_observation_strata": len(obs),
        "n_genomes_in_observation_model": int(y_genome.sum()),
        "n_projects": len(projects),
        "use_project_effects": not args.disable_project_effects,
        "n_phylogeographic_event_lineage_rows": len(event_table),
        "n_unique_phylogeographic_events": int(event_table["event_id"].nunique()),
        "case_counts_used_for_lineage_definition": False,
        "exact_tmrca_used": False,
    }
    (args.output_dir / "joint_model_data_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
