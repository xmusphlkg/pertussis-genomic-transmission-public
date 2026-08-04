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
    p.add_argument(
        "--transition-threshold",
        type=float,
        default=0.5,
        help="Minimum strongest post-transition support used to build exposure.",
    )
    p.add_argument(
        "--exposure-time-rule",
        choices=("lower", "midpoint", "interval-uniform"),
        default="lower",
        help=(
            "Map the interval-censored earliest attributed sample to its lower "
            "bound, midpoint, or uniformly across all intersecting months."
        ),
    )
    p.add_argument(
        "--initial-prior-mode",
        choices=("historical", "symmetric"),
        default="historical",
        help=(
            "Construct the initial-lineage prior from samples strictly before "
            "the model period, or use a weak symmetric Dirichlet(0.5) prior."
        ),
    )
    p.add_argument(
        "--initial-prior-end-year",
        type=int,
        default=START.year - 1,
        help=(
            "Last sampling year eligible for the historical initial-lineage "
            "prior. The default (2018) is strictly before the 2019 model start."
        ),
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not 0 <= args.transition_threshold <= 1:
        raise ValueError("--transition-threshold must be between 0 and 1")
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
    lineages["joint_model_lineage_id"] = lineages[
        "primary_model_lineage_id"
    ].where(lineages["primary_model_lineage_id"].isin(LINEAGES), "Other")
    lineages["model_month"] = (
        pd.to_datetime(lineages["date_lower"], errors="coerce")
        .dt.to_period("M")
        .dt.to_timestamp()
    )

    if args.initial_prior_end_year >= START.year:
        raise ValueError(
            "--initial-prior-end-year must be strictly before the model start year"
        )

    # The initial state is separated temporally from every monthly model
    # observation. Historical priors may include interval-censored year-only
    # samples, but only when their reported year precedes the model period.
    pre = lineages.iloc[0:0].copy()
    if args.initial_prior_mode == "historical":
        pre = lineages[
            lineages["country_iso3"].isin(COUNTRIES)
            & lineages["year"].le(args.initial_prior_end_year)
        ].copy()
    initial_alpha = np.full((len(COUNTRIES), len(LINEAGES)), 0.5)
    for (country, lineage), n in pre.groupby(
        ["country_iso3", "joint_model_lineage_id"]
    ).size().items():
        initial_alpha[country_index[country] - 1, lineage_index[lineage] - 1] += int(n)

    initial_prior_rows = []
    for country in COUNTRIES:
        for lineage in LINEAGES:
            count = int(
                (
                    pre["country_iso3"].eq(country)
                    & pre["joint_model_lineage_id"].eq(lineage)
                ).sum()
            )
            initial_prior_rows.append(
                {
                    "country_iso3": country,
                    "primary_model_lineage_id": lineage,
                    "initial_prior_mode": args.initial_prior_mode,
                    "initial_prior_end_year": args.initial_prior_end_year,
                    "n_historical_tips": count,
                    "dirichlet_alpha": 0.5 + count,
                }
            )
    initial_prior_table = pd.DataFrame(initial_prior_rows)

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
            "date_upper",
            "date_resolution",
            "country_iso3",
            "primary_model_lineage_id",
            "epidemic_period",
        ]
    ]
    attribution = attribution.drop(
        columns=[
            x
            for x in [
                "date_lower",
                "date_upper",
                "date_resolution",
                "country_iso3",
                "primary_model_lineage_id",
                "epidemic_period",
            ]
            if x in attribution.columns
        ]
    ).merge(tip_fields, on="tree_sample_id", how="left", validate="one_to_one")
    attribution["date_lower"] = pd.to_datetime(attribution["date_lower"], errors="coerce")
    attribution["date_upper"] = pd.to_datetime(attribution["date_upper"], errors="coerce")
    event_tips = attribution[
        attribution["country_iso3"].isin(COUNTRIES)
        & attribution["epidemic_period"].eq("resurgence")
        & attribution["strongest_post_event_id"].fillna("").ne("")
        & attribution["strongest_post_transition_support"].ge(
            args.transition_threshold
        )
    ].copy()

    # One event contributes its transition support once, split across model
    # lineages according to attributed descendants. The event-time proxy is the
    # interval for the earliest attributed sample: its lower and upper bounds
    # are min(date_lower) and min(date_upper), respectively. Sensitivity inputs
    # place this proxy at the lower bound, at the interval midpoint, or spread
    # the same total event weight uniformly across all intersecting months.
    event_rows = []
    for (event_id, country), group in event_tips.groupby(
        ["strongest_post_event_id", "country_iso3"]
    ):
        if group["date_lower"].isna().any() or group["date_upper"].isna().any():
            raise ValueError(f"Missing date bounds among tips attributed to {event_id}")
        event_date_lower = group["date_lower"].min()
        event_date_upper = group["date_upper"].min()
        if event_date_upper < event_date_lower:
            raise ValueError(f"Inconsistent earliest-sample interval for {event_id}")

        if args.exposure_time_rule == "lower":
            allocated_months = [
                event_date_lower.to_period("M").to_timestamp()
            ]
        elif args.exposure_time_rule == "midpoint":
            midpoint = event_date_lower + (event_date_upper - event_date_lower) / 2
            allocated_months = [midpoint.to_period("M").to_timestamp()]
        else:
            allocated_months = list(
                pd.date_range(
                    event_date_lower.to_period("M").to_timestamp(),
                    event_date_upper.to_period("M").to_timestamp(),
                    freq="MS",
                )
            )
        allocated_months = [
            month for month in allocated_months if START <= month <= END
        ]
        if not allocated_months:
            continue

        support = float(group["strongest_post_transition_support"].max())
        counts = group["primary_model_lineage_id"].fillna("Other").value_counts()
        total = int(counts.sum())
        for lineage, n in counts.items():
            lineage = lineage if lineage in lineage_index else "Other"
            total_lineage_weight = support * int(n) / total
            allocated_lineage_weight = total_lineage_weight / len(allocated_months)
            for month in allocated_months:
                event_rows.append(
                    {
                        "event_id": event_id,
                        "country_iso3": country,
                        "model_month": month,
                        "event_date_lower": event_date_lower,
                        "event_date_upper": event_date_upper,
                        "exposure_time_rule": args.exposure_time_rule,
                        "n_allocated_months": len(allocated_months),
                        "primary_model_lineage_id": lineage,
                        "transition_support": support,
                        "total_lineage_event_weight": total_lineage_weight,
                        "lineage_event_weight": allocated_lineage_weight,
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
        "import_exposure": exposure.tolist(),
        "obs_country": obs["country_id"].astype(int).tolist(),
        "obs_month": obs["month_id"].astype(int).tolist(),
        "obs_project": obs["project_id_numeric"].astype(int).tolist(),
        "y_genome": y_genome.tolist(),
    }
    (args.output_dir / "joint_model_data.json").write_text(
        json.dumps(stan_data, separators=(",", ":"))
    )
    initial_prior_table.to_csv(
        args.output_dir / "initial_lineage_prior.tsv", sep="\t", index=False
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
        "initial_prior_mode": args.initial_prior_mode,
        "initial_prior_end_year": args.initial_prior_end_year,
        "initial_prior_strictly_precedes_model": (
            args.initial_prior_end_year < START.year
        ),
        "n_historical_tips_in_initial_prior": int(len(pre)),
        "n_case_months_per_country": len(months),
        "n_genome_observation_strata": len(obs),
        "n_genomes_in_observation_model": int(y_genome.sum()),
        "n_projects": len(projects),
        "use_project_effects": not args.disable_project_effects,
        "transition_threshold": args.transition_threshold,
        "exposure_time_rule": args.exposure_time_rule,
        "n_phylogeographic_event_lineage_rows": len(event_table),
        "n_unique_phylogeographic_events": int(event_table["event_id"].nunique()),
        "n_unique_phylogeographic_event_lineages": int(
            event_table[["event_id", "primary_model_lineage_id"]]
            .drop_duplicates()
            .shape[0]
        ),
        "total_phylogeographic_exposure_weight": float(
            event_table["lineage_event_weight"].sum()
        ),
        "case_counts_used_for_lineage_definition": False,
        "exact_tmrca_used": False,
    }
    (args.output_dir / "joint_model_data_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
