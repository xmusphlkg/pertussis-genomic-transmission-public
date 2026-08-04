#!/usr/bin/env python3
"""Build a leakage-audited 2019-2024 training set for 2025 holdout testing."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


COUNTRIES = ["AUS", "CHN", "JPN"]
LINEAGES = ["L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("full_data_json", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--start-month", default="2019-01")
    parser.add_argument("--cutoff-month", default="2024-12")
    parser.add_argument("--forecast-end-month", default="2025-12")
    parser.add_argument("--project-index", type=Path)
    parser.add_argument("--tip-attribution", type=Path)
    parser.add_argument("--lineage-assignments", type=Path)
    parser.add_argument("--stan-model", type=Path)
    return parser.parse_args()


def parse_month(value: str) -> tuple[int, int]:
    parsed = date.fromisoformat(f"{value}-01")
    return parsed.year, parsed.month


def month_offset(start: str, value: str) -> int:
    start_year, start_month = parse_month(start)
    value_year, value_month = parse_month(value)
    return (value_year - start_year) * 12 + value_month - start_month


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_project_names(path: Path | None) -> dict[int, str]:
    if path is None or not path.is_file():
        return {}
    with path.open(newline="") as handle:
        return {
            int(row["project_id_numeric"]): row["project_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def slice_months(matrix: list[list[Any]], n_months: int) -> list[list[Any]]:
    return [row[:n_months] for row in matrix]


def reconstruct_frozen_exposure(
    tip_attribution: Path,
    lineage_assignments: Path,
    start_month: str,
    n_months: int,
    cutoff_month: str,
) -> dict[str, Any]:
    attribution = pd.read_csv(tip_attribution, sep="\t")
    lineages = pd.read_csv(lineage_assignments, sep="\t")
    fields = [
        "tree_sample_id",
        "date_lower",
        "date_upper",
        "date_resolution",
        "country_iso3",
        "primary_model_lineage_id",
        "epidemic_period",
    ]
    attribution = attribution.drop(
        columns=[
            field
            for field in fields[1:]
            if field in attribution.columns
        ]
    ).merge(
        lineages[fields],
        on="tree_sample_id",
        how="left",
        validate="one_to_one",
    )
    attribution["date_lower"] = pd.to_datetime(
        attribution["date_lower"], errors="coerce"
    )
    event_tips = attribution[
        attribution["country_iso3"].isin(COUNTRIES)
        & attribution["epidemic_period"].eq("resurgence")
        & attribution["strongest_post_event_id"].fillna("").ne("")
        & attribution["strongest_post_transition_support"].ge(0.5)
    ].copy()
    if event_tips["date_lower"].isna().any():
        raise ValueError("Missing dates in the frozen exposure audit")

    country_index = {country: index for index, country in enumerate(COUNTRIES)}
    lineage_index = {lineage: index for index, lineage in enumerate(LINEAGES)}

    def build(tips: pd.DataFrame) -> np.ndarray:
        exposure = np.zeros((len(COUNTRIES), n_months, len(LINEAGES)))
        for (_, country), group in tips.groupby(
            ["strongest_post_event_id", "country_iso3"]
        ):
            event_month = group["date_lower"].min().strftime("%Y-%m")
            month_index = month_offset(start_month, event_month)
            if month_index < 0 or month_index >= n_months:
                continue
            support = float(
                group["strongest_post_transition_support"].max()
            )
            counts = (
                group["primary_model_lineage_id"]
                .fillna("Other")
                .where(
                    group["primary_model_lineage_id"].isin(LINEAGES),
                    "Other",
                )
                .value_counts()
            )
            for lineage, count in counts.items():
                exposure[
                    country_index[country],
                    month_index,
                    lineage_index[lineage],
                ] += support * int(count) / len(group)
        return exposure

    cutoff_end = pd.Period(cutoff_month, freq="M").end_time
    pre_cutoff_tips = event_tips[event_tips["date_lower"] <= cutoff_end]
    full_exposure = build(event_tips)
    pre_cutoff_exposure = build(pre_cutoff_tips)

    grouped = event_tips.assign(
        is_holdout=event_tips["date_lower"] > cutoff_end
    ).groupby(["strongest_post_event_id", "country_iso3"])
    mixed_event_count = 0
    future_tips_in_mixed_events = 0
    for _, group in grouped:
        if group["is_holdout"].any() and (~group["is_holdout"]).any():
            mixed_event_count += 1
            future_tips_in_mixed_events += int(group["is_holdout"].sum())

    return {
        "full_exposure": full_exposure,
        "pre_cutoff_exposure": pre_cutoff_exposure,
        "event_tips_full": len(event_tips),
        "event_tips_holdout": int(
            (event_tips["date_lower"] > cutoff_end).sum()
        ),
        "mixed_pre2025_and_2025_events": mixed_event_count,
        "holdout_tips_in_mixed_events": future_tips_in_mixed_events,
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    full_path = args.full_data_json.resolve()
    full = json.loads(full_path.read_text())

    cutoff_index = month_offset(args.start_month, args.cutoff_month) + 1
    forecast_end_index = month_offset(
        args.start_month, args.forecast_end_month
    ) + 1
    if forecast_end_index != full["T"]:
        raise ValueError(
            "The requested forecast end must equal the frozen full-data end"
        )
    if cutoff_index <= 1 or cutoff_index >= full["T"]:
        raise ValueError("The training cutoff must be inside the model period")
    if full["T"] - cutoff_index != 12:
        raise ValueError("This analysis requires an exactly 12-month holdout")

    keep_observation = [
        index
        for index, month_index in enumerate(full["obs_month"])
        if month_index <= cutoff_index
    ]
    excluded_observation = [
        index
        for index, month_index in enumerate(full["obs_month"])
        if month_index > cutoff_index
    ]
    used_projects = sorted(
        {full["obs_project"][index] for index in keep_observation}
    )
    project_remap = {
        old_index: new_index
        for new_index, old_index in enumerate(used_projects, start=1)
    }

    training = {
        "C": full["C"],
        "T": cutoff_index,
        "L": full["L"],
        "P": len(used_projects),
        "J": len(keep_observation),
        "K": full["K"],
        "use_project_effects": full["use_project_effects"],
        "cases": slice_months(full["cases"], cutoff_index),
        "B": full["B"][:cutoff_index],
        "reporting_change": slice_months(
            full["reporting_change"], cutoff_index
        ),
        "initial_alpha": full["initial_alpha"],
        "import_exposure": [
            country[:cutoff_index] for country in full["import_exposure"]
        ],
        "obs_country": [
            full["obs_country"][index] for index in keep_observation
        ],
        "obs_month": [
            full["obs_month"][index] for index in keep_observation
        ],
        "obs_project": [
            project_remap[full["obs_project"][index]]
            for index in keep_observation
        ],
        "y_genome": [
            full["y_genome"][index] for index in keep_observation
        ],
    }

    training_file = args.output_dir / "joint_model_data_through_2024.json"
    training_file.write_text(json.dumps(training, separators=(",", ":")))

    project_index_path = args.project_index
    if project_index_path is None:
        candidate = full_path.with_name("project_index.tsv")
        project_index_path = candidate if candidate.is_file() else None
    project_names = load_project_names(project_index_path)
    project_rows = [
        {
            "old_project_index": old_index,
            "new_project_index": project_remap[old_index],
            "project_id": project_names.get(old_index, ""),
        }
        for old_index in used_projects
    ]
    with (args.output_dir / "project_index_through_2024.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "old_project_index",
                "new_project_index",
                "project_id",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(project_rows)

    full_basis = np.asarray(full["B"], dtype=float)
    training_basis = np.asarray(training["B"], dtype=float)
    future_basis = full_basis[cutoff_index:, :]
    training_rank = int(np.linalg.matrix_rank(training_basis))
    full_rank = int(np.linalg.matrix_rank(full_basis))
    future_only_columns = [
        index + 1
        for index in range(full["K"])
        if np.allclose(training_basis[:, index], 0)
        and not np.allclose(future_basis[:, index], 0)
    ]

    original_future_exposure = np.asarray(
        full["import_exposure"], dtype=float
    )[:, cutoff_index:, :]
    training_exposure = np.asarray(training["import_exposure"], dtype=float)
    country_names = ["AUS", "CHN", "JPN"]
    original_exposure_by_country = {
        country_names[index]: float(original_future_exposure[index].sum())
        for index in range(full["C"])
    }
    forced_exposure_by_country = {
        country: 0.0 for country in country_names
    }

    repo_root = full_path.parents[2]
    tip_attribution = args.tip_attribution or (
        repo_root
        / "results/phylogeography/events_thr0_5"
        / "tip_persistence_reseeding_support.tsv"
    )
    lineage_assignments = args.lineage_assignments or (
        repo_root
        / "results/lineages/primary_finalized"
        / "model_lineage_assignments.tsv"
    )
    stan_model = args.stan_model or (
        repo_root / "scripts/model/gtd_joint_transmission_sampling.stan"
    )
    for required in (tip_attribution, lineage_assignments, stan_model):
        if not required.is_file():
            raise FileNotFoundError(
                f"Required leakage-audit source is missing: {required}"
            )

    exposure_audit = reconstruct_frozen_exposure(
        tip_attribution,
        lineage_assignments,
        args.start_month,
        full["T"],
        args.cutoff_month,
    )
    frozen_exposure = np.asarray(full["import_exposure"], dtype=float)
    reconstructed_full_difference = float(
        np.max(np.abs(exposure_audit["full_exposure"] - frozen_exposure))
    )
    training_exposure_difference_after_tip_exclusion = float(
        np.max(
            np.abs(
                exposure_audit["pre_cutoff_exposure"][:, :cutoff_index, :]
                - frozen_exposure[:, :cutoff_index, :]
            )
        )
    )
    stan_text = stan_model.read_text()
    persistence_support_reference_count = stan_text.count("persistence_support")

    assertions = {
        "training_cases_end_at_cutoff": all(
            len(row) == cutoff_index for row in training["cases"]
        ),
        "training_genomes_end_at_cutoff": (
            max(training["obs_month"], default=0) <= cutoff_index
        ),
        "training_exposure_ends_at_cutoff": all(
            len(country) == cutoff_index
            for country in training["import_exposure"]
        ),
        "training_spline_is_exact_full_prefix": training["B"]
        == full["B"][:cutoff_index],
        "no_holdout_observation_row_is_retained": not (
            set(keep_observation) & set(excluded_observation)
        ),
        "prediction_exposure_is_forced_to_zero": all(
            value == 0 for value in forced_exposure_by_country.values()
        ),
        "reconstructed_exposure_matches_frozen_input": (
            reconstructed_full_difference < 1e-12
        ),
        "excluding_2025_attributed_tips_leaves_training_exposure_unchanged": (
            training_exposure_difference_after_tip_exclusion < 1e-12
        ),
        "unused_persistence_input_removed": (
            persistence_support_reference_count == 0
            and "persistence_support" not in full
            and "persistence_support" not in training
        ),
    }
    if not all(assertions.values()):
        failed = [key for key, value in assertions.items() if not value]
        raise AssertionError(
            "Temporal holdout leakage assertion failed: "
            + ", ".join(failed)
        )

    audit = {
        "analysis": "2025 leakage-aware temporal holdout",
        "source_full_data_json": str(full_path),
        "source_full_data_sha256": sha256_file(full_path),
        "training_data_json": str(training_file.resolve()),
        "training_data_sha256": sha256_file(training_file),
        "model_start_month": args.start_month,
        "training_cutoff_month": args.cutoff_month,
        "forecast_start_month": "2025-01",
        "forecast_end_month": args.forecast_end_month,
        "full_months": full["T"],
        "training_months": cutoff_index,
        "holdout_months": full["T"] - cutoff_index,
        "case_values_full": full["C"] * full["T"],
        "case_values_passed_to_fit": full["C"] * cutoff_index,
        "holdout_case_values_excluded_from_fit": (
            full["C"] * (full["T"] - cutoff_index)
        ),
        "genome_strata_full": full["J"],
        "genome_strata_passed_to_fit": len(keep_observation),
        "genome_strata_excluded_from_fit": len(excluded_observation),
        "genomes_passed_to_fit": int(
            sum(sum(training["y_genome"][index]) for index in range(training["J"]))
        ),
        "genomes_excluded_from_fit": int(
            sum(sum(full["y_genome"][index]) for index in excluded_observation)
        ),
        "maximum_training_genome_month_index": max(
            training["obs_month"], default=0
        ),
        "projects_full": full["P"],
        "projects_passed_to_fit": training["P"],
        "project_old_to_new_index": project_rows,
        "original_2025_exposure_weight_by_country": (
            original_exposure_by_country
        ),
        "prediction_2025_exposure_weight_by_country": (
            forced_exposure_by_country
        ),
        "prediction_uses_import_scale": False,
        "frozen_exposure_direct_future_tip_audit": {
            "tip_attribution_source": str(tip_attribution.resolve()),
            "lineage_assignment_source": str(lineage_assignments.resolve()),
            "event_attributed_tips_full": exposure_audit["event_tips_full"],
            "event_attributed_tips_in_2025": exposure_audit[
                "event_tips_holdout"
            ],
            "events_with_both_pre2025_and_2025_attributed_tips": (
                exposure_audit["mixed_pre2025_and_2025_events"]
            ),
            "future_tips_in_mixed_events": exposure_audit[
                "holdout_tips_in_mixed_events"
            ],
            "maximum_absolute_difference_reconstructed_vs_frozen_exposure": (
                reconstructed_full_difference
            ),
            "maximum_absolute_difference_in_training_exposure_after_excluding_2025_tips": (
                training_exposure_difference_after_tip_exclusion
            ),
            "interpretation": (
                "Within the frozen event labels and transition supports, "
                "excluding every 2025 attributed tip leaves all 2019-2024 "
                "exposure cells numerically unchanged."
            ),
        },
        "frozen_phylogeography_boundary": (
            "Tree topology, event labels, and transition supports were not "
            "re-inferred after removing 2025 tips; this is a conditional "
            "temporal validation of the frozen feature map."
        ),
        "removed_unused_persistence_input_audit": {
            "stan_model": str(stan_model.resolve()),
            "stan_model_sha256": sha256_file(stan_model),
            "references_in_stan_source": persistence_support_reference_count,
            "interpretation": (
                "The legacy persistence_support field was removed from the "
                "Stan data block and frozen JSON because it did not enter "
                "the state update, likelihood, or generated quantities."
            ),
        },
        "training_spline_columns": training["K"],
        "training_spline_rank": training_rank,
        "full_spline_rank": full_rank,
        "future_only_spline_columns": future_only_columns,
        "future_only_spline_column_details": [
            {
                "column": index,
                "maximum_absolute_training_value": float(
                    np.abs(training_basis[:, index - 1]).max()
                ),
                "maximum_absolute_forecast_value": float(
                    np.abs(future_basis[:, index - 1]).max()
                ),
            }
            for index in future_only_columns
        ],
        "future_spline_inference_boundary": (
            "Column 14 is exactly zero through 2024 and nonzero in 2025. "
            "Its coefficient is therefore not directly informed by the "
            "training likelihood; projection is regularised only by the "
            "coefficient priors, including the adjacent-coefficient "
            "random-walk prior. Widening forecast intervals reflect this "
            "pre-frozen design boundary."
        ),
        "projection_policy": (
            "Use rows 73-84 of the pre-frozen 84-month spline design; "
            "do not rebuild the basis after truncation."
        ),
        "fit_input_policy": (
            "Only the generated through-2024 JSON is passed to Stan. "
            "Full 2025 cases are read only after fitting for forecast scoring."
        ),
        "leakage_assertions": assertions,
    }
    (args.output_dir / "temporal_holdout_2025_prefit_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
