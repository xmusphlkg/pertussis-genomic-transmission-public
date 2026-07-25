#!/usr/bin/env python3
"""Calibrate Ct-dependent sequencing/profile success in the Australian cohort."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("sampling_records", type=Path)
    p.add_argument("output_dir", type=Path)
    return p.parse_args()


def auc_rank(y: np.ndarray, score: np.ndarray) -> float:
    order = pd.Series(score).rank(method="average").to_numpy()
    n1 = y.sum()
    n0 = len(y) - n1
    return float((order[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.sampling_records, sep="\t")
    data["ct"] = pd.to_numeric(data["ct_is481"], errors="coerce")
    data["success"] = data["profile_complete"].astype(int)
    data["specimen_group"] = np.where(
        data["specimen_type"].str.contains("swab|nose|throat|np|op", case=False, na=False),
        "upper_respiratory_swab",
        "other_or_unknown",
    )
    fit_data = data.dropna(subset=["ct", "success"]).copy()
    fit_data["ct_centered"] = fit_data["ct"] - fit_data["ct"].mean()

    model = smf.glm(
        "success ~ ct_centered + C(specimen_group)",
        data=fit_data,
        family=sm.families.Binomial(),
    ).fit()
    null_model = smf.glm(
        "success ~ 1", data=fit_data, family=sm.families.Binomial()
    ).fit()
    pred = model.predict(fit_data)
    brier = float(np.mean((fit_data["success"].to_numpy() - pred) ** 2))
    auc = auc_rank(fit_data["success"].to_numpy(), pred.to_numpy())

    params = pd.DataFrame(
        {
            "term": model.params.index,
            "estimate_log_odds": model.params.values,
            "standard_error": model.bse.values,
            "odds_ratio": np.exp(model.params.values),
            "ci_lower_odds_ratio": np.exp(model.conf_int()[0].values),
            "ci_upper_odds_ratio": np.exp(model.conf_int()[1].values),
            "p_value": model.pvalues.values,
        }
    )
    params.to_csv(args.output_dir / "australia_sampling_model_parameters.tsv", sep="\t", index=False)

    grid = pd.DataFrame(
        {
            "ct": np.linspace(10, 40, 121),
            "specimen_group": "upper_respiratory_swab",
        }
    )
    grid["ct_centered"] = grid["ct"] - fit_data["ct"].mean()
    prediction = model.get_prediction(grid).summary_frame(alpha=0.05)
    grid["success_probability"] = prediction["mean"].to_numpy()
    grid["ci_lower"] = prediction["mean_ci_lower"].to_numpy()
    grid["ci_upper"] = prediction["mean_ci_upper"].to_numpy()
    grid.to_csv(args.output_dir / "australia_ct_success_curve.tsv", sep="\t", index=False)

    bins = pd.cut(fit_data["ct"], bins=[0, 20, 25, 30, 35, 40, 100], include_lowest=True)
    calibration = (
        fit_data.assign(ct_bin=bins)
        .groupby("ct_bin", observed=True)
        .agg(
            n_specimens=("success", "size"),
            n_success=("success", "sum"),
            observed_success_probability=("success", "mean"),
            mean_predicted_probability=("ct_centered", lambda x: float(model.predict(fit_data.loc[x.index]).mean())),
        )
        .reset_index()
    )
    calibration["ct_bin"] = calibration["ct_bin"].astype(str)
    calibration.to_csv(args.output_dir / "australia_sampling_calibration.tsv", sep="\t", index=False)

    report = {
        "n_direct_specimens": int(len(data)),
        "n_with_numeric_ct": int(len(fit_data)),
        "n_profile_complete": int(fit_data["success"].sum()),
        "observed_success_fraction": float(fit_data["success"].mean()),
        "ct_mean": float(fit_data["ct"].mean()),
        "ct_sd": float(fit_data["ct"].std()),
        "auc": auc,
        "brier_score": brier,
        "likelihood_ratio_vs_intercept_only": float(2 * (model.llf - null_model.llf)),
        "ct_odds_ratio_per_cycle": float(np.exp(model.params["ct_centered"])),
        "ct_odds_ratio_ci95": [
            float(np.exp(model.conf_int().loc["ct_centered", 0])),
            float(np.exp(model.conf_int().loc["ct_centered", 1])),
        ],
        "use_in_joint_model": (
            "External calibration for Ct-dependent sequencing success and "
            "sampling-sensitivity analyses; not assumed to transfer exactly "
            "to other countries"
        ),
    }
    (args.output_dir / "australia_sampling_model_validation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
