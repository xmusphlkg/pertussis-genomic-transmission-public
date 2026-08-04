from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/model_temporal_holdout_2025"


def test_temporal_holdout_training_input_is_leakage_free() -> None:
    full = json.loads(
        (ROOT / "data/model_inputs/joint_model_data.json").read_text()
    )
    training = json.loads(
        (RESULTS / "joint_model_data_through_2024.json").read_text()
    )
    audit = json.loads(
        (
            RESULTS / "temporal_holdout_2025_prefit_audit.json"
        ).read_text()
    )

    assert training["T"] == 72
    assert full["T"] == 84
    assert all(len(row) == 72 for row in training["cases"])
    assert training["cases"] == [row[:72] for row in full["cases"]]
    assert training["B"] == full["B"][:72]
    assert all(len(country) == 72 for country in training["import_exposure"])
    assert max(training["obs_month"]) <= 72
    assert training["J"] == 60
    assert sum(map(sum, training["y_genome"])) == 220
    assert training["P"] == 9

    assert audit["holdout_case_values_excluded_from_fit"] == 36
    assert audit["genome_strata_excluded_from_fit"] == 12
    assert audit["genomes_excluded_from_fit"] == 41
    assert audit["training_spline_rank"] == 13
    assert audit["full_spline_rank"] == 14
    assert audit["future_only_spline_columns"] == [14]
    assert all(audit["leakage_assertions"].values())
    assert audit["original_2025_exposure_weight_by_country"]["JPN"] > 0
    assert set(
        audit["prediction_2025_exposure_weight_by_country"].values()
    ) == {0.0}
    direct_exposure_audit = audit[
        "frozen_exposure_direct_future_tip_audit"
    ]
    assert direct_exposure_audit["event_attributed_tips_in_2025"] == 43
    assert (
        direct_exposure_audit[
            "maximum_absolute_difference_reconstructed_vs_frozen_exposure"
        ]
        == 0
    )
    assert (
        direct_exposure_audit[
            "maximum_absolute_difference_in_training_exposure_after_excluding_2025_tips"
        ]
        == 0
    )
    assert audit["removed_unused_persistence_input_audit"][
        "references_in_stan_source"
    ] == 0
    assert "not re-inferred" in audit["frozen_phylogeography_boundary"]


def test_temporal_holdout_forecasts_and_scoring_contract() -> None:
    monthly = pd.read_csv(
        RESULTS / "temporal_holdout_2025_monthly_forecasts.tsv",
        sep="\t",
    )
    metrics = pd.read_csv(
        RESULTS / "temporal_holdout_2025_metrics.tsv",
        sep="\t",
    )
    diagnostics = json.loads(
        (
            RESULTS / "temporal_holdout_2025_sampling_diagnostics.json"
        ).read_text()
    )
    audit = json.loads(
        (
            RESULTS / "temporal_holdout_2025_leakage_audit.json"
        ).read_text()
    )

    methods = {
        "joint_model_expected",
        "joint_model_negative_binomial",
        "seasonal_naive",
    }
    assert len(monthly) == 108
    assert set(monthly["forecast_method"]) == methods
    assert set(monthly["country_iso3"]) == {"AUS", "CHN", "JPN"}
    assert set(monthly["forecast_month"]) == {
        f"2025-{month:02d}" for month in range(1, 13)
    }
    assert monthly.groupby(
        ["forecast_method", "country_iso3"]
    ).size().eq(12).all()

    ordered_quantiles = [
        "lower_95",
        "lower_90",
        "lower_80",
        "lower_50",
        "median",
        "upper_50",
        "upper_80",
        "upper_90",
        "upper_95",
    ]
    quantile_values = monthly[ordered_quantiles]
    assert (
        quantile_values.diff(axis=1).iloc[:, 1:].ge(-1e-10).all().all()
    )
    assert monthly["point_forecast"].ge(0).all()

    observed_by_method = monthly.pivot(
        index=["country_iso3", "forecast_month"],
        columns="forecast_method",
        values="observed_cases",
    )
    assert observed_by_method.nunique(axis=1).eq(1).all()

    full = json.loads(
        (ROOT / "data/model_inputs/joint_model_data.json").read_text()
    )
    seasonal = monthly[
        monthly["forecast_method"].eq("seasonal_naive")
    ].sort_values(["country_iso3", "horizon_month"])
    country_index = {"AUS": 0, "CHN": 1, "JPN": 2}
    expected_naive = [
        full["cases"][country_index[row.country_iso3]][60 + row.horizon_month - 1]
        for row in seasonal.itertuples(index=False)
    ]
    assert seasonal["point_forecast"].tolist() == expected_naive

    assert len(metrics) == 12
    assert set(metrics["forecast_method"]) == methods
    assert set(metrics["country_iso3"]) == {"AUS", "CHN", "JPN", "ALL"}
    assert metrics["coverage_95"].between(0, 1).all()
    assert metrics[
        [
            "log_rmse",
            "mean_absolute_error",
            "median_95_interval_width",
            "mean_wis",
            "log_rmse_ratio_to_seasonal_naive",
            "mean_wis_ratio_to_seasonal_naive",
        ]
    ].notna().all().all()

    assert diagnostics["chains"] == 4
    assert diagnostics["iterations_per_chain"] == 2000
    assert diagnostics["post_warmup_draws"] == 4000
    assert diagnostics["divergent_transitions"] == 0
    assert diagnostics["maximum_treedepth_hits"] == 0
    assert diagnostics["maximum_rhat"] < 1.01
    assert diagnostics["minimum_neff"] >= 400

    model_predictive = metrics[
        metrics["forecast_method"].eq(
            "joint_model_negative_binomial"
        )
    ]
    assert model_predictive[
        "log_rmse_ratio_to_seasonal_naive"
    ].lt(1).all()
    assert model_predictive[
        "mean_wis_ratio_to_seasonal_naive"
    ].lt(1).all()

    postfit = audit["postfit_prediction_audit"]
    assert postfit["full_outcome_data_loaded_only_after_fit_or_cache_validation"]
    assert postfit["forecast_exposure_forced_zero"]
    assert set(postfit["forecast_exposure_weight_used"].values()) == {0}
    assert postfit["holdout_cases_used_only_for_scoring"]
    assert not postfit["holdout_genomes_used_for_fit"]
    assert not postfit["spline_basis_rebuilt_after_truncation"]
