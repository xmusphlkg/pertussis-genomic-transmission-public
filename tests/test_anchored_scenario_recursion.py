from __future__ import annotations

import subprocess
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_anchored_scenario_recursion_contract() -> None:
    subprocess.run(
        ["Rscript", "tests/test_anchored_scenario_recursion.R"],
        cwd=ROOT,
        check=True,
    )


def test_anchor_sensitivity_output_contract() -> None:
    output = pd.read_csv(
        ROOT
        / "results/model_anchor_sensitivity/anchor_scenario_sensitivity.tsv",
        sep="\t",
    )
    assert len(output) == 48
    assert set(output["model_specification"]) == {"main", "no_project"}
    assert set(output["anchor_month"]) == {
        "2022-06",
        "2022-12",
        "2023-06",
        "2023-12",
    }
    assert set(output["scenario"]) == {
        "no_post_2022_exposure",
        "l10207_neutral_growth",
    }
    assert set(output["country_iso3"]) == {"AUS", "CHN", "JPN"}
    assert output["n_posterior_draws"].eq(4000).all()
    assert (
        output["first_month_baseline_max_abs_error_from_stan_mu"] < 1e-8
    ).all()

    china_no_exposure = output[
        output["country_iso3"].eq("CHN")
        & output["scenario"].eq("no_post_2022_exposure")
    ]
    assert china_no_exposure["scenario_equals_baseline_all_draws"].all()
    for column in (
        "cumulative_case_difference_mean",
        "cumulative_case_difference_median",
        "cumulative_case_difference_lower_95",
        "cumulative_case_difference_upper_95",
        "cumulative_difference_fraction_mean",
        "cumulative_difference_fraction_median",
        "cumulative_difference_fraction_lower_95",
        "cumulative_difference_fraction_upper_95",
    ):
        assert china_no_exposure[column].eq(0).all()
