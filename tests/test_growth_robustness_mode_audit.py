from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/model_growth_robustness"
N_ANALYSES = 15
AUDITED_ANALYSIS_ID = "omit_country_JPN"


def _read_tsv(filename: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / filename, sep="\t")


def _assert_unique_analysis_ids(table: pd.DataFrame) -> set[str]:
    assert len(table) == N_ANALYSES
    assert table["analysis_id"].notna().all()
    assert table["analysis_id"].is_unique
    return set(table["analysis_id"])


def test_formal_growth_outputs_pass_and_include_the_audited_japan_omission() -> None:
    target = _read_tsv("l1_02_07_growth_robustness.tsv")
    diagnostics = _read_tsv("fit_diagnostics.tsv")

    target_ids = _assert_unique_analysis_ids(target)
    diagnostic_ids = _assert_unique_analysis_ids(diagnostics)
    assert target_ids == diagnostic_ids
    assert target["diagnostic_pass"].eq(True).all()  # noqa: E712
    assert diagnostics["diagnostic_pass"].eq(True).all()  # noqa: E712
    assert diagnostics["maximum_rhat"].le(1.01).all()
    assert diagnostics["divergent_transitions"].eq(0).all()
    assert diagnostics["maximum_treedepth_hits"].eq(0).all()

    audited_target = target.set_index("analysis_id").loc[AUDITED_ANALYSIS_ID]
    audited_diagnostics = diagnostics.set_index("analysis_id").loc[
        AUDITED_ANALYSIS_ID
    ]

    assert audited_target["median"] == pytest.approx(
        1.1052518502578, abs=1e-12
    )
    assert audited_target["lower_95"] == pytest.approx(
        1.06947234569106, abs=1e-12
    )
    assert audited_target["upper_95"] == pytest.approx(
        1.15385870616661, abs=1e-12
    )
    assert audited_target["diagnostic_pass"] == True  # noqa: E712
    assert audited_diagnostics["maximum_rhat"] <= 1.01
    assert audited_diagnostics["divergent_transitions"] == 0
    assert audited_diagnostics["maximum_treedepth_hits"] == 0
    assert audited_diagnostics["diagnostic_pass"] == True  # noqa: E712


def test_detailed_growth_outputs_are_complete_and_have_unique_composite_keys() -> None:
    target = _read_tsv("l1_02_07_growth_robustness.tsv")
    formal_ids = set(target["analysis_id"])
    lineage = _read_tsv("lineage_relative_growth_by_analysis.tsv")
    pairwise = _read_tsv("l1_02_07_pairwise_probabilities.tsv")

    assert len(lineage) == N_ANALYSES * 5
    assert set(lineage["analysis_id"]) == formal_ids
    assert lineage.groupby("analysis_id", observed=True).size().eq(5).all()
    assert not lineage.duplicated(["analysis_id", "lineage"]).any()

    assert len(pairwise) == N_ANALYSES * 4
    assert set(pairwise["analysis_id"]) == formal_ids
    assert pairwise.groupby("analysis_id", observed=True).size().eq(4).all()
    assert not pairwise.duplicated(
        ["analysis_id", "target_lineage", "comparator_lineage"]
    ).any()


def test_failed_legacy_japan_omission_is_preserved_outside_formal_tables() -> None:
    legacy_files = {
        "legacy_omit_country_JPN_failed_summary.tsv",
        "legacy_omit_country_JPN_failed_diagnostics.tsv",
        "legacy_omit_country_JPN_chain_medians.tsv",
        "legacy_omit_country_JPN_parameter_diagnostics.tsv",
    }
    assert all((RESULTS / filename).is_file() for filename in legacy_files)

    formal_target = _read_tsv("l1_02_07_growth_robustness.tsv").set_index(
        "analysis_id"
    )
    formal_diagnostics = _read_tsv("fit_diagnostics.tsv").set_index("analysis_id")
    legacy_target = _read_tsv(
        "legacy_omit_country_JPN_failed_summary.tsv"
    ).set_index("analysis_id")
    legacy_diagnostics = _read_tsv(
        "legacy_omit_country_JPN_failed_diagnostics.tsv"
    ).set_index("analysis_id")

    assert list(legacy_target.index) == [AUDITED_ANALYSIS_ID]
    assert list(legacy_diagnostics.index) == [AUDITED_ANALYSIS_ID]
    assert legacy_target.loc[AUDITED_ANALYSIS_ID, "diagnostic_pass"] == False  # noqa: E712
    assert legacy_diagnostics.loc[AUDITED_ANALYSIS_ID, "diagnostic_pass"] == False  # noqa: E712
    assert legacy_diagnostics.loc[AUDITED_ANALYSIS_ID, "maximum_rhat"] > 1.01

    assert formal_target.loc[AUDITED_ANALYSIS_ID, "diagnostic_pass"] == True  # noqa: E712
    assert formal_diagnostics.loc[AUDITED_ANALYSIS_ID, "diagnostic_pass"] == True  # noqa: E712
    assert formal_target.loc[AUDITED_ANALYSIS_ID, "median"] != pytest.approx(
        legacy_target.loc[AUDITED_ANALYSIS_ID, "median"], abs=1e-3
    )


def test_formal_run_configuration_identifies_all_three_posterior_components() -> None:
    configuration = _read_tsv("run_configuration.tsv")
    assert len(configuration) == 3
    assert set(configuration["run_component"]) == {
        "base_refits",
        "no_project_country_omissions",
        "omit_japan_mode_audit",
    }
    assert int(configuration["n_analyses"].sum()) == N_ANALYSES

    audit = configuration.set_index("run_component").loc[
        "omit_japan_mode_audit"
    ]
    assert audit["chains"] == 4
    assert audit["iterations_per_chain"] == 4000
    assert audit["warmup_per_chain"] == 2000
    assert audit["adapt_delta"] == pytest.approx(0.99)
    assert audit["max_treedepth"] == 15
    assert audit["init_mode"] == "mode_audit"
