from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/molecular_characterisation"


def test_molecular_annotations_reconnect_without_identifier_conflicts() -> None:
    validation = json.loads(
        (RESULTS / "molecular_characterisation_validation.json").read_text()
    )
    assert validation["status"] == "PASS"
    assert validation["n_focal_tree_tips"] == 774
    assert validation["n_focal_l10207"] == 271
    assert validation["n_focal_source_record_matches"] == 774
    assert validation["n_focal_frozen_archive_matches"] == 495
    assert validation["n_identifier_conflicts"] == 0


def test_l10207_typed_antigen_profile_and_molecular_heterogeneity() -> None:
    findings = pd.read_csv(RESULTS / "key_findings.tsv", sep="\t").set_index(
        "finding_id"
    )
    assert findings.loc["l10207_ptxP3_among_typed", "estimate"] == "93/93"
    assert (
        findings.loc[
            "l10207_ptxP3_fim3_1_among_jointly_typed", "estimate"
        ]
        == "83/83"
    )
    assert (
        findings.loc["l10207_frozen_23S_MR_A2047G", "estimate"]
        == "33 MR_A2047G; 47 MS"
    )
    assert (
        findings.loc["l10207_frozen_PRN_disrupted", "estimate"]
        == "17 disrupted; 61 intact"
    )
    assert (
        findings.loc["l10207_frozen_PRN_disruption_mechanisms", "estimate"]
        == "11 IS481; 1 inversion/rearrangement; 5 other"
    )


def test_l10207_period_shift_persists_in_observed_resistance_and_prn_strata() -> None:
    shifts = pd.read_csv(
        RESULTS / "l10207_period_shift_within_molecular_strata.tsv",
        sep="\t",
    ).set_index(["country_iso3", "feature", "molecular_stratum"])

    expected = {
        ("CHN", "frozen_23S_status", "MR_A2047G"): (1, 79, 27, 1),
        ("CHN", "frozen_23S_status", "MS"): (2, 16, 18, 0),
        ("JPN", "frozen_23S_status", "MS"): (0, 19, 7, 1),
        ("CHN", "frozen_PRN_outcome", "intact"): (4, 36, 27, 3),
        ("JPN", "frozen_PRN_outcome", "intact"): (0, 18, 14, 2),
        ("CHN", "source_antigen_profile", "ptxP3/fim3-1"): (2, 10, 45, 0),
    }
    for key, counts in expected.items():
        row = shifts.loc[key]
        observed = tuple(
            int(row[column])
            for column in (
                "prepandemic_l10207",
                "prepandemic_other",
                "resurgence_l10207",
                "resurgence_other",
            )
        )
        assert observed == counts
        assert row["analysis_status"] == "fisher_exact_prespecified"
        assert row["fisher_p_bh_exploratory"] < 0.001


def test_recent_china_japan_molecular_missingness_is_explicit() -> None:
    coverage = pd.read_csv(
        RESULTS / "molecular_annotation_coverage.tsv", sep="\t"
    )
    china = coverage[
        coverage["country_iso3"].eq("CHN")
        & coverage["epidemic_period"].eq("resurgence")
    ]
    japan = coverage[
        coverage["country_iso3"].eq("JPN")
        & coverage["epidemic_period"].eq("resurgence")
    ]
    assert china["n_tree_tips"].sum() == 153
    assert china["n_frozen_23S_annotated"].sum() == 46
    assert china["n_source_antigen_profile_annotated"].sum() == 46
    assert japan["n_tree_tips"].sum() == 60
    assert japan["n_frozen_23S_annotated"].sum() == 8
    assert japan["n_source_antigen_profile_annotated"].sum() == 0
