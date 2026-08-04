from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/lineage_nomenclature"
DISPLAY_MAP = ROOT / "data/derived/model_lineage_display_map.tsv"


def test_mt28_nomenclature_crosswalk_passes() -> None:
    validation = json.loads(
        (RESULTS / "mt28_lineage_nomenclature_validation.json").read_text()
    )
    assert validation["status"] == "PASS"
    assert validation["analysis_identifier"] == "L1_02.07"
    assert (
        validation["recommended_display_name"]
        == "MT28-associated genomic lineage"
    )
    assert validation["n_tree_tips"] == 989
    assert validation["n_l10207_tree_tips"] == 288
    assert validation["n_l10207_focal_tips"] == 271
    assert validation["n_l10207_global_background_tips"] == 17


def test_published_mt28_labels_and_mrca_map_exactly_to_l10207() -> None:
    validation = json.loads(
        (RESULTS / "mt28_lineage_nomenclature_validation.json").read_text()
    )
    assert validation["n_published_mt28_labelled_tips"] == 99
    assert validation["n_published_mt28_labels_in_l10207"] == 99
    assert validation["n_published_mt28_labels_outside_l10207"] == 0
    assert validation["n_published_mt28_mrca_descendants"] == 288
    assert validation["n_non_l10207_inside_published_mt28_mrca"] == 0
    assert validation["n_l10207_outside_published_mt28_mrca"] == 0
    assert validation["l10207_monophyletic"] is True
    assert (
        validation["published_mt28_mrca_exactly_matches_l10207"] is True
    )

    tips = pd.read_csv(
        RESULTS / "mt28_lineage_tip_crosswalk.tsv", sep="\t"
    )
    labelled = tips[tips["has_published_mt28_label"]]
    assert len(labelled) == 99
    assert labelled["is_l10207"].all()
    assert labelled["inside_published_mt28_mrca"].all()


def test_mt28_alias_does_not_impute_missing_per_tip_typing() -> None:
    tips = pd.read_csv(
        RESULTS / "mt28_lineage_tip_crosswalk.tsv", sep="\t"
    )
    target = tips[tips["is_l10207"]]
    assert len(target) == 288
    assert target["published_mt28_label"].notna().sum() == 99
    assert target["published_mt28_label"].isna().sum() == 189


def test_reader_facing_lineage_map_is_complete_and_unique() -> None:
    display = pd.read_csv(DISPLAY_MAP, sep="\t")
    assert set(display["internal_id"]) == {
        "L1_01.02",
        "L1_02.05",
        "L1_02.06",
        "L1_02.07",
        "Other",
    }
    assert display["internal_id"].is_unique
    assert display["display_name"].is_unique
    assert display["short_label"].is_unique
    target = display.set_index("internal_id").loc["L1_02.07"]
    assert target["display_name"] == "MT28-associated genomic lineage"
    assert target["short_label"] == "MT28-associated"
