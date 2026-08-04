from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/cohort_gate_sensitivity"
SCRIPT = ROOT / "scripts/pipeline/gtd_36_build_cohort_gate_sensitivity.py"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_frozen_pre_rescue_inputs_and_sha_manifest() -> None:
    manifest = read_tsv(RESULTS / "input_sha256_manifest.tsv")
    assert len(manifest) == 4
    assert {
        (row["cohort_version"], row["analysis_domain"]) for row in manifest
    } == {
        ("pre_rescue", "ancestry"),
        ("pre_rescue", "growth"),
        ("final", "ancestry"),
        ("final", "growth"),
    }

    for row in manifest:
        input_path = ROOT / row["input_path"]
        assert input_path.is_file()
        assert row["sha256"] == sha256_file(input_path)
        assert int(row["size_bytes"]) == input_path.stat().st_size
        assert int(row["data_rows"]) == 5

    expected_frozen_hashes = {
        "ancestry": (
            "5f280a1e57632a1eba7df8f810d0ff24"
            "9ea7a5f6599d3fa85ee58cccd18f2640"
        ),
        "growth": (
            "269eeafe041adee1f894da07d3411c73d"
            "2c34899ea75b5109a61ef25a50a9464"
        ),
    }
    pre_rescue = {
        row["analysis_domain"]: row["sha256"]
        for row in manifest
        if row["cohort_version"] == "pre_rescue"
    }
    assert pre_rescue == expected_frozen_hashes


def test_long_comparison_preserves_all_values_used_for_ranking() -> None:
    rows = read_tsv(RESULTS / "cohort_gate_sensitivity_long.tsv")
    assert len(rows) == 22
    assert {row["cohort_version"] for row in rows} == {
        "pre_rescue",
        "final",
    }
    assert {
        row["metric"] for row in rows
    } == {
        "relative_transmission_median",
        "mean_local_persistence_support",
        "mean_post_import_support",
    }
    for version in ("pre_rescue", "final"):
        version_rows = [
            row for row in rows if row["cohort_version"] == version
        ]
        assert len(version_rows) == 11
        for metric in (
            "mean_local_persistence_support",
            "mean_post_import_support",
        ):
            metric_rows = [
                row for row in version_rows if row["metric"] == metric
            ]
            assert len(metric_rows) == 5
            assert {int(row["rank_descending"]) for row in metric_rows} == {
                1,
                2,
                3,
                4,
                5,
            }

    keyed = {
        (row["cohort_version"], row["metric"], row["entity"]): row
        for row in rows
    }
    expected_growth = {
        "pre_rescue": (
            1.17515629339598,
            1.11043594801603,
            1.26846394470403,
        ),
        "final": (
            1.11204828279221,
            1.07829085892382,
            1.1555086421104499,
        ),
    }
    for version, expected in expected_growth.items():
        row = keyed[
            (version, "relative_transmission_median", "L1_02.07")
        ]
        assert (
            float(row["estimate"]),
            float(row["lower_95"]),
            float(row["upper_95"]),
        ) == expected
        assert row["prespecified_gate_member"] == "true"

    assert keyed[
        ("pre_rescue", "mean_local_persistence_support", "CHN")
    ]["rank_descending"] == "1"
    assert keyed[
        ("final", "mean_local_persistence_support", "CHN")
    ]["rank_descending"] == "1"
    assert {
        entity
        for (version, metric, entity), row in keyed.items()
        if version == "pre_rescue"
        and metric == "mean_post_import_support"
        and int(row["rank_descending"]) <= 2
    } == {"FRA", "JPN"}
    assert {
        entity
        for (version, metric, entity), row in keyed.items()
        if version == "final"
        and metric == "mean_post_import_support"
        and int(row["rank_descending"]) <= 2
    } == {"FRA", "JPN"}


def test_prespecified_cohort_gate_conclusions_pass_in_both_versions() -> None:
    validation = json.loads(
        (
            RESULTS / "cohort_gate_sensitivity_validation.json"
        ).read_text(encoding="utf-8")
    )
    assert validation["status"] == "PASS"
    assert validation["all_versions_all_prespecified_gates_pass"]

    expected = {
        "pre_rescue": {
            "lower": 1.11043594801603,
            "local_support": 0.998828373617346,
            "reseeding_order": ["JPN", "FRA"],
            "reseeding_support": {
                "JPN": 0.891006356848033,
                "FRA": 0.8760571161439243,
            },
        },
        "final": {
            "lower": 1.07829085892382,
            "local_support": 0.9986657453233722,
            "reseeding_order": ["FRA", "JPN"],
            "reseeding_support": {
                "FRA": 0.8721294777604437,
                "JPN": 0.8628041492705159,
            },
        },
    }
    for version, values in expected.items():
        version_result = validation["versions"][version]
        assert version_result["all_prespecified_gates_pass"]
        gates = version_result["gates"]
        assert all(gate["pass"] for gate in gates.values())

        growth_gate = gates["l10207_lower_95_above_1"]
        assert growth_gate["observed_lower_95"] == values["lower"]
        assert growth_gate["observed_lower_95"] > 1

        local_gate = gates["china_highest_local_persistence"]
        assert local_gate["observed_top_country"] == "CHN"
        assert local_gate["observed_top_support"] == values["local_support"]

        reseeding_gate = gates["france_japan_top2_reseeding"]
        assert reseeding_gate["observed_top2_ordered"] == values[
            "reseeding_order"
        ]
        assert reseeding_gate["observed_top2_set_sorted"] == ["FRA", "JPN"]
        assert reseeding_gate[
            "observed_support_by_top2_order"
        ] == values["reseeding_support"]


def test_pre_rescue_and_final_fits_use_strict_pre_model_prior_rule() -> None:
    audit_paths = {
        "pre_rescue": (
            ROOT
            / "data/model_sensitivity/cohort_gate_pre_rescue"
            / "joint_model_data_validation.json",
            ROOT
            / "data/model_sensitivity/cohort_gate_pre_rescue"
            / "initial_lineage_prior.tsv",
            ROOT
            / "data/model_sensitivity/cohort_gate_pre_rescue"
            / "sampling_diagnostics.json",
        ),
        "final": (
            ROOT / "data/model_inputs/joint_model_data_validation.json",
            ROOT / "data/model_inputs/initial_lineage_prior.tsv",
            ROOT / "results/model_main/sampling_diagnostics.json",
        ),
    }
    expected_historical_tips = {"pre_rescue": 99, "final": 164}
    for version, (
        validation_path,
        prior_path,
        diagnostics_path,
    ) in audit_paths.items():
        validation = json.loads(validation_path.read_text(encoding="utf-8"))
        assert validation["start_month"] == "2019-01-01"
        assert validation["initial_prior_mode"] == "historical"
        assert validation["initial_prior_end_year"] == 2018
        assert validation["initial_prior_strictly_precedes_model"]
        assert (
            validation["n_historical_tips_in_initial_prior"]
            == expected_historical_tips[version]
        )

        prior = read_tsv(prior_path)
        assert len(prior) == 15
        assert {row["initial_prior_mode"] for row in prior} == {
            "historical"
        }
        assert {int(row["initial_prior_end_year"]) for row in prior} == {2018}
        assert (
            sum(int(row["n_historical_tips"]) for row in prior)
            == expected_historical_tips[version]
        )

        diagnostics = json.loads(
            diagnostics_path.read_text(encoding="utf-8")
        )
        assert diagnostics["chains"] == 4
        assert diagnostics["divergent_transitions"] == 0
        assert diagnostics["maximum_treedepth_hits"] == 0


def test_builder_reproduces_committed_outputs(tmp_path: Path) -> None:
    regenerated = tmp_path / "cohort_gate_sensitivity"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--repo-root",
            str(ROOT),
            "--output-dir",
            str(regenerated),
        ],
        check=True,
    )
    for filename in (
        "cohort_gate_sensitivity_long.tsv",
        "input_sha256_manifest.tsv",
        "cohort_gate_sensitivity_validation.json",
    ):
        assert (regenerated / filename).read_bytes() == (
            RESULTS / filename
        ).read_bytes()
