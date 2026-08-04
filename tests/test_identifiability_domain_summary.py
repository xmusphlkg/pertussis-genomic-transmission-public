from __future__ import annotations

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/pipeline/gtd_39_summarise_identifiability_domains.R"
)
RECOVERY = ROOT / "results/model_main/all_recovery.tsv"
DIAGNOSTICS = ROOT / "results/model_main/all_diagnostics.tsv"
GATE = ROOT / "results/model_main/identifiability_gate.json"


pytestmark = pytest.mark.skipif(
    shutil.which("Rscript") is None,
    reason="Rscript is required for the domain-summary integration tests",
)


def _run_summary(
    outdir: Path,
    gate: Path = GATE,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            "Rscript",
            str(SCRIPT),
            str(RECOVERY),
            str(DIAGNOSTICS),
            str(gate),
            str(outdir),
        ],
        text=True,
        capture_output=True,
        check=False,
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_real_recovery_outputs_are_summarised_without_replacing_joint_gate(
    tmp_path: Path,
) -> None:
    outdir = tmp_path / "domain-summary"
    completed = _run_summary(outdir)
    assert completed.returncode == 0, completed.stderr

    rows = {
        row["estimand_group"]: row
        for row in _read_tsv(outdir / "estimand_recovery_summary.tsv")
    }
    report = json.loads(
        (outdir / "identifiability_domain_interpretation.json").read_text()
    )
    source_gate = json.loads(GATE.read_text())
    source_diagnostics = _read_tsv(DIAGNOSTICS)

    assert set(rows) == {
        "all_lineage_growth",
        "target_lineage_growth",
        "import_scale",
    }
    assert report["joint_all_parameters_pass"] is source_gate["pass"]
    assert report["joint_all_parameters_pass"] is False
    assert report["lineage_growth_domain_pass"] is True
    assert report["import_scale_domain_pass"] is False
    assert report["provenance"]["joint_gate_value_copied_not_recomputed"] is True
    assert report["provenance"]["domain_pass_does_not_replace_joint_gate"] is True
    assert report["provenance"]["new_inferential_thresholds_added"] is False

    lineage = rows["all_lineage_growth"]
    assert lineage["n_replicates"] == "6"
    assert lineage["n_parameters"] == "5"
    assert lineage["n_parameter_replicate_pairs"] == "30"
    assert lineage["domain_pass"] == "TRUE"
    assert float(lineage["minimum_95_interval_coverage_threshold"]) == (
        source_gate["thresholds_defined_before_recovery_results"][
            "minimum_95_interval_coverage"
        ]
    )
    assert float(lineage["maximum_median_absolute_log_error_threshold"]) == (
        source_gate["thresholds_defined_before_recovery_results"][
            "maximum_median_absolute_log_error_lineage_growth"
        ]
    )
    assert float(lineage["highest_lineage_rank_recovery"]) == (
        source_gate["observed"]["highest_lineage_rank_recovery"]
    )

    target = rows["target_lineage_growth"]
    assert target["analysis_identifier"] == "L1_02.07"
    assert target["n_replicates"] == "6"
    assert target["n_parameters"] == "1"
    assert target["n_parameter_replicate_pairs"] == "6"
    assert target["recovery_status"] == "descriptive_only_n_6_replicates"
    assert target["domain_pass"] == ""
    assert target["minimum_95_interval_coverage_threshold"] == ""
    assert (
        report["domains"]["target_lineage_growth"]["threshold_evaluated"]
        is False
    )

    import_scale = rows["import_scale"]
    assert import_scale["n_replicates"] == "6"
    assert import_scale["n_parameters"] == "3"
    assert import_scale["n_parameter_replicate_pairs"] == "18"
    assert import_scale["domain_pass"] == "FALSE"
    assert float(
        import_scale["maximum_median_absolute_log_error_threshold"]
    ) == source_gate["thresholds_defined_before_recovery_results"][
        "maximum_median_absolute_log_error_import_scale"
    ]

    assert int(lineage["divergent_transitions"]) == sum(
        int(row["divergent_transitions"]) for row in source_diagnostics
    )
    assert int(lineage["maximum_treedepth_hits"]) == sum(
        int(row["maximum_treedepth_hits"]) for row in source_diagnostics
    )
    assert float(lineage["maximum_rhat"]) == max(
        float(row["maximum_rhat"]) for row in source_diagnostics
    )
    assert float(lineage["minimum_neff"]) == min(
        float(row["minimum_neff"]) for row in source_diagnostics
    )
    assert "structurally less recoverable" in (
        report["interpretation"]["import_scale"]
    )
    assert "must remain conditional" in (
        report["interpretation"]["import_scale"]
    )


def test_summary_rejects_gate_from_a_different_recovery_run(
    tmp_path: Path,
) -> None:
    mismatched_gate = json.loads(GATE.read_text())
    mismatched_gate["observed"]["lineage_growth_coverage"] += 0.1
    mismatched_gate_path = tmp_path / "mismatched-gate.json"
    mismatched_gate_path.write_text(json.dumps(mismatched_gate))

    completed = _run_summary(
        tmp_path / "should-not-be-written",
        gate=mismatched_gate_path,
    )
    assert completed.returncode != 0
    assert "does not match existing gate field" in completed.stderr
