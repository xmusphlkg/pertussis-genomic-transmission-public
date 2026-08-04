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
    / "scripts/pipeline/gtd_38_summarise_initial_state_prior_sensitivity.R"
)
LINEAGES = {"L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07", "Other"}


pytestmark = pytest.mark.skipif(
    shutil.which("Rscript") is None,
    reason="Rscript is required for the summary-script integration test",
)


def _write_posterior(path: Path, n_lineages: int = 5) -> None:
    expression = (
        "set.seed(38); "
        f"x <- matrix(exp(rnorm(40 * {n_lineages}, 0, 0.08)), "
        f"nrow=40, ncol={n_lineages}); "
        f"saveRDS(list(lineage_relative_transmission=x), {str(path)!r})"
    )
    subprocess.run(["Rscript", "-e", expression], check=True)


def _write_diagnostics(path: Path, minimum_field: str) -> None:
    payload = {
        "chains": 4,
        "iterations_per_chain": 20,
        "post_warmup_draws": 40,
        "seed": 38,
        "divergent_transitions": 0,
        "maximum_treedepth_hits": 0,
        "maximum_rhat": 1.001,
        minimum_field: 120.5,
    }
    path.write_text(json.dumps(payload))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _run_summary(tmp_path: Path) -> tuple[subprocess.CompletedProcess[str], Path]:
    strict_posterior = tmp_path / "strict.rds"
    strict_diagnostics = tmp_path / "strict.json"
    symmetric_posterior = tmp_path / "symmetric.rds"
    symmetric_diagnostics = tmp_path / "symmetric.json"
    outdir = tmp_path / "summary"
    _write_posterior(strict_posterior)
    _write_posterior(symmetric_posterior)
    _write_diagnostics(strict_diagnostics, "minimum_bulk_neff")
    _write_diagnostics(symmetric_diagnostics, "minimum_neff")

    completed = subprocess.run(
        [
            "Rscript",
            str(SCRIPT),
            str(strict_posterior),
            str(strict_diagnostics),
            str(symmetric_posterior),
            str(symmetric_diagnostics),
            str(outdir),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    return completed, outdir


def test_summary_outputs_are_complete_and_preserve_machine_lineage_ids(
    tmp_path: Path,
) -> None:
    completed, outdir = _run_summary(tmp_path)
    assert completed.returncode == 0, completed.stderr

    lineage = _read_tsv(
        outdir / "lineage_relative_transmission_by_initial_state_prior.tsv"
    )
    pairwise = _read_tsv(
        outdir / "l10207_pairwise_growth_by_initial_state_prior.tsv"
    )
    diagnostics = _read_tsv(
        outdir / "initial_state_prior_sampling_diagnostics.tsv"
    )
    manifest = _read_tsv(
        outdir / "initial_state_prior_sensitivity_manifest.tsv"
    )

    analysis_ids = {"historical_pre2019", "symmetric_dirichlet"}
    assert len(lineage) == 10
    assert {row["analysis_id"] for row in lineage} == analysis_ids
    assert {row["lineage"] for row in lineage} == LINEAGES
    assert len(pairwise) == 8
    assert {row["analysis_id"] for row in pairwise} == analysis_ids
    assert {row["numerator"] for row in pairwise} == {"L1_02.07"}
    assert {row["denominator"] for row in pairwise} == (
        LINEAGES - {"L1_02.07"}
    )
    assert len(diagnostics) == 2
    assert {
        row["minimum_effective_sample_size_field"] for row in diagnostics
    } == {"minimum_bulk_neff", "minimum_neff"}
    assert len(manifest) == 2
    assert {row["n_lineages"] for row in manifest} == {"5"}
    assert {row["n_posterior_draws"] for row in manifest} == {"40"}
    assert all(
        row["lineage_order"]
        == "L1_01.02;L1_02.05;L1_02.06;L1_02.07;Other"
        for row in manifest
    )


def test_summary_rejects_a_non_five_lineage_posterior(tmp_path: Path) -> None:
    completed, _ = _run_summary(tmp_path)
    assert completed.returncode == 0, completed.stderr

    invalid_posterior = tmp_path / "invalid.rds"
    _write_posterior(invalid_posterior, n_lineages=4)
    invalid = subprocess.run(
        [
            "Rscript",
            str(SCRIPT),
            str(invalid_posterior),
            str(tmp_path / "strict.json"),
            str(tmp_path / "symmetric.rds"),
            str(tmp_path / "symmetric.json"),
            str(tmp_path / "invalid-summary"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    assert invalid.returncode != 0
    assert "must have exactly 5 columns" in invalid.stderr
