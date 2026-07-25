from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/pipeline/gtd_01_build_transmission_model_tables.py"
SPEC = importlib.util.spec_from_file_location("gtd_model_tables", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_date_intervals_preserve_uncertainty() -> None:
    day = MODULE.parse_date_interval("2024-08-17")
    month = MODULE.parse_date_interval("2024-08")
    year = MODULE.parse_date_interval("2024")

    assert day.lower == day.upper == "2024-08-17"
    assert day.resolution == "day"
    assert month.lower == "2024-08-01"
    assert month.upper == "2024-08-31"
    assert month.month_or_quarter
    assert year.lower == "2024-01-01"
    assert year.upper == "2024-12-31"
    assert not year.month_or_quarter


def test_marker_lineage_normalization_is_consistent() -> None:
    assert MODULE.marker_lineage("ptxP_3", "fim3_1") == "ptxP3/fim3-1"
    assert MODULE.marker_lineage("ptxP(3)", "fim3(1)") == "ptxP3/fim3-1"


def test_stage0_gate_outputs_keep_clock_and_lineages_locked() -> None:
    gate_path = ROOT / "data/derived/country_data_gate_summary.tsv"
    if not gate_path.exists():
        return
    gates = pd.read_csv(gate_path, sep="\t")
    assert set(gates["country_iso3"]) == {"CHN", "AUS", "JPN"}
    assert gates["temporal_metadata_gate_pass"].all()
    assert gates["phylogenetic_clock_gate_status"].str.startswith("PENDING_").all()
    assert gates["formal_lineage_gate_status"].str.startswith("PENDING_").all()
