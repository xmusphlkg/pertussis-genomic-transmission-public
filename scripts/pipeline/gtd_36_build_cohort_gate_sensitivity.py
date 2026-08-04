#!/usr/bin/env python3
"""Compare prespecified ancestry and growth conclusions before and after cohort rescue."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any


EXPECTED_COUNTRIES = {"AUS", "BEL", "CHN", "FRA", "JPN"}
EXPECTED_RESEEDING_TOP2 = {"FRA", "JPN"}
TARGET_LINEAGE = "L1_02.07"


def parse_args() -> argparse.Namespace:
    root_default = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(
        description=(
            "Audit whether prespecified ancestry and L1_02.07 growth "
            "conclusions survive the cohort gate-rescue step."
        )
    )
    parser.add_argument("--repo-root", type=Path, default=root_default)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader)


def require_fields(
    rows: list[dict[str, str]],
    required: set[str],
    path: Path,
) -> None:
    if not rows:
        raise ValueError(f"Input has no data rows: {path}")
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(
            f"Input {path} is missing required fields: {sorted(missing)}"
        )


def parse_float(value: str, field: str, path: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"Non-numeric value for {field} in {path}: {value!r}"
        ) from error
    if result != result or result in (float("inf"), float("-inf")):
        raise ValueError(f"Non-finite value for {field} in {path}: {value!r}")
    return result


def load_ancestry(path: Path) -> dict[str, dict[str, float]]:
    rows = read_tsv(path)
    required = {
        "country_iso3",
        "mean_post_import_support",
        "mean_local_persistence_support",
    }
    require_fields(rows, required, path)
    countries = [row["country_iso3"] for row in rows]
    if len(countries) != len(set(countries)):
        raise ValueError(f"Duplicate countries in {path}")
    if set(countries) != EXPECTED_COUNTRIES:
        raise ValueError(
            f"Unexpected country set in {path}: {sorted(set(countries))}"
        )

    parsed: dict[str, dict[str, float]] = {}
    for row in rows:
        country = row["country_iso3"]
        parsed[country] = {}
        for field in (
            "mean_post_import_support",
            "mean_local_persistence_support",
        ):
            value = parse_float(row[field], field, path)
            if not 0.0 <= value <= 1.0:
                raise ValueError(
                    f"{field} must be in [0, 1] in {path}; observed {value}"
                )
            parsed[country][field] = value
    return parsed


def load_growth(path: Path) -> dict[str, float]:
    rows = read_tsv(path)
    required = {"lineage", "mean", "median", "lower_95", "upper_95"}
    require_fields(rows, required, path)
    lineages = [row["lineage"] for row in rows]
    if len(lineages) != len(set(lineages)):
        raise ValueError(f"Duplicate lineages in {path}")
    target = [row for row in rows if row["lineage"] == TARGET_LINEAGE]
    if len(target) != 1:
        raise ValueError(
            f"Expected exactly one {TARGET_LINEAGE} row in {path}"
        )
    parsed = {
        field: parse_float(target[0][field], field, path)
        for field in ("mean", "median", "lower_95", "upper_95")
    }
    if not (
        parsed["lower_95"]
        <= parsed["median"]
        <= parsed["upper_95"]
    ):
        raise ValueError(f"Invalid uncertainty interval ordering in {path}")
    return parsed


def descending_ranks(values: dict[str, float]) -> dict[str, int]:
    ordered = sorted(values, key=lambda item: (-values[item], item))
    return {item: index for index, item in enumerate(ordered, start=1)}


def relative_path(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def write_tsv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, Any]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir
        else repo_root / "results" / "cohort_gate_sensitivity"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs = {
        ("pre_rescue", "ancestry"): (
            repo_root
            / "data/model_sensitivity/cohort_gate_pre_rescue/"
            "country_phylogeographic_summary.tsv"
        ),
        ("pre_rescue", "growth"): (
            repo_root
            / "data/model_sensitivity/cohort_gate_pre_rescue/"
            "lineage_relative_transmission.tsv"
        ),
        ("final", "ancestry"): (
            repo_root
            / "results/phylogeography/events_thr0_5/"
            "country_phylogeographic_summary.tsv"
        ),
        ("final", "growth"): (
            repo_root
            / "results/model_main/lineage_relative_transmission.tsv"
        ),
    }
    missing = [path for path in inputs.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing required inputs:\n" + "\n".join(map(str, missing))
        )

    ancestry = {
        version: load_ancestry(inputs[(version, "ancestry")])
        for version in ("pre_rescue", "final")
    }
    growth = {
        version: load_growth(inputs[(version, "growth")])
        for version in ("pre_rescue", "final")
    }

    comparison_rows: list[dict[str, Any]] = []
    validation_versions: dict[str, Any] = {}
    for version in ("pre_rescue", "final"):
        comparison_rows.append(
            {
                "cohort_version": version,
                "analysis_domain": "growth",
                "metric": "relative_transmission_median",
                "entity": TARGET_LINEAGE,
                "estimate": format(growth[version]["median"], ".17g"),
                "lower_95": format(growth[version]["lower_95"], ".17g"),
                "upper_95": format(growth[version]["upper_95"], ".17g"),
                "rank_descending": "",
                "prespecified_gate_member": "true",
            }
        )

        local_values = {
            country: values["mean_local_persistence_support"]
            for country, values in ancestry[version].items()
        }
        reseeding_values = {
            country: values["mean_post_import_support"]
            for country, values in ancestry[version].items()
        }
        local_ranks = descending_ranks(local_values)
        reseeding_ranks = descending_ranks(reseeding_values)
        for metric, values, ranks, gate_members in (
            (
                "mean_local_persistence_support",
                local_values,
                local_ranks,
                {"CHN"},
            ),
            (
                "mean_post_import_support",
                reseeding_values,
                reseeding_ranks,
                EXPECTED_RESEEDING_TOP2,
            ),
        ):
            for country in sorted(values):
                comparison_rows.append(
                    {
                        "cohort_version": version,
                        "analysis_domain": "ancestry",
                        "metric": metric,
                        "entity": country,
                        "estimate": format(values[country], ".17g"),
                        "lower_95": "",
                        "upper_95": "",
                        "rank_descending": ranks[country],
                        "prespecified_gate_member": str(
                            country in gate_members
                        ).lower(),
                    }
                )

        observed_local_top = min(
            local_ranks, key=lambda country: local_ranks[country]
        )
        ordered_reseeding = sorted(
            reseeding_ranks, key=lambda country: reseeding_ranks[country]
        )
        observed_reseeding_top2 = ordered_reseeding[:2]
        version_gates = {
            "l10207_lower_95_above_1": {
                "threshold": 1.0,
                "observed_lower_95": growth[version]["lower_95"],
                "pass": growth[version]["lower_95"] > 1.0,
            },
            "china_highest_local_persistence": {
                "expected_top_country": "CHN",
                "observed_top_country": observed_local_top,
                "observed_top_support": local_values[observed_local_top],
                "pass": observed_local_top == "CHN",
            },
            "france_japan_top2_reseeding": {
                "expected_top2_set_sorted": sorted(
                    EXPECTED_RESEEDING_TOP2
                ),
                "observed_top2_ordered": observed_reseeding_top2,
                "observed_top2_set_sorted": sorted(
                    observed_reseeding_top2
                ),
                "observed_support_by_top2_order": {
                    country: reseeding_values[country]
                    for country in observed_reseeding_top2
                },
                "pass": set(observed_reseeding_top2)
                == EXPECTED_RESEEDING_TOP2,
            },
        }
        validation_versions[version] = {
            "gates": version_gates,
            "all_prespecified_gates_pass": all(
                gate["pass"] for gate in version_gates.values()
            ),
        }

    comparison_path = output_dir / "cohort_gate_sensitivity_long.tsv"
    comparison_fields = [
        "cohort_version",
        "analysis_domain",
        "metric",
        "entity",
        "estimate",
        "lower_95",
        "upper_95",
        "rank_descending",
        "prespecified_gate_member",
    ]
    write_tsv(comparison_path, comparison_fields, comparison_rows)

    manifest_rows = []
    for version in ("pre_rescue", "final"):
        for domain in ("ancestry", "growth"):
            path = inputs[(version, domain)]
            manifest_rows.append(
                {
                    "cohort_version": version,
                    "analysis_domain": domain,
                    "input_path": relative_path(path, repo_root),
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                    "data_rows": len(read_tsv(path)),
                }
            )
    write_tsv(
        output_dir / "input_sha256_manifest.tsv",
        [
            "cohort_version",
            "analysis_domain",
            "input_path",
            "sha256",
            "size_bytes",
            "data_rows",
        ],
        manifest_rows,
    )

    all_pass = all(
        result["all_prespecified_gates_pass"]
        for result in validation_versions.values()
    )
    validation = {
        "analysis_id": "cohort_gate_pre_rescue_vs_final",
        "comparison_table": "cohort_gate_sensitivity_long.tsv",
        "input_manifest": "input_sha256_manifest.tsv",
        "prespecified_gate_definition": {
            "growth": f"{TARGET_LINEAGE} lower_95 > 1 in each cohort",
            "local_persistence": (
                "CHN has rank 1 mean_local_persistence_support "
                "in each cohort"
            ),
            "reseeding": (
                "rank 1-2 mean_post_import_support country set is "
                "{FRA, JPN} in each cohort"
            ),
        },
        "versions": validation_versions,
        "all_versions_all_prespecified_gates_pass": all_pass,
        "status": "PASS" if all_pass else "FAIL",
    }
    (output_dir / "cohort_gate_sensitivity_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not all_pass:
        raise RuntimeError("One or more prespecified cohort-gate checks failed")


if __name__ == "__main__":
    main()
