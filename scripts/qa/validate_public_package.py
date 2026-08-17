#!/usr/bin/env python3
"""Validate the Letter-first public reproducibility package.

The checks intentionally cover the reader-facing information architecture as
well as the small set of numerical claims surfaced on the repository landing
page.  The script uses only the Python standard library.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[2]
FILE_MANIFEST = ROOT / "provenance/FILE_MANIFEST_SHA256.tsv"
IGNORED_PARTS = {
    ".git",
    ".pytest_cache",
    ".venv",
    ".venv-tests",
    "__pycache__",
}
TEXT_SUFFIXES = {
    "",
    ".cff",
    ".csv",
    ".iqtree",
    ".json",
    ".md",
    ".py",
    ".r",
    ".sh",
    ".stan",
    ".tsv",
    ".txt",
}
PRIVATE_PATH_MARKERS = (
    "/mnt/nas/",
    "/home/ctm/",
    "pertussis_genomic_transmission_joi",
)
FONG_SOURCE_DATASET = "australia_fong_2026_direct_specimens"
RESTRICTED_FONG_COLUMNS = {
    "sample_id",
    "genome_qc_status",
    "sampling_process_observed",
    "sequencing_success",
    "ct_value",
    "specimen_type",
    "preliminary_lineage_id",
    "lineage_definition_status",
    "lineage_stratum",
    "published_branch",
    "published_lineage",
    "published_sublineage",
    "ptxP_label",
    "fim3_label",
    "marker_23s_status",
}
FONG_SANITISED_TABLES = (
    "data/derived/transmission_genome_records.tsv",
    "results/phylogeny/focal_phylogeny_selection.tsv",
    "results/phylogeny/primary_phylogeny_manifest.tsv",
    "results/phylogeny/tree_tip_metadata.tsv",
    "results/phylogeny/uniform_sequence_qc.tsv",
    "figures/source_data/figure2a_tree_tip_metadata.tsv",
)


class Validation:
    def __init__(self) -> None:
        self.failures: list[str] = []

    def require(self, condition: bool, message: str) -> None:
        if not condition:
            self.failures.append(message)

    def finish(self) -> None:
        if self.failures:
            print(f"FAIL: {len(self.failures)} public-package check(s) failed")
            for failure in self.failures:
                print(f"  - {failure}")
            raise SystemExit(1)
        print("PASS: public package is structurally and numerically valid")


def read_tsv(relative_path: str) -> list[dict[str, str]]:
    with (ROOT / relative_path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def released_files() -> list[Path]:
    paths = []
    for path in sorted(ROOT.rglob("*")):
        if (
            not path.is_file()
            or IGNORED_PARTS.intersection(path.parts)
            or path.suffix in {".pyc", ".pyo"}
            or path == FILE_MANIFEST
        ):
            continue
        paths.append(path)
    return paths


def validate_structure(check: Validation) -> set[str]:
    required = {
        "README.md",
        "DATA_DICTIONARY.md",
        "LETTER_EVIDENCE_MAP.md",
        "CHANGELOG.md",
        "CITATION.cff",
        "figures/letter/Figure_1_observation_structure.png",
        "figures/letter/Figure_1_observation_structure.pdf",
        "figures/letter/Figure_2_growth_robustness.png",
        "figures/letter/Figure_2_growth_robustness.pdf",
        "figures/letter/Supplementary_Figure_S1_genomic_validation.png",
        "figures/letter/Supplementary_Figure_S2_posterior_predictive_checks.png",
        "figures/letter/Supplementary_Figure_S3_ancestry_input_sensitivity.png",
        "figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv",
        "figures/letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv",
        "figures/letter/RENDER_MANIFEST.tsv",
        "figures/letter/SUPPLEMENTARY_RENDER_MANIFEST.tsv",
        "scripts/figures/render_letter_figures.R",
        "scripts/figures/render_letter_supplementary_figures.R",
        "EID_README.md",
        "EID_DATA_DICTIONARY.md",
        "data/derived/public_genome_availability.tsv",
        "data/derived/public_availability_ena_runs.tsv",
        "data/derived/public_availability_ena_project_audit.tsv",
        "results/public_availability/eid_detection_clock_shift.tsv",
        "results/public_availability/eid_milestone_visibility.tsv",
        "results/public_availability/eid_threshold_sensitivity.tsv",
        "results/public_availability/eid_case_clock_sensitivity.tsv",
        "results/public_availability/eid_project_batch_release.tsv",
        "results/public_availability/eid_geography_audit.tsv",
        "results/public_availability/eid_external_candidate_summary.tsv",
        "figures/eid/Figure_1_release_clock_pertussis_eid.png",
        "figures/eid/Figure_1_release_clock_pertussis_eid.pdf",
        "figures/eid/Figure_1_release_clock_pertussis_eid.svg",
        "figures/eid/Figure_1_release_clock_pertussis_eid.tiff",
        "manuscript/eid_dispatch_appendix.md",
        "provenance/EID_NAS_SNAPSHOT_MANIFEST.tsv",
        "scripts/pipeline/gtd_40_build_public_availability.py",
        "scripts/pipeline/gtd_43_build_eid_dispatch_tables.py",
        "scripts/figures/render_eid_dispatch_release_clock_figure.R",
    }
    for relative_path in sorted(required):
        check.require((ROOT / relative_path).is_file(), f"missing {relative_path}")

    check.require(not (ROOT / "figures/main").exists(), "obsolete figures/main remains")
    check.require(
        not (ROOT / "figures/supplementary").exists(),
        "obsolete figures/supplementary remains",
    )
    check.require(
        not (ROOT / "data/derived/australia_sampling_process_records.tsv").exists(),
        "non-redistributable Australian specimen-level table remains",
    )
    check.require(
        not list((ROOT / "figures/letter").glob("*.svg")),
        "unreleased SVG figure export remains in figures/letter",
    )
    check.require(
        not list((ROOT / "figures/letter").glob("*.tiff")),
        "unreleased TIFF figure export remains in figures/letter",
    )

    mapped_sources: set[str] = set()
    for manifest in (
        "figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv",
        "figures/letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv",
    ):
        if not (ROOT / manifest).is_file():
            continue
        rows = read_tsv(manifest)
        check.require(bool(rows), f"empty mapping manifest: {manifest}")
        for row in rows:
            source_file = row.get("source_file", "")
            mapped_sources.add(source_file)
            check.require(
                bool(source_file) and (ROOT / source_file).is_file(),
                f"mapped source does not exist: {source_file or '<blank>'}",
            )

    expected_panel_sources = {
        "figures/source_data/figure1a_cases.tsv",
        "figures/source_data/figure1b_annual_genomes.tsv",
        "figures/source_data/figure2_tree_manifest.tsv",
        "figures/source_data/figure2a_tree_tip_metadata.tsv",
        "figures/source_data/figure2b_tip_ancestry_support.tsv",
        "figures/source_data/figure3a_lineage_growth_main.tsv",
        "figures/source_data/figure3b_l10207_pairwise_growth.tsv",
        "figures/source_data/figure3c_l10207_growth_robustness.tsv",
        "figures/source_data/figure3d_selection_cap_weighted_l10207_shares.tsv",
        "figures/source_data/figure4abc_monthly_counterfactuals.tsv",
        "figures/source_data/figure4d_counterfactual_summary.tsv",
        "figures/source_data/figure4e_australia_ct_curve.tsv",
        "figures/source_data/figure4f_identifiability_recovery.tsv",
        "figures/source_data/eid_figure1a_cases.tsv",
        "figures/source_data/eid_figure1a_selected_detection.tsv",
        "figures/source_data/eid_figure1b_release_lags.tsv",
        "figures/source_data/eid_figure1c_clock_shift.tsv",
    }
    observed_panel_sources = {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "figures/source_data").glob("*.tsv")
    }
    check.require(
        observed_panel_sources == expected_panel_sources,
        "figures/source_data does not exactly match the released panel set",
    )
    check.require(
        expected_panel_sources.difference(
            {
                "figures/source_data/eid_figure1a_cases.tsv",
                "figures/source_data/eid_figure1a_selected_detection.tsv",
                "figures/source_data/eid_figure1b_release_lags.tsv",
                "figures/source_data/eid_figure1c_clock_shift.tsv",
            }
        ).issubset(mapped_sources),
        "one or more Letter panel sources are absent from the mapping manifests",
    )
    return mapped_sources


def validate_links(check: Validation) -> None:
    link_pattern = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
    for document in (
        "README.md",
        "EID_README.md",
        "EID_DATA_DICTIONARY.md",
        "DATA_DICTIONARY.md",
        "LETTER_EVIDENCE_MAP.md",
        "docs/EID_DISPATCH_CLAIM_EVIDENCE_MAP.md",
        "docs/EID_FIGURE_CONTRACT.md",
        "docs/PRJNA1071282_BOUNDARY_AUDIT.md",
        "figures/source_data/README.md",
    ):
        path = ROOT / document
        if not path.is_file():
            continue
        for raw_target in link_pattern.findall(path.read_text(encoding="utf-8")):
            target = raw_target.strip().split()[0].strip("<>")
            if target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target = unquote(target.split("#", 1)[0])
            resolved = (path.parent / target).resolve()
            check.require(
                resolved == ROOT or ROOT in resolved.parents,
                f"link escapes repository in {document}: {raw_target}",
            )
            check.require(resolved.exists(), f"broken link in {document}: {raw_target}")


def validate_citation(check: Validation) -> None:
    path = ROOT / "CITATION.cff"
    if not path.is_file():
        return
    text = path.read_text(encoding="utf-8")
    check.require("version: 2.0.0" in text, "CITATION.cff version is not 2.0.0")
    check.require(
        "date-released: 2026-08-04" in text,
        "CITATION.cff release date is not 2026-08-04",
    )
    check.require(
        len(re.findall(r"^  - family-names:", text, flags=re.MULTILINE)) == 15,
        "CITATION.cff does not contain the 15 Letter authors",
    )


def validate_core_numbers(check: Validation) -> None:
    validation_path = (
        ROOT
        / "results/lineage_nomenclature/mt28_lineage_nomenclature_validation.json"
    )
    if validation_path.is_file():
        payload = json.loads(validation_path.read_text(encoding="utf-8"))
        check.require(payload.get("status") == "PASS", "MT28 crosswalk does not pass")
        check.require(payload.get("n_tree_tips") == 989, "tree-tip count is not 989")
        check.require(
            payload.get("n_l10207_tree_tips") == 288,
            "MT28-associated population does not contain 288 tips",
        )
        check.require(
            payload.get("n_published_mt28_labelled_tips") == 99
            and payload.get("n_published_mt28_labels_in_l10207") == 99,
            "published MT28 label crosswalk is not 99 of 99",
        )
        check.require(
            payload.get("published_mt28_mrca_exactly_matches_l10207") is True,
            "published MT28 MRCA does not exactly match L1_02.07",
        )
    else:
        check.require(False, "missing MT28 nomenclature validation JSON")

    growth_path = "results/model_main/lineage_relative_transmission.tsv"
    if (ROOT / growth_path).is_file():
        rows = [row for row in read_tsv(growth_path) if row["lineage"] == "L1_02.07"]
        check.require(len(rows) == 1, "primary growth table lacks one L1_02.07 row")
        if len(rows) == 1:
            check.require(
                math.isclose(as_float(rows[0], "median"), 1.11204828279221, rel_tol=1e-10),
                "primary L1_02.07 median changed from the frozen result",
            )
            check.require(
                as_float(rows[0], "lower_95") > 1,
                "primary L1_02.07 credible interval is not above reference",
            )

    robustness_path = "results/model_growth_robustness/l1_02_07_growth_robustness.tsv"
    if (ROOT / robustness_path).is_file():
        rows = read_tsv(robustness_path)
        check.require(all(row["diagnostic_pass"] == "TRUE" for row in rows), "a reported robustness fit fails diagnostics")
        by_id = {row["analysis_id"]: row for row in rows}
        for analysis_id in ("country_CHN_project_adjusted", "country_JPN_project_adjusted"):
            check.require(
                analysis_id in by_id and as_float(by_id[analysis_id], "lower_95") > 1,
                f"{analysis_id} is not completely above reference",
            )
        aus = by_id.get("country_AUS_project_adjusted")
        check.require(
            aus is not None
            and as_float(aus, "lower_95") < 1 < as_float(aus, "upper_95"),
            "Australian country-only interval is not uninformative around reference",
        )
        omission_rows = [
            row
            for row in rows
            if row["analysis_type"] in {"country_omission", "dominant_project_omission"}
        ]
        check.require(
            bool(omission_rows) and all(as_float(row, "lower_95") > 1 for row in omission_rows),
            "one or more reported omission intervals are not above reference",
        )

    input_path = "results/model_input_sensitivity_summary/l10207_input_sensitivity.tsv"
    if (ROOT / input_path).is_file():
        rows = read_tsv(input_path)
        check.require(len(rows) == 6, "input-sensitivity table does not contain six refits")
        check.require(
            all(as_float(row, "lower_95") > 1 for row in rows),
            "one or more input-refit intervals are not above reference",
        )

    recovery_path = "results/model_main/recovery_summary.tsv"
    if (ROOT / recovery_path).is_file():
        rows = {row["parameter_type"]: row for row in read_tsv(recovery_path)}
        growth = rows.get("lineage_growth")
        imports = rows.get("import_scale")
        check.require(
            growth is not None
            and as_float(growth, "coverage_95") >= 0.8
            and as_float(growth, "median_absolute_log_error") <= 0.2,
            "lineage-growth recovery gate is not retained",
        )
        check.require(
            imports is not None
            and as_float(imports, "median_absolute_log_error") > 0.5
            and as_float(imports, "correlation_truth_posterior_median") < 0.7,
            "import-scale non-identifiability boundary is not retained",
        )

    annual_path = "figures/source_data/figure1b_annual_genomes.tsv"
    if (ROOT / annual_path).is_file():
        total = sum(int(row["n_sampled_genomes"]) for row in read_tsv(annual_path))
        check.require(total == 774, "annual focal-genome table does not sum to 774")


def validate_eid_archive_clock(check: Validation) -> None:
    public_rows = read_tsv("data/derived/public_genome_availability.tsv")
    required_fields = {
        "collection_upper_effective",
        "ena_first_public_date",
        "assembly_release_date",
        "public_route",
        "subnational_location",
        "location_source",
        "location_resolution",
        "nas_snapshot_record_id",
        "lag_min_days",
        "lag_max_days",
    }
    check.require(
        bool(public_rows) and required_fields.issubset(public_rows[0]),
        "EID accession-level availability schema is incomplete",
    )
    check.require(
        len(public_rows) == 774
        and len({row["tree_sample_id"] for row in public_rows}) == 774,
        "EID accession-level table is not the unique frozen 774-record focal cohort",
    )

    expected_shifts = {
        "AUS": ("3", "2024-08-01", "2024-08-31", "2025-01-22", "144", "174"),
        "CHN": ("5", "2023-01-01", "2023-03-27", "2024-10-04", "557", "642"),
        "JPN": ("5", "2024-01-01", "2024-12-09", "2025-05-12", "154", "497"),
    }
    shift_rows = {
        row["country_iso3"]: row
        for row in read_tsv("results/public_availability/eid_detection_clock_shift.tsv")
    }
    shift_fields = (
        "cumulative_genome_threshold",
        "collection_detection_lower",
        "collection_detection_upper",
        "public_detection_date",
        "clock_shift_min_days",
        "clock_shift_max_days",
    )
    for country, expected in expected_shifts.items():
        row = shift_rows.get(country, {})
        observed = tuple(row.get(field, "") for field in shift_fields)
        check.require(observed == expected, f"EID {country} clock displacement changed")

    expected_lags = {
        "AUS": ("138", "168"),
        "CHN": ("398.0", "446.0"),
        "JPN": ("249.5", "261.0"),
    }
    lag_rows = {
        row["country_iso3"]: row
        for row in read_tsv("results/public_availability/eid_country_lineage_lag_summary.tsv")
        if row.get("lineage") == "L1_02.07"
    }
    for country, expected in expected_lags.items():
        row = lag_rows.get(country, {})
        observed = (
            row.get("median_lag_min_days", ""),
            row.get("median_lag_max_days", ""),
        )
        check.require(observed == expected, f"EID {country} median lag interval changed")

    candidates = {
        row["project_id"]: row
        for row in read_tsv("results/public_availability/eid_external_candidate_summary.tsv")
    }
    prj = candidates.get("PRJNA1071282", {})
    observed_boundary = tuple(
        prj.get(field, "")
        for field in (
            "n_explicit_pertussis_runs",
            "n_frozen_tree_tips",
            "n_frozen_target_lineage",
            "n_frozen_resurgence_target",
        )
    )
    check.require(
        observed_boundary == ("734", "16", "6", "3"),
        "PRJNA1071282 public boundary changed from 734/16/6/3",
    )

    snapshot_rows = read_tsv("provenance/EID_NAS_SNAPSHOT_MANIFEST.tsv")
    check.require(
        len(snapshot_rows) == 5
        and all(row["record_type"] == "repository_snapshot" for row in snapshot_rows),
        "public EID snapshot manifest must contain only five released snapshots",
    )


def is_blank(value: str) -> bool:
    return value.strip() in {"", "NA", "N/A", "na", "null", "None"}


def fong_row(row: dict[str, str], biosamples: set[str]) -> bool:
    return (
        row.get("source_dataset") == FONG_SOURCE_DATASET
        or row.get("biosample_accession") in biosamples
        or row.get("sample_accession") in biosamples
        or row.get("sample_id", "").startswith("24-BPE-")
    )


def validate_redistribution_boundary(check: Validation) -> None:
    transmission_path = ROOT / "data/derived/transmission_genome_records.tsv"
    if not transmission_path.is_file():
        check.require(False, "missing data/derived/transmission_genome_records.tsv")
        return

    transmission_rows = read_tsv("data/derived/transmission_genome_records.tsv")
    fong_biosamples = {
        row["biosample_accession"]
        for row in transmission_rows
        if row.get("source_dataset") == FONG_SOURCE_DATASET
        and row.get("biosample_accession")
    }
    check.require(
        len(fong_biosamples) == 154,
        "Fong-linked BioSample count changed from the frozen public boundary",
    )
    check.require(
        (ROOT / "figures/source_data/figure4e_australia_ct_curve.tsv").is_file(),
        "aggregate Australian Ct recovery curve is missing",
    )
    check.require(
        not (ROOT / "data/derived/australia_sampling_process_records.tsv").exists(),
        "Australian specimen-level sampling-process table is still present",
    )

    observed_fong_rows = 0
    for relative_path in FONG_SANITISED_TABLES:
        path = ROOT / relative_path
        if not path.is_file():
            check.require(False, f"missing sanitised table: {relative_path}")
            continue
        rows = read_tsv(relative_path)
        for row_index, row in enumerate(rows, start=2):
            if not fong_row(row, fong_biosamples):
                continue
            observed_fong_rows += 1
            for column in RESTRICTED_FONG_COLUMNS.intersection(row):
                check.require(
                    is_blank(row[column]),
                    f"{relative_path}:{row_index} retains restricted {column}",
                )
            if "provenance_note" in row:
                check.require(
                    row["provenance_note"]
                    == "Third-party specimen metadata not redistributed",
                    f"{relative_path}:{row_index} has unsanitised provenance_note",
                )

    check.require(
        observed_fong_rows > 0,
        "no Fong-linked rows were found during redistribution-boundary checks",
    )


def validate_private_paths(check: Validation) -> None:
    for path in released_files():
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for marker in PRIVATE_PATH_MARKERS:
            if marker in text:
                check.require(
                    False,
                    f"private workspace path remains in {path.relative_to(ROOT)}",
                )
                break


def validate_manifest(check: Validation) -> None:
    if not FILE_MANIFEST.is_file():
        check.require(False, "missing provenance/FILE_MANIFEST_SHA256.tsv")
        return
    rows = read_tsv("provenance/FILE_MANIFEST_SHA256.tsv")
    observed = {
        path.relative_to(ROOT).as_posix(): (str(path.stat().st_size), sha256(path))
        for path in released_files()
    }
    recorded = {
        row["relative_path"]: (row["bytes"], row["sha256"])
        for row in rows
    }
    check.require(
        recorded == observed,
        "FILE_MANIFEST_SHA256.tsv is stale; run scripts/qa/build_file_manifest.py",
    )


def main() -> None:
    check = Validation()
    validate_structure(check)
    validate_links(check)
    validate_citation(check)
    validate_core_numbers(check)
    validate_eid_archive_clock(check)
    validate_redistribution_boundary(check)
    validate_private_paths(check)
    validate_manifest(check)
    check.finish()


if __name__ == "__main__":
    main()
