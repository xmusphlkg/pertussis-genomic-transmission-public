#!/usr/bin/env python3
"""Build panel-keyed source tables for four main figures and linked supplements."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/source_data"


def copy(source: Path, name: str, provenance: list[dict[str, str]]) -> None:
    destination = OUT / name
    shutil.copy2(source, destination)
    provenance.append({"figure_source_file": name, "derived_from": str(source)})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    provenance: list[dict[str, str]] = []

    copy(ROOT / "data/derived/country_month_cases.tsv", "figure1a_cases.tsv", provenance)
    copy(
        ROOT / "data/derived/figure1a_native_resolution_surveillance.tsv",
        "figure1a_native_resolution_surveillance.tsv",
        provenance,
    )
    copy(
        ROOT / "data/derived/figure1a_belgium_france_annual_cases.tsv",
        "figure1a_belgium_france_annual_cases.tsv",
        provenance,
    )
    lineage_path = ROOT / "results/lineages/primary_finalized/model_lineage_assignments.tsv"
    lineage = pd.read_csv(lineage_path, sep="\t")
    dated = lineage[lineage["tree_role"].eq("focal")].copy()
    dated["sampling_year"] = pd.to_datetime(
        dated["date_lower"], errors="raise"
    ).dt.year
    annual = (
        dated.groupby(["country_iso3", "sampling_year", "epidemic_period"])
        .size()
        .rename("n_sampled_genomes")
        .reset_index()
    )
    if int(annual["n_sampled_genomes"].sum()) != len(dated):
        raise ValueError("Figure 1B annual bins do not contain every focal genome")
    annual.to_csv(OUT / "figure1b_annual_genomes.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure1b_annual_genomes.tsv",
            "derived_from": str(lineage_path),
        }
    )

    qc = pd.read_csv(
        ROOT / "results/phylogeny/uniform_sequence_qc.tsv", sep="\t"
    )
    alignment = pd.read_csv(
        ROOT / "results/phylogeny/alignment_missingness_final.tsv",
        sep="\t",
    )
    flow = pd.DataFrame(
        [
            ("selected", len(qc)),
            ("uniform_qc_pass", int(qc["tree_include_after_uniform_qc"].sum())),
            ("final_alignment_pass", int(alignment["alignment_qc_pass"].sum())),
        ],
        columns=["stage", "n_genomes"],
    )
    flow.to_csv(OUT / "figure1c_cohort_flow.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure1c_cohort_flow.tsv",
            "derived_from": "uniform_sequence_qc.tsv;alignment_missingness_final.tsv",
        }
    )

    copy(
        ROOT / "results/phylogeny/tree_tip_metadata.tsv",
        "figure2a_tree_tip_metadata.tsv",
        provenance,
    )
    copy(
        ROOT
        / "results/phylogeography/events_thr0_5/tip_persistence_reseeding_support.tsv",
        "figure2b_tip_ancestry_support.tsv",
        provenance,
    )
    copy(
        ROOT
        / "results/phylogeography/events_thr0_5/country_lineage_phylogeographic_summary.tsv",
        "figure2c_country_lineage_ancestry.tsv",
        provenance,
    )
    tree_manifest = pd.DataFrame(
        [
            {
                "tree_file": str(
                    ROOT / "results/phylogeny/gtd_gate_rescue4_depth_core_snp.treefile"
                ),
                "support_tree_file": str(
                    ROOT / "results/phylogeny/gtd_gate_rescue4_depth_core_snp.contree"
                ),
                "rooted_tree_file": str(
                    ROOT / "results/phylogeography/gtd_primary_final_midpoint_rooted.tree"
                ),
            }
        ]
    )
    tree_manifest.to_csv(OUT / "figure2_tree_manifest.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure2_tree_manifest.tsv",
            "derived_from": "final IQ-TREE and midpoint-rooted outputs",
        }
    )

    event_path = (
        ROOT / "results/phylogeography/events_thr0_5/introduction_events.tsv"
    )
    events = pd.read_csv(event_path, sep="\t")
    post = events[events["high_support_post_reseeding"]].copy()
    post.to_csv(OUT / "figureS2de_post_reseeding_edges.tsv", sep="\t", index=False)
    post[post["interval_compatible_success"]].to_csv(
        OUT / "figure2d_successful_sampled_clusters.tsv", sep="\t", index=False
    )
    provenance.extend(
        [
            {
                "figure_source_file": "figureS2de_post_reseeding_edges.tsv",
                "derived_from": str(event_path),
            },
            {
                "figure_source_file": "figure2d_successful_sampled_clusters.tsv",
                "derived_from": str(event_path),
            },
        ]
    )
    threshold_rows = []
    for threshold in ("0_5", "0_7", "0_9"):
        path = (
            ROOT
            / f"results/phylogeography/events_thr{threshold}/phylogeography_validation.json"
        )
        payload = json.loads(path.read_text())
        threshold_rows.append(payload)
    pd.DataFrame(threshold_rows).to_csv(
        OUT / "figureS2f_threshold_sensitivity.tsv", sep="\t", index=False
    )
    provenance.append(
        {
            "figure_source_file": "figureS2f_threshold_sensitivity.tsv",
            "derived_from": "events_thr0_5;events_thr0_7;events_thr0_9 validation JSON",
        }
    )

    main_results = ROOT / "results/model_main"
    no_project = ROOT / "results/model_no_project"
    copy(
        main_results / "lineage_relative_transmission.tsv",
        "figure3a_lineage_growth_main.tsv",
        provenance,
    )
    main_growth = pd.read_csv(main_results / "lineage_relative_transmission.tsv", sep="\t")
    no_growth = pd.read_csv(no_project / "lineage_relative_transmission.tsv", sep="\t")
    growth = main_growth.merge(
        no_growth, on="lineage", suffixes=("_main", "_no_project"), validate="one_to_one"
    )
    growth.to_csv(OUT / "figureS4a_lineage_growth_sensitivity.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figureS4a_lineage_growth_sensitivity.tsv",
            "derived_from": "main and no-project lineage_relative_transmission.tsv",
        }
    )
    main_pairwise = pd.read_csv(
        main_results / "l10207_pairwise_growth.tsv", sep="\t"
    )
    no_pairwise = pd.read_csv(
        no_project / "l10207_pairwise_growth.tsv", sep="\t"
    )
    pairwise = main_pairwise.merge(
        no_pairwise,
        on=["numerator", "denominator"],
        suffixes=("_main", "_no_project"),
        validate="one_to_one",
    )
    pairwise.to_csv(
        OUT / "figure3b_l10207_pairwise_growth.tsv",
        sep="\t",
        index=False,
    )
    provenance.append(
        {
            "figure_source_file": "figure3b_l10207_pairwise_growth.tsv",
            "derived_from": "main and no-project l10207_pairwise_growth.tsv",
        }
    )
    copy(
        ROOT
        / "results/model_growth_robustness/l1_02_07_growth_robustness.tsv",
        "figure3c_l10207_growth_robustness.tsv",
        provenance,
    )
    copy(
        ROOT
        / "results/model_growth_robustness/"
        "selection_cap_weighted_l1_02_07_shares.tsv",
        "figure3d_selection_cap_weighted_l10207_shares.tsv",
        provenance,
    )
    copy(
        main_results / "sampling_corrected_lineage_shares.tsv",
        "figureS4b_raw_vs_corrected_shares.tsv",
        provenance,
    )

    copy(
        main_results / "monthly_model_and_counterfactuals.tsv",
        "figure4abc_monthly_counterfactuals.tsv",
        provenance,
    )
    main_cf = pd.read_csv(main_results / "counterfactual_summary.tsv", sep="\t")
    no_cf = pd.read_csv(no_project / "counterfactual_summary.tsv", sep="\t")
    cf = main_cf.merge(
        no_cf,
        on=["country_iso3", "scenario"],
        suffixes=("_main", "_no_project"),
        validate="one_to_one",
    )
    cf.to_csv(
        OUT / "figure4d_counterfactual_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    provenance.append(
        {
            "figure_source_file": "figure4d_counterfactual_summary.tsv",
            "derived_from": "main and no-project counterfactual_summary.tsv",
        }
    )
    copy(
        ROOT / "results/australia_sampling/australia_ct_success_curve.tsv",
        "figure4e_australia_ct_curve.tsv",
        provenance,
    )
    copy(
        ROOT / "results/model_main/recovery_summary.tsv",
        "figure4f_identifiability_recovery.tsv",
        provenance,
    )
    pd.DataFrame(provenance).to_csv(
        OUT / "FIGURE_SOURCE_PROVENANCE.tsv", sep="\t", index=False
    )
    print(f"wrote {len(provenance)} source-data artifacts to {OUT}")


if __name__ == "__main__":
    main()
