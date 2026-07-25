#!/usr/bin/env python3
"""Build panel-keyed source tables for the five final manuscript figures."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
ANALYSIS = ROOT / "analysis/genomic_transmission_dynamics"
OUT = ANALYSIS / "figures/source_data"


def copy(source: Path, name: str, provenance: list[dict[str, str]]) -> None:
    destination = OUT / name
    shutil.copy2(source, destination)
    provenance.append({"figure_source_file": name, "derived_from": str(source)})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    provenance: list[dict[str, str]] = []

    copy(ANALYSIS / "derived/country_month_cases.tsv", "figure1a_cases.tsv", provenance)
    lineage_path = (
        ANALYSIS
        / "lineages_gate_rescue4_depth/primary_finalized/model_lineage_assignments.tsv"
    )
    lineage = pd.read_csv(lineage_path, sep="\t")
    dated = lineage[
        lineage["tree_role"].eq("focal")
        & lineage["date_resolution"].isin(["day", "month"])
    ].copy()
    dated["model_month"] = (
        pd.to_datetime(dated["date_lower"]).dt.to_period("M").dt.to_timestamp()
    )
    monthly = (
        dated.groupby(["country_iso3", "model_month", "epidemic_period"])
        .size()
        .rename("n_sampled_genomes")
        .reset_index()
    )
    monthly.to_csv(OUT / "figure1b_monthly_genomes.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure1b_monthly_genomes.tsv",
            "derived_from": str(lineage_path),
        }
    )

    qc = pd.read_csv(
        ANALYSIS / "phylogeny/qc_gate_rescue4/uniform_sequence_qc.tsv", sep="\t"
    )
    alignment = pd.read_csv(
        ANALYSIS
        / "phylogeny/results/gate_rescue4_depth/alignment_missingness_final.tsv",
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
        ANALYSIS
        / "phylogeny/results/gate_rescue4_depth/tree_tip_metadata.tsv",
        "figure2a_tree_tip_metadata.tsv",
        provenance,
    )
    copy(
        ANALYSIS
        / "phylogeography_gate_rescue4_depth/events_thr0_5/tip_persistence_reseeding_support.tsv",
        "figure2b_tip_ancestry_support.tsv",
        provenance,
    )
    copy(
        ANALYSIS
        / "phylogeography_gate_rescue4_depth/events_thr0_5/country_lineage_phylogeographic_summary.tsv",
        "figure2c_country_lineage_ancestry.tsv",
        provenance,
    )
    tree_manifest = pd.DataFrame(
        [
            {
                "tree_file": str(
                    ANALYSIS
                    / "phylogeny/results/gate_rescue4_depth/gtd_gate_rescue4_depth_core_snp.treefile"
                ),
                "support_tree_file": str(
                    ANALYSIS
                    / "phylogeny/results/gate_rescue4_depth/gtd_gate_rescue4_depth_core_snp.contree"
                ),
                "rooted_tree_file": str(
                    ANALYSIS
                    / "phylogeography_gate_rescue4_depth/gtd_primary_final_midpoint_rooted.tree"
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
        ANALYSIS
        / "phylogeography_gate_rescue4_depth/events_thr0_5/introduction_events.tsv"
    )
    events = pd.read_csv(event_path, sep="\t")
    post = events[events["high_support_post_reseeding"]].copy()
    post.to_csv(OUT / "figure3a_post_reseeding_edges.tsv", sep="\t", index=False)
    post[post["successful_sampled_cluster"]].to_csv(
        OUT / "figure3b_successful_sampled_clusters.tsv", sep="\t", index=False
    )
    provenance.extend(
        [
            {
                "figure_source_file": "figure3a_post_reseeding_edges.tsv",
                "derived_from": str(event_path),
            },
            {
                "figure_source_file": "figure3b_successful_sampled_clusters.tsv",
                "derived_from": str(event_path),
            },
        ]
    )
    copy(
        ANALYSIS
        / "results_gate_rescue4_depth/main/successful_sampled_cluster_probability.tsv",
        "figure3c_cluster_probability.tsv",
        provenance,
    )
    threshold_rows = []
    for threshold in ("0_5", "0_7", "0_9"):
        path = (
            ANALYSIS
            / f"phylogeography_gate_rescue4_depth/events_thr{threshold}/phylogeography_validation.json"
        )
        payload = json.loads(path.read_text())
        threshold_rows.append(payload)
    pd.DataFrame(threshold_rows).to_csv(
        OUT / "figure3d_threshold_sensitivity.tsv", sep="\t", index=False
    )
    provenance.append(
        {
            "figure_source_file": "figure3d_threshold_sensitivity.tsv",
            "derived_from": "events_thr0_5;events_thr0_7;events_thr0_9 validation JSON",
        }
    )

    main_results = ANALYSIS / "results_gate_rescue4_depth/main"
    no_project = ANALYSIS / "results_gate_rescue4_depth/no_project"
    copy(
        main_results / "lineage_relative_transmission.tsv",
        "figure4a_lineage_growth_main.tsv",
        provenance,
    )
    main_growth = pd.read_csv(main_results / "lineage_relative_transmission.tsv", sep="\t")
    no_growth = pd.read_csv(no_project / "lineage_relative_transmission.tsv", sep="\t")
    growth = main_growth.merge(
        no_growth, on="lineage", suffixes=("_main", "_no_project"), validate="one_to_one"
    )
    growth.to_csv(OUT / "figure4b_lineage_growth_sensitivity.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure4b_lineage_growth_sensitivity.tsv",
            "derived_from": "main and no-project lineage_relative_transmission.tsv",
        }
    )
    copy(
        main_results / "monthly_sampling_corrected_lineage_shares.tsv",
        "figure4c_monthly_latent_lineage_shares.tsv",
        provenance,
    )
    copy(
        main_results / "sampling_corrected_lineage_shares.tsv",
        "figure4d_raw_vs_corrected_shares.tsv",
        provenance,
    )

    copy(
        main_results / "monthly_model_and_counterfactuals.tsv",
        "figure5abc_monthly_counterfactuals.tsv",
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
    cf.to_csv(OUT / "figure5d_counterfactual_summary.tsv", sep="\t", index=False)
    provenance.append(
        {
            "figure_source_file": "figure5d_counterfactual_summary.tsv",
            "derived_from": "main and no-project counterfactual_summary.tsv",
        }
    )
    copy(
        ANALYSIS / "model/australia_sampling/australia_ct_success_curve.tsv",
        "figure5e_australia_ct_curve.tsv",
        provenance,
    )
    copy(
        ANALYSIS
        / "model_gate_rescue4_depth/identifiability_recovery/recovery_summary.tsv",
        "figure5f_identifiability_recovery.tsv",
        provenance,
    )
    pd.DataFrame(provenance).to_csv(
        OUT / "FIGURE_SOURCE_PROVENANCE.tsv", sep="\t", index=False
    )
    print(f"wrote {len(provenance)} source-data artifacts to {OUT}")


if __name__ == "__main__":
    main()
