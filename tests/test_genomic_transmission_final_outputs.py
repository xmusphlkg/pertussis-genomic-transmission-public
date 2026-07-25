from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_lineages_were_frozen_without_case_counts() -> None:
    validation = json.loads(
        (ROOT / "results/lineages/primary_finalized/model_lineage_validation.json").read_text()
    )
    assert validation["case_data_used_for_lineage_definition"] is False
    assert validation["n_primary_model_lineages_excluding_other"] >= 3
    assert validation["gate_pass"] is True


def test_lineage_partition_sensitivity_preserves_high_growth_lineage() -> None:
    validation = json.loads(
        (ROOT / "results/lineages/sensitivity/lineage_partition_sensitivity.json").read_text()
    )
    assert validation["adjusted_rand_index_level2"] >= 0.95
    assert validation["minimum_primary_member_recovery_among_model_lineages"] >= 0.95
    assert validation["high_growth_lineage_L1_02_07_jaccard"] >= 0.95


def test_primary_stan_sampling_passed_diagnostics() -> None:
    diagnostics = json.loads(
        (ROOT / "results/model_main/sampling_diagnostics.json").read_text()
    )
    assert diagnostics["divergent_transitions"] == 0
    assert diagnostics["maximum_treedepth_hits"] == 0
    assert diagnostics["maximum_rhat"] < 1.01
    assert diagnostics["minimum_bulk_neff"] >= 400


def test_background_genomes_do_not_enter_joint_observation_denominator() -> None:
    validation = json.loads(
        (ROOT / "data/model_inputs/joint_model_data_validation.json").read_text()
    )
    strata = pd.read_csv(ROOT / "data/model_inputs/genome_observation_strata.tsv", sep="\t")
    assert validation["n_genomes_in_observation_model"] == int(
        strata.filter(regex=r"^n_").to_numpy().sum()
    )
    assert set(strata["country_iso3"]) == {"AUS", "CHN", "JPN"}
    assert validation["case_counts_used_for_lineage_definition"] is False
    assert validation["exact_tmrca_used"] is False


def test_no_new_introduction_case_reduction_is_bounded() -> None:
    counterfactual = pd.read_csv(
        ROOT / "results/model_main/counterfactual_summary.tsv", sep="\t"
    )
    rows = counterfactual[
        counterfactual["scenario"].eq("no_new_introduction_case_reduction_fraction")
    ]
    for column in ("mean", "median", "lower_95", "upper_95"):
        assert rows[column].between(0, 1).all()


def test_only_l10207_has_robust_positive_relative_transmission() -> None:
    primary = pd.read_csv(
        ROOT / "results/model_main/lineage_relative_transmission.tsv", sep="\t"
    ).set_index("lineage")
    sensitivity = pd.read_csv(
        ROOT / "results/model_no_project/lineage_relative_transmission.tsv", sep="\t"
    ).set_index("lineage")
    robust = set(
        primary.index[
            primary["lower_95"].gt(1) & sensitivity["lower_95"].gt(1)
        ]
    )
    assert robust == {"L1_02.07"}


def test_post_reseeding_events_exclude_pre_2023_country_descendants() -> None:
    events = pd.read_csv(
        ROOT / "results/phylogeography/events_thr0_5/introduction_events.tsv",
        sep="\t",
    )
    post = events[events["high_support_post_reseeding"]]
    assert len(post) == 75
    assert post["n_descendant_prepandemic"].eq(0).all()
    assert post["n_descendant_pandemic"].eq(0).all()
    assert post["n_descendant_resurgence"].gt(0).all()
    assert int(post["successful_sampled_cluster"].sum()) == 3


def test_country_ancestry_is_robust_to_prespecified_root() -> None:
    primary = pd.read_csv(
        ROOT / "results/phylogeography/events_thr0_5/country_phylogeographic_summary.tsv",
        sep="\t",
    ).set_index("country_iso3")
    alternate = pd.read_csv(
        ROOT
        / "results/phylogeography_alternative_root/events_thr0_5/country_phylogeographic_summary.tsv",
        sep="\t",
    ).set_index("country_iso3")
    columns = ["mean_post_import_support", "mean_local_persistence_support"]
    assert primary.index.equals(alternate.index)
    assert (primary[columns] - alternate[columns]).abs().to_numpy().max() < 0.001


def test_identifiability_supports_growth_but_not_import_scale() -> None:
    gate = json.loads((ROOT / "results/model_main/identifiability_gate.json").read_text())
    observed = gate["observed"]
    assert gate["pass"] is False
    assert observed["lineage_growth_coverage"] >= 0.8
    assert observed["lineage_growth_median_absolute_log_error"] <= 0.2
    assert observed["lineage_growth_truth_median_correlation"] >= 0.7
    assert observed["highest_lineage_rank_recovery"] >= 0.75
    assert (
        observed["import_scale_median_absolute_log_error"] > 0.5
        or observed["import_scale_truth_median_correlation"] < 0.7
    )


def test_cgmlst_sensitivity_is_complete_and_concordant() -> None:
    concordance = json.loads(
        (ROOT / "results/cgmlst/cgmlst_core_snp_concordance.json").read_text()
    )
    assert concordance["n_profiles_with_frozen_core_snp_lineage"] == 582
    assert concordance["n_pairwise_core_snp_comparisons"] == 169071
    assert concordance["core_snp_cgmlst_spearman_rho"] > 0.8
    assert concordance["nearest_neighbour_lineage_agreement_formal_lineages"] > 0.99


def test_five_main_figures_have_frozen_source_data() -> None:
    source = ROOT / "figures/source_data"
    provenance = pd.read_csv(source / "FIGURE_SOURCE_PROVENANCE.tsv", sep="\t")
    expected = set(provenance["figure_source_file"])
    assert len(expected) == 19
    assert all((source / name).is_file() for name in expected)
    flow = pd.read_csv(source / "figure1c_cohort_flow.tsv", sep="\t")
    assert dict(zip(flow["stage"], flow["n_genomes"])) == {
        "selected": 1188,
        "uniform_qc_pass": 1078,
        "final_alignment_pass": 989,
    }


def test_figure_file_contract() -> None:
    main = ROOT / "figures/main"
    supplementary = ROOT / "figures/supplementary"
    main_manifest = pd.read_csv(main / "RENDER_MANIFEST.tsv", sep="\t")
    supp_manifest = pd.read_csv(supplementary / "RENDER_MANIFEST.tsv", sep="\t")
    assert len(main_manifest) == 5
    assert len(supp_manifest) == 4
    assert set(main_manifest["formats"]) == {"PDF-vector;PNG-600dpi"}
    assert set(supp_manifest["format"]) == {"PNG-600dpi"}
    assert len(list(main.glob("*.png"))) == 5
    assert len(list(main.glob("*.pdf"))) == 5
    assert not list(main.glob("*.svg"))
    assert not list(main.glob("*.tif*"))
    assert len(list(supplementary.glob("*.png"))) == 4
    assert not list(supplementary.glob("*.pdf"))
    assert not list(supplementary.glob("*.svg"))
    assert not list(supplementary.glob("*.tif*"))
    assert all(path.stat().st_size > 50_000 for path in main.glob("*.png"))
    assert all(path.stat().st_size > 10_000 for path in main.glob("*.pdf"))
    assert all(path.stat().st_size > 50_000 for path in supplementary.glob("*.png"))
