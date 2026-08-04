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


def test_initial_lineage_prior_precedes_model_and_unused_input_was_removed() -> None:
    validation = json.loads(
        (ROOT / "data/model_inputs/joint_model_data_validation.json").read_text()
    )
    model_data = json.loads(
        (ROOT / "data/model_inputs/joint_model_data.json").read_text()
    )
    prior = pd.read_csv(
        ROOT / "data/model_inputs/initial_lineage_prior.tsv", sep="\t"
    )
    symmetric = json.loads(
        (
            ROOT
            / "data/model_sensitivity/initial_state_symmetric/joint_model_data.json"
        ).read_text()
    )
    stan_source = (
        ROOT / "scripts/model/gtd_joint_transmission_sampling.stan"
    ).read_text()

    assert validation["initial_prior_mode"] == "historical"
    assert validation["initial_prior_end_year"] == 2018
    assert validation["initial_prior_strictly_precedes_model"] is True
    assert validation["n_historical_tips_in_initial_prior"] == 164
    assert int(prior["n_historical_tips"].sum()) == 164
    assert set(prior.groupby("country_iso3")["n_historical_tips"].sum()) == {
        49,
        51,
        64,
    }
    assert "persistence_support" not in model_data
    assert "persistence_support" not in stan_source
    assert all(
        alpha == 0.5
        for country_alpha in symmetric["initial_alpha"]
        for alpha in country_alpha
    )


def test_reader_facing_narrative_obeys_inference_boundaries() -> None:
    reader_text = "\n".join(
        (ROOT / path).read_text()
        for path in ("README.md", "DATA_DICTIONARY.md", "LETTER_EVIDENCE_MAP.md")
    ).lower()

    for banned in (
        "sampling adjustment",
        "sampling-aware",
        "independently supported",
        "growth neutralisation",
        "growth-neutralised",
        "country-specific resurgence ancestry",
    ):
        assert banned not in reader_text
    assert "sampled ancestry" in reader_text
    assert "project adjustment" in reader_text
    assert "national lineage prevalence" in reader_text
    assert "individual transmission links" in reader_text
    assert "causal contribution" in reader_text
    assert "biological fitness" in reader_text
    assert "identifiable import scale" in reader_text


def test_no_new_exposure_case_difference_is_bounded() -> None:
    counterfactual = pd.read_csv(
        ROOT / "results/model_main/counterfactual_summary.tsv", sep="\t"
    )
    rows = counterfactual[
        counterfactual["scenario"].eq("no_new_exposure_case_difference_fraction")
    ]
    assert set(rows["country_iso3"]) == {"AUS", "CHN", "JPN"}
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


def test_l10207_exceeds_every_other_lineage_and_anchor_sensitivity_is_positive() -> None:
    for model_dir in ("model_main", "model_no_project"):
        pairwise = pd.read_csv(
            ROOT / f"results/{model_dir}/l10207_pairwise_growth.tsv",
            sep="\t",
        )
        assert set(pairwise["denominator"]) == {
            "L1_01.02",
            "L1_02.05",
            "L1_02.06",
            "Other",
        }
        assert pairwise["lower_95"].gt(1).all()
        assert pairwise["posterior_probability_above_one"].eq(1).all()

    anchors = pd.read_csv(
        ROOT / "results/model_anchor_sensitivity/anchor_scenario_sensitivity.tsv",
        sep="\t",
    )
    neutral = anchors[anchors["scenario"].eq("l10207_neutral_growth")]
    assert len(neutral) == 24
    assert neutral["cumulative_difference_fraction_lower_95"].gt(0).all()
    assert neutral["posterior_probability_difference_above_zero"].eq(1).all()


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
    assert int(post["robust_success"].sum()) == 1
    assert int(post["interval_compatible_success"].sum()) == 3
    assert post["successful_sampled_cluster"].equals(post["robust_success"])
    qualifying = post[post["interval_compatible_success"]]
    assert set(qualifying["destination_country"]) == {"FRA", "JPN"}
    japan = qualifying[qualifying["destination_country"].eq("JPN")]
    assert japan["assigned_span_min_days"].tolist() == [189]
    assert japan["assigned_span_max_days"].tolist() == [189]


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


def test_country_ancestry_order_survives_topology_and_reference_perturbations() -> None:
    ranks = pd.read_csv(
        ROOT
        / "results/phylogeography_sensitivity_summary/country_ancestry_rank_stability.tsv",
        sep="\t",
    )
    assert ranks["analysis_id"].nunique() == 6
    for analysis_id, group in ranks.groupby("analysis_id"):
        local_order = (
            group.sort_values("local_rank_median")["country_iso3"].tolist()
        )
        reseeding_order = (
            group.sort_values("reseeding_rank_median")["country_iso3"].tolist()
        )
        assert local_order == ["CHN", "BEL", "AUS", "JPN", "FRA"], analysis_id
        assert reseeding_order == ["FRA", "JPN", "AUS", "BEL", "CHN"], analysis_id
    jackknife = ranks[ranks["analysis_class"].isin(
        ["reference_sampling", "background_sampling"]
    )]
    assert jackknife.loc[
        jackknife["country_iso3"].eq("CHN"), "probability_highest_local"
    ].eq(1).all()
    assert jackknife.loc[
        jackknife["country_iso3"].eq("FRA"), "probability_highest_reseeding"
    ].ge(0.90).all()


def test_phylogeographic_exposure_sensitivity_inputs_change_only_intended_rules() -> None:
    base = json.loads(
        (ROOT / "data/model_inputs/joint_model_data_validation.json").read_text()
    )
    threshold_07 = json.loads(
        (
            ROOT
            / "data/model_sensitivity/threshold_0_7/joint_model_data_validation.json"
        ).read_text()
    )
    threshold_09 = json.loads(
        (
            ROOT
            / "data/model_sensitivity/threshold_0_9/joint_model_data_validation.json"
        ).read_text()
    )
    midpoint = json.loads(
        (
            ROOT
            / "data/model_sensitivity/time_midpoint/joint_model_data_validation.json"
        ).read_text()
    )
    uniform = json.loads(
        (
            ROOT
            / "data/model_sensitivity/time_interval_uniform/"
            "joint_model_data_validation.json"
        ).read_text()
    )
    assert [
        base["n_unique_phylogeographic_events"],
        threshold_07["n_unique_phylogeographic_events"],
        threshold_09["n_unique_phylogeographic_events"],
    ] == [31, 30, 28]
    assert (
        base["total_phylogeographic_exposure_weight"]
        > threshold_07["total_phylogeographic_exposure_weight"]
        > threshold_09["total_phylogeographic_exposure_weight"]
    )
    assert midpoint["total_phylogeographic_exposure_weight"] == (
        base["total_phylogeographic_exposure_weight"]
    )
    assert uniform["total_phylogeographic_exposure_weight"] == (
        base["total_phylogeographic_exposure_weight"]
    )
    assert uniform["n_phylogeographic_event_lineage_rows"] > (
        base["n_phylogeographic_event_lineage_rows"]
    )


def test_model_refits_preserve_l10207_signal_across_input_definitions() -> None:
    summary = ROOT / "results/model_input_sensitivity_summary"
    growth = pd.read_csv(summary / "l10207_input_sensitivity.tsv", sep="\t")
    diagnostics = pd.read_csv(
        summary / "model_input_sensitivity_diagnostics.tsv", sep="\t"
    )
    scenarios = pd.read_csv(summary / "scenario_input_sensitivity.tsv", sep="\t")
    assert len(growth) == 6
    assert growth["lower_95"].gt(1).all()
    assert diagnostics["divergent_transitions"].eq(0).all()
    assert diagnostics["maximum_treedepth_hits"].eq(0).all()
    assert diagnostics["maximum_rhat"].lt(1.01).all()
    neutral = scenarios[
        scenarios["scenario"].eq("l10207_growth_scenario_difference_fraction")
    ]
    assert len(neutral) == 18
    assert neutral["lower_95"].gt(0).all()
    china_exposure = scenarios[
        scenarios["country_iso3"].eq("CHN")
        & scenarios["scenario"].eq("no_new_exposure_case_difference_fraction")
    ]
    assert china_exposure[["mean", "median", "lower_95", "upper_95"]].eq(0).all().all()


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


def test_letter_figures_have_frozen_source_data() -> None:
    source = ROOT / "figures/source_data"
    manifests = [
        ROOT / "figures/letter/LETTER_FIGURE_SOURCE_FILES.tsv",
        ROOT / "figures/letter/LETTER_SUPPLEMENTARY_FIGURE_SOURCE_FILES.tsv",
    ]
    expected = set()
    for manifest in manifests:
        mapping = pd.read_csv(manifest, sep="\t")
        expected.update(mapping["source_file"])
    assert len(expected) == 19
    assert all((ROOT / name).is_file() for name in expected)

    letter_source_files = {
        path.relative_to(ROOT).as_posix()
        for path in source.glob("*.tsv")
    }
    assert letter_source_files == {
        name for name in expected if name.startswith("figures/source_data/")
    }

    annual = pd.read_csv(source / "figure1b_annual_genomes.tsv", sep="\t")
    assert int(annual["n_sampled_genomes"].sum()) == 774
    assert (
        annual.groupby("country_iso3")["n_sampled_genomes"].sum().to_dict()
        == {"AUS": 97, "BEL": 119, "CHN": 388, "FRA": 59, "JPN": 111}
    )
    france = annual.loc[annual["country_iso3"].eq("FRA")]
    assert france.set_index("sampling_year")["n_sampled_genomes"].to_dict() == {
        2000: 1,
        2007: 1,
        2008: 1,
        2009: 1,
        2014: 2,
        2023: 7,
        2024: 46,
    }


def test_figure1a_belgium_france_annual_cases_are_not_disaggregated() -> None:
    path = ROOT / "data/derived/figure1a_belgium_france_annual_cases.tsv"
    cases = pd.read_csv(path, sep="\t")
    assert len(cases) == 22
    assert set(cases["country_iso3"]) == {"BEL", "FRA"}
    assert set(cases["year"]) == set(range(2015, 2026))
    assert set(cases["time_resolution"]) == {"annual"}
    assert cases["comparability_note"].str.contains(
        "must not be disaggregated into monthly counts", regex=False
    ).all()
    observed_2024 = (
        cases.loc[cases["year"].eq(2024)]
        .set_index("country_iso3")["reported_cases"]
        .to_dict()
    )
    assert observed_2024 == {"BEL": 3078.0, "FRA": 464.0}


def test_figure1a_preserves_native_surveillance_resolution_and_scope() -> None:
    source = ROOT / "data/derived"
    native = pd.read_csv(
        source / "figure1a_native_resolution_surveillance.tsv",
        sep="\t",
    )
    assert len(native) == 664
    assert native.groupby("country_iso3").size().to_dict() == {
        "AUS": 84,
        "BEL": 36,
        "CHN": 84,
        "FRA": 96,
        "JPN": 364,
    }
    resolution = (
        native.groupby("country_iso3")["time_resolution"].unique().map(list)
    )
    assert resolution.to_dict() == {
        "AUS": ["monthly"],
        "BEL": ["monthly"],
        "CHN": ["monthly"],
        "FRA": ["monthly"],
        "JPN": ["weekly"],
    }
    observed_titles = (
        native.groupby("country_iso3")["y_axis_title"]
        .apply(lambda values: set(values))
        .to_dict()
    )
    assert observed_titles == {
        "AUS": {"Notifications per month"},
        "BEL": {"NRC-confirmed cases per month"},
        "CHN": {"Reported cases per month"},
        "FRA": {
            "PCR positivity per month (%)",
            "PCR-positive tests per month",
        },
        "JPN": {"Reported cases per week"},
    }
    assert set(native.loc[native["country_iso3"].eq("BEL"), "year"]) == {
        2019,
        2023,
        2024,
    }
    assert native.loc[
        native["country_iso3"].eq("BEL"),
        "value_status",
    ].eq("approximate").all()
    france_counts = native.loc[
        native["metric"].eq("pcr_positive_tests")
    ]
    assert len(france_counts) == 24
    assert france_counts["value_status"].eq("reported").all()
    assert france_counts.groupby("year")["value"].sum().to_dict() == {
        2023: 518,
        2024: 38847,
    }

    observed_japan = (
        native.loc[
            native["country_iso3"].eq("JPN"),
            ["observation_date", "value"],
        ]
        .assign(
            observation_date=lambda frame: pd.to_datetime(
                frame["observation_date"]
            ),
        )
    )
    assert observed_japan["observation_date"].is_monotonic_increasing


def test_figure_file_contract() -> None:
    letter = ROOT / "figures/letter"
    main_manifest = pd.read_csv(letter / "RENDER_MANIFEST.tsv", sep="\t")
    supp_manifest = pd.read_csv(
        letter / "SUPPLEMENTARY_RENDER_MANIFEST.tsv", sep="\t"
    )
    assert len(main_manifest) == 2
    assert len(supp_manifest) == 3
    assert set(main_manifest["panels"]) == {"A-E", "A-F"}
    assert set(supp_manifest["panels"]) == {"A-C", "A-D", "A-B"}
    assert len(list(letter.glob("Figure_*.png"))) == 2
    assert len(list(letter.glob("Figure_*.pdf"))) == 2
    assert len(list(letter.glob("Supplementary_Figure_*.png"))) == 3
    assert not list(letter.glob("Supplementary_Figure_*.pdf"))
    assert not list(letter.glob("*.svg"))
    assert not list(letter.glob("*.tif*"))
    assert all(path.stat().st_size > 50_000 for path in letter.glob("*.png"))
    assert all(path.stat().st_size > 10_000 for path in letter.glob("*.pdf"))
