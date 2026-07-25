#!/usr/bin/env python3
"""Convert PastML marginal states into sampling-aware introduction evidence.

This is a modular ("cut posterior") phylogeographic layer. It deliberately
reports support scores and date intervals rather than exact tMRCAs or
individual transmission links.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

FOCAL = ("AUS", "BEL", "CHN", "FRA", "JPN")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", type=Path)
    parser.add_argument("marginals", type=Path)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("lineages", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--transition-threshold", type=float, default=0.5)
    parser.add_argument("--success-size", type=int, default=5)
    parser.add_argument("--success-span-days", type=int, default=180)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    tree = Phylo.read(args.tree, "newick")
    metadata = pd.read_csv(args.metadata, sep="\t")
    lineage = pd.read_csv(args.lineages, sep="\t")[
        ["tree_sample_id", "primary_model_lineage_id"]
    ]
    metadata = metadata.merge(lineage, on="tree_sample_id", how="left", validate="one_to_one")
    metadata = metadata.set_index("tree_sample_id")

    marg = pd.read_csv(args.marginals, sep="\t").set_index("node")
    states = list(marg.columns)
    for clade in tree.find_clades():
        if not clade.name:
            raise ValueError("All internal nodes must be named before event inference")
        if clade.name not in marg.index:
            if clade.is_terminal() and clade.name in metadata.index:
                state = metadata.at[clade.name, "country_iso3"]
                if state not in states:
                    continent = metadata.at[clade.name, "continent"]
                    state = (
                        f"{str(continent).replace(' ', '_')}_other"
                        if pd.notna(continent)
                        else "Unknown_other"
                    )
                marg.loc[clade.name, :] = 0.0
                if state in marg.columns:
                    marg.loc[clade.name, state] = 1.0
            else:
                raise ValueError(f"Missing marginal state for node {clade.name}")

    parent_of = {}
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            parent_of[child.name] = parent.name

    descendant_tips: dict[str, list[str]] = {}
    for clade in tree.find_clades(order="postorder"):
        descendant_tips[clade.name] = [tip.name for tip in clade.get_terminals()]

    edge_rows = []
    source_rows = []
    for child_name, parent_name in parent_of.items():
        tips = descendant_tips[child_name]
        tip_meta = metadata.loc[[x for x in tips if x in metadata.index]]
        for country in FOCAL:
            country_tips = tip_meta.index[tip_meta["country_iso3"].eq(country)].tolist()
            if not country_tips:
                continue
            p_child = float(marg.at[child_name, country])
            p_parent = float(marg.at[parent_name, country])
            transition_support = (1.0 - p_parent) * p_child
            source_prob = marg.loc[parent_name].drop(labels=[country]).astype(float)
            denom = float(source_prob.sum())
            if denom > 0:
                source_prob = source_prob / denom
                top_source = str(source_prob.idxmax())
                top_source_probability = float(source_prob.max())
            else:
                top_source = "unresolved"
                top_source_probability = np.nan
            dates = pd.to_datetime(metadata.loc[country_tips, "date_lower"], errors="coerce")
            periods = metadata.loc[country_tips, "epidemic_period"]
            edge_id = f"{parent_name}->{child_name}:{country}"
            edge_rows.append(
                {
                    "event_id": edge_id,
                    "parent_node": parent_name,
                    "child_node": child_name,
                    "destination_country": country,
                    "transition_support": transition_support,
                    "parent_country_probability": p_parent,
                    "child_country_probability": p_child,
                    "top_source_state": top_source,
                    "top_source_probability_conditional": top_source_probability,
                    "n_descendant_country_tips": len(country_tips),
                    "n_descendant_prepandemic": int(periods.eq("prepandemic").sum()),
                    "n_descendant_pandemic": int(periods.eq("pandemic").sum()),
                    "n_descendant_resurgence": int(periods.eq("resurgence").sum()),
                    "descendant_date_lower": dates.min(),
                    "descendant_date_upper": pd.to_datetime(
                        metadata.loc[country_tips, "date_upper"], errors="coerce"
                    ).max(),
                }
            )
            for state, probability in source_prob.items():
                source_rows.append(
                    {
                        "event_id": edge_id,
                        "source_state": state,
                        "conditional_source_probability": float(probability),
                    }
                )
    edges = pd.DataFrame(edge_rows)

    # Attribute each focal tip to its strongest country-entry edge for
    # diagnostics. A separate resurgence-reseeding score excludes edges whose
    # descendant country clade contains any pre-2023 sample. This prevents a
    # historical country-entry edge from being relabelled with the first
    # resurgence sampling month and passed to the dynamic model as a new input.
    path_cache = {}
    for tip in tree.get_terminals():
        path_cache[tip.name] = [x.name for x in tree.get_path(tip)]

    attribution = []
    edge_index = edges.set_index(["child_node", "destination_country"])
    for tip_id, row in metadata.loc[metadata["country_iso3"].isin(FOCAL)].iterrows():
        country = row["country_iso3"]
        path_nodes = path_cache[tip_id]
        candidates = []
        post_candidates = []
        persistence_candidates = []
        for node in path_nodes:
            key = (node, country)
            if key not in edge_index.index:
                continue
            item = edge_index.loc[key]
            if isinstance(item, pd.DataFrame):
                item = item.iloc[0]
            candidates.append(item)
            if (
                int(item["n_descendant_prepandemic"]) == 0
                and int(item["n_descendant_pandemic"]) == 0
                and int(item["n_descendant_resurgence"]) > 0
            ):
                post_candidates.append(item)
            else:
                persistence_candidates.append(
                    float(marg.at[node, country])
                )
        best = max(candidates, key=lambda x: float(x["transition_support"]))
        best_post = max(
            post_candidates,
            key=lambda x: float(x["transition_support"]),
            default=None,
        )
        best_post_support = (
            float(best_post["transition_support"]) if best_post is not None else 0.0
        )
        persistence_support = max(persistence_candidates, default=0.0)
        denom = best_post_support + persistence_support
        normalised_post_import = best_post_support / denom if denom > 0 else np.nan
        attribution.append(
            {
                "tree_sample_id": tip_id,
                "country_iso3": country,
                "epidemic_period": row["epidemic_period"],
                "primary_model_lineage_id": row["primary_model_lineage_id"],
                "strongest_event_id": best["event_id"],
                "strongest_transition_support": float(best["transition_support"]),
                "strongest_post_event_id": (
                    best_post["event_id"] if best_post is not None else ""
                ),
                "strongest_post_transition_support": best_post_support,
                "post_reseeding_support": best_post_support,
                "local_persistence_support": persistence_support,
                "normalised_post_import_support": normalised_post_import,
            }
        )
    attribution = pd.DataFrame(attribution)

    # Successful sampled transmission clusters are defined before looking at
    # cases: at least five attributed genomes spanning at least six months.
    event_assignments = attribution[
        attribution["epidemic_period"].eq("resurgence")
        & attribution["strongest_post_event_id"].ne("")
        & attribution["strongest_post_transition_support"].ge(
            args.transition_threshold
        )
    ].copy()
    assigned = (
        event_assignments.groupby("strongest_post_event_id")
        .agg(
            n_attributed_tips=("tree_sample_id", "size"),
            n_resurgence=("epidemic_period", lambda x: int((x == "resurgence").sum())),
            n_lineages=("primary_model_lineage_id", "nunique"),
        )
        .reset_index()
        .rename(columns={"strongest_post_event_id": "event_id"})
    )
    events = edges.merge(assigned, on="event_id", how="left")
    events[["n_attributed_tips", "n_resurgence", "n_lineages"]] = events[
        ["n_attributed_tips", "n_resurgence", "n_lineages"]
    ].fillna(0).astype(int)
    events["sample_span_days"] = (
        pd.to_datetime(events["descendant_date_upper"])
        - pd.to_datetime(events["descendant_date_lower"])
    ).dt.days
    events["high_support_introduction"] = events["transition_support"].ge(
        args.transition_threshold
    )
    events["post_resurgence_candidate"] = (
        events["n_descendant_prepandemic"].eq(0)
        & events["n_descendant_pandemic"].eq(0)
        & events["n_descendant_resurgence"].gt(0)
    )
    events["high_support_post_reseeding"] = (
        events["post_resurgence_candidate"]
        & events["high_support_introduction"]
    )
    events["successful_sampled_cluster"] = (
        events["high_support_post_reseeding"]
        & events["n_attributed_tips"].ge(args.success_size)
        & events["sample_span_days"].ge(args.success_span_days)
    )

    post = attribution[attribution["epidemic_period"].eq("resurgence")].copy()
    summary = (
        post.groupby(["country_iso3", "primary_model_lineage_id"], dropna=False)
        .agg(
            n_resurgence_tips=("tree_sample_id", "size"),
            mean_post_import_support=("normalised_post_import_support", "mean"),
            mean_local_persistence_support=("local_persistence_support", "mean"),
            mean_raw_reseeding_support=("post_reseeding_support", "mean"),
        )
        .reset_index()
    )
    country_summary = (
        post.groupby("country_iso3")
        .agg(
            n_resurgence_tips=("tree_sample_id", "size"),
            mean_post_import_support=("normalised_post_import_support", "mean"),
            mean_local_persistence_support=("local_persistence_support", "mean"),
        )
        .reset_index()
    )

    date_cols = ["descendant_date_lower", "descendant_date_upper"]
    for col in date_cols:
        events[col] = pd.to_datetime(events[col]).dt.strftime("%Y-%m-%d")
    edges.to_csv(args.output_dir / "all_country_entry_edges.tsv", sep="\t", index=False)
    events.to_csv(args.output_dir / "introduction_events.tsv", sep="\t", index=False)
    pd.DataFrame(source_rows).to_csv(
        args.output_dir / "introduction_source_probabilities.tsv", sep="\t", index=False
    )
    attribution.to_csv(
        args.output_dir / "tip_persistence_reseeding_support.tsv", sep="\t", index=False
    )
    summary.to_csv(
        args.output_dir / "country_lineage_phylogeographic_summary.tsv",
        sep="\t",
        index=False,
    )
    country_summary.to_csv(
        args.output_dir / "country_phylogeographic_summary.tsv", sep="\t", index=False
    )

    validation = {
        "n_tree_tips": len(tree.get_terminals()),
        "n_focal_tip_attributions": len(attribution),
        "transition_threshold": args.transition_threshold,
        "successful_cluster_minimum_genomes": args.success_size,
        "successful_cluster_minimum_span_days": args.success_span_days,
        "n_high_support_introduction_edges": int(events["high_support_introduction"].sum()),
        "n_high_support_post_reseeding_edges": int(
            events["high_support_post_reseeding"].sum()
        ),
        "n_successful_sampled_clusters": int(events["successful_sampled_cluster"].sum()),
        "uncertainty_semantics": (
            "PastML marginal-state support propagated through a modular "
            "phylogeographic layer; not a full joint phylogenetic posterior"
        ),
    }
    (args.output_dir / "phylogeography_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
