#!/usr/bin/env python3
"""Compare Pasteur cgMLST allelic structure with the frozen core-SNP analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--phylip-map", type=Path, required=True)
    parser.add_argument("--cg-tree", type=Path, required=True)
    parser.add_argument("--lineages", type=Path, required=True)
    parser.add_argument("--core-distances", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    profiles = pd.read_csv(args.profiles, sep="\t", dtype=str).fillna("")
    ids = pd.read_csv(args.phylip_map, sep="\t", dtype=str)
    mapping = dict(zip(ids["phylip_id"], ids["tree_sample_id"]))
    named_tree = args.cg_tree.read_text(encoding="utf-8")
    for phylip_id in sorted(mapping, key=len, reverse=True):
        named_tree = named_tree.replace(phylip_id, mapping[phylip_id])
    (args.output_dir / "cgmlst_rapidnj_named.tree").write_text(
        named_tree, encoding="utf-8"
    )

    sample_names = profiles["tree_sample_id"].tolist()
    allele = profiles.drop(columns="tree_sample_id").to_numpy(dtype=str)
    profile_index = {sample: index for index, sample in enumerate(sample_names)}
    core = pd.read_csv(args.core_distances, sep="\t")
    common = core[
        core["Sample1"].isin(profile_index) & core["Sample2"].isin(profile_index)
    ].copy()
    cg_fraction: list[float] = []
    cg_mismatches: list[int] = []
    jointly_called: list[int] = []
    for row in common.itertuples(index=False):
        left = allele[profile_index[row.Sample1]]
        right = allele[profile_index[row.Sample2]]
        called = (left != "") & (right != "")
        jointly_called.append(int(called.sum()))
        mismatches = int((left[called] != right[called]).sum())
        cg_mismatches.append(mismatches)
        cg_fraction.append(mismatches / called.sum())
    common["cgmlst_jointly_called_loci"] = jointly_called
    common["cgmlst_allelic_mismatches"] = cg_mismatches
    common["cgmlst_mismatch_fraction"] = cg_fraction
    common.to_csv(
        args.output_dir / "cgmlst_core_snp_pairwise_comparison.tsv",
        sep="\t",
        index=False,
    )
    rho, pvalue = spearmanr(common["Distance"], common["cgmlst_allelic_mismatches"])

    lineages = pd.read_csv(args.lineages, sep="\t", dtype=str).set_index(
        "tree_sample_id"
    )
    eligible = [sample for sample in sample_names if sample in lineages.index]
    eligible_idx = np.array([profile_index[sample] for sample in eligible])
    eligible_profiles = allele[eligible_idx]
    eligible_lineages = lineages.loc[eligible, "primary_model_lineage_id"].to_numpy()
    nearest_agreement: list[dict[str, object]] = []
    for i, sample in enumerate(eligible):
        called = (eligible_profiles != "") & (eligible_profiles[i] != "")
        denominator = called.sum(axis=1)
        mismatch = (
            ((eligible_profiles != eligible_profiles[i]) & called).sum(axis=1)
            / np.maximum(denominator, 1)
        )
        mismatch[i] = np.inf
        neighbours = np.argsort(mismatch)[:5]
        same = eligible_lineages[neighbours] == eligible_lineages[i]
        nearest_agreement.append(
            {
                "tree_sample_id": sample,
                "primary_model_lineage_id": eligible_lineages[i],
                "nearest_cgmlst_same_lineage": bool(same[0]),
                "five_neighbour_same_lineage_fraction": float(same.mean()),
            }
        )
    neighbour_table = pd.DataFrame(nearest_agreement)
    neighbour_table.to_csv(
        args.output_dir / "cgmlst_nearest_neighbour_lineage_concordance.tsv",
        sep="\t",
        index=False,
    )
    formal = neighbour_table[
        neighbour_table["primary_model_lineage_id"].ne("Other")
    ]
    summary = {
        "n_cgmlst_profiles": len(profiles),
        "n_profiles_with_frozen_core_snp_lineage": len(eligible),
        "n_pairwise_core_snp_comparisons": len(common),
        "core_snp_cgmlst_spearman_rho": float(rho),
        "core_snp_cgmlst_spearman_p": float(pvalue),
        "nearest_neighbour_lineage_agreement_all": float(
            neighbour_table["nearest_cgmlst_same_lineage"].mean()
        ),
        "nearest_neighbour_lineage_agreement_formal_lineages": float(
            formal["nearest_cgmlst_same_lineage"].mean()
        ),
        "mean_five_neighbour_same_lineage_fraction_formal_lineages": float(
            formal["five_neighbour_same_lineage_fraction"].mean()
        ),
        "boundary": (
            "Assembly-available focal subset; cgMLST validates cluster concordance "
            "and is not used for temporal or branch-length inference."
        ),
    }
    (args.output_dir / "cgmlst_core_snp_concordance.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
