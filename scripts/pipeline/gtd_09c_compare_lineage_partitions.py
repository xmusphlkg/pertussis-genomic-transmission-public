#!/usr/bin/env python3
"""Compare primary and alternate pre-effect hierBAPS partitions."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from math import comb
from pathlib import Path


def read_partition(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {r["tree_sample_id"]: r["model_sublineage_id"] for r in rows}


def adjusted_rand_index(pairs: list[tuple[str, str]]) -> float:
    table: Counter[tuple[str, str]] = Counter(pairs)
    primary = Counter(a for a, _ in pairs)
    alternate = Counter(b for _, b in pairs)
    n = len(pairs)
    if n < 2:
        return 1.0
    nij = sum(comb(x, 2) for x in table.values())
    ai = sum(comb(x, 2) for x in primary.values())
    bj = sum(comb(x, 2) for x in alternate.values())
    total = comb(n, 2)
    expected = ai * bj / total
    maximum = 0.5 * (ai + bj)
    return 1.0 if maximum == expected else (nij - expected) / (maximum - expected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary", type=Path)
    parser.add_argument("alternate", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    primary = read_partition(args.primary)
    alternate = read_partition(args.alternate)
    if primary.keys() != alternate.keys():
        raise ValueError("Primary and alternate partitions do not contain identical tips")

    pairs = [(primary[k], alternate[k]) for k in primary]
    primary_members: dict[str, Counter[str]] = defaultdict(Counter)
    alternate_sizes = Counter(alternate.values())
    contingency = Counter(pairs)
    for p, a in pairs:
        primary_members[p][a] += 1

    rows = []
    for cluster, candidates in primary_members.items():
        best, overlap = candidates.most_common(1)[0]
        n_primary = sum(candidates.values())
        n_alternate = alternate_sizes[best]
        rows.append(
            {
                "primary_cluster": cluster,
                "best_alternate_cluster": best,
                "overlap_n": overlap,
                "primary_n": n_primary,
                "alternate_n": n_alternate,
                "primary_member_recovery": overlap / n_primary,
                "jaccard": overlap / (n_primary + n_alternate - overlap),
            }
        )
    rows.sort(key=lambda r: (-r["primary_n"], r["primary_cluster"]))
    fields = list(rows[0])
    with (args.output_dir / "lineage_partition_sensitivity.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    alt_labels = sorted(set(alternate.values()))
    with (args.output_dir / "lineage_partition_contingency.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["primary_cluster", *alt_labels])
        for p in sorted(set(primary.values())):
            writer.writerow([p, *(contingency[(p, a)] for a in alt_labels)])

    primary_model = {"L1_01.02", "L1_02.05", "L1_02.06", "L1_02.07"}
    key_rows = [r for r in rows if r["primary_cluster"] in primary_model]
    growth = next(r for r in key_rows if r["primary_cluster"] == "L1_02.07")
    result = {
        "n_tips": len(pairs),
        "adjusted_rand_index_level2": adjusted_rand_index(pairs),
        "minimum_primary_member_recovery_among_model_lineages": min(
            r["primary_member_recovery"] for r in key_rows
        ),
        "minimum_jaccard_among_model_lineages": min(
            r["jaccard"] for r in key_rows
        ),
        "high_growth_lineage_L1_02_07_member_recovery": growth[
            "primary_member_recovery"
        ],
        "high_growth_lineage_L1_02_07_jaccard": growth["jaccard"],
    }
    (args.output_dir / "lineage_partition_sensitivity.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
