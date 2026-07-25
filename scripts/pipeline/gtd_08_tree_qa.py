#!/usr/bin/env python3
"""Audit final tree coverage and temporal signal without asserting a molecular clock."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import linregress


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tree",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/results/primary/gtd_primary_core_snp.treefile",
    )
    parser.add_argument(
        "--qc",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/qc/uniform_sequence_qc.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "analysis/genomic_transmission_dynamics/phylogeny/results/primary",
    )
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260724)
    return parser.parse_args()


def decimal_date(frame: pd.DataFrame) -> pd.Series:
    lower = pd.to_datetime(frame["date_lower"], errors="coerce")
    upper = pd.to_datetime(frame["date_upper"], errors="coerce")
    midpoint = lower + (upper - lower) / 2
    start = pd.to_datetime(
        midpoint.dt.year.astype("Int64").astype("string") + "-01-01",
        errors="coerce",
    )
    end = pd.to_datetime(
        (midpoint.dt.year + 1).astype("Int64").astype("string") + "-01-01",
        errors="coerce",
    )
    return midpoint.dt.year + (midpoint - start).dt.total_seconds() / (
        end - start
    ).dt.total_seconds()


def regression_payload(frame: pd.DataFrame, rng: np.random.Generator, permutations: int) -> dict:
    dated = frame.dropna(subset=["decimal_date", "root_to_tip"]).copy()
    if len(dated) < 10 or dated["decimal_date"].max() - dated["decimal_date"].min() < 3:
        return {"n": int(len(dated)), "status": "INSUFFICIENT"}
    fit = linregress(dated["decimal_date"], dated["root_to_tip"])
    permuted_slopes = []
    dates = dated["decimal_date"].to_numpy(copy=True)
    distances = dated["root_to_tip"].to_numpy()
    for _ in range(permutations):
        permuted_slopes.append(linregress(rng.permutation(dates), distances).slope)
    slope_p = (1 + sum(value >= fit.slope for value in permuted_slopes)) / (
        permutations + 1
    )
    return {
        "n": int(len(dated)),
        "status": "ESTIMATED",
        "year_span": round(float(dates.max() - dates.min()), 3),
        "slope_substitutions_per_site_per_year": float(fit.slope),
        "intercept": float(fit.intercept),
        "r_squared": float(fit.rvalue**2),
        "parametric_slope_p": float(fit.pvalue),
        "date_randomisation_one_sided_p": float(slope_p),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tree = Phylo.read(args.tree, "newick")
    tree.root_at_midpoint()
    root_distances = {
        terminal.name: float(tree.distance(tree.root, terminal))
        for terminal in tree.get_terminals()
    }
    qc = pd.read_csv(args.qc, sep="\t", dtype=str).fillna("")
    tips = qc[qc["tree_sample_id"].isin(root_distances)].copy()
    tips["root_to_tip"] = tips["tree_sample_id"].map(root_distances)
    tips["decimal_date"] = decimal_date(tips)
    tips["epidemic_period"] = np.select(
        [
            pd.to_numeric(tips["year"], errors="coerce").le(2019),
            pd.to_numeric(tips["year"], errors="coerce").between(2020, 2022),
            pd.to_numeric(tips["year"], errors="coerce").ge(2023),
        ],
        ["prepandemic", "pandemic", "resurgence"],
        default="unknown",
    )
    tips.to_csv(args.output_dir / "tree_tip_metadata.tsv", sep="\t", index=False)

    rng = np.random.default_rng(args.seed)
    regressions = {"all": regression_payload(tips, rng, args.permutations)}
    for country, group in tips.groupby("country_iso3"):
        if len(group) >= 20:
            regressions[f"country:{country}"] = regression_payload(
                group, rng, args.permutations
            )
    coverage = {
        f"{role}:{country}:{period}": int(value)
        for (role, country, period), value in tips.groupby(
            ["tree_role", "country_iso3", "epidemic_period"]
        ).size().items()
    }
    report = {
        "n_tree_tips": int(len(tree.get_terminals())),
        "n_metadata_matched": int(len(tips)),
        "n_dated_tips": int(tips["decimal_date"].notna().sum()),
        "tree_coverage": coverage,
        "root_to_tip_regressions": regressions,
        "interpretation_rule": (
            "Use calendar-time tMRCA only if temporal signal remains positive under "
            "date randomisation and clock-model sensitivity; otherwise use time-stratified "
            "SNP clusters and interval-valued introduction times."
        ),
    }
    (args.output_dir / "tree_qa_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# 核心 SNP 传播树 QA",
        "",
        f"- 最终树 tips：{report['n_tree_tips']}",
        f"- 成功匹配元数据：{report['n_metadata_matched']}",
        f"- 有日期区间的 tips：{report['n_dated_tips']}",
        "",
        "## Root-to-tip 与日期随机化",
        "",
        "| 分层 | n | 年跨度 | 斜率 | R² | 日期随机化 p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for label, payload in regressions.items():
        lines.append(
            f"| {label} | {payload.get('n', '')} | {payload.get('year_span', '')} | "
            f"{payload.get('slope_substitutions_per_site_per_year', '')} | "
            f"{payload.get('r_squared', '')} | "
            f"{payload.get('date_randomisation_one_sided_p', '')} |"
        )
    lines.extend(
        [
            "",
            "时间信号结论必须结合日期随机化及后续严格钟/松弛钟敏感性；"
            "若不通过，不报告精确 tMRCA，改用时间分层 SNP 簇和区间输入概率。",
            "",
        ]
    )
    (args.output_dir / "TREE_QA_REPORT.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
