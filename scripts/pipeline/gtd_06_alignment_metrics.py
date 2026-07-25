#!/usr/bin/env python3
"""Measure per-tip missingness in a FASTA alignment."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def iter_fasta(path: Path):
    name = ""
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name:
            yield name, "".join(chunks)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--max-missing", type=float, default=0.20)
    parser.add_argument("--exclude-list", type=Path)
    args = parser.parse_args()
    rows = []
    for name, sequence in iter_fasta(args.alignment):
        length = len(sequence)
        upper = sequence.upper()
        called = sum(upper.count(base) for base in ("A", "C", "G", "T"))
        missing = length - called
        rows.append(
            {
                "tree_sample_id": name,
                "alignment_length": length,
                "missing_sites": missing,
                "missing_fraction": missing / length if length else 1.0,
            }
        )
    frame = pd.DataFrame(rows)
    frame["alignment_qc_pass"] = frame["missing_fraction"].le(args.max_missing)
    frame.to_csv(args.output, sep="\t", index=False)
    if args.exclude_list:
        excluded = frame.loc[~frame["alignment_qc_pass"], "tree_sample_id"].tolist()
        args.exclude_list.write_text(
            "\n".join(excluded) + ("\n" if excluded else ""),
            encoding="utf-8",
        )
    print(frame["alignment_qc_pass"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
