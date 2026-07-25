#!/usr/bin/env python3
"""Compute final-alignment SNP distances only for cgMLST-profile samples."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", type=Path)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    requested = set(
        pd.read_csv(args.profiles, sep="\t", usecols=["tree_sample_id"])[
            "tree_sample_id"
        ]
    )
    records = [
        record for record in SeqIO.parse(args.alignment, "fasta") if record.id in requested
    ]
    if len(records) < 2:
        raise RuntimeError("Fewer than two profile samples occur in the final alignment")
    length = {len(record.seq) for record in records}
    if len(length) != 1:
        raise RuntimeError("Alignment records have unequal lengths")

    names = [record.id for record in records]
    sequence = np.frombuffer(
        "".join(str(record.seq).upper() for record in records).encode("ascii"),
        dtype=np.uint8,
    ).reshape(len(records), length.pop())
    valid = np.isin(sequence, np.frombuffer(b"ACGT", dtype=np.uint8))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "Sample1",
                "Sample2",
                "Distance",
                "Mismatches (proportion)",
                "Match count",
                "Mismatch count",
            ]
        )
        for left in range(len(names) - 1):
            called = valid[left + 1 :] & valid[left]
            jointly_called = called.sum(axis=1)
            mismatches = ((sequence[left + 1 :] != sequence[left]) & called).sum(axis=1)
            for offset, (called_n, mismatch_n) in enumerate(
                zip(jointly_called, mismatches), start=left + 1
            ):
                writer.writerow(
                    [
                        names[left],
                        names[offset],
                        int(mismatch_n),
                        float(mismatch_n / called_n) if called_n else "",
                        int(called_n - mismatch_n),
                        int(mismatch_n),
                    ]
                )

    print(
        f"profiles_requested={len(requested)} profiles_in_final_alignment={len(records)} "
        f"pairs={len(records) * (len(records) - 1) // 2}"
    )


if __name__ == "__main__":
    main()
