#!/usr/bin/env python3
"""Call the frozen Pasteur B. pertussis 2038-locus cgMLST scheme via BIGSdb."""

from __future__ import annotations

import argparse
import base64
import concurrent.futures
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

SCHEME_URL = (
    "https://bigsdb.pasteur.fr/api/db/"
    "pubmlst_bordetella_seqdef/schemes/4/sequence"
)
LOCI_URL = (
    "https://bigsdb.pasteur.fr/api/db/"
    "pubmlst_bordetella_seqdef/schemes/4/loci"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--minimum-call-fraction", type=float, default=0.90)
    return parser.parse_args()


def get_loci(output_dir: Path) -> list[str]:
    freeze = output_dir / "pasteur_scheme4_loci.json"
    if freeze.is_file():
        payload = json.loads(freeze.read_text(encoding="utf-8"))
    else:
        response = requests.get(LOCI_URL, timeout=120)
        response.raise_for_status()
        payload = response.json()
        freeze.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return [uri.rsplit("/", 1)[-1] for uri in payload["loci"]]


def call_one(row: dict[str, str], cache_dir: Path) -> dict[str, object]:
    sample = row["tree_sample_id"]
    cache = cache_dir / f"{sample}.json"
    if cache.is_file() and cache.stat().st_size:
        payload = json.loads(cache.read_text(encoding="utf-8"))
        return {"tree_sample_id": sample, "status": "PASS_CACHE", "payload": payload}
    assembly = Path(row["cgmlst_assembly_path"])
    encoded = base64.b64encode(assembly.read_bytes()).decode("ascii")
    errors: list[str] = []
    for attempt in range(1, 4):
        try:
            response = requests.post(
                SCHEME_URL,
                json={"base64": True, "details": False, "sequence": encoded},
                timeout=300,
            )
            response.raise_for_status()
            payload = response.json()
            temporary = cache.with_suffix(".json.partial")
            temporary.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            temporary.replace(cache)
            return {"tree_sample_id": sample, "status": "PASS", "payload": payload}
        except Exception as error:  # noqa: BLE001
            errors.append(repr(error))
            if attempt < 3:
                time.sleep(2**attempt)
    return {
        "tree_sample_id": sample,
        "status": "FAIL_API",
        "payload": {},
        "error": " | ".join(errors)[-1500:],
    }


def phylip_distance(matrix: np.ndarray, names: list[str], output: Path) -> None:
    with output.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(names)}\n")
        for i, name in enumerate(names):
            row: list[str] = []
            for j in range(len(names)):
                both = (matrix[i] != "") & (matrix[j] != "")
                distance = (
                    float(np.mean(matrix[i, both] != matrix[j, both]))
                    if np.any(both)
                    else 1.0
                )
                row.append(f"{distance:.8f}")
            handle.write(f"{name[:10]:<10} {' '.join(row)}\n")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.output_dir / "api_cache"
    cache_dir.mkdir(exist_ok=True)
    loci = get_loci(args.output_dir)
    manifest = pd.read_csv(args.assembly_manifest, sep="\t", dtype=str).fillna("")
    manifest = manifest[manifest["assembly_status"].str.startswith("PASS")].copy()
    rows = manifest.to_dict(orient="records")
    calls: list[dict[str, object]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(call_one, row, cache_dir): row for row in rows}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            calls.append(result)
            print(
                f"{index}/{len(rows)} {result['tree_sample_id']} {result['status']}",
                flush=True,
            )
    call_status = pd.DataFrame(
        {
            "tree_sample_id": result["tree_sample_id"],
            "cgmlst_api_status": result["status"],
            "cgST": result["payload"].get("fields", {}).get("cgST", ""),
            "n_exact_loci": len(result["payload"].get("exact_matches", {})),
            "api_error": result.get("error", ""),
        }
        for result in calls
    )
    call_status.to_csv(args.output_dir / "cgmlst_call_status.tsv", sep="\t", index=False)
    passed = [
        result
        for result in calls
        if result["status"].startswith("PASS")
        and len(result["payload"].get("exact_matches", {}))
        >= args.minimum_call_fraction * len(loci)
    ]
    sample_names = [result["tree_sample_id"] for result in passed]
    profiles = np.full((len(passed), len(loci)), "", dtype=object)
    locus_index = {locus: index for index, locus in enumerate(loci)}
    for i, result in enumerate(passed):
        for locus, matches in result["payload"].get("exact_matches", {}).items():
            if locus in locus_index and matches:
                profiles[i, locus_index[locus]] = str(matches[0]["allele_id"])
    profile_table = pd.DataFrame(profiles, columns=loci)
    profile_table.insert(0, "tree_sample_id", sample_names)
    profile_table.to_csv(args.output_dir / "cgmlst_allele_profiles.tsv", sep="\t", index=False)
    phylip_distance(
        profiles,
        [f"S{i:08d}" for i in range(len(sample_names))],
        args.output_dir / "cgmlst_normalised_allelic_distance.phy",
    )
    pd.DataFrame(
        {"phylip_id": [f"S{i:08d}" for i in range(len(sample_names))],
         "tree_sample_id": sample_names}
    ).to_csv(args.output_dir / "cgmlst_phylip_id_map.tsv", sep="\t", index=False)
    report = {
        "scheme": "Pasteur cgMLST_pertussis",
        "scheme_id": 4,
        "scheme_loci": len(loci),
        "minimum_call_fraction": args.minimum_call_fraction,
        "n_assembly_inputs": len(manifest),
        "n_api_pass": int(call_status["cgmlst_api_status"].str.startswith("PASS").sum()),
        "n_profiles_passing_call_fraction": len(passed),
    }
    (args.output_dir / "cgmlst_call_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
