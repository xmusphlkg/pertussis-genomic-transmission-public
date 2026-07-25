#!/usr/bin/env python3
"""Stream a fixed paired-read prefix from ENA and record FASTQ QC metrics.

The fixed-pair primary cap removes extreme depth differences while avoiding
full downloads of multi-gigabyte isolate runs. Samples failing SKA coverage
can later be rescued with a larger cap.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import os
import time
from pathlib import Path

import pandas as pd
import requests


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stream_one_read_file(
    url: str,
    output_path: Path,
    max_records: int,
) -> dict[str, object]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial_path = output_path.with_suffix(output_path.suffix + ".partial")
    records = 0
    bases = 0
    q30_bases = 0
    started = time.time()
    with requests.get(url, stream=True, timeout=(30, 300)) as response:
        response.raise_for_status()
        response.raw.decode_content = False
        with gzip.GzipFile(fileobj=response.raw, mode="rb") as input_handle:
            with gzip.open(partial_path, "wb", compresslevel=3) as output_handle:
                while records < max_records:
                    header = input_handle.readline()
                    if not header:
                        break
                    sequence = input_handle.readline()
                    plus = input_handle.readline()
                    quality = input_handle.readline()
                    if not sequence or not plus or not quality:
                        raise ValueError(f"truncated FASTQ record from {url}")
                    if not header.startswith(b"@") or not plus.startswith(b"+"):
                        raise ValueError(f"invalid FASTQ structure from {url}")
                    sequence_stripped = sequence.rstrip(b"\r\n")
                    quality_stripped = quality.rstrip(b"\r\n")
                    if len(sequence_stripped) != len(quality_stripped):
                        raise ValueError(f"sequence/quality length mismatch from {url}")
                    output_handle.write(header)
                    output_handle.write(sequence)
                    output_handle.write(plus)
                    output_handle.write(quality)
                    records += 1
                    bases += len(sequence_stripped)
                    q30_bases += sum(value >= 63 for value in quality_stripped)
    os.replace(partial_path, output_path)
    return {
        "records": records,
        "bases": bases,
        "q30_fraction": q30_bases / bases if bases else 0.0,
        "output_bytes": output_path.stat().st_size,
        "sha256": sha256(output_path),
        "elapsed_seconds": round(time.time() - started, 2),
    }


def process_run(row: dict[str, str], output_root: Path, max_records: int) -> dict[str, object]:
    run = row["run_accession"]
    urls = [value for value in row["fastq_ftp"].split(";") if value]
    paired_urls = [
        value for value in urls if value.endswith("_1.fastq.gz") or value.endswith("_2.fastq.gz")
    ]
    if len(paired_urls) == 2:
        urls = sorted(paired_urls, key=lambda value: 1 if value.endswith("_1.fastq.gz") else 2)
    if len(urls) != 2:
        return {
            "run_accession": run,
            "status": "FAIL_NOT_TWO_PAIRED_FILES",
            "error": f"observed {len(urls)} candidate files",
        }
    outputs = [
        output_root / f"{run}_1.subset{max_records}.fastq.gz",
        output_root / f"{run}_2.subset{max_records}.fastq.gz",
    ]
    done_path = output_root / f"{run}.subset{max_records}.done.json"
    if done_path.exists() and all(path.exists() for path in outputs):
        done = json.loads(done_path.read_text(encoding="utf-8"))
        return done
    errors = []
    for attempt in range(1, 4):
        try:
            metrics = []
            for url, output in zip(urls, outputs):
                # A small subset of ENA paths returns an HTML proxy page over HTTPS
                # while the identical archive object is valid over HTTP. Preserve
                # HTTPS for the first two attempts and use the archive fallback only
                # on the final attempt; output hashes still make the result auditable.
                scheme = "https" if attempt < 3 else "http"
                metrics.append(
                    stream_one_read_file(f"{scheme}://{url}", output, max_records)
                )
            records = [int(metric["records"]) for metric in metrics]
            status = (
                "PASS"
                if records[0] == records[1] and records[0] >= 50_000
                else "FAIL_READ_COUNT"
            )
            result = {
                "run_accession": run,
                "status": status,
                "r1_path": str(outputs[0]),
                "r2_path": str(outputs[1]),
                "paired_records": min(records),
                "r1_bases": metrics[0]["bases"],
                "r2_bases": metrics[1]["bases"],
                "r1_q30_fraction": round(float(metrics[0]["q30_fraction"]), 6),
                "r2_q30_fraction": round(float(metrics[1]["q30_fraction"]), 6),
                "r1_output_bytes": metrics[0]["output_bytes"],
                "r2_output_bytes": metrics[1]["output_bytes"],
                "r1_sha256": metrics[0]["sha256"],
                "r2_sha256": metrics[1]["sha256"],
                "elapsed_seconds": round(
                    float(metrics[0]["elapsed_seconds"])
                    + float(metrics[1]["elapsed_seconds"]),
                    2,
                ),
                "attempts": attempt,
                "error": "",
            }
            done_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
            return result
        except Exception as error:  # noqa: BLE001
            errors.append(repr(error))
            for output in outputs:
                partial = output.with_suffix(output.suffix + ".partial")
                partial.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(2**attempt)
    return {
        "run_accession": run,
        "status": "FAIL_DOWNLOAD_OR_FASTQ",
        "attempts": 3,
        "error": " | ".join(errors)[-2000:],
    }


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/acquisition/primary_phylogeny_manifest_resolved.tsv",
    )
    parser.add_argument(
        "--ena-resolution",
        type=Path,
        default=root
        / "analysis/genomic_transmission_dynamics/phylogeny/acquisition/ena_run_resolution.tsv",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=root
        / "pertussis_data/pertussis_gene/genomic_transmission_dynamics/reads/subset",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=root / "analysis/genomic_transmission_dynamics/phylogeny/qc",
    )
    parser.add_argument("--max-records", type=int, default=500_000)
    parser.add_argument("--workers", type=int, default=16)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    ena = pd.read_csv(args.ena_resolution, sep="\t", dtype=str).fillna("")
    required = manifest.loc[
        manifest["sequence_acquisition"].eq("ena_fastq_download"),
        ["tree_sample_id", "run_accession"],
    ].merge(ena[["run_accession", "fastq_ftp"]], on="run_accession", how="left")
    rows = required.to_dict("records")
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_map = {
            executor.submit(process_run, row, args.output_root, args.max_records): row
            for row in rows
        }
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            print(
                f"{len(results)}/{len(rows)} {result.get('run_accession', '')} "
                f"{result.get('status', '')}",
                flush=True,
            )
    qc = pd.DataFrame(results).sort_values("run_accession")
    qc_path = args.report_dir / f"fastq_subset_{args.max_records}_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)

    manifest = manifest.merge(
        qc[["run_accession", "status", "r1_path", "r2_path", "paired_records"]],
        on="run_accession",
        how="left",
        suffixes=("", "_subset"),
    ).fillna("")
    is_fastq = manifest["sequence_acquisition"].eq("ena_fastq_download")
    passed = is_fastq & manifest["status"].eq("PASS")
    manifest.loc[passed, "fastq_r1"] = manifest.loc[passed, "r1_path"]
    manifest.loc[passed, "fastq_r2"] = manifest.loc[passed, "r2_path"]
    manifest.loc[passed, "sequence_input_path"] = (
        manifest.loc[passed, "r1_path"] + ";" + manifest.loc[passed, "r2_path"]
    )
    manifest.loc[passed, "source_resolution_status"] = "resolved_streamed_fastq_subset"
    manifest.loc[is_fastq & ~passed, "source_resolution_status"] = "fastq_subset_failed"
    output_manifest = args.report_dir / "primary_phylogeny_manifest_with_fastq.tsv"
    manifest.to_csv(output_manifest, sep="\t", index=False)
    report = {
        "n_required_runs": len(rows),
        "n_pass": int(qc["status"].eq("PASS").sum()),
        "n_fail": int((~qc["status"].eq("PASS")).sum()),
        "max_records_per_read_file": args.max_records,
        "total_subset_gib": round(
            float(
                pd.to_numeric(qc.get("r1_output_bytes", 0), errors="coerce").fillna(0).sum()
                + pd.to_numeric(qc.get("r2_output_bytes", 0), errors="coerce").fillna(0).sum()
            )
            / (1024**3),
            2,
        ),
    }
    (args.report_dir / f"fastq_subset_{args.max_records}_report.json").write_text(
        json.dumps(report, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
