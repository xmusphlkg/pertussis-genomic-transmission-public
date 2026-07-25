#!/usr/bin/env python3
"""Create a stable SHA-256 manifest for repository files."""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "provenance/FILE_MANIFEST_SHA256.tsv"
IGNORED_PARTS = {
    ".git",
    ".pytest_cache",
    ".venv",
    ".venv-tests",
    "__pycache__",
}


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


rows = []
for path in sorted(ROOT.rglob("*")):
    if (
        not path.is_file()
        or IGNORED_PARTS.intersection(path.parts)
        or path.suffix in {".pyc", ".pyo"}
        or path == OUTPUT
    ):
        continue
    rows.append((path.relative_to(ROOT).as_posix(), path.stat().st_size, digest(path)))

with OUTPUT.open("w", encoding="utf-8") as handle:
    handle.write("relative_path\tbytes\tsha256\n")
    for row in rows:
        handle.write("\t".join(map(str, row)) + "\n")
