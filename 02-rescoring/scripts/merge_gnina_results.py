#!/usr/bin/env python3
"""
Parse per-ligand GNINA text outputs into a single CSV.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


NUMERIC_LINE = re.compile(r"^\s*\d+\s+[-+]?\d")


def parse_gnina_file(path: Path) -> dict[str, float | str]:
    conformer_id = path.stem.replace("_gnina", "")
    with path.open("r") as handle:
        for line in handle:
            if not NUMERIC_LINE.match(line):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            return {
                "conformer_id": conformer_id,
                "gnina_affinity": float(parts[1]),
                "gnina_cnnscore": float(parts[2]),
                "gnina_cnnaffinity": float(parts[3]),
            }
    raise ValueError(f"No GNINA score row found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge GNINA per-ligand outputs into a CSV.")
    parser.add_argument("--input-dir", required=True, help="Directory containing *_gnina.txt files")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = [parse_gnina_file(path) for path in sorted(input_dir.glob("*_gnina.txt"))]
    if not rows:
        raise FileNotFoundError(f"No *_gnina.txt files found in {input_dir}")

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["conformer_id", "gnina_affinity", "gnina_cnnscore", "gnina_cnnaffinity"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} GNINA rows to {output_path}")


if __name__ == "__main__":
    main()
