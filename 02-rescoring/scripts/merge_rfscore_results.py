#!/usr/bin/env python3
"""
Combine per-ligand RF-Score outputs into a single CSV.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_score_file(path: Path) -> tuple[str, float]:
    conformer_id = path.stem.replace("_out", "")
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    for line in lines:
        if "," in line:
            parts = [part.strip() for part in line.split(",")]
        else:
            parts = [line]
        for part in reversed(parts):
            try:
                return conformer_id, float(part)
            except ValueError:
                continue
    raise ValueError(f"No RF-Score value found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge per-ligand RF-Score outputs into a CSV.")
    parser.add_argument("--input-dir", required=True, help="Directory containing *.score files")
    parser.add_argument("--output", required=True, help="Output CSV path")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in sorted(input_dir.glob("*.score")):
        conformer_id, score = parse_score_file(path)
        rows.append({"conformer_id": conformer_id, "rfscore_v3": score})

    if not rows:
        raise FileNotFoundError(f"No .score files found in {input_dir}")

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["conformer_id", "rfscore_v3"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} RF-Score rows to {output_path}")


if __name__ == "__main__":
    main()
