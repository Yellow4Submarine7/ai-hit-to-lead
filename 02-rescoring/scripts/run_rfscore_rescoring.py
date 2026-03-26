#!/usr/bin/env python3
"""
Run an external RF-Score command template over a shortlist CSV.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run an RF-Score command template over shortlist poses.")
    parser.add_argument("--shortlist-csv", required=True, help="Shortlist CSV from stage 01")
    parser.add_argument("--output-dir", required=True, help="Directory for *.score outputs")
    parser.add_argument("--command-template", required=True, help="Shell template with {ligand}, {output}, {conformer_id}, and optionally {receptor}")
    parser.add_argument("--pose-column", default="pose_path", help="Pose path column in shortlist CSV")
    parser.add_argument("--receptor", default="", help="Optional receptor path for command templates")
    parser.add_argument("--limit", type=int, default=0, help="Optional limit for smoke runs")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(args.shortlist_csv, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader, start=1):
            if args.limit and index > args.limit:
                break
            conformer_id = row.get("zinc_id") or row.get("conformer_id")
            pose_path = row.get(args.pose_column, "")
            if not conformer_id or not pose_path:
                raise ValueError(f"Missing conformer_id/pose_path in row {index}")
            output_path = output_dir / f"{conformer_id}_out.score"
            command = args.command_template.format(
                ligand=pose_path,
                output=str(output_path),
                conformer_id=conformer_id,
                receptor=args.receptor,
            )
            result = subprocess.run(command, shell=True, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                raise RuntimeError(f"RF-Score command failed for {conformer_id}: {result.stderr}")


if __name__ == "__main__":
    main()
