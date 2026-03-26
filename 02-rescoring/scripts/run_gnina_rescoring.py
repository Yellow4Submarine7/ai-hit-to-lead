#!/usr/bin/env python3
"""
Run GNINA `--score_only` over a shortlist CSV.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Run GNINA score-only rescoring on shortlist poses.")
    parser.add_argument("--shortlist-csv", required=True, help="Shortlist CSV from stage 01")
    parser.add_argument("--receptor", required=True, help="Receptor PDBQT file")
    parser.add_argument("--output-dir", required=True, help="Directory for *_gnina.txt outputs")
    parser.add_argument("--gnina-path", default="gnina", help="Path to GNINA binary")
    parser.add_argument("--pose-column", default="pose_path", help="Pose path column in shortlist CSV")
    parser.add_argument("--device", default="0", help="CUDA device id")
    parser.add_argument("--cnn", default="crossdock_default2018", help="GNINA CNN model")
    parser.add_argument("--cnn-scoring", default="rescore", help="GNINA CNN scoring mode")
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

            output_file = output_dir / f"{conformer_id}_gnina.txt"
            command = [
                args.gnina_path,
                "--receptor",
                args.receptor,
                "--ligand",
                pose_path,
                "--score_only",
                "--cnn",
                args.cnn,
                "--cnn_scoring",
                args.cnn_scoring,
                "--device",
                args.device,
            ]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            output_file.write_text((result.stdout or "") + (result.stderr or ""))
            if result.returncode != 0:
                raise RuntimeError(f"GNINA failed for {conformer_id}: {result.stderr}")


if __name__ == "__main__":
    main()
