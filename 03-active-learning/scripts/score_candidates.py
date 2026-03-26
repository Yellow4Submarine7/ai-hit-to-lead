#!/usr/bin/env python3
"""
Rank candidate molecules from precomputed predictions or gninatorch inference.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

try:  # pragma: no cover - gninatorch is optional in tests
    from gninatorch.inference import inference, options as inference_options
except Exception:  # pragma: no cover
    inference = None  # type: ignore
    inference_options = None  # type: ignore


def load_predictions(predictions_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(predictions_csv)
    if "student_score" in df.columns:
        return df[["ligand_id", "student_score"]]
    if "affinity_pred" in df.columns:
        renamed = df.rename(columns={"affinity_pred": "student_score"})
        return renamed[["ligand_id", "student_score"]]
    raise ValueError("Predictions CSV must contain student_score or affinity_pred")


def run_gninatorch(candidates: pd.DataFrame, receptor_path: Path, checkpoint: Path, output_dir: Path) -> pd.DataFrame:
    if inference is None or inference_options is None:
        raise ImportError("gninatorch is required unless --predictions-csv is supplied")

    types_path = output_dir / "candidates.types"
    with types_path.open("w") as handle:
        for _, row in candidates.iterrows():
            handle.write(f"0 0.0 {receptor_path} {row['pdbqt_path']}\n")

    args = inference_options(
        [
            str(types_path),
            "default2018",
            str(checkpoint),
            "--affinity_pos",
            "1",
            "--label_pos",
            "0",
            "--out_dir",
            str(output_dir),
            "--log_file",
            "candidates.log",
            "--no_roc_auc",
        ]
    )
    inference(args)
    return load_predictions(output_dir / "candidates_results.csv")


def main() -> None:
    parser = argparse.ArgumentParser(description="Rank candidate molecules for active learning.")
    parser.add_argument("--candidates-csv", required=True, help="CSV with ligand_id, smiles_canonical, and pdbqt_path")
    parser.add_argument("--output", required=True, help="Output ranked CSV")
    parser.add_argument("--predictions-csv", help="Precomputed prediction CSV for offline scoring")
    parser.add_argument("--checkpoint", help="Optional gninatorch checkpoint")
    parser.add_argument("--receptor-path", help="Required when using --checkpoint")
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates_csv)
    if args.predictions_csv:
        predictions = load_predictions(Path(args.predictions_csv))
    elif args.checkpoint and args.receptor_path:
        output_dir = Path(args.output).resolve().parent / "student_inference"
        output_dir.mkdir(parents=True, exist_ok=True)
        predictions = run_gninatorch(candidates, Path(args.receptor_path), Path(args.checkpoint), output_dir)
    else:
        raise ValueError("Provide either --predictions-csv or both --checkpoint and --receptor-path")

    merged = candidates.merge(predictions, on="ligand_id", how="inner")
    ranked = merged[["ligand_id", "smiles_canonical", "student_score"]].sort_values(
        "student_score", ascending=True
    )
    ranked.insert(0, "rank", range(1, len(ranked) + 1))

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(output_path, index=False)
    print(f"Wrote ranked candidates to {output_path}")


if __name__ == "__main__":
    main()
