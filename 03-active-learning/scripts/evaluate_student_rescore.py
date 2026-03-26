#!/usr/bin/env python3
"""
Evaluate student checkpoints on a validation .types file.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

try:  # pragma: no cover
    from gninatorch.inference import inference, options as inference_options
except Exception:  # pragma: no cover
    inference = None  # type: ignore
    inference_options = None  # type: ignore


def load_targets(val_types: Path) -> pd.Series:
    values = []
    for line in val_types.read_text().splitlines():
        if not line.strip():
            continue
        values.append(float(line.split()[1]))
    return pd.Series(values, dtype=float)


def run_inference(types_file: Path, checkpoint: Path, output_dir: Path) -> Path:
    if inference is None or inference_options is None:
        raise ImportError("gninatorch is required to evaluate checkpoints")
    output_dir.mkdir(parents=True, exist_ok=True)
    args = inference_options(
        [
            str(types_file),
            "default2018",
            str(checkpoint),
            "--affinity_pos",
            "1",
            "--label_pos",
            "0",
            "--out_dir",
            str(output_dir),
            "--log_file",
            "val.log",
            "--no_roc_auc",
        ]
    )
    inference(args)
    return output_dir / "val_results.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate student checkpoints on validation types.")
    parser.add_argument("--val-types", required=True, help="Validation .types file")
    parser.add_argument("--checkpoints-root", required=True, help="Root directory containing per-K checkpoint folders")
    parser.add_argument("--k-values", default="0,300,1000,4000", help="Comma-separated K values to evaluate")
    parser.add_argument("--output", required=True, help="Output summary CSV")
    args = parser.parse_args()

    val_types = Path(args.val_types)
    checkpoints_root = Path(args.checkpoints_root)
    y_true = load_targets(val_types)
    rows = []

    for raw_k in [item.strip() for item in args.k_values.split(",") if item.strip()]:
        k = int(raw_k)
        checkpoint_dir = checkpoints_root / f"student_rescore_v4_K{k}"
        candidates = sorted(checkpoint_dir.rglob("*.pt"), key=lambda path: path.stat().st_mtime)
        if not candidates:
            raise FileNotFoundError(f"No checkpoints found in {checkpoint_dir}")
        checkpoint = candidates[-1]
        output_dir = Path(args.output).resolve().parent / f"student_rescore_v4_K{k}_val"
        csv_path = run_inference(val_types, checkpoint, output_dir)
        df = pd.read_csv(csv_path)
        correlation, pvalue = spearmanr(df["affinity_pred"], y_true.to_numpy())
        if correlation != correlation:
            correlation, pvalue = 0.0, 1.0
        rows.append(
            {
                "K": k,
                "val_spearman_signed": float(correlation),
                "val_pvalue": float(pvalue),
                "n_samples": int(len(df)),
                "checkpoint_path": str(checkpoint),
            }
        )

    output = pd.DataFrame(rows).sort_values("K")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    print(f"Wrote evaluation summary to {output_path}")


if __name__ == "__main__":
    main()
