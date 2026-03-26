#!/usr/bin/env python3
"""
Train the public Teacher ensemble with multi-seed scaffold validation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import joblib
import pandas as pd
from scipy.stats import spearmanr
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_ROOT))

from ai_hit_to_lead_active_learning.features import extract_features
from ai_hit_to_lead_active_learning.teacher import TeacherEnsemble


def parse_seeds(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def split_dataframe(df: pd.DataFrame, test_size: float, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    if "scaffold" in df.columns and df["scaffold"].nunique() > 1:
        scaffolds = df["scaffold"].dropna().unique().tolist()
        train_scaffolds, val_scaffolds = train_test_split(scaffolds, test_size=test_size, random_state=seed)
        train_df = df[df["scaffold"].isin(set(train_scaffolds))].copy()
        val_df = df[df["scaffold"].isin(set(val_scaffolds))].copy()
        if not train_df.empty and not val_df.empty:
            return train_df, val_df
    train_df, val_df = train_test_split(df, test_size=test_size, random_state=seed)
    return train_df.copy(), val_df.copy()


def safe_spearman(y_true, y_pred) -> float:
    corr, _ = spearmanr(y_true, y_pred)
    if corr != corr:  # NaN
        return 0.0
    return float(corr)


def main() -> None:
    parser = argparse.ArgumentParser(description="Train the public Teacher ensemble.")
    parser.add_argument("--experimental-csv", required=True, help="Experimental labels CSV")
    parser.add_argument("--checkpoint-output", required=True, help="Output teacher checkpoint")
    parser.add_argument("--metrics-output", required=True, help="Output metrics CSV")
    parser.add_argument("--test-size", type=float, default=0.2, help="Validation fraction")
    parser.add_argument("--seeds", default="42,43,44,45,46", help="Comma-separated random seeds")
    args = parser.parse_args()

    df = pd.read_csv(args.experimental_csv)
    seeds = parse_seeds(args.seeds)
    metrics_rows = []

    final_teacher = None
    final_scaler = None
    final_seed = seeds[0]

    for seed in seeds:
        train_df, val_df = split_dataframe(df, test_size=args.test_size, seed=seed)
        X_train, scaler = extract_features(train_df["smiles_canonical"].tolist(), fit_scaler=True)
        X_val, _ = extract_features(val_df["smiles_canonical"].tolist(), scaler=scaler)

        teacher = TeacherEnsemble().fit(X_train, train_df["binding_level_nm"].to_numpy())
        y_pred = teacher.predict(X_val)
        metrics_rows.append(
            {
                "seed": seed,
                "spearman": safe_spearman(val_df["binding_level_nm"].to_numpy(), y_pred),
                "train_size": len(train_df),
                "val_size": len(val_df),
            }
        )

        if seed == final_seed:
            final_teacher = teacher
            final_scaler = scaler

    metrics_df = pd.DataFrame(metrics_rows)
    metrics_output = Path(args.metrics_output)
    metrics_output.parent.mkdir(parents=True, exist_ok=True)
    metrics_df.to_csv(metrics_output, index=False)

    checkpoint_output = Path(args.checkpoint_output)
    checkpoint_output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": final_teacher,
            "scaler": final_scaler,
            "split_seed": final_seed,
            "test_size": args.test_size,
            "feature_dim": 2056,
        },
        checkpoint_output,
    )
    print(f"Wrote teacher checkpoint to {checkpoint_output}")
    print(f"Wrote metrics to {metrics_output}")


if __name__ == "__main__":
    main()
