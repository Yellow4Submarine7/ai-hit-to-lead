#!/usr/bin/env python3
"""
Generate teacher predictions and `.types` files for student training.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import joblib
import pandas as pd
from sklearn.model_selection import train_test_split


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_ROOT))

from ai_hit_to_lead_active_learning.features import extract_features


def parse_k_values(raw: str) -> list[int]:
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


def write_types_file(output_path: Path, rows: pd.DataFrame, receptor_path: str, score_column: str, threshold: float) -> None:
    with output_path.open("w") as handle:
        for _, row in rows.iterrows():
            affinity = float(row[score_column])
            label = 1 if affinity >= threshold else 0
            handle.write(f"{label} {affinity:.4f} {receptor_path} {row['pdbqt_path']}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate pseudo labels and .types files.")
    parser.add_argument("--teacher-checkpoint", required=True, help="Teacher checkpoint from train_teacher.py")
    parser.add_argument("--experimental-csv", required=True, help="Experimental labels CSV")
    parser.add_argument("--candidates-csv", required=True, help="Candidate CSV to pseudo-label")
    parser.add_argument("--receptor-path", required=True, help="Receptor path to place in generated .types files")
    parser.add_argument("--output-dir", required=True, help="Directory for predictions and .types files")
    parser.add_argument("--k-values", default="1000,2000,4000", help="Comma-separated pseudo-label sizes")
    args = parser.parse_args()

    teacher_bundle = joblib.load(args.teacher_checkpoint)
    teacher = teacher_bundle["model"]
    scaler = teacher_bundle["scaler"]
    seed = int(teacher_bundle.get("split_seed", 42))
    test_size = float(teacher_bundle.get("test_size", 0.2))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    experimental_df = pd.read_csv(args.experimental_csv)
    candidates_df = pd.read_csv(args.candidates_csv)

    train_df, val_df = split_dataframe(experimental_df, test_size=test_size, seed=seed)
    threshold = float(train_df["binding_level_nm"].median())

    candidate_features, _ = extract_features(candidates_df["smiles_canonical"].tolist(), scaler=scaler)
    candidate_preds = teacher.predict(candidate_features)
    candidate_predictions = candidates_df.copy()
    candidate_predictions["teacher_pred"] = candidate_preds
    candidate_predictions = candidate_predictions.sort_values("teacher_pred", ascending=True).reset_index(drop=True)
    candidate_predictions.to_csv(output_dir / "teacher_predictions.csv", index=False)

    val_features, _ = extract_features(val_df["smiles_canonical"].tolist(), scaler=scaler)
    val_predictions = val_df.copy()
    val_predictions["teacher_pred"] = teacher.predict(val_features)
    val_predictions.to_csv(output_dir / "val_teacher_preds.csv", index=False)

    train_types = train_df.copy()
    train_types["score_for_types"] = train_types["binding_level_nm"]
    val_types = val_df.copy()
    val_types["score_for_types"] = val_types["binding_level_nm"]
    write_types_file(output_dir / "val.types", val_types, args.receptor_path, "score_for_types", threshold)

    for k in parse_k_values(args.k_values):
        pseudo_subset = candidate_predictions.head(min(k, len(candidate_predictions))).copy()
        pseudo_subset["score_for_types"] = pseudo_subset["teacher_pred"]
        combined = pd.concat([train_types, pseudo_subset], ignore_index=True)
        write_types_file(output_dir / f"train_K{k}.types", combined, args.receptor_path, "score_for_types", threshold)

    print(f"Wrote pseudo-label outputs to {output_dir}")


if __name__ == "__main__":
    main()
