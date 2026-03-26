#!/usr/bin/env python3
"""
Merge Uni-Dock shortlist, GNINA scores, and RF-Score results into one ranking table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def normalise_shortlist(df: pd.DataFrame) -> pd.DataFrame:
    renamed = df.rename(
        columns={
            "zinc_id": "conformer_id",
            "score": "unidock_affinity",
            "rank": "unidock_rank",
        }
    ).copy()
    if "conformer_id" not in renamed.columns:
        raise ValueError("Shortlist CSV must contain zinc_id or conformer_id")
    if "unidock_affinity" not in renamed.columns:
        raise ValueError("Shortlist CSV must contain unidock_affinity or score")
    if "unidock_rank" not in renamed.columns:
        renamed["unidock_rank"] = renamed["unidock_affinity"].rank(method="min", ascending=True)
    return renamed[["conformer_id", "unidock_affinity", "unidock_rank"]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Integrate Uni-Dock, GNINA, and RF-Score results.")
    parser.add_argument("--shortlist-csv", required=True, help="Shortlist CSV from stage 01")
    parser.add_argument("--gnina-csv", required=True, help="Merged GNINA CSV")
    parser.add_argument("--rfscore-csv", required=True, help="Merged RF-Score CSV")
    parser.add_argument("--output", required=True, help="Output integrated CSV")
    args = parser.parse_args()

    shortlist = normalise_shortlist(pd.read_csv(args.shortlist_csv))
    gnina = pd.read_csv(args.gnina_csv)
    rfscore = pd.read_csv(args.rfscore_csv)

    merged = shortlist.merge(gnina, on="conformer_id", how="inner").merge(
        rfscore, on="conformer_id", how="inner"
    )
    if merged.empty:
        raise ValueError("No overlapping conformers across shortlist, GNINA, and RF-Score inputs")

    merged["rfscore_rank"] = merged["rfscore_v3"].rank(method="min", ascending=False)
    merged["gnina_rank"] = merged["gnina_cnnscore"].rank(method="min", ascending=False)
    merged["composite_rank"] = merged[["unidock_rank", "rfscore_rank", "gnina_rank"]].mean(axis=1)

    output = merged[
        [
            "conformer_id",
            "unidock_affinity",
            "rfscore_v3",
            "gnina_cnnscore",
            "gnina_cnnaffinity",
            "unidock_rank",
            "rfscore_rank",
            "gnina_rank",
            "composite_rank",
        ]
    ].sort_values(["composite_rank", "conformer_id"], ascending=[True, True])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    print(f"Wrote integrated scores to {output_path}")


if __name__ == "__main__":
    main()
