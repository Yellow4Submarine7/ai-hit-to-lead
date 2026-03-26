#!/usr/bin/env python3
"""
Export Uni-Dock results into a public shortlist contract.

This script parses `*_out.pdbqt` pose files, merges the scores with metadata from
the extraction stage, ranks the poses by Uni-Dock affinity, and optionally writes
a manifest of shortlisted pose files for downstream rescoring.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict


def parse_pdbqt_score(filepath: Path) -> tuple[str, float | None, str | None]:
    """Extract the Uni-Dock/Vina score from a PDBQT file."""
    try:
        zinc_id = filepath.stem.replace("_out", "")
        with filepath.open("r") as handle:
            for line in handle:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        return zinc_id, float(parts[3]), None
        return zinc_id, None, "No VINA RESULT found"
    except Exception as exc:  # pragma: no cover - defensive parsing
        return filepath.stem, None, str(exc)


def process_batch(filepaths: list[Path]) -> list[tuple[str, float, str]]:
    """Parse a batch of PDBQT files in a worker process."""
    rows: list[tuple[str, float, str]] = []
    for filepath in filepaths:
        zinc_id, score, error = parse_pdbqt_score(filepath)
        if error is None and score is not None:
            rows.append((zinc_id, score, str(filepath)))
    return rows


def _candidate_keys(zinc_id: str) -> list[str]:
    keys = [zinc_id]
    if "." in zinc_id:
        without_variant = ".".join(zinc_id.split(".")[:-1])
        if without_variant and without_variant not in keys:
            keys.append(without_variant)
        base_id = zinc_id.split(".")[0]
        if base_id not in keys:
            keys.append(base_id)
    return keys


def load_metadata_lookup(csv_path: Path) -> Dict[str, dict[str, str]]:
    """Load smiles and source metadata, storing fallback keys for conformer variants."""
    lookup: Dict[str, dict[str, str]] = {}
    with csv_path.open("r", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            zinc_id = row.get("zinc_id") or row.get("conformer_id") or row.get("ligand_id")
            if not zinc_id:
                continue

            payload = {
                "smiles": row.get("smiles")
                or row.get("smiles_canonical")
                or row.get("canonical_smiles")
                or "",
                "subset": row.get("subset", ""),
                "source_pdbqt_path": row.get("pdbqt_path", ""),
            }

            for key in _candidate_keys(zinc_id):
                lookup.setdefault(key, payload)
    return lookup


def resolve_metadata(zinc_id: str, lookup: Dict[str, dict[str, str]]) -> dict[str, str]:
    for key in _candidate_keys(zinc_id):
        if key in lookup:
            return lookup[key]
    return {"smiles": "", "subset": "", "source_pdbqt_path": ""}


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export Uni-Dock results to a shortlist CSV for downstream rescoring."
    )
    parser.add_argument("--results-dir", required=True, help="Directory containing *_out.pdbqt files")
    parser.add_argument("--smiles-csv", required=True, help="CSV containing smiles metadata")
    parser.add_argument(
        "--input-csv",
        help="Optional extraction-stage CSV with subset and pdbqt_path columns. Defaults to --smiles-csv.",
    )
    parser.add_argument("--output", required=True, help="Output shortlist CSV")
    parser.add_argument("--manifest-output", help="Optional output manifest of shortlisted pose paths")
    parser.add_argument("--top-n", type=int, default=0, help="Limit to the top N poses (0 = all)")
    parser.add_argument("--ncpus", type=int, default=os.cpu_count() or 1, help="Number of workers")
    parser.add_argument("--batch-size", type=int, default=10000, help="Batch size for parallel parsing")
    return parser


def main() -> None:
    parser = build_argument_parser()
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    output_path = Path(args.output)
    metadata_path = Path(args.input_csv or args.smiles_csv)

    if not results_dir.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_files = sorted(results_dir.glob("*_out.pdbqt"))
    if not result_files:
        raise FileNotFoundError(f"No *_out.pdbqt files found in {results_dir}")

    print(f"Scanning {len(result_files):,} Uni-Dock poses...", file=sys.stderr)
    metadata_lookup = load_metadata_lookup(metadata_path)
    print(f"Loaded {len(metadata_lookup):,} metadata keys", file=sys.stderr)

    batches = [result_files[i : i + args.batch_size] for i in range(0, len(result_files), args.batch_size)]
    parsed_rows: list[tuple[str, float, str]] = []

    with ProcessPoolExecutor(max_workers=args.ncpus) as executor:
        futures = [executor.submit(process_batch, batch) for batch in batches]
        for future in as_completed(futures):
            parsed_rows.extend(future.result())

    parsed_rows.sort(key=lambda row: row[1])
    if args.top_n > 0:
        parsed_rows = parsed_rows[: args.top_n]

    fieldnames = [
        "rank",
        "zinc_id",
        "unidock_affinity",
        "smiles",
        "subset",
        "pose_path",
        "source_pdbqt_path",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, (zinc_id, score, pose_path) in enumerate(parsed_rows, start=1):
            metadata = resolve_metadata(zinc_id, metadata_lookup)
            writer.writerow(
                {
                    "rank": index,
                    "zinc_id": zinc_id,
                    "unidock_affinity": f"{score:.4f}",
                    "smiles": metadata["smiles"],
                    "subset": metadata["subset"],
                    "pose_path": pose_path,
                    "source_pdbqt_path": metadata["source_pdbqt_path"],
                }
            )

    if args.manifest_output:
        manifest_path = Path(args.manifest_output)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with manifest_path.open("w") as handle:
            for _, _, pose_path in parsed_rows:
                handle.write(f"{pose_path}\n")

    print(f"Wrote shortlist CSV: {output_path}", file=sys.stderr)
    if args.manifest_output:
        print(f"Wrote shortlist manifest: {args.manifest_output}", file=sys.stderr)


if __name__ == "__main__":
    main()
