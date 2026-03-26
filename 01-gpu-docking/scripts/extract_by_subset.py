#!/usr/bin/env python3
"""
Extract ligand PDBQT files for a specific subset from a stage-00 CSV.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from extract_ligands import load_csv, extract_from_tgz  # type: ignore


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract only one subset from a stage-00 CSV.")
    parser.add_argument("--input-csv", required=True, help="CSV with pdbqt_path and optional subset columns")
    parser.add_argument("--subset", required=True, help="Subset name to filter, e.g. zinc-22g")
    parser.add_argument("--data-root", required=True, help="Directory containing TGZ archives")
    parser.add_argument("--output-dir", required=True, help="Output directory for extracted ligands")
    parser.add_argument("--ligand-list", help="Optional manifest of extracted ligand paths")
    parser.add_argument("--ncpus", type=int, default=1, help="Reserved for interface compatibility")
    args = parser.parse_args()

    input_csv = Path(args.input_csv)
    data_root = Path(args.data_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = load_csv(input_csv, args.subset)
    extracted_paths: list[str] = []
    for row in rows:
        output_path = output_dir / f"{row['zinc_id']}.pdbqt"
        success, _, error = extract_from_tgz(
            Path(row["tgz_path"]),
            row["internal_path"],
            output_path,
            data_root,
        )
        if not success:
            raise RuntimeError(error)
        extracted_paths.append(str(output_path))

    if args.ligand_list:
        ligand_list = Path(args.ligand_list)
        ligand_list.parent.mkdir(parents=True, exist_ok=True)
        ligand_list.write_text("\n".join(sorted(extracted_paths)) + ("\n" if extracted_paths else ""))

    print(f"Extracted {len(extracted_paths)} ligands for subset {args.subset}")


if __name__ == "__main__":
    main()
