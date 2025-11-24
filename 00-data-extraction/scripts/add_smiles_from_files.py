#!/usr/bin/env python3
"""
Add SMILES to ZINC22 catalog using downloaded 2D .smi.gz files
Builds lookup dictionary from .smi.gz files and merges with catalog CSV.

Usage: uv run python add_smiles_from_files.py
"""

import csv
import gzip
import os
from pathlib import Path
from collections import defaultdict
import time

# Configuration
BASE_DIR = "/data/petretto/home/tiesunlong/docking_zinc22_H27-29_P300-400"
SMILES_DIR = Path(BASE_DIR) / "smiles_2d"
INPUT_CSV = Path(BASE_DIR) / "zinc22_catalog_complete.csv"
OUTPUT_CSV = Path(BASE_DIR) / "zinc22_catalog_with_smiles.csv"

def build_smiles_lookup(smiles_dir):
    """
    Build ZINC ID → SMILES lookup dictionary from .smi.gz files.

    .smi file format (tab-separated):
        SMILES\tZINC_ID\n

    Returns:
        dict: {zinc_id: smiles}
    """
    print("=" * 60)
    print("Building SMILES Lookup Dictionary")
    print("=" * 60)

    smiles_lookup = {}
    smi_files = list(smiles_dir.glob("**/*.smi.gz"))

    print(f"Found {len(smi_files)} .smi.gz files")
    print()

    start_time = time.time()

    for i, smi_file in enumerate(smi_files, 1):
        # Skip empty files (failed downloads)
        if smi_file.stat().st_size == 0:
            print(f"[{i}/{len(smi_files)}] Skipping empty file: {smi_file.name}")
            continue

        file_size_gb = smi_file.stat().st_size / 1e9
        print(f"[{i}/{len(smi_files)}] Processing: {smi_file.name} ({file_size_gb:.1f} GB)")

        try:
            with gzip.open(smi_file, 'rt') as f:
                count = 0
                for line in f:
                    line = line.strip()
                    if not line:
                        continue

                    # Parse: SMILES\tZINC_ID
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        smiles = parts[0]
                        zinc_id = parts[1]
                        smiles_lookup[zinc_id] = smiles
                        count += 1

                print(f"  Added {count:,} SMILES entries")

        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    elapsed = time.time() - start_time
    print()
    print(f"Total SMILES entries: {len(smiles_lookup):,}")
    print(f"Time elapsed: {elapsed:.1f}s")
    print()

    return smiles_lookup

def merge_smiles_to_catalog(input_csv, output_csv, smiles_lookup):
    """
    Merge SMILES column into catalog CSV.

    Args:
        input_csv: Input catalog file
        output_csv: Output catalog with SMILES
        smiles_lookup: Dict mapping ZINC ID to SMILES
    """
    print("=" * 60)
    print("Merging SMILES to Catalog CSV")
    print("=" * 60)
    print(f"Input: {input_csv}")
    print(f"Output: {output_csv}")
    print()

    # Count total rows
    with open(input_csv, 'r') as f:
        total_rows = sum(1 for _ in f) - 1  # Exclude header

    print(f"Total molecules: {total_rows:,}")
    print()

    # Process
    found = 0
    not_found = 0
    start_time = time.time()

    with open(input_csv, 'r') as infile, open(output_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['smiles']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for i, row in enumerate(reader, 1):
            zinc_id = row['zinc_id']

            # Lookup SMILES
            smiles = smiles_lookup.get(zinc_id, '')

            if smiles:
                found += 1
            else:
                not_found += 1

            row['smiles'] = smiles
            writer.writerow(row)

            # Progress update
            if i % 10000 == 0:
                elapsed = time.time() - start_time
                rate = i / elapsed if elapsed > 0 else 0
                remaining = (total_rows - i) / rate if rate > 0 else 0

                print(f"Progress: {i:,}/{total_rows:,} ({i*100/total_rows:.1f}%)")
                print(f"  Found: {found:,}, Not found: {not_found:,}")
                print(f"  Rate: {rate:.0f} rows/sec, ETA: {remaining:.0f}s")
                print()

    elapsed = time.time() - start_time

    # Summary
    print()
    print("=" * 60)
    print("MERGE COMPLETE")
    print("=" * 60)
    print(f"Total processed: {total_rows:,}")
    print(f"SMILES found: {found:,} ({found*100/total_rows:.2f}%)")
    print(f"SMILES not found: {not_found:,} ({not_found*100/total_rows:.2f}%)")
    print(f"Time elapsed: {elapsed/60:.1f} minutes")
    print(f"Output saved to: {output_csv}")
    print("=" * 60)

def main():
    print("=" * 60)
    print("ZINC22 SMILES Addition from .smi.gz Files")
    print("=" * 60)
    print()

    # Check input files
    if not INPUT_CSV.exists():
        print(f"ERROR: Input catalog not found: {INPUT_CSV}")
        return

    if not SMILES_DIR.exists():
        print(f"ERROR: SMILES directory not found: {SMILES_DIR}")
        return

    # Step 1: Build lookup dictionary
    smiles_lookup = build_smiles_lookup(SMILES_DIR)

    if not smiles_lookup:
        print("ERROR: No SMILES entries found!")
        return

    # Step 2: Merge to catalog
    merge_smiles_to_catalog(INPUT_CSV, OUTPUT_CSV, smiles_lookup)

if __name__ == "__main__":
    main()
