#!/usr/bin/env python3
"""
Export Docking Results to Sorted CSV

Parses Uni-Dock output files (*_out.pdbqt) to extract docking scores,
merges with SMILES data, and exports a sorted CSV ranked by binding affinity.

Usage:
    uv run python export_results.py \
        --results-dir ./results \
        --smiles-csv ./smiles.csv \
        --output ./docking_results.csv

Features:
    - Multi-key SMILES lookup (handles format variations)
    - Parallel processing for large result sets
    - Memory-efficient streaming for millions of results
    - Automatic sorting by docking score (best first)
"""

import argparse
import csv
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional


def parse_pdbqt_score(filepath: Path) -> tuple:
    """
    Extract docking score from a PDBQT file.

    Looks for: REMARK VINA RESULT: -8.5 0.000 0.000

    Returns: (zinc_id, score, error_msg)
    """
    try:
        zinc_id = filepath.stem.replace('_out', '')

        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT:'):
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        score = float(parts[3])
                        return zinc_id, score, None

        return zinc_id, None, "No VINA RESULT found"

    except Exception as e:
        return filepath.stem, None, str(e)


def load_smiles_lookup(csv_path: Path) -> dict:
    """
    Load SMILES data into a lookup dictionary.

    Handles multiple zinc_id column names and creates fallback keys
    for format variations.
    """
    lookup = {}
    zinc_id_columns = ['zinc_id', 'ZINC_ID', 'zincid', 'id', 'ID']

    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)

        # Find the zinc_id column
        zinc_col = None
        for col in zinc_id_columns:
            if col in reader.fieldnames:
                zinc_col = col
                break

        if not zinc_col:
            raise ValueError(f"No zinc_id column found. Available: {reader.fieldnames}")

        # Find the SMILES column
        smiles_columns = ['smiles', 'SMILES', 'Smiles', 'canonical_smiles']
        smiles_col = None
        for col in smiles_columns:
            if col in reader.fieldnames:
                smiles_col = col
                break

        if not smiles_col:
            raise ValueError(f"No SMILES column found. Available: {reader.fieldnames}")

        for row in reader:
            zinc_id = row[zinc_col]
            smiles = row[smiles_col]

            if zinc_id and smiles:
                # Store original key
                lookup[zinc_id] = smiles

                # Create fallback keys for format variations
                # ZINCxxx.0 -> ZINCxxx (remove conformer)
                if '.' in zinc_id:
                    base_id = zinc_id.rsplit('.', 1)[0]
                    if base_id not in lookup:
                        lookup[base_id] = smiles

    return lookup


def lookup_smiles(zinc_id: str, smiles_lookup: dict) -> Optional[str]:
    """
    Multi-key SMILES lookup with format fallbacks.

    Tries multiple key formats to handle variations:
    1. Full ID: ZINCxxx.0.N
    2. Without variant: ZINCxxx.0
    3. Base only: ZINCxxx
    """
    # Try exact match first
    if zinc_id in smiles_lookup:
        return smiles_lookup[zinc_id]

    # Try removing variant letter (ZINCxxx.0.N -> ZINCxxx.0)
    if '.' in zinc_id:
        parts = zinc_id.rsplit('.', 1)
        if len(parts[1]) == 1 and parts[1].isalpha():
            without_variant = parts[0]
            if without_variant in smiles_lookup:
                return smiles_lookup[without_variant]

            # Try base only (ZINCxxx.0 -> ZINCxxx)
            if '.' in without_variant:
                base_id = without_variant.rsplit('.', 1)[0]
                if base_id in smiles_lookup:
                    return smiles_lookup[base_id]

    # Try base ID directly
    if '.' in zinc_id:
        base_id = zinc_id.split('.')[0]
        if base_id in smiles_lookup:
            return smiles_lookup[base_id]

    return None


def process_batch(filepaths: list) -> list:
    """Process a batch of PDBQT files."""
    results = []
    for fp in filepaths:
        zinc_id, score, error = parse_pdbqt_score(fp)
        if score is not None:
            results.append((zinc_id, score))
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Export docking results to sorted CSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic export
  uv run python export_results.py \\
      --results-dir results \\
      --smiles-csv data/smiles.csv \\
      --output docking_results.csv

  # Export top 10000 hits only
  uv run python export_results.py \\
      --results-dir results \\
      --smiles-csv data/smiles.csv \\
      --output top_hits.csv \\
      --top-n 10000

  # High-performance parallel processing
  uv run python export_results.py \\
      --results-dir results \\
      --smiles-csv data/smiles.csv \\
      --output results.csv \\
      --ncpus 48
"""
    )
    parser.add_argument('--results-dir', required=True,
                        help='Directory containing *_out.pdbqt files')
    parser.add_argument('--smiles-csv', required=True,
                        help='CSV file with zinc_id and smiles columns')
    parser.add_argument('--output', required=True,
                        help='Output CSV path')
    parser.add_argument('--top-n', type=int, default=0,
                        help='Limit output to top N results (0 = all)')
    parser.add_argument('--ncpus', type=int, default=os.cpu_count(),
                        help='Number of parallel workers')
    parser.add_argument('--batch-size', type=int, default=10000,
                        help='Files per batch for parallel processing')

    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    smiles_csv = Path(args.smiles_csv)
    output_path = Path(args.output)

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)

    if not smiles_csv.exists():
        print(f"Error: SMILES CSV not found: {smiles_csv}", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Find all result files
    print(f"Scanning results directory: {results_dir}", file=sys.stderr)
    result_files = list(results_dir.glob('*_out.pdbqt'))
    print(f"Found {len(result_files):,} result files", file=sys.stderr)

    if not result_files:
        print("Error: No *_out.pdbqt files found", file=sys.stderr)
        sys.exit(1)

    # Load SMILES lookup
    print(f"Loading SMILES data: {smiles_csv}", file=sys.stderr)
    smiles_lookup = load_smiles_lookup(smiles_csv)
    print(f"Loaded {len(smiles_lookup):,} SMILES entries", file=sys.stderr)

    # Process results in parallel
    print(f"Processing with {args.ncpus} workers...", file=sys.stderr)

    all_results = []
    batches = [result_files[i:i+args.batch_size]
               for i in range(0, len(result_files), args.batch_size)]

    with ProcessPoolExecutor(max_workers=args.ncpus) as executor:
        futures = {executor.submit(process_batch, batch): i
                   for i, batch in enumerate(batches)}

        for future in as_completed(futures):
            batch_results = future.result()
            all_results.extend(batch_results)

            if len(all_results) % 100000 < args.batch_size:
                print(f"  Processed: {len(all_results):,} results", file=sys.stderr)

    print(f"Total results parsed: {len(all_results):,}", file=sys.stderr)

    # Sort by score (ascending = best first, more negative is better)
    print("Sorting by docking score...", file=sys.stderr)
    all_results.sort(key=lambda x: x[1])

    # Limit if requested
    if args.top_n > 0:
        all_results = all_results[:args.top_n]
        print(f"Limited to top {args.top_n:,} results", file=sys.stderr)

    # Write output with SMILES
    print(f"Writing output: {output_path}", file=sys.stderr)
    smiles_found = 0
    smiles_missing = 0

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['zinc_id', 'score', 'smiles'])

        for zinc_id, score in all_results:
            smiles = lookup_smiles(zinc_id, smiles_lookup)

            if smiles:
                smiles_found += 1
            else:
                smiles_missing += 1
                smiles = ''

            writer.writerow([zinc_id, f"{score:.2f}", smiles])

    # Summary
    print(f"\nExport complete:", file=sys.stderr)
    print(f"  Total rows: {len(all_results):,}", file=sys.stderr)
    print(f"  SMILES found: {smiles_found:,}", file=sys.stderr)
    print(f"  SMILES missing: {smiles_missing:,}", file=sys.stderr)
    print(f"  Output: {output_path}", file=sys.stderr)

    if smiles_missing > 0:
        pct_missing = 100.0 * smiles_missing / len(all_results)
        print(f"\nWarning: {pct_missing:.1f}% of entries missing SMILES", file=sys.stderr)


if __name__ == '__main__':
    main()
