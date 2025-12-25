#!/usr/bin/env python3
"""
Extract Ligand PDBQT Files from TGZ Archives

Extracts individual PDBQT files from ZINC22 tar.gz archives in parallel.
Supports filtering by subset and resuming interrupted extractions.

Usage:
    uv run python extract_ligands.py \
        --input-csv sampled.csv \
        --output-dir ./ligands \
        --data-root ./data/zinc22 \
        --ncpus 24
"""

import argparse
import csv
import os
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional


def extract_from_tgz(
    tgz_path: Path,
    internal_path: str,
    output_path: Path,
    data_root: Path
) -> tuple:
    """
    Extract a single file from a tgz archive.

    Returns: (success: bool, output_path: str, error_msg: str)
    """
    try:
        full_tgz_path = data_root / tgz_path

        if not full_tgz_path.exists():
            return False, str(output_path), f"Archive not found: {full_tgz_path}"

        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with tarfile.open(full_tgz_path, 'r:gz') as tar:
            try:
                member = tar.getmember(internal_path)
                # Extract to temp file first, then move
                with tempfile.NamedTemporaryFile(delete=False) as tmp:
                    extracted = tar.extractfile(member)
                    if extracted:
                        tmp.write(extracted.read())
                        tmp.flush()
                        os.rename(tmp.name, output_path)
                        return True, str(output_path), ""
                    else:
                        return False, str(output_path), f"Could not extract: {internal_path}"
            except KeyError:
                return False, str(output_path), f"File not in archive: {internal_path}"

    except Exception as e:
        return False, str(output_path), str(e)


def parse_pdbqt_path(pdbqt_path: str) -> tuple:
    """
    Parse pdbqt_path column into tgz_path and internal_path.

    Format: data/zinc22/H27P310-K.tgz:H27/H27P310/FZ/to/ZINCrB000001FZto.0.K.pdbqt
    Returns: (tgz_path, internal_path, zinc_id)
    """
    if ':' not in pdbqt_path:
        return None, None, None

    parts = pdbqt_path.split(':', 1)
    tgz_path = parts[0]
    internal_path = parts[1]

    # Extract zinc_id from internal path (filename without extension)
    filename = os.path.basename(internal_path)
    zinc_id = filename.replace('.pdbqt', '')

    return tgz_path, internal_path, zinc_id


def load_csv(csv_path: Path, subset_filter: Optional[str] = None) -> list:
    """Load CSV and optionally filter by subset."""
    rows = []
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdbqt_path = row.get('pdbqt_path', '')
            if not pdbqt_path:
                continue

            # Apply subset filter if specified
            if subset_filter and subset_filter not in pdbqt_path:
                continue

            tgz_path, internal_path, zinc_id = parse_pdbqt_path(pdbqt_path)
            if tgz_path and internal_path and zinc_id:
                rows.append({
                    'tgz_path': tgz_path,
                    'internal_path': internal_path,
                    'zinc_id': zinc_id
                })

    return rows


def main():
    parser = argparse.ArgumentParser(
        description='Extract ligand PDBQT files from TGZ archives',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract all ligands from CSV
  uv run python extract_ligands.py \\
      --input-csv outputs/sampled.csv \\
      --output-dir ligands \\
      --data-root data/zinc22

  # Extract only zinc-22g subset
  uv run python extract_ligands.py \\
      --input-csv outputs/sampled.csv \\
      --output-dir ligands/zinc-22g \\
      --data-root data/zinc22 \\
      --subset-filter zinc-22g

  # Resume interrupted extraction
  uv run python extract_ligands.py \\
      --input-csv outputs/sampled.csv \\
      --output-dir ligands \\
      --data-root data/zinc22 \\
      --resume
"""
    )
    parser.add_argument('--input-csv', required=True, help='CSV with pdbqt_path column')
    parser.add_argument('--output-dir', required=True, help='Directory for extracted PDBQT files')
    parser.add_argument('--data-root', required=True, help='Root directory containing TGZ archives')
    parser.add_argument('--subset-filter', help='Filter by subset pattern (e.g., zinc-22g)')
    parser.add_argument('--ligand-list', help='Output file listing all extracted ligands')
    parser.add_argument('--ncpus', type=int, default=os.cpu_count(), help='Number of workers')
    parser.add_argument('--resume', action='store_true', help='Skip existing files')

    args = parser.parse_args()

    input_csv = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    data_root = Path(args.data_root)

    if not input_csv.exists():
        print(f"Error: Input CSV not found: {input_csv}", file=sys.stderr)
        sys.exit(1)

    if not data_root.exists():
        print(f"Error: Data root not found: {data_root}", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading CSV: {input_csv}", file=sys.stderr)
    rows = load_csv(input_csv, args.subset_filter)
    print(f"Found {len(rows):,} ligands to extract", file=sys.stderr)

    if args.subset_filter:
        print(f"Subset filter: {args.subset_filter}", file=sys.stderr)

    # Prepare extraction tasks
    tasks = []
    skipped = 0
    for row in rows:
        output_path = output_dir / f"{row['zinc_id']}.pdbqt"

        if args.resume and output_path.exists():
            skipped += 1
            continue

        tasks.append({
            'tgz_path': Path(row['tgz_path']),
            'internal_path': row['internal_path'],
            'output_path': output_path
        })

    if skipped > 0:
        print(f"Skipping {skipped:,} existing files (--resume)", file=sys.stderr)

    print(f"Extracting {len(tasks):,} files with {args.ncpus} workers...", file=sys.stderr)

    success_count = 0
    error_count = 0
    extracted_files = []

    with ThreadPoolExecutor(max_workers=args.ncpus) as executor:
        futures = {
            executor.submit(
                extract_from_tgz,
                task['tgz_path'],
                task['internal_path'],
                task['output_path'],
                data_root
            ): task
            for task in tasks
        }

        for i, future in enumerate(as_completed(futures), 1):
            success, output_path, error_msg = future.result()

            if success:
                success_count += 1
                extracted_files.append(output_path)
            else:
                error_count += 1
                if error_count <= 10:
                    print(f"Error: {error_msg}", file=sys.stderr)

            if i % 10000 == 0:
                print(f"  Progress: {i:,}/{len(tasks):,} ({100*i/len(tasks):.1f}%)", file=sys.stderr)

    print(f"\nExtraction complete:", file=sys.stderr)
    print(f"  Success: {success_count:,}", file=sys.stderr)
    print(f"  Errors: {error_count:,}", file=sys.stderr)
    print(f"  Skipped: {skipped:,}", file=sys.stderr)

    # Write ligand list if requested
    if args.ligand_list:
        ligand_list_path = Path(args.ligand_list)
        ligand_list_path.parent.mkdir(parents=True, exist_ok=True)

        # Include both new and existing files
        all_files = list(output_dir.glob('*.pdbqt'))
        with open(ligand_list_path, 'w') as f:
            for filepath in sorted(all_files):
                f.write(f"{filepath}\n")

        print(f"\nLigand list written: {ligand_list_path}", file=sys.stderr)
        print(f"  Total files: {len(all_files):,}", file=sys.stderr)


if __name__ == '__main__':
    main()
