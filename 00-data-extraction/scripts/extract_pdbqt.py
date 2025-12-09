#!/usr/bin/env python3
"""
Extract PDBQT files from tgz archives using multiprocessing.
Parallel extraction for high throughput (~350 molecules/second with 48 cores).
"""
import pandas as pd
import subprocess
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
import time
import argparse


def extract_one(args):
    """Extract a single PDBQT file from tgz archive."""
    zinc_id, pdbqt_path, data_root, output_dir = args
    try:
        # Format: tgz_path:internal_path
        tgz_path, internal_path = pdbqt_path.split(":")
        full_tgz = os.path.join(data_root, tgz_path)
        out_file = os.path.join(output_dir, f"{zinc_id}.pdbqt")

        cmd = f"tar -xzOf {full_tgz} {internal_path} > {out_file}"
        result = subprocess.run(cmd, shell=True, capture_output=True)

        if result.returncode == 0 and os.path.exists(out_file):
            return (zinc_id, out_file, True)
        else:
            return (zinc_id, None, False)
    except Exception as e:
        return (zinc_id, None, False)


def main():
    parser = argparse.ArgumentParser(
        description='Extract PDBQT files from tgz archives in parallel'
    )
    parser.add_argument(
        '--input', required=True,
        help='Input CSV file with zinc_id, smiles, pdbqt_path columns'
    )
    parser.add_argument(
        '--output_dir', default='./ligands',
        help='Output directory for extracted PDBQT files'
    )
    parser.add_argument(
        '--data_root', required=True,
        help='Root directory containing the tgz files'
    )
    parser.add_argument(
        '--ligand_list', default='./outputs/ligand_files.txt',
        help='Output file listing all extracted PDBQT paths'
    )
    parser.add_argument(
        '--ncpus', type=int, default=0,
        help='Number of CPU cores (0 = auto-detect)'
    )

    args = parser.parse_args()

    # Auto-detect CPU count
    n_workers = args.ncpus if args.ncpus > 0 else cpu_count()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.ligand_list), exist_ok=True)

    print("=" * 60)
    print("PDBQT Extraction - PARALLEL VERSION")
    print("=" * 60)
    print(f"Input CSV: {args.input}")
    print(f"Data root: {args.data_root}")
    print(f"Output dir: {args.output_dir}")
    print(f"Ligand list: {args.ligand_list}")
    print(f"CPU cores: {n_workers}")
    print()

    # Load data
    print("Loading compounds...")
    df = pd.read_csv(args.input)
    print(f"Total molecules: {len(df):,}")

    # Prepare arguments
    args_list = [
        (row["zinc_id"], row["pdbqt_path"], args.data_root, args.output_dir)
        for _, row in df.iterrows()
    ]

    print()
    start_time = time.time()
    success = 0

    with open(args.ligand_list, "w") as f_list:
        with Pool(n_workers) as pool:
            for i, result in enumerate(pool.imap_unordered(extract_one, args_list, chunksize=100)):
                zinc_id, out_file, ok = result
                if ok:
                    f_list.write(f"{out_file}\n")
                    success += 1

                if (i + 1) % 10000 == 0:
                    elapsed = time.time() - start_time
                    rate = (i + 1) / elapsed
                    eta = (len(df) - i - 1) / rate / 60
                    print(f"Progress: {i+1:,}/{len(df):,} ({100*(i+1)/len(df):.1f}%) - {rate:.0f} mol/s - ETA: {eta:.1f} min")

    elapsed = time.time() - start_time
    print()
    print("=" * 60)
    print(f"Done! Extracted {success:,}/{len(df):,} in {elapsed/60:.1f} minutes")
    print(f"Rate: {success/elapsed:.0f} molecules/second")
    print(f"Ligand list: {args.ligand_list}")
    print("=" * 60)


if __name__ == "__main__":
    main()
