#!/usr/bin/env python3
"""
ZINC22 Full-scale PDBQT to SMILES Conversion using Open Babel
Processes all 2,141 .tgz files to extract ~10.7 million SMILES
"""

import tarfile
import tempfile
import os
from pathlib import Path
import csv
from openbabel import pybel
import sys
import argparse
from multiprocessing import Pool, Manager
from tqdm import tqdm
import json
from datetime import datetime

def extract_zinc_id_from_pdbqt(content):
    """Extract ZINC ID from REMARK line in PDBQT file"""
    for line in content.split('\n'):
        if line.startswith('REMARK') and 'Name' in line:
            parts = line.split('=')
            if len(parts) > 1:
                return parts[-1].strip()
    return None

def pdbqt_to_smiles(pdbqt_path):
    """
    Convert PDBQT to SMILES using Open Babel
    Returns: (smiles, zinc_id, error_msg)
    """
    try:
        with open(pdbqt_path, 'r') as f:
            content = f.read()
        zinc_id = extract_zinc_id_from_pdbqt(content)
        if not zinc_id:
            zinc_id = Path(pdbqt_path).stem

        mol = next(pybel.readfile("pdbqt", str(pdbqt_path)))
        smiles_output = mol.write("smi").strip()
        smiles = smiles_output.split()[0] if smiles_output else None

        if not smiles:
            return None, zinc_id, "SMILES generation failed"

        return smiles, zinc_id, None

    except StopIteration:
        zinc_id = Path(pdbqt_path).stem
        return None, zinc_id, "No molecule found in PDBQT"
    except Exception as e:
        zinc_id = Path(pdbqt_path).stem
        return None, zinc_id, f"Exception: {str(e)}"

def process_single_tgz(args):
    """
    Process a single .tgz file (called by multiprocessing workers)
    Returns: (tgz_filename, results_list, errors_list)
    """
    tgz_path, base_output_dir = args
    results = []
    errors = []

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Extract
            with tarfile.open(tgz_path, 'r:gz') as tar:
                tar.extractall(tmpdir)

            # Find all PDBQT files
            pdbqt_files = list(Path(tmpdir).rglob('*.pdbqt'))

            # Process each PDBQT
            for pdbqt_file in pdbqt_files:
                smiles, zinc_id, error = pdbqt_to_smiles(str(pdbqt_file))

                if error:
                    errors.append({
                        'zinc_id': zinc_id or 'UNKNOWN',
                        'file': str(pdbqt_file.name),
                        'tgz': Path(tgz_path).name,
                        'error': error
                    })
                else:
                    # Store full path relative to zinc22_H27-29_P300-400
                    # Extract relative path from tgz_path
                    tgz_rel = str(Path(tgz_path).relative_to(base_output_dir.parent))
                    pdbqt_rel = pdbqt_file.relative_to(tmpdir)
                    full_path = f"{tgz_rel}:{pdbqt_rel}"

                    results.append({
                        'zinc_id': zinc_id,
                        'smiles': smiles,
                        'pdbqt_path': full_path
                    })

        return (Path(tgz_path).name, results, errors)

    except Exception as e:
        return (Path(tgz_path).name, [], [{'error': f"TGZ processing failed: {str(e)}"}])

def main():
    parser = argparse.ArgumentParser(description='ZINC22 PDBQT to SMILES - Full Processing')
    parser.add_argument('--input_dir', required=True, help='Base directory containing .tgz files')
    parser.add_argument('--output_csv', required=True, help='Output CSV file path')
    parser.add_argument('--error_log', default='errors.log', help='Error log file')
    parser.add_argument('--ncpus', type=int, default=16, help='Number of CPU cores')
    parser.add_argument('--resume', action='store_true', help='Resume from progress file')
    parser.add_argument('--progress_file', default='progress.json', help='Progress tracking file')

    args = parser.parse_args()

    print("="*70)
    print("ZINC22 PDBQT to SMILES - Full Processing (Open Babel)")
    print("="*70)
    print(f"Input directory: {args.input_dir}")
    print(f"Output CSV: {args.output_csv}")
    print(f"CPU cores: {args.ncpus}")
    print(f"Resume mode: {args.resume}")
    print()

    # Find all .tgz files
    base_dir = Path(args.input_dir)
    all_tgz_files = sorted(base_dir.rglob("*.tgz"))
    print(f"Found {len(all_tgz_files)} .tgz files")

    # Load progress if resuming
    processed_files = set()
    if args.resume and os.path.exists(args.progress_file):
        with open(args.progress_file, 'r') as f:
            progress_data = json.load(f)
            processed_files = set(progress_data.get('processed_tgz', []))
        print(f"Resuming: {len(processed_files)} files already processed")

    # Filter out already processed files
    tgz_to_process = [(str(f), base_dir) for f in all_tgz_files
                      if f.name not in processed_files]
    print(f"Files to process: {len(tgz_to_process)}")

    if not tgz_to_process:
        print("No files to process!")
        return

    # Prepare output
    csv_mode = 'a' if args.resume and os.path.exists(args.output_csv) else 'w'
    csv_file = open(args.output_csv, csv_mode, newline='')
    csv_writer = csv.DictWriter(csv_file, fieldnames=['zinc_id', 'smiles', 'pdbqt_path'])
    if csv_mode == 'w':
        csv_writer.writeheader()

    error_file = open(args.error_log, 'a')

    # Statistics
    total_molecules = 0
    total_errors = 0
    processed_tgz_count = len(processed_files)

    # Process with multiprocessing
    print(f"\nProcessing with {args.ncpus} workers...")
    print()

    with Pool(args.ncpus) as pool:
        with tqdm(total=len(tgz_to_process), desc="Processing .tgz files") as pbar:
            for tgz_name, results, errors in pool.imap_unordered(process_single_tgz, tgz_to_process):
                # Write results
                if results:
                    csv_writer.writerows(results)
                    csv_file.flush()
                    total_molecules += len(results)

                # Write errors
                if errors:
                    for err in errors:
                        error_file.write(f"{datetime.now()} | {err}\n")
                    error_file.flush()
                    total_errors += len(errors)

                # Update progress
                processed_files.add(tgz_name)
                processed_tgz_count += 1

                # Save progress
                if processed_tgz_count % 10 == 0:
                    progress_data = {
                        'processed_tgz': list(processed_files),
                        'total_molecules': total_molecules,
                        'total_errors': total_errors,
                        'last_update': datetime.now().isoformat()
                    }
                    with open(args.progress_file, 'w') as f:
                        json.dump(progress_data, f, indent=2)

                pbar.update(1)
                pbar.set_postfix({
                    'molecules': total_molecules,
                    'errors': total_errors,
                    'success_rate': f'{100*(total_molecules/(total_molecules+total_errors) if total_molecules+total_errors>0 else 0):.1f}%'
                })

    # Close files
    csv_file.close()
    error_file.close()

    # Final progress save
    progress_data = {
        'processed_tgz': list(processed_files),
        'total_molecules': total_molecules,
        'total_errors': total_errors,
        'last_update': datetime.now().isoformat(),
        'completed': True
    }
    with open(args.progress_file, 'w') as f:
        json.dump(progress_data, f, indent=2)

    # Summary
    print("\n" + "="*70)
    print("Processing Complete!")
    print("="*70)
    print(f"Total .tgz files processed: {len(all_tgz_files)}")
    print(f"Total molecules extracted: {total_molecules:,}")
    print(f"Total errors: {total_errors:,}")
    if total_molecules + total_errors > 0:
        success_rate = 100 * total_molecules / (total_molecules + total_errors)
        print(f"Success rate: {success_rate:.2f}%")
    print(f"\nOutput CSV: {args.output_csv}")
    print(f"Error log: {args.error_log}")
    print(f"Progress file: {args.progress_file}")

    # Write summary
    summary_file = args.output_csv.replace('.csv', '_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"ZINC22 PDBQT to SMILES Processing Summary\n")
        f.write(f"Generated: {datetime.now()}\n")
        f.write(f"="*70 + "\n\n")
        f.write(f"Input directory: {args.input_dir}\n")
        f.write(f"Total .tgz files: {len(all_tgz_files)}\n")
        f.write(f"Total molecules: {total_molecules:,}\n")
        f.write(f"Total errors: {total_errors:,}\n")
        if total_molecules + total_errors > 0:
            f.write(f"Success rate: {success_rate:.2f}%\n")
        f.write(f"\nOutput CSV: {args.output_csv}\n")
        f.write(f"Error log: {args.error_log}\n")

    print(f"Summary written to: {summary_file}")

if __name__ == "__main__":
    main()
