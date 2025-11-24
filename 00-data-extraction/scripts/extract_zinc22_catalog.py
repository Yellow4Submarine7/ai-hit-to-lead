#!/usr/bin/env python3
"""
ZINC22 Catalog Extraction Script
Extracts file information from TGZ archives and generates CSV mapping table.

Usage: uv run python extract_zinc22_catalog.py
"""

import os
import tarfile
import re
import csv
from collections import defaultdict
from pathlib import Path
import time

# Configuration
BASE_DIR = "/data/petretto/home/tiesunlong/docking_zinc22_H27-29_P300-400"
DOWNLOADS_DIR = os.path.join(BASE_DIR, "downloads")
OUTPUT_DIR = BASE_DIR

def parse_pdbqt_filename(filepath):
    """
    Parse PDBQT filename to extract metadata.

    Example: H29/H29P380/Ai/37/ZINCtI000009Ai37.0.N.pdbqt
    Returns: (zinc_id, h_class, p_class, conformer, protonation, rel_path)
    """
    filename = os.path.basename(filepath)

    # Pattern: ZINCxxxxxx.0.N.pdbqt
    match = re.match(r'^(ZINC[a-zA-Z0-9]+)\.(\d+)\.([MNLOP])\.pdbqt$', filename)
    if not match:
        return None

    zinc_id = match.group(1)
    conformer = match.group(2)
    protonation = match.group(3)

    # Extract H and P class from path
    # Example path: H29/H29P380/Ai/37/filename.pdbqt
    parts = filepath.split('/')
    h_class = None
    p_class = None

    for part in parts:
        if re.match(r'^H\d+$', part):
            h_class = part
        elif re.match(r'^H\d+P\d+$', part):
            p_class = part

    return {
        'zinc_id': zinc_id,
        'h_class': h_class,
        'p_class': p_class,
        'conformer': conformer,
        'protonation': protonation,
        'file_path': filepath
    }

def extract_from_tgz(tgz_path):
    """Extract file list from a TGZ archive."""
    results = []
    try:
        with tarfile.open(tgz_path, 'r:gz') as tar:
            for member in tar.getmembers():
                if member.name.endswith('.pdbqt'):
                    parsed = parse_pdbqt_filename(member.name)
                    if parsed:
                        results.append(parsed)
    except Exception as e:
        print(f"Error reading {tgz_path}: {e}")
    return results

def find_all_tgz_files(downloads_dir):
    """Find all TGZ files recursively."""
    tgz_files = []
    for root, dirs, files in os.walk(downloads_dir):
        for f in files:
            if f.endswith('.tgz'):
                tgz_files.append(os.path.join(root, f))
    return sorted(tgz_files)

def aggregate_by_zinc_id(all_files):
    """Aggregate files by ZINC ID."""
    aggregated = defaultdict(lambda: {
        'h_class': set(),
        'p_class': set(),
        'conformers': set(),
        'protonations': set(),
        'file_paths': []
    })

    for f in all_files:
        zinc_id = f['zinc_id']
        aggregated[zinc_id]['h_class'].add(f['h_class'])
        aggregated[zinc_id]['p_class'].add(f['p_class'])
        aggregated[zinc_id]['conformers'].add(f['conformer'])
        aggregated[zinc_id]['protonations'].add(f['protonation'])
        aggregated[zinc_id]['file_paths'].append(f['file_path'])

    return aggregated

def main():
    print("=" * 60)
    print("ZINC22 Catalog Extraction")
    print("=" * 60)
    print(f"Downloads directory: {DOWNLOADS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()

    # Find all TGZ files
    print("Step 1: Finding TGZ files...")
    tgz_files = find_all_tgz_files(DOWNLOADS_DIR)
    print(f"  Found {len(tgz_files)} TGZ files")
    print()

    # Extract file information from all TGZ
    print("Step 2: Extracting file information from TGZ archives...")
    all_files = []
    start_time = time.time()

    for i, tgz_path in enumerate(tgz_files):
        files = extract_from_tgz(tgz_path)
        all_files.extend(files)

        if (i + 1) % 50 == 0:
            elapsed = time.time() - start_time
            remaining = (elapsed / (i + 1)) * (len(tgz_files) - i - 1)
            print(f"  Processed {i + 1}/{len(tgz_files)} TGZ files... "
                  f"(ETA: {remaining:.0f}s)")

    print(f"  Total PDBQT files found: {len(all_files)}")
    print()

    # Aggregate by ZINC ID
    print("Step 3: Aggregating by ZINC ID...")
    aggregated = aggregate_by_zinc_id(all_files)
    print(f"  Unique ZINC IDs: {len(aggregated)}")
    print()

    # Save unique ZINC IDs
    unique_ids_file = os.path.join(OUTPUT_DIR, "unique_zinc_ids.txt")
    with open(unique_ids_file, 'w') as f:
        for zinc_id in sorted(aggregated.keys()):
            f.write(zinc_id + '\n')
    print(f"  Saved unique IDs to: {unique_ids_file}")

    # Generate CSV
    print("Step 4: Generating CSV catalog...")
    csv_file = os.path.join(OUTPUT_DIR, "zinc22_catalog_complete.csv")

    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'zinc_id', 'h_class', 'p_class', 'num_conformers',
            'protonation_states', 'file_count', 'file_paths'
        ])

        for zinc_id in sorted(aggregated.keys()):
            data = aggregated[zinc_id]
            writer.writerow([
                zinc_id,
                ';'.join(sorted(data['h_class'])),
                ';'.join(sorted(data['p_class'])),
                len(data['conformers']),
                ';'.join(sorted(data['protonations'])),
                len(data['file_paths']),
                ';'.join(sorted(data['file_paths']))
            ])

    print(f"  Saved CSV to: {csv_file}")
    print()

    # Generate statistics
    print("Step 5: Generating statistics...")
    stats_file = os.path.join(OUTPUT_DIR, "extraction_stats_complete.txt")

    # Calculate statistics
    h_class_counts = defaultdict(int)
    p_class_counts = defaultdict(int)
    conformer_counts = defaultdict(int)
    protonation_counts = defaultdict(int)

    for zinc_id, data in aggregated.items():
        for h in data['h_class']:
            h_class_counts[h] += 1
        for p in data['p_class']:
            p_class_counts[p] += 1
        conformer_counts[len(data['conformers'])] += 1
        for prot in data['protonations']:
            protonation_counts[prot] += 1

    with open(stats_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("ZINC22 Extraction Statistics (COMPLETE)\n")
        f.write("=" * 60 + "\n\n")

        f.write(f"Total TGZ files: {len(tgz_files)}\n")
        f.write(f"Total PDBQT files: {len(all_files)}\n")
        f.write(f"Unique ZINC IDs: {len(aggregated)}\n\n")

        f.write("H-Class Distribution:\n")
        for h_class in sorted(h_class_counts.keys()):
            f.write(f"  {h_class}: {h_class_counts[h_class]}\n")
        f.write("\n")

        f.write("P-Class Distribution:\n")
        for p_class in sorted(p_class_counts.keys()):
            f.write(f"  {p_class}: {p_class_counts[p_class]}\n")
        f.write("\n")

        f.write("Conformer Count Distribution:\n")
        for num_conf in sorted(conformer_counts.keys()):
            f.write(f"  {num_conf} conformer(s): {conformer_counts[num_conf]} molecules\n")
        f.write("\n")

        f.write("Protonation State Distribution:\n")
        for prot in sorted(protonation_counts.keys()):
            f.write(f"  {prot}: {protonation_counts[prot]}\n")
        f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("Extraction completed successfully!\n")

    print(f"  Saved statistics to: {stats_file}")
    print()

    # Print summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total TGZ files: {len(tgz_files)}")
    print(f"Total PDBQT files: {len(all_files)}")
    print(f"Unique ZINC IDs: {len(aggregated)}")
    print()
    print("H-Class Distribution:")
    for h_class in sorted(h_class_counts.keys()):
        print(f"  {h_class}: {h_class_counts[h_class]}")
    print()
    print("Conformer Count Distribution:")
    for num_conf in sorted(conformer_counts.keys()):
        print(f"  {num_conf} conformer(s): {conformer_counts[num_conf]} molecules")
    print()
    print(f"Output files saved to: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
