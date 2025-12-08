#!/usr/bin/env python3
"""
Add SMILES to ZINC22 catalog using ZINC API
Queries ZINC API with batch processing, retry logic, and checkpoint support.

Usage: uv run python add_smiles_via_api.py
"""

import csv
import json
import time
import requests
from pathlib import Path
from collections import defaultdict

# Configuration - use relative paths
# Run from the 00-data-extraction directory
BASE_DIR = "./outputs"
INPUT_CSV = Path(BASE_DIR) / "zinc22_catalog_complete.csv"
OUTPUT_CSV = Path(BASE_DIR) / "zinc22_catalog_with_smiles.csv"
CHECKPOINT_FILE = Path(BASE_DIR) / "smiles_api_checkpoint.json"
FAILED_IDS_FILE = Path(BASE_DIR) / "failed_zinc_ids.txt"

# API Configuration
ZINC_API_BASE = "https://zinc.docking.org/substances"
BATCH_SIZE = 100  # Process 100 IDs at a time
REQUEST_DELAY = 0.5  # 0.5 seconds between requests (2 req/sec)
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds

def load_checkpoint():
    """Load progress from checkpoint file."""
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE, 'r') as f:
            data = json.load(f)
            return data.get('processed_count', 0), data.get('smiles_cache', {})
    return 0, {}

def save_checkpoint(processed_count, smiles_cache):
    """Save progress to checkpoint file."""
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump({
            'processed_count': processed_count,
            'smiles_cache': smiles_cache,
            'timestamp': time.time()
        }, f)

def query_zinc_api(zinc_id, retries=MAX_RETRIES):
    """
    Query ZINC API for a single ZINC ID.

    Returns:
        str: SMILES string or empty string if failed
    """
    for attempt in range(retries):
        try:
            # ZINC API endpoint: https://zinc.docking.org/substances/{zinc_id}.json
            url = f"{ZINC_API_BASE}/{zinc_id}.json"
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                data = response.json()
                # Extract SMILES from response
                smiles = data.get('smiles', '')
                return smiles
            elif response.status_code == 404:
                # ID not found in ZINC
                return ''
            else:
                print(f"  Warning: HTTP {response.status_code} for {zinc_id}")
                if attempt < retries - 1:
                    time.sleep(RETRY_DELAY)
                    continue
                return ''

        except requests.exceptions.Timeout:
            print(f"  Timeout for {zinc_id} (attempt {attempt + 1}/{retries})")
            if attempt < retries - 1:
                time.sleep(RETRY_DELAY)
                continue
            return ''
        except Exception as e:
            print(f"  Error querying {zinc_id}: {e}")
            if attempt < retries - 1:
                time.sleep(RETRY_DELAY)
                continue
            return ''

    return ''

def main():
    print("=" * 60)
    print("ZINC22 SMILES Addition via API")
    print("=" * 60)
    print()

    # Load checkpoint
    start_from, smiles_cache = load_checkpoint()
    if start_from > 0:
        print(f"Resuming from checkpoint: {start_from:,} molecules processed")
        print(f"SMILES cache size: {len(smiles_cache):,}")
        print()

    # Read all ZINC IDs from catalog
    print("Step 1: Reading ZINC IDs from catalog...")
    zinc_ids = []
    with open(INPUT_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            zinc_ids.append(row['zinc_id'])

    total_molecules = len(zinc_ids)
    print(f"  Total molecules: {total_molecules:,}")
    print()

    # Skip already processed IDs
    remaining_ids = zinc_ids[start_from:]
    print(f"Step 2: Querying ZINC API for SMILES...")
    print(f"  Remaining IDs to process: {len(remaining_ids):,}")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Request rate: ~{1/REQUEST_DELAY:.1f} req/sec")
    print()

    start_time = time.time()
    found_count = len([s for s in smiles_cache.values() if s])
    not_found_count = len(smiles_cache) - found_count
    failed_ids = []

    for i, zinc_id in enumerate(remaining_ids, 1):
        current_index = start_from + i

        # Query API
        smiles = query_zinc_api(zinc_id)
        smiles_cache[zinc_id] = smiles

        if smiles:
            found_count += 1
        else:
            not_found_count += 1
            failed_ids.append(zinc_id)

        # Progress update every 100 molecules
        if i % 100 == 0:
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            remaining = len(remaining_ids) - i
            eta = remaining / rate if rate > 0 else 0

            print(f"Progress: {current_index:,}/{total_molecules:,} ({current_index*100/total_molecules:.1f}%)")
            print(f"  Found: {found_count:,}, Not found: {not_found_count:,}")
            print(f"  Rate: {rate:.1f} IDs/sec, ETA: {eta/3600:.1f} hours")
            print()

            # Save checkpoint every 100 molecules
            save_checkpoint(current_index, smiles_cache)

        # Rate limiting
        time.sleep(REQUEST_DELAY)

    # Final checkpoint
    save_checkpoint(total_molecules, smiles_cache)

    # Save failed IDs
    if failed_ids:
        with open(FAILED_IDS_FILE, 'w') as f:
            for zinc_id in failed_ids:
                f.write(zinc_id + '\n')
        print(f"Saved {len(failed_ids):,} failed IDs to: {FAILED_IDS_FILE}")
        print()

    print("Step 3: Merging SMILES to catalog CSV...")
    print(f"  Input: {INPUT_CSV}")
    print(f"  Output: {OUTPUT_CSV}")
    print()

    # Merge SMILES to catalog
    with open(INPUT_CSV, 'r') as infile, open(OUTPUT_CSV, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['smiles']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            zinc_id = row['zinc_id']
            row['smiles'] = smiles_cache.get(zinc_id, '')
            writer.writerow(row)

    elapsed_total = time.time() - start_time

    # Summary
    print()
    print("=" * 60)
    print("API QUERY COMPLETE")
    print("=" * 60)
    print(f"Total molecules: {total_molecules:,}")
    print(f"SMILES found: {found_count:,} ({found_count*100/total_molecules:.2f}%)")
    print(f"SMILES not found: {not_found_count:,} ({not_found_count*100/total_molecules:.2f}%)")
    print(f"Time elapsed: {elapsed_total/3600:.1f} hours")
    print(f"Output saved to: {OUTPUT_CSV}")
    print("=" * 60)

if __name__ == "__main__":
    main()
