#!/usr/bin/env python3
"""
ZINC22 PDBQT to SMILES Conversion Test Script using Open Babel
Uses Open Babel's native PDBQT support
"""

import tarfile
import tempfile
import os
from pathlib import Path
import csv
from openbabel import pybel
import sys

def extract_zinc_id_from_pdbqt(content):
    """Extract ZINC ID from REMARK line in PDBQT file"""
    for line in content.split('\n'):
        if line.startswith('REMARK') and 'Name' in line:
            # Format: REMARK  Name = ZINCtI000009ADZy
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
        # Extract ZINC ID
        with open(pdbqt_path, 'r') as f:
            content = f.read()
        zinc_id = extract_zinc_id_from_pdbqt(content)
        if not zinc_id:
            zinc_id = Path(pdbqt_path).stem

        # Read PDBQT using Open Babel (native support!)
        mol = next(pybel.readfile("pdbqt", str(pdbqt_path)))

        # Generate SMILES (hydrogen atoms removed by default)
        smiles_output = mol.write("smi").strip()
        # SMILES output format: "SMILES molecule_name", extract first part only
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

def process_tgz_file(tgz_path, output_csv, max_molecules=None):
    """
    Extract and process all PDBQT files from a .tgz archive
    """
    results = []
    errors = []

    print(f"Processing {tgz_path}...")

    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract .tgz file
        try:
            with tarfile.open(tgz_path, 'r:gz') as tar:
                tar.extractall(tmpdir)
        except Exception as e:
            print(f"ERROR: Failed to extract {tgz_path}: {e}")
            return results, errors

        # Find all .pdbqt files
        pdbqt_files = list(Path(tmpdir).rglob('*.pdbqt'))
        print(f"Found {len(pdbqt_files)} PDBQT files")

        # Process each PDBQT file
        processed = 0
        for pdbqt_file in pdbqt_files:
            if max_molecules and processed >= max_molecules:
                break

            smiles, zinc_id, error = pdbqt_to_smiles(str(pdbqt_file))

            if error:
                errors.append({
                    'zinc_id': zinc_id or 'UNKNOWN',
                    'file': str(pdbqt_file.name),
                    'error': error
                })
            else:
                # Store relative path from the original .tgz structure
                rel_path = pdbqt_file.relative_to(tmpdir)
                results.append({
                    'zinc_id': zinc_id,
                    'smiles': smiles,
                    'pdbqt_path': str(rel_path)
                })

            processed += 1
            if processed % 100 == 0:
                print(f"  Processed {processed}/{len(pdbqt_files)} molecules...")

    # Write results to CSV
    if results:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['zinc_id', 'smiles', 'pdbqt_path'])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSUCCESS: Wrote {len(results)} molecules to {output_csv}")
    else:
        print("\nWARNING: No valid molecules found!")

    # Report errors
    if errors:
        print(f"\nERRORS: {len(errors)} molecules failed to process")
        print("First 10 errors:")
        for err in errors[:10]:
            print(f"  {err['zinc_id']}: {err['error']}")

    return results, errors

def main():
    # Configuration
    test_tgz = "/data/petretto/home/tiesunlong/zinc22_H27-29_P300-400/zinc22/zinc-22s/H29/H29P380/a/H29P380-N-saj.pdbqt.tgz"
    output_csv = "/data/petretto/home/tiesunlong/zinc22_H27-29_P300-400/outputs/test_openbabel.csv"

    print("="*60)
    print("ZINC22 PDBQT to SMILES - Open Babel Version")
    print("="*60)
    print(f"Input: {test_tgz}")
    print(f"Output: {output_csv}")
    print()

    # Check if test file exists
    if not os.path.exists(test_tgz):
        print(f"ERROR: Test file not found: {test_tgz}")
        print("\nSearching for alternative test file...")
        # Try to find any .tgz file
        base_dir = Path("/data/petretto/home/tiesunlong/zinc22_H27-29_P300-400")
        tgz_files = list(base_dir.rglob("*.tgz"))
        if tgz_files:
            test_tgz = str(tgz_files[0])
            print(f"Using: {test_tgz}")
        else:
            print("ERROR: No .tgz files found!")
            sys.exit(1)

    # Process the file (test first 100 molecules)
    results, errors = process_tgz_file(test_tgz, output_csv, max_molecules=100)

    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    print(f"Success: {len(results)}")
    print(f"Failed: {len(errors)}")
    if results:
        success_rate = len(results) / (len(results) + len(errors)) * 100
        print(f"Success rate: {success_rate:.1f}%")

        # Show first few results
        print("\nFirst 5 results:")
        for i, r in enumerate(results[:5]):
            print(f"  {i+1}. {r['zinc_id']}: {r['smiles'][:50]}...")

    # Exit code
    if len(results) > 0:
        print("\n✅ TEST PASSED: Successfully extracted SMILES with Open Babel")
        sys.exit(0)
    else:
        print("\n❌ TEST FAILED: No SMILES extracted")
        sys.exit(1)

if __name__ == "__main__":
    main()
