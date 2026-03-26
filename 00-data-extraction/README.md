# Stage 00: ZINC22 Data Extraction

This stage prepares ZINC22 ligands for downstream docking.

## Responsibilities

- download ZINC22 AutoDock archives on the login node
- extract SMILES from PDBQT structures
- preserve source `pdbqt_path` references
- optionally sample or filter rows
- extract ligand PDBQT files and build a ligand manifest

## Public interface

Key scripts:

- `scripts/extract_smiles_openbabel.py`
- `scripts/extract_pdbqt.py`
- `scripts/test_openbabel_parser.py`

Key outputs:

- stage CSV with `zinc_id`, `smiles`, `pdbqt_path`
- extracted ligand directory
- ligand manifest for Uni-Dock

## Typical usage

```bash
uv run python scripts/extract_smiles_openbabel.py \
  --input_dir ./data \
  --output_csv ./outputs/zinc22_smiles.csv \
  --error_log ./logs/errors.log \
  --progress_file ./logs/progress.json \
  --ncpus 16 \
  --resume
```

```bash
uv run python scripts/extract_pdbqt.py \
  --input ./outputs/zinc22_smiles.csv \
  --output_dir ./ligands \
  --data_root ./data \
  --ligand_list ./outputs/ligand_files.txt \
  --ncpus 24
```

## Notes

- All paths are supplied through CLI flags or PBS environment variables.
- The public repo does not ship raw downloaded ZINC22 archives.
- Vendor enrichment and tranche-specific filtering are optional follow-up steps, not hardcoded defaults.
