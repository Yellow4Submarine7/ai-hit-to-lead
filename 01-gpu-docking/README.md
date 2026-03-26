# Stage 01: Uni-Dock Shortlist Docking

This stage performs GPU docking with Uni-Dock and exports a shortlist contract for rescoring.

## Responsibilities

- extract ligands from stage-00 CSVs
- optionally extract a specific subset
- split large ligand manifests for parallel jobs
- run Uni-Dock with a box configuration and receptor
- export ranked shortlist tables and pose manifests

## Public interface

Key scripts:

- `scripts/extract_ligands.py`
- `scripts/extract_by_subset.py`
- `scripts/split_ligand_list.py`
- `scripts/run_docking.py`
- `scripts/export_results.py`

Shortlist output contract:

- `rank`
- `zinc_id`
- `unidock_affinity`
- `smiles`
- `subset`
- `pose_path`
- `source_pdbqt_path`

## Typical usage

```bash
uv run python scripts/run_docking.py \
  --receptor receptor.pdbqt \
  --ligand-index outputs/ligand_files.txt \
  --output-dir results \
  --box-config box_config.txt
```

```bash
uv run python scripts/export_results.py \
  --results-dir results \
  --smiles-csv ../00-data-extraction/outputs/zinc22_smiles.csv \
  --input-csv ../00-data-extraction/outputs/zinc22_smiles.csv \
  --output outputs/shortlist.csv \
  --manifest-output outputs/shortlist_manifest.txt
```

## Notes

- The public repo uses config files, CLI arguments, and PBS environment variables instead of source edits.
- The canonical public story for this stage is the 9M Uni-Dock shortlist workflow.
