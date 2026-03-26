# 9M Reference Workflow

This document captures the canonical public workflow for the WWP2-C2 program.

## Stage 00: ZINC22 acquisition and preprocessing

- Download the required ZINC22 AutoDock tranches.
- Extract SMILES and preserve source `pdbqt_path` references.
- Optionally sample or filter rows before ligand extraction.
- Extract ligand PDBQT files and build a ligand manifest.

Output contract:
- a CSV with `zinc_id`, `smiles`, and `pdbqt_path`
- an extracted ligand directory
- a ligand manifest text file

## Stage 01: Uni-Dock shortlist generation

- Dock extracted ligands against the configured receptor and box.
- Export ranked docking results from `*_out.pdbqt` outputs.
- Produce a shortlist CSV and optional pose manifest for rescoring.

Output contract:
- `rank`
- `zinc_id`
- `unidock_affinity`
- `smiles`
- `subset`
- `pose_path`
- `source_pdbqt_path`

## Stage 02: Rescoring and integration

- Run GNINA score-only rescoring on shortlisted poses.
- Run RF-Score through a configurable wrapper.
- Merge Uni-Dock, GNINA, and RF-Score into one integrated rescoring table.

Output contract:
- `conformer_id`
- `unidock_affinity`
- `rfscore_v3`
- `gnina_cnnscore`
- `gnina_cnnaffinity`
- `unidock_rank`
- `rfscore_rank`
- `gnina_rank`
- `composite_rank`

## Stage 03: Active learning

- Train the Teacher ensemble on wet-lab labels.
- Generate pseudo-labels for shortlisted candidates.
- Produce `.types` files for student training.
- Rank candidates with either precomputed predictions or gninatorch inference.

Output contract:
- teacher checkpoint
- teacher metrics CSV
- pseudo-label CSVs
- `train_K*.types` and `val.types`
- ranked candidate CSVs
