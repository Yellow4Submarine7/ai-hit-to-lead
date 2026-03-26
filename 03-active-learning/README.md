# Active Learning

This stage packages the reusable active-learning components from the WWP2-C2 workflow.

## Public contract

Inputs:
- an experimental dataset with `mcule_id`, `smiles_canonical`, `binding_level_nm`, `scaffold`, and `pdbqt_path`
- a ranked candidate CSV with `ligand_id`, `smiles_canonical`, and `pdbqt_path`
- an optional integrated rescoring table from `02-rescoring`

Outputs:
- a trained Teacher ensemble checkpoint
- teacher evaluation metrics
- pseudo-label prediction tables
- `.types` files for student training
- ranked candidate tables for downstream wet-lab selection

## Scripts

- `scripts/train_teacher.py`: train the Teacher ensemble with multi-seed scaffold splits
- `scripts/generate_pseudo_labels.py`: build teacher predictions and `.types` files
- `scripts/score_candidates.py`: rank candidates from precomputed predictions or gninatorch inference
- `scripts/evaluate_student_rescore.py`: evaluate student checkpoints on a validation `.types` file
