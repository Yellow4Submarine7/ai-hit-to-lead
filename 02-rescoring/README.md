# Rescoring and Score Integration

This stage turns a Uni-Dock shortlist into a ranked rescoring table using:

- GNINA score-only rescoring
- RF-Score result collection
- score integration across Uni-Dock, GNINA, and RF-Score

## Public contract

Input:
- a shortlist CSV from `01-gpu-docking/scripts/export_results.py`
- pose files referenced by `pose_path`
- a receptor PDBQT file

Output:
- per-ligand GNINA text outputs
- RF-Score score files or a merged RF-Score CSV
- an integrated rescoring CSV with these columns:
  - `conformer_id`
  - `unidock_affinity`
  - `rfscore_v3`
  - `gnina_cnnscore`
  - `gnina_cnnaffinity`
  - `unidock_rank`
  - `rfscore_rank`
  - `gnina_rank`
  - `composite_rank`

## Scripts

- `scripts/run_gnina_rescoring.py`: run GNINA `--score_only` on shortlisted poses
- `scripts/run_rfscore_rescoring.py`: run an external RF-Score command template on shortlisted poses
- `scripts/merge_gnina_results.py`: parse GNINA text outputs into a CSV
- `scripts/merge_rfscore_results.py`: combine RF-Score outputs into a CSV
- `scripts/integrate_scores.py`: merge Uni-Dock, GNINA, and RF-Score into a single ranking table
