# AI Hit-to-Lead

Public, config-driven hit-to-lead workflow for the WWP2-C2 program.

The canonical public pipeline is:

`ZINC22 acquisition -> 9M Uni-Dock shortlist docking -> GNINA/RF rescoring -> active learning`

This repository deliberately excludes raw HPC outputs, checkpoints, experiment artifacts,
and unrelated `aidd_public_docking_v1-v4` benchmark infrastructure.

## Repository layout

```text
00-data-extraction/   ZINC22 download, SMILES extraction, ligand extraction
01-gpu-docking/       Uni-Dock shortlist generation and export
02-rescoring/         GNINA parsing, RF-Score merging, score integration
03-active-learning/   Teacher training, pseudo-labeling, candidate ranking
docs/                 Decisions, pipeline notes, provenance
tests/                Smoke tests for the public stage contracts
```

## Quick start

1. Create a Python environment and install dependencies:

   ```bash
   uv venv
   uv pip install -r requirements.txt
   ```

2. Run the smoke-test suite:

   ```bash
   uv run pytest -q
   ```

3. Follow the stage READMEs in order:

- [00-data-extraction](./00-data-extraction/README.md)
- [01-gpu-docking](./01-gpu-docking/README.md)
- [02-rescoring](./02-rescoring/README.md)
- [03-active-learning](./03-active-learning/README.md)

## Key public decisions

- The public target provenance is `C2_shifted` / WWP2-C2.
- `6M0J` is treated as stale documentation metadata, not as the runtime target used for the 9M workflow.
- All runnable stages must be configured through CLI arguments, environment variables, or example config files rather than source edits.

See:

- [Public workflow decision](./docs/decisions/0001-public-wwp2-c2-workflow.md)
- [9M reference workflow](./docs/pipeline/9m-reference-workflow.md)
- [WWP2-C2 provenance note](./docs/provenance/wwp2-c2-target-note.md)

## Requirements

- Python 3.9+
- UV
- Open Babel
- RDKit
- scikit-learn
- XGBoost
- optional: GNINA and gninatorch for GPU rescoring and student inference

## License

MIT License
