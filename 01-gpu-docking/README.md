# GPU-Accelerated Molecular Docking

High-throughput virtual screening using Uni-Dock on NVIDIA GPUs, processing millions of molecules efficiently.

## Overview

This pipeline performs GPU-accelerated molecular docking using Uni-Dock:
1. Extract ligand PDBQT files from compressed archives
2. Configure docking box parameters
3. Run GPU-accelerated docking
4. Export results as sorted CSV with SMILES

### Key Features

- **GPU Acceleration**: ~16 molecules/second on NVIDIA H200
- **Parallel Processing**: Multi-core extraction and export
- **Progress Tracking**: Real-time progress with ETA estimates
- **Resume Capability**: Continue interrupted jobs seamlessly
- **Multi-key SMILES Lookup**: Handles format variations automatically

### Configuration

All PBS scripts use a configuration section at the top:
```bash
# ===== CONFIGURATION - MODIFY THESE =====
PROJECT_DIR="/path/to/your/project"
RECEPTOR="${PROJECT_DIR}/receptor.pdbqt"
# ========================================
```

Modify these paths to match your project structure.

### Performance Benchmarks

| GPU Model | Speed | Throughput |
|-----------|-------|------------|
| NVIDIA H200-141GB | ~16 mol/s | ~59K mol/hour |
| NVIDIA A100-40GB | ~12 mol/s | ~43K mol/hour |
| NVIDIA V100-32GB | ~8 mol/s | ~29K mol/hour |

**Example**: 9.2 million molecules on H200:
- Docking time: ~156 GPU hours (~6.5 days with 1 GPU)
- With 7 parallel jobs: ~1 day

## Prerequisites

### HPC Environment
- PBS/Torque job scheduler
- NVIDIA GPU with CUDA support
- 24-48 CPU cores for extraction/export

### Software Requirements
```bash
# Python 3.10+
python --version

# UV package manager
curl -LsSf https://astral.sh/uv/install.sh | sh

# Verify UV installation
uv --version
```

### Uni-Dock Installation
```bash
# Via conda (recommended)
conda create -n unidock -c conda-forge unidock
conda activate unidock

# Verify installation
unidock --version
```

## Installation

1. **Clone this repository**
```bash
git clone https://github.com/Yellow4Submarine7/ai-hit-to-lead.git
cd ai-hit-to-lead/01-gpu-docking
```

2. **Verify GPU and Uni-Dock**
```bash
bash utils/check_unidock.sh
```

## Usage

### Step 0: Prepare Input Data

Ensure you have outputs from Stage 00:
- `sampled.csv` - CSV with pdbqt_path column
- ZINC22 TGZ archives in data directory

### Step 1: Configure Docking Box

Create your box configuration file based on your receptor's binding site:

```bash
cp templates/box_config.txt.example box_config.txt
# Edit box_config.txt with your binding site coordinates
```

**Box configuration format**:
```
# Box center (from binding site analysis)
CENTER_X=-26.60
CENTER_Y=0.48
CENTER_Z=8.06

# Box dimensions (binding pocket + 4A margin)
SIZE_X=24.60
SIZE_Y=22.96
SIZE_Z=17.83
```

**How to determine box parameters**:
1. Open receptor in PyMOL/Chimera
2. Identify binding site residues
3. Calculate centroid as box center
4. Set dimensions to encompass pocket with ~4 Angstrom margin

### Step 2: Extract Ligand Files

1. **Configure the job**:
```bash
# Edit pbs/extract_ligands.pbs
# Set PROJECT_DIR, INPUT_CSV, DATA_ROOT, OUTPUT_DIR
```

2. **Submit extraction job**:
```bash
qsub pbs/extract_ligands.pbs
```

3. **Monitor progress**:
```bash
qstat -u $USER
tail -f extract_ligands.o*
```

4. **Job parameters**:
   - **CPU cores**: 48
   - **Memory**: 64 GB
   - **Walltime**: 8 hours
   - **Queue**: long

**Performance**: ~350 molecules/second with 48 cores

**Output**:
- `ligands/*.pdbqt` - Extracted ligand files
- `data/ligand_list.txt` - List of all ligand paths

### Step 3: Run GPU Docking

1. **Estimate runtime**:
```bash
bash utils/estimate_runtime.sh 1000000 h200
```

2. **Configure the job**:
```bash
# Edit pbs/dock_batch.pbs
# Set PROJECT_DIR, RECEPTOR, LIGAND_INDEX, BOX_CONFIG, OUTPUT_DIR
```

3. **Submit docking job**:
```bash
qsub pbs/dock_batch.pbs
```

4. **Monitor progress**:
```bash
qstat -u $USER
tail -f dock_batch.o*

# Check completed results
ls results/*_out.pdbqt | wc -l
```

5. **Job parameters**:
   - **GPU**: 1x H200 (or A100)
   - **CPU cores**: 24
   - **Memory**: 250 GB
   - **Walltime**: 240 hours (10 days)
   - **Queue**: gpu-h200

**Output**: `results/*_out.pdbqt` - Docked poses with scores

### Step 4: Export Results

1. **Configure the job**:
```bash
# Edit pbs/export_results.pbs
# Set PROJECT_DIR, RESULTS_DIR, SMILES_CSV, OUTPUT_CSV
```

2. **Submit export job**:
```bash
qsub pbs/export_results.pbs
```

3. **Job parameters**:
   - **CPU cores**: 48
   - **Memory**: 128 GB
   - **Walltime**: 4 hours
   - **Queue**: long

**Performance**: ~5,000 molecules/second with 48 cores

**Output**: `outputs/docking_results.csv` - Sorted by binding affinity

### Step 5: Verify Results

```bash
# Check result count
wc -l outputs/docking_results.csv

# Preview top hits
head -20 outputs/docking_results.csv

# Check SMILES coverage
grep -c '""$' outputs/docking_results.csv  # Empty SMILES count
```

## Output Format

### Docking Results CSV
```csv
zinc_id,score,smiles
ZINCrA0000001DcV.0.N,-10.25,CC(C)C1=CC=C(C=C1)...
ZINCsB000002XYZa.1.K,-9.87,COC1=CC(=CC(=C1O)...
```

**Columns**:
- `zinc_id`: ZINC identifier with conformer and variant
- `score`: Uni-Dock binding affinity (kcal/mol, lower is better)
- `smiles`: Canonical SMILES string

### PDBQT Output Files
Docking results are stored in `*_out.pdbqt` files with embedded scores:
```
REMARK VINA RESULT:    -8.5      0.000      0.000
REMARK  Name = ZINCrA0000001DcV.0.N
```

## File Organization

```
01-gpu-docking/
├── README.md                       # This documentation
├── scripts/
│   ├── extract_ligands.py          # Extract PDBQT from TGZ archives
│   ├── run_docking.py              # Uni-Dock wrapper with progress
│   ├── export_results.py           # Parse results to sorted CSV
│   └── split_ligand_list.py        # Split for parallel jobs
├── pbs/
│   ├── extract_ligands.pbs         # CPU job: ligand extraction
│   ├── dock_batch.pbs              # GPU job: molecular docking
│   └── export_results.pbs          # CPU job: result export
├── utils/
│   ├── check_unidock.sh            # Verify Uni-Dock installation
│   └── estimate_runtime.sh         # Estimate docking time
├── templates/
│   └── box_config.txt.example      # Example box configuration
└── .gitignore                      # Stage-specific ignores
```

## Advanced Usage

### Parallel Docking (Multiple GPUs)

For large datasets exceeding single-GPU walltime limits:

1. **Split ligand list**:
```bash
uv run python scripts/split_ligand_list.py \
    --input data/ligand_list.txt \
    --output-prefix data/ligand_batch \
    --chunks 7
```

2. **Submit multiple jobs**:
```bash
# Edit dock_batch.pbs for each batch
# Change LIGAND_INDEX to data/ligand_batch_01.txt, etc.
for i in 01 02 03 04 05 06 07; do
    qsub -v BATCH=${i} pbs/dock_batch.pbs
done
```

### Resume Interrupted Extraction

```bash
uv run python scripts/extract_ligands.py \
    --input-csv data/sampled.csv \
    --output-dir ligands \
    --data-root data/zinc22 \
    --resume
```

### Filter by Subset

Extract ligands from specific ZINC22 subsets:

```bash
uv run python scripts/extract_ligands.py \
    --input-csv data/sampled.csv \
    --output-dir ligands/zinc-22g \
    --data-root data/zinc22 \
    --subset-filter zinc-22g
```

### Export Top Hits Only

```bash
uv run python scripts/export_results.py \
    --results-dir results \
    --smiles-csv data/smiles.csv \
    --output top_10000.csv \
    --top-n 10000
```

### Custom Docking Parameters

```bash
uv run python scripts/run_docking.py \
    --receptor receptor.pdbqt \
    --ligand-index ligands.txt \
    --output-dir results \
    --box-config box_config.txt \
    --exhaustiveness 32 \
    --search-mode detail \
    --num-modes 20
```

**Search modes**:
- `fast`: Quick screening (~2x faster, less accurate)
- `balance`: Default, good balance of speed/accuracy
- `detail`: Thorough search (~2x slower, more accurate)

## Troubleshooting

### Issue: Uni-Dock not found
```bash
# Check installation
which unidock

# If using conda, ensure environment is activated
conda activate unidock

# Verify with check script
bash utils/check_unidock.sh
```

### Issue: GPU not detected
```bash
# Check NVIDIA driver
nvidia-smi

# Check CUDA
nvcc --version

# Ensure GPU queue is selected
#PBS -q gpu-h200
```

### Issue: Missing SMILES in output
- Check SMILES CSV has correct column names (`zinc_id`, `smiles`)
- Verify zinc_id format matches between files
- Export script uses multi-key fallback for format variations

### Issue: Low docking speed
- Verify GPU is being used: `nvidia-smi` during job
- Check search_mode setting (use `fast` for screening)
- Ensure sufficient memory for batch processing

### Issue: PDBQT parse errors
```bash
# Check file format
head -20 ligands/ZINCxxx.pdbqt

# Verify extraction worked correctly
ls ligands/*.pdbqt | wc -l
```

### Issue: Job exceeds walltime
- Split into multiple parallel jobs
- Use `estimate_runtime.sh` to calculate required time
- Consider using `fast` search mode for initial screening

## Queue Information

| Queue | GPU | Max Walltime | Max Jobs |
|-------|-----|--------------|----------|
| gpu-h200 | H200-141GB | 240h | 7 queued |
| gpu-h200-int | H200-141GB | 24h | 1 queued |
| gpu | A100-40GB | 240h | 1 queued |

**CPU Queues for extraction/export**:
| Queue | Max CPUs | Max Walltime |
|-------|----------|--------------|
| short | 24 | 48h |
| long | 50 | 240h |

## Next Steps

After completing docking, you will have:
- `outputs/docking_results.csv` - Sorted results with SMILES

These outputs are ready for the next stage: **02-hit-selection** (Active learning and hit prioritization).

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{gpu_docking_pipeline,
  title = {GPU-Accelerated Molecular Docking Pipeline},
  author = {Tensorlong},
  year = {2025},
  url = {https://github.com/Yellow4Submarine7/ai-hit-to-lead}
}
```

## References

- Uni-Dock: https://github.com/dptech-corp/Uni-Dock
- ZINC22 Database: https://zinc.docking.org/
- AutoDock Vina: https://vina.scripps.edu/
- UV Package Manager: https://docs.astral.sh/uv/

## License

MIT License - See LICENSE file in repository root
