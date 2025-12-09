# ZINC22 Data Extraction

Download ZINC22 3D PDBQT files, extract SMILES strings, sample molecules, and prepare ligand files for docking.

## Overview

This pipeline downloads ZINC22 3D structures in PDBQT format and processes them through multiple steps:
1. Download .tgz archives from ZINC22
2. Test extraction on sample data
3. Extract SMILES strings from all PDBQT files
4. Sample a subset of molecules for docking
5. Extract PDBQT files for sampled molecules

### Configuration

All PBS scripts use a configuration variable at the top:
```bash
# ===== CONFIGURATION - MODIFY THESE =====
DATASET_NAME="zinc22"
# ========================================
```

Change `DATASET_NAME` to match your dataset. Output files will be named `${DATASET_NAME}_smiles.csv`, `${DATASET_NAME}_sampled.csv`, etc.

### Example Results (zinc22 H27-H29 P300-400)

- **2,141** .tgz archive files downloaded
- **9,297,704** molecules successfully extracted
- **100%** success rate
- **Output size**: 2.0 GB CSV file

## Prerequisites

### HPC Environment
- PBS/Torque job scheduler
- 16-48 CPU cores recommended
- 32-64 GB RAM minimum

### Software Requirements
```bash
# Python 3.10+
python --version

# UV package manager (https://docs.astral.sh/uv/)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Open Babel (for SMILES extraction)
uv pip install openbabel-wheel

# Pandas (for PDBQT extraction)
uv pip install pandas
```

## Installation

1. **Clone this repository**
```bash
git clone https://github.com/Yellow4Submarine7/ai-hit-to-lead.git
cd ai-hit-to-lead/00-data-extraction
```

2. **Install Python dependencies**
```bash
uv pip install -r ../requirements.txt
```

## Usage

### Step 1: Download ZINC22 Files

1. **Obtain wget script** from ZINC22 3D Tranches:
   - Visit https://cartblanche.docking.org/tranches/3d
   - Select your target tranches (e.g., H27-H29, P300-400)
   - Format: **AutoDock**
   - Method: **WGET**
   - Click "Download" to get the wget script

2. **Run download on login node** (compute nodes cannot access internet):
```bash
# Create data directory
mkdir -p data/zinc22

# Run wget script directly on login node
cd data/zinc22
bash /path/to/your/ZINC22-downloader.wget
```

**Note**: Download must be executed on the login node as PBS compute nodes do not have internet access.

Expected download time: ~2 hours
Expected download size: ~6.8 GB (2,141 files)

### Step 2: Test SMILES Extraction

Before processing all files, test on a sample dataset:

```bash
# Edit scripts/test_openbabel_parser.py to set paths
uv run python scripts/test_openbabel_parser.py
```

Expected output:
- Process 100 molecules from one .tgz file
- Success rate should be >90%
- Verify SMILES format is correct

### Step 3: Full SMILES Extraction

1. **Configure the job**:
```bash
# Edit pbs/extract_smiles.pbs
# Set DATASET_NAME at the top of the file
```

2. **Submit extraction job**:
```bash
qsub pbs/extract_smiles.pbs
```

3. **Monitor progress**:
```bash
qstat -u $USER
tail -f smiles_extraction.o*

# Check progress file (updates every 10 files)
cat logs/progress.json
```

4. **Job parameters**:
   - **CPU cores**: 16
   - **Memory**: 32 GB
   - **Walltime**: 12 hours (typically finishes in <10 minutes)
   - **Queue**: short

**Output**: `outputs/${DATASET_NAME}_smiles.csv`

### Step 4: Sample Molecules

Sample a subset of molecules for docking (e.g., 1 million):

1. **Configure the job**:
```bash
# Edit pbs/sample_molecules.pbs
# Set DATASET_NAME and SAMPLE_SIZE
```

2. **Submit sampling job**:
```bash
qsub pbs/sample_molecules.pbs
```

3. **Job parameters**:
   - **CPU cores**: 4
   - **Memory**: 16 GB
   - **Walltime**: 30 minutes
   - **Queue**: short

**Output**: `outputs/${DATASET_NAME}_sampled.csv`

### Step 5: Extract PDBQT Files

Extract 3D structure files for sampled molecules:

1. **Configure the job**:
```bash
# Edit pbs/extract_pdbqt.pbs
# Set DATASET_NAME and DATA_ROOT
```

2. **Submit extraction job**:
```bash
qsub pbs/extract_pdbqt.pbs
```

3. **Monitor progress**:
```bash
qstat -u $USER
tail -f logs/extract_pdbqt.log
```

4. **Job parameters**:
   - **CPU cores**: 48
   - **Memory**: 64 GB
   - **Walltime**: 3 hours
   - **Queue**: long

**Performance**: ~350 molecules/second with 48 cores (~47 minutes for 1M molecules)

**Output**:
- `ligands/*.pdbqt` - Individual PDBQT files
- `outputs/ligand_files.txt` - List of all extracted files (for docking)

### Step 6: Verify Results

```bash
# Check SMILES output
ls -lh outputs/${DATASET_NAME}_smiles.csv
head -5 outputs/${DATASET_NAME}_smiles.csv

# Check sampled output
wc -l outputs/${DATASET_NAME}_sampled.csv

# Check extracted PDBQT files
ls ligands/ | wc -l
wc -l outputs/ligand_files.txt
```

## Output Format

### SMILES CSV Structure
```csv
zinc_id,smiles,pdbqt_path
ZINCrB000001FZto.0,C1(=[C]C(=[C][C]=C1OC(F)(F)F)...,data/zinc22/H27P310-K-aaaaaa.pdbqt.tgz:H27/...
```

**Columns**:
- `zinc_id`: ZINC database identifier with conformer number
- `smiles`: Canonical SMILES string
- `pdbqt_path`: Path to 3D structure file (tgz_path:internal_path)

### Ligand Files List
```
./ligands/ZINCrB000001FZto.pdbqt
./ligands/ZINCrB000002ABcd.pdbqt
...
```

This file is used as input for Uni-Dock's `--ligand_index` parameter.

## File Organization

```
00-data-extraction/
├── scripts/
│   ├── extract_smiles_openbabel.py    # SMILES extraction (Step 3)
│   ├── extract_pdbqt.py               # PDBQT extraction (Step 5)
│   ├── test_openbabel_parser.py       # Testing/validation (Step 2)
│   ├── extract_zinc22_catalog.py      # Catalog generation
│   ├── add_smiles_via_api.py          # API alternative
│   └── add_smiles_from_files.py       # File-based alternative
├── pbs/
│   ├── extract_smiles.pbs             # SMILES extraction job (Step 3)
│   ├── sample_molecules.pbs           # Sampling job (Step 4)
│   └── extract_pdbqt.pbs              # PDBQT extraction job (Step 5)
└── utils/
    ├── run_wget_fixed.sh              # Wget execution helper
    └── zinc_retry_download.sh         # Retry failed downloads
```

## Advanced Usage

### Resume Interrupted SMILES Extraction

```bash
# The --resume flag will skip already processed files
uv run python scripts/extract_smiles_openbabel.py \
  --input_dir ./data \
  --output_csv ./outputs/${DATASET_NAME}_smiles.csv \
  --error_log ./logs/errors.log \
  --progress_file ./logs/progress.json \
  --ncpus 16 \
  --resume
```

### Retry Failed Downloads

```bash
# Use the retry download script on login node
./utils/zinc_retry_download.sh /path/to/failed_urls.txt
```

### Alternative SMILES Methods

#### API-based SMILES Extraction
```bash
uv run python scripts/add_smiles_via_api.py
```
Note: Slower but doesn't require downloading 3D structures

#### Extract from 2D .smi Files
```bash
uv run python scripts/add_smiles_from_files.py
```
Note: Requires downloading .smi.gz files separately

## Troubleshooting

### Issue: Open Babel not found
```bash
# Solution: Install via UV
uv pip install openbabel-wheel

# Verify installation
uv run python -c "from openbabel import pybel; print('OK')"
```

### Issue: PBS job pending
```bash
# Check queue status
qstat -Q

# Check your job details
qstat -f <job_id>
```

### Issue: Low success rate (<90%)
- Check error log: `cat logs/errors.log`
- Common causes:
  - Corrupted .tgz files: Re-download affected files
  - Insufficient memory: Increase memory in PBS script
  - Invalid PDBQT format: Report to ZINC database

### Issue: Slow PDBQT extraction
- Increase CPU cores (up to 50 in long queue)
- Use faster storage (SSD vs HDD)
- Check disk I/O bottleneck

## Performance Benchmarks

| Step | Speed | Resources |
|------|-------|-----------|
| Download | ~50 MB/min | Network dependent |
| SMILES extraction | ~1M mol/min | 16 cores |
| Sampling | ~1M mol/min | 4 cores |
| PDBQT extraction | ~350 mol/s | 48 cores |

## Next Steps

After completing all steps, you will have:
- `outputs/${DATASET_NAME}_sampled.csv` - Sampled molecules with SMILES
- `ligands/*.pdbqt` - 3D structure files
- `outputs/ligand_files.txt` - Ligand list for docking

These outputs are ready for the next stage: **01-docking** (Uni-Dock molecular docking).

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{zinc22_extraction,
  title = {ZINC22 Data Extraction Pipeline},
  author = {Tensorlong},
  year = {2025},
  url = {https://github.com/Yellow4Submarine7/ai-hit-to-lead}
}
```

## References

- ZINC22 Database: https://zinc.docking.org/
- Open Babel: https://openbabel.org/
- UV Package Manager: https://docs.astral.sh/uv/

## License

MIT License - See LICENSE file in repository root
