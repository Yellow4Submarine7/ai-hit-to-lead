# ZINC22 Data Extraction

Download ZINC22 3D PDBQT files and extract SMILES strings for drug discovery applications.

## Overview

This pipeline downloads ZINC22 3D structures in PDBQT format and converts them to SMILES strings using Open Babel. The dataset covers molecules with:
- **Heavy atoms**: 27-29 (H27-H29)
- **LogP**: 300-400 (P300-400)

### Results

- **2,141** .tgz archive files downloaded
- **9,297,704** molecules successfully extracted
- **100%** success rate
- **Output size**: 2.0 GB CSV file

## Prerequisites

### HPC Environment
- PBS/Torque job scheduler
- 16 CPU cores recommended
- 32 GB RAM minimum

### Software Requirements
```bash
# Python 3.10+
python --version

# UV package manager (https://docs.astral.sh/uv/)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Open Babel
uv pip install openbabel-wheel
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

1. **Obtain wget script** from ZINC22 database:
   - Visit https://zinc.docking.org/
   - Navigate to 3D structures section
   - Generate wget script for H27-H29, P300-400 range

2. **Submit download job**:
```bash
# Edit pbs/download_zinc22.pbs to set correct paths
qsub pbs/download_zinc22.pbs
```

3. **Monitor download**:
```bash
qstat -u $USER
tail -f download_zinc22.o*
```

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

1. **Submit extraction job**:
```bash
# Edit pbs/extract_smiles.pbs to set correct paths
qsub pbs/extract_smiles.pbs
```

2. **Monitor progress**:
```bash
qstat -u $USER
tail -f extract_smiles.o*

# Check progress file (updates every 10 files)
cat logs/progress.json
```

3. **Job parameters**:
   - **CPU cores**: 16
   - **Memory**: 32 GB
   - **Walltime**: 12 hours (typically finishes in <10 minutes)
   - **Queue**: short

### Step 4: Verify Results

```bash
# Check output file
ls -lh outputs/zinc22_H27-H29_P300-400_smiles.csv

# View summary
cat outputs/zinc22_H27-H29_P300-400_smiles_summary.txt

# Check first few rows
head -10 outputs/zinc22_H27-H29_P300-400_smiles.csv
```

## Output Format

### CSV Structure
```csv
zinc_id,smiles,pdbqt_path
ZINCrB000001FZto.0,C1(=[C]C(=[C][C]=C1OC(F)(F)F)...,zinc22_H27-29.../H27P310-K-aaaaaa.pdbqt.tgz:H27/...
```

**Columns**:
- `zinc_id`: ZINC database identifier with conformer number
- `smiles`: Canonical SMILES string
- `pdbqt_path`: Path to 3D structure file (.tgz:internal_path)

### Statistics File
The summary file includes:
- Total .tgz files processed
- Total molecules extracted
- Success rate percentage
- Error count
- Processing timestamp

## Advanced Usage

### Resume Interrupted Job

If the extraction job is interrupted, you can resume:

```bash
# The --resume flag will skip already processed files
uv run python scripts/extract_smiles_openbabel.py \
  --input_dir /path/to/zinc22_H27-29_P300-400 \
  --output_csv outputs/zinc22_smiles.csv \
  --error_log logs/errors.log \
  --progress_file logs/progress.json \
  --ncpus 16 \
  --resume
```

### Retry Failed Downloads

If some files failed to download:

```bash
# Generate retry script
./utils/zinc_retry_download.sh

# Submit retry job
qsub pbs/retry_download.pbs
```

### Alternative Methods

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

### Issue: Slow processing
- Increase CPU cores in PBS script (up to 32)
- Use faster storage (SSD vs HDD)
- Split processing across multiple nodes

## File Organization

```
00-data-extraction/
├── scripts/
│   ├── extract_smiles_openbabel.py    # Main extraction script
│   ├── test_openbabel_parser.py        # Testing/validation
│   ├── extract_zinc22_catalog.py       # Catalog generation
│   ├── add_smiles_via_api.py           # API alternative
│   └── add_smiles_from_files.py        # File-based alternative
├── pbs/
│   ├── download_zinc22.pbs             # Download job script
│   ├── extract_smiles.pbs              # Extraction job script
│   └── retry_download.pbs              # Retry failed downloads
└── utils/
    ├── run_wget_fixed.sh               # Wget execution helper
    └── zinc_retry_download.sh          # Generate retry commands
```

## Performance Benchmarks

- **Download speed**: ~50 MB/min (depends on network)
- **Extraction speed**: ~1,000,000 molecules per minute (16 cores)
- **Memory usage**: ~2 GB per worker process
- **Disk I/O**: Sequential read, minimal random access

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{zinc22_extraction,
  title = {ZINC22 SMILES Extraction Pipeline},
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
