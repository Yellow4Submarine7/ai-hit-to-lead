#!/bin/bash
# Check Uni-Dock installation and GPU availability

echo "=========================================="
echo "Uni-Dock Installation Check"
echo "=========================================="

# Check if unidock is in PATH
if command -v unidock &> /dev/null; then
    UNIDOCK_PATH=$(which unidock)
    echo "[OK] Uni-Dock found: $UNIDOCK_PATH"
    echo ""
    echo "Version info:"
    unidock --version 2>&1 || echo "(version info not available)"
else
    echo "[WARNING] Uni-Dock not found in PATH"
    echo ""
    echo "Common installation locations:"
    echo "  - ~/.conda/envs/unidock/bin/unidock"
    echo "  - /usr/local/bin/unidock"
    echo ""
    echo "To install via conda:"
    echo "  conda create -n unidock -c conda-forge unidock"
    echo "  conda activate unidock"
fi

echo ""
echo "=========================================="
echo "GPU Check"
echo "=========================================="

# Check for NVIDIA GPU
if command -v nvidia-smi &> /dev/null; then
    echo "[OK] NVIDIA driver found"
    echo ""
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv
else
    echo "[WARNING] nvidia-smi not found"
    echo "GPU docking requires NVIDIA GPU with CUDA support"
fi

echo ""
echo "=========================================="
echo "CUDA Check"
echo "=========================================="

if command -v nvcc &> /dev/null; then
    echo "[OK] CUDA compiler found"
    nvcc --version | head -4
elif [ -d "/usr/local/cuda" ]; then
    echo "[OK] CUDA installation found at /usr/local/cuda"
    ls /usr/local/cuda/version.txt 2>/dev/null && cat /usr/local/cuda/version.txt
else
    echo "[WARNING] CUDA not found in standard locations"
fi

echo ""
echo "=========================================="
