#!/bin/bash
# Estimate docking runtime based on molecule count and GPU type

usage() {
    echo "Usage: $0 <molecule_count> [gpu_type]"
    echo ""
    echo "Arguments:"
    echo "  molecule_count  Number of molecules to dock"
    echo "  gpu_type        GPU model (h200|a100|v100) [default: h200]"
    echo ""
    echo "Example:"
    echo "  $0 1000000        # 1M molecules on H200"
    echo "  $0 500000 a100    # 500K molecules on A100"
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

MOLECULE_COUNT=$1
GPU_TYPE=${2:-h200}

# Speed estimates (molecules per second)
case $GPU_TYPE in
    h200)
        SPEED=16
        GPU_NAME="NVIDIA H200-141GB"
        ;;
    a100)
        SPEED=12
        GPU_NAME="NVIDIA A100-40GB"
        ;;
    v100)
        SPEED=8
        GPU_NAME="NVIDIA V100-32GB"
        ;;
    *)
        echo "Unknown GPU type: $GPU_TYPE"
        echo "Supported: h200, a100, v100"
        exit 1
        ;;
esac

# Calculate times
SECONDS_TOTAL=$((MOLECULE_COUNT / SPEED))
MINUTES_TOTAL=$((SECONDS_TOTAL / 60))
HOURS_TOTAL=$((MINUTES_TOTAL / 60))
DAYS_TOTAL=$(echo "scale=2; $HOURS_TOTAL / 24" | bc)

# Format readable time
if [ $HOURS_TOTAL -ge 24 ]; then
    TIME_READABLE="${DAYS_TOTAL} days"
elif [ $HOURS_TOTAL -ge 1 ]; then
    TIME_READABLE="${HOURS_TOTAL} hours"
else
    TIME_READABLE="${MINUTES_TOTAL} minutes"
fi

echo "=========================================="
echo "Docking Runtime Estimate"
echo "=========================================="
echo "Molecules:     $(printf "%'d" $MOLECULE_COUNT)"
echo "GPU:           $GPU_NAME"
echo "Speed:         ~${SPEED} mol/sec"
echo "=========================================="
echo "Estimated time: $TIME_READABLE"
echo "  - Seconds:    $(printf "%'d" $SECONDS_TOTAL)"
echo "  - Minutes:    $(printf "%'d" $MINUTES_TOTAL)"
echo "  - Hours:      $(printf "%'d" $HOURS_TOTAL)"
echo "=========================================="
echo ""
echo "Parallel GPU recommendations:"
if [ $HOURS_TOTAL -gt 240 ]; then
    GPUS_NEEDED=$((HOURS_TOTAL / 240 + 1))
    echo "  - Exceeds 10-day walltime limit"
    echo "  - Recommend splitting into $GPUS_NEEDED parallel jobs"
    BATCH_SIZE=$((MOLECULE_COUNT / GPUS_NEEDED))
    echo "  - Batch size: ~$(printf "%'d" $BATCH_SIZE) per job"
elif [ $HOURS_TOTAL -gt 48 ]; then
    echo "  - Consider 2 parallel jobs for faster completion"
else
    echo "  - Single GPU job is sufficient"
fi
echo "=========================================="
