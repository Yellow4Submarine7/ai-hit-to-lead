#!/bin/bash

# ZINC22 Fixed Download Script - Execute wget commands line by line
# Date: 2025-11-21

# Configuration - use relative paths
# Run from the 00-data-extraction directory
BASE_DIR="."
DOWNLOAD_DIR="$BASE_DIR/data/zinc22"
LOG_DIR="$BASE_DIR/logs"
WGET_SCRIPT="${1:-./ZINC22-downloader.wget}"

mkdir -p "$DOWNLOAD_DIR"
mkdir -p "$LOG_DIR"
cd "$DOWNLOAD_DIR"

# Initialize logs
echo "=== ZINC22 Fixed Download Started at $(date) ===" > "$LOG_DIR/wget_fixed_progress.log"
echo "Total commands to execute: $(wc -l < "$WGET_SCRIPT")" >> "$LOG_DIR/wget_fixed_progress.log"
echo "" > "$LOG_DIR/wget_fixed_success.log"
echo "" > "$LOG_DIR/wget_fixed_failed.log"

# Counters
TOTAL_LINES=$(wc -l < "$WGET_SCRIPT")
CURRENT=0
SUCCESS_COUNT=0
FAILED_COUNT=0

echo "============================================"
echo "ZINC22 Fixed Download - Retry missing files"
echo "============================================"
echo "Total commands: $TOTAL_LINES"
echo "Download directory: $DOWNLOAD_DIR"
echo "Log directory: $LOG_DIR"
echo ""
echo "Starting download..."
echo ""

# Read and execute each wget command
while IFS= read -r wget_cmd; do
    CURRENT=$((CURRENT + 1))

    # Skip empty lines
    if [ -z "$wget_cmd" ]; then
        continue
    fi

    echo "[$CURRENT/$TOTAL_LINES] Executing..."

    # Execute wget command directly
    if eval "$wget_cmd" 2>&1 | tee -a "$LOG_DIR/wget_fixed_detailed.log"; then
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        echo "[SUCCESS] $wget_cmd" >> "$LOG_DIR/wget_fixed_success.log"
        echo "[$CURRENT/$TOTAL_LINES] SUCCESS" >> "$LOG_DIR/wget_fixed_progress.log"
    else
        FAILED_COUNT=$((FAILED_COUNT + 1))
        echo "[FAILED $(date)]" >> "$LOG_DIR/wget_fixed_failed.log"
        echo "Command: $wget_cmd" >> "$LOG_DIR/wget_fixed_failed.log"
        echo "---" >> "$LOG_DIR/wget_fixed_failed.log"
        echo "[$CURRENT/$TOTAL_LINES] FAILED" >> "$LOG_DIR/wget_fixed_progress.log"
        echo "  FAILED - Check logs for details"
    fi

    # Print progress every 50 commands
    if [ $((CURRENT % 50)) -eq 0 ]; then
        echo ""
        echo "=========================================="
        echo "Progress: $CURRENT/$TOTAL_LINES ($((CURRENT * 100 / TOTAL_LINES))%)"
        echo "Success: $SUCCESS_COUNT, Failed: $FAILED_COUNT"
        echo "=========================================="
        echo ""
    fi

done < "$WGET_SCRIPT"

# Generate final summary
echo "" | tee -a "$LOG_DIR/wget_fixed_progress.log"
echo "============================================" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "=== Download Complete at $(date) ===" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "============================================" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
echo "Total commands: $TOTAL_LINES" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "Successful: $SUCCESS_COUNT" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "Failed: $FAILED_COUNT" | tee -a "$LOG_DIR/wget_fixed_progress.log" "$LOG_DIR/wget_fixed_summary.txt"
echo "" | tee -a "$LOG_DIR/wget_fixed_summary.txt"

# File statistics
echo "=== Downloaded Files ===" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
TGZ_COUNT=$(find "$DOWNLOAD_DIR" -name "*.tgz" -type f | wc -l)
echo "TGZ files: $TGZ_COUNT" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
du -sh "$DOWNLOAD_DIR" | awk '{print "TGZ total size: " $1}' | tee -a "$LOG_DIR/wget_fixed_summary.txt"
echo "" | tee -a "$LOG_DIR/wget_fixed_summary.txt"

if [ $FAILED_COUNT -gt 0 ]; then
    echo "=== Failed Downloads ===" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
    echo "See wget_fixed_failed.log for retry commands" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
    echo "" | tee -a "$LOG_DIR/wget_fixed_summary.txt"
fi

echo "Download process completed!"
echo "Check logs at: $LOG_DIR/"
