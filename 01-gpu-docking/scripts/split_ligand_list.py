#!/usr/bin/env python3
"""
Split Ligand List for Parallel Docking Jobs

Splits a large ligand list file into smaller chunks for parallel GPU jobs.
This is useful when docking millions of molecules across multiple GPUs.

Usage:
    uv run python split_ligand_list.py \
        --input ligand_list.txt \
        --output-prefix ligand_batch \
        --chunks 7
"""

import argparse
import math
import sys
from pathlib import Path


def count_lines(filepath: Path) -> int:
    """Count lines in a file efficiently."""
    count = 0
    with open(filepath, 'r') as f:
        for _ in f:
            count += 1
    return count


def split_file(input_path: Path, output_prefix: str, num_chunks: int) -> list:
    """Split input file into chunks."""
    total_lines = count_lines(input_path)
    lines_per_chunk = math.ceil(total_lines / num_chunks)

    print(f"Total lines: {total_lines:,}", file=sys.stderr)
    print(f"Chunks: {num_chunks}", file=sys.stderr)
    print(f"Lines per chunk: ~{lines_per_chunk:,}", file=sys.stderr)

    output_files = []
    current_chunk = 1
    current_lines = 0
    current_file = None

    with open(input_path, 'r') as infile:
        for line in infile:
            if current_file is None or current_lines >= lines_per_chunk:
                if current_file is not None:
                    current_file.close()
                    print(f"  Wrote {output_files[-1]}: {current_lines:,} lines", file=sys.stderr)

                if current_chunk <= num_chunks:
                    output_path = f"{output_prefix}_{current_chunk:02d}.txt"
                    output_files.append(output_path)
                    current_file = open(output_path, 'w')
                    current_chunk += 1
                    current_lines = 0

            current_file.write(line)
            current_lines += 1

    if current_file is not None:
        current_file.close()
        print(f"  Wrote {output_files[-1]}: {current_lines:,} lines", file=sys.stderr)

    return output_files


def split_by_size(input_path: Path, output_prefix: str, chunk_size: int) -> list:
    """Split input file by chunk size."""
    output_files = []
    current_chunk = 1
    current_lines = 0
    current_file = None

    with open(input_path, 'r') as infile:
        for line in infile:
            if current_file is None or current_lines >= chunk_size:
                if current_file is not None:
                    current_file.close()
                    print(f"  Wrote {output_files[-1]}: {current_lines:,} lines", file=sys.stderr)

                output_path = f"{output_prefix}_{current_chunk:02d}.txt"
                output_files.append(output_path)
                current_file = open(output_path, 'w')
                current_chunk += 1
                current_lines = 0

            current_file.write(line)
            current_lines += 1

    if current_file is not None:
        current_file.close()
        print(f"  Wrote {output_files[-1]}: {current_lines:,} lines", file=sys.stderr)

    return output_files


def main():
    parser = argparse.ArgumentParser(
        description='Split ligand list for parallel docking jobs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Split into 7 chunks (for gpu-h200 queue with 7 slots)
  uv run python split_ligand_list.py --input ligands.txt --output-prefix batch --chunks 7

  # Split into 100K molecule chunks
  uv run python split_ligand_list.py --input ligands.txt --output-prefix batch --chunk-size 100000
"""
    )
    parser.add_argument('--input', required=True, help='Input ligand list file')
    parser.add_argument('--output-prefix', required=True, help='Prefix for output files')
    parser.add_argument('--chunks', type=int, default=7, help='Number of chunks (default: 7)')
    parser.add_argument('--chunk-size', type=int, help='Lines per chunk (overrides --chunks)')

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Splitting {input_path}...", file=sys.stderr)

    if args.chunk_size:
        output_files = split_by_size(input_path, args.output_prefix, args.chunk_size)
    else:
        output_files = split_file(input_path, args.output_prefix, args.chunks)

    print(f"\nCreated {len(output_files)} files:", file=sys.stderr)
    for f in output_files:
        print(f"  {f}")


if __name__ == '__main__':
    main()
