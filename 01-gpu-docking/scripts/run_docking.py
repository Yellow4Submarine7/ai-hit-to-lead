#!/usr/bin/env python3
"""
Uni-Dock Wrapper with Progress Tracking

Python wrapper for GPU-accelerated molecular docking with Uni-Dock.
Provides progress monitoring, performance logging, and graceful error handling.

Usage:
    uv run python run_docking.py \
        --receptor receptor.pdbqt \
        --ligand-index ligands.txt \
        --output-dir results \
        --box-config box_config.txt

Features:
    - Real-time progress tracking by counting output files
    - Performance metrics (molecules/second, ETA)
    - Input validation before docking
    - Detailed logging with timestamps
    - Graceful interrupt handling
"""

import argparse
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional


class DockingRunner:
    """Manages Uni-Dock execution with progress tracking."""

    def __init__(
        self,
        receptor: Path,
        ligand_index: Path,
        output_dir: Path,
        box_config: Path,
        unidock_path: str = 'unidock',
        num_modes: int = 9,
        exhaustiveness: int = 8,
        search_mode: str = 'balance',
        log_file: Optional[Path] = None,
        progress_interval: int = 60
    ):
        self.receptor = receptor
        self.ligand_index = ligand_index
        self.output_dir = output_dir
        self.box_config = box_config
        self.unidock_path = unidock_path
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness
        self.search_mode = search_mode
        self.log_file = log_file
        self.progress_interval = progress_interval

        self.box_params = {}
        self.total_ligands = 0
        self.start_time = None
        self.process = None
        self.interrupted = False

    def log(self, message: str):
        """Log message to stderr and optionally to file."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_line = f"[{timestamp}] {message}"
        print(log_line, file=sys.stderr)

        if self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(log_line + '\n')

    def parse_box_config(self) -> bool:
        """Parse box configuration file."""
        required_keys = ['CENTER_X', 'CENTER_Y', 'CENTER_Z',
                         'SIZE_X', 'SIZE_Y', 'SIZE_Z']

        try:
            with open(self.box_config, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue

                    if '=' in line:
                        key, value = line.split('=', 1)
                        key = key.strip()
                        value = value.strip()

                        if key in required_keys:
                            self.box_params[key] = float(value)

            missing = set(required_keys) - set(self.box_params.keys())
            if missing:
                self.log(f"Error: Missing box parameters: {missing}")
                return False

            return True

        except Exception as e:
            self.log(f"Error parsing box config: {e}")
            return False

    def count_ligands(self) -> int:
        """Count total ligands in index file."""
        count = 0
        with open(self.ligand_index, 'r') as f:
            for line in f:
                if line.strip():
                    count += 1
        return count

    def count_output_files(self) -> int:
        """Count completed docking output files."""
        return len(list(self.output_dir.glob('*_out.pdbqt')))

    def validate_inputs(self) -> bool:
        """Validate all inputs before running."""
        errors = []

        if not self.receptor.exists():
            errors.append(f"Receptor not found: {self.receptor}")

        if not self.ligand_index.exists():
            errors.append(f"Ligand index not found: {self.ligand_index}")

        if not self.box_config.exists():
            errors.append(f"Box config not found: {self.box_config}")

        # Check Uni-Dock installation
        try:
            result = subprocess.run(
                [self.unidock_path, '--version'],
                capture_output=True,
                text=True
            )
            if result.returncode != 0:
                errors.append(f"Uni-Dock check failed: {result.stderr}")
        except FileNotFoundError:
            errors.append(f"Uni-Dock not found: {self.unidock_path}")

        if errors:
            for error in errors:
                self.log(f"Validation Error: {error}")
            return False

        return True

    def build_command(self) -> list:
        """Build Uni-Dock command line."""
        cmd = [
            self.unidock_path,
            '--receptor', str(self.receptor),
            '--ligand_index', str(self.ligand_index),
            '--dir', str(self.output_dir),
            '--center_x', str(self.box_params['CENTER_X']),
            '--center_y', str(self.box_params['CENTER_Y']),
            '--center_z', str(self.box_params['CENTER_Z']),
            '--size_x', str(self.box_params['SIZE_X']),
            '--size_y', str(self.box_params['SIZE_Y']),
            '--size_z', str(self.box_params['SIZE_Z']),
            '--num_modes', str(self.num_modes),
            '--exhaustiveness', str(self.exhaustiveness),
            '--search_mode', self.search_mode
        ]
        return cmd

    def format_duration(self, seconds: float) -> str:
        """Format duration as human-readable string."""
        td = timedelta(seconds=int(seconds))
        hours, remainder = divmod(td.seconds, 3600)
        minutes, seconds = divmod(remainder, 60)

        if td.days > 0:
            return f"{td.days}d {hours:02d}:{minutes:02d}:{seconds:02d}"
        else:
            return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

    def monitor_progress(self):
        """Monitor docking progress in background."""
        last_count = self.count_output_files()
        last_time = time.time()

        while self.process and self.process.poll() is None and not self.interrupted:
            time.sleep(self.progress_interval)

            current_count = self.count_output_files()
            current_time = time.time()

            # Calculate metrics
            elapsed = current_time - self.start_time
            delta_count = current_count - last_count
            delta_time = current_time - last_time

            if delta_time > 0:
                rate = delta_count / delta_time
            else:
                rate = 0

            # Overall rate
            if elapsed > 0:
                overall_rate = current_count / elapsed
            else:
                overall_rate = 0

            # ETA
            remaining = self.total_ligands - current_count
            if overall_rate > 0:
                eta_seconds = remaining / overall_rate
                eta_str = self.format_duration(eta_seconds)
            else:
                eta_str = "calculating..."

            # Progress percentage
            if self.total_ligands > 0:
                pct = 100.0 * current_count / self.total_ligands
            else:
                pct = 0

            self.log(
                f"Progress: {current_count:,}/{self.total_ligands:,} ({pct:.1f}%) | "
                f"Rate: {rate:.1f} mol/s | ETA: {eta_str}"
            )

            last_count = current_count
            last_time = current_time

    def handle_interrupt(self, signum, frame):
        """Handle interrupt signal gracefully."""
        self.log("Interrupt received, stopping docking...")
        self.interrupted = True
        if self.process:
            self.process.terminate()

    def run(self) -> int:
        """Run the docking workflow."""
        # Validate inputs
        self.log("Validating inputs...")
        if not self.validate_inputs():
            return 1

        # Parse box config
        self.log(f"Parsing box config: {self.box_config}")
        if not self.parse_box_config():
            return 1

        self.log(f"  Center: ({self.box_params['CENTER_X']}, "
                 f"{self.box_params['CENTER_Y']}, {self.box_params['CENTER_Z']})")
        self.log(f"  Size: ({self.box_params['SIZE_X']}, "
                 f"{self.box_params['SIZE_Y']}, {self.box_params['SIZE_Z']})")

        # Count ligands
        self.total_ligands = self.count_ligands()
        self.log(f"Total ligands: {self.total_ligands:,}")

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Check for existing results (resume capability)
        existing = self.count_output_files()
        if existing > 0:
            self.log(f"Found {existing:,} existing results (will be overwritten by Uni-Dock)")

        # Build command
        cmd = self.build_command()
        self.log(f"Command: {' '.join(cmd)}")

        # Set up signal handlers
        signal.signal(signal.SIGINT, self.handle_interrupt)
        signal.signal(signal.SIGTERM, self.handle_interrupt)

        # Start docking
        self.log("Starting Uni-Dock...")
        self.start_time = time.time()

        try:
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            # Start progress monitoring in main thread
            # (simplified - in production could use threading)
            import threading
            monitor_thread = threading.Thread(target=self.monitor_progress, daemon=True)
            monitor_thread.start()

            # Stream output
            for line in self.process.stdout:
                print(line, end='', file=sys.stderr)

            self.process.wait()
            return_code = self.process.returncode

        except Exception as e:
            self.log(f"Error running Uni-Dock: {e}")
            return 1

        # Final statistics
        end_time = time.time()
        total_time = end_time - self.start_time
        final_count = self.count_output_files()

        self.log("=" * 60)
        self.log("Docking Complete")
        self.log("=" * 60)
        self.log(f"Total time: {self.format_duration(total_time)}")
        self.log(f"Results: {final_count:,} / {self.total_ligands:,}")

        if total_time > 0:
            rate = final_count / total_time
            self.log(f"Average rate: {rate:.2f} mol/s ({rate * 3600:.0f} mol/hour)")

        self.log(f"Output directory: {self.output_dir}")
        self.log(f"Exit code: {return_code}")

        return return_code


def main():
    parser = argparse.ArgumentParser(
        description='Uni-Dock wrapper with progress tracking',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic docking run
  uv run python run_docking.py \\
      --receptor receptor.pdbqt \\
      --ligand-index ligands.txt \\
      --output-dir results \\
      --box-config box_config.txt

  # High-exhaustiveness docking with logging
  uv run python run_docking.py \\
      --receptor receptor.pdbqt \\
      --ligand-index ligands.txt \\
      --output-dir results \\
      --box-config box_config.txt \\
      --exhaustiveness 32 \\
      --search-mode detail \\
      --log-file docking.log

  # Fast screening mode
  uv run python run_docking.py \\
      --receptor receptor.pdbqt \\
      --ligand-index ligands.txt \\
      --output-dir results \\
      --box-config box_config.txt \\
      --exhaustiveness 4 \\
      --search-mode fast
"""
    )

    parser.add_argument('--receptor', required=True,
                        help='Receptor PDBQT file')
    parser.add_argument('--ligand-index', required=True,
                        help='File listing ligand PDBQT paths (one per line)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for results')
    parser.add_argument('--box-config', required=True,
                        help='Box configuration file (CENTER/SIZE)')
    parser.add_argument('--unidock-path', default='unidock',
                        help='Path to Uni-Dock binary')
    parser.add_argument('--num-modes', type=int, default=9,
                        help='Number of binding modes (default: 9)')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                        help='Search exhaustiveness (default: 8)')
    parser.add_argument('--search-mode', default='balance',
                        choices=['balance', 'fast', 'detail'],
                        help='Search mode (default: balance)')
    parser.add_argument('--log-file',
                        help='Log file path (optional)')
    parser.add_argument('--progress-interval', type=int, default=60,
                        help='Progress update interval in seconds (default: 60)')

    args = parser.parse_args()

    runner = DockingRunner(
        receptor=Path(args.receptor),
        ligand_index=Path(args.ligand_index),
        output_dir=Path(args.output_dir),
        box_config=Path(args.box_config),
        unidock_path=args.unidock_path,
        num_modes=args.num_modes,
        exhaustiveness=args.exhaustiveness,
        search_mode=args.search_mode,
        log_file=Path(args.log_file) if args.log_file else None,
        progress_interval=args.progress_interval
    )

    sys.exit(runner.run())


if __name__ == '__main__':
    main()
