#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import traceback
from pathlib import Path

# Fix module imports when executed inside Nextflow's bin/ directory
sys.path.append(str(Path(__file__).resolve().parents[1]))

# Import your summary logic module
from src.summary_statistics.summary import SummaryStats

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the summary process."""
    parser = argparse.ArgumentParser(description="Processes an initial motif analysis CSV and computes summary statistics.")
    parser.add_argument("-i", "--initial_csv", 
                        type=Path, required=True, 
                        help="Path to the initial motif analysis results CSV file.")
    return parser.parse_args()

def main() -> int:
    try:
        args = parse_args()

        if not args.initial_csv.exists():
            raise FileNotFoundError(f"Initial analysis CSV file not found: {args.initial_csv}")

        sys.stdout.write(f"\nGenerating summary statistics for: {args.initial_csv}\n")

        # Instantiate and run your summary logic
        summary = SummaryStats(args.initial_csv)
        summary.summary_findings()

        sys.stdout.write(f"Summary statistics compilation finished successfully. Output: {summary.output_file}\n")
        return 0

    except Exception as e:
        traceback.print_exc(file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())