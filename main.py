#!usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse, sys
from pathlib import Path

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the program.

    Returns
    -------
    argparse.Namespace
    Object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(description="File paths for input and output files.")

    parser.add_argument(
        "-i", "--input_blast_files",
        required=True,
        default=["example_data/blast_files/",
                    "data/Binsularis_BLAST_Hepatozoon.tsv",
                    "data/Binsularis_BLAST_Mitochondrion.tsv",
                    "data/Binsularis_BLAST_SexualChromosomes.tsv"],
        type=str,
        help="Path to input files containing BLAST results."
    )
    parser.add_argument(
        "--coverage_threshold",
        required=False,
        default=0.9,
        type=float,
        help="Threshold for coverage to classify contigs. Default is 0.9."
    )
    parser.add_argument(
        "--contig_size_threshold",
        required=False,
        default=3000,
        type=int,
        help="Threshold for contig size to classify contigs. Default is 3000."
    )

    return parser.parse_args()

def validate_args(args: argparse.Namespace):
    """
    Validate the parsed command-line arguments.

    Parameters
    ----------
    args : argparse.Namespace
    Object containing parsed arguments.

    Raises
    ------
    ValueError
    If any of the arguments are invalid.
    """

# Checking blast directory path and contained file paths
    directory = Path(args.input_blast_files)

    if not directory.is_dir():
        raise ValueError("Input BLAST files path must be a directory containing TSV files.")

    for file in directory.iterdir():
        if not file.name.endswith(".tsv"):
            raise ValueError("Input files must be TSV files with .tsv extension.")
        filepath = Path(file)
        if not filepath.is_file():
            raise ValueError(f"Input file {file} does not exist.")

# Checking the files for appropriate extensions and existence.
    if not args.contig_file.endswith(".txt"):
        raise ValueError("Contig file must be a text file with .txt extension.")
    
# Checking the thresholds.
    if not 0 <= args.coverage_threshold <= 1:
        raise ValueError("Coverage threshold must be a float between 0 and 1.")
    if args.contig_size_threshold <= 0:
        raise ValueError("Contig size threshold must be a positive integer.")
    if not isinstance(args.contig_size_threshold, int):
        raise ValueError("Contig size threshold must be an integer.")

def main() -> int: 
    args = parse_args()

    sys.stdout.write("""
    Arguments received:
        Input BLAST File Paths: {input_blast_files}
        Coverage Threshold: {coverage_threshold}
        Contig Size Threshold: {contig_size_threshold}
    """.format(
        input_blast_files=args.input_blast_files,
        coverage_threshold=args.coverage_threshold,
        contig_size_threshold=args.contig_size_threshold,
    ))

    validate_args(args)


    # Placeholder for the main logic of the program



    sys.stdout.write("""
    Program executed successfully.
        Output Data Frame File Path: {data_frame_file}
        Output Summary Statistics File Path: {summary_stats_file}
        Contigs Bar Plot generated at: {contigs_barplot}
        Base Pairs Bar Plot generated at: {bps_barplot}
    """.format(
        data_frame_file=args.data_frame_file,
        summary_stats_file=args.summary_stats_file,
        contigs_barplot=args.contigs_barplot,
        bps_barplot=args.bps_barplot
    ))

    return 0

if __name__ == "__main__":
    sys.exit(main())