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
    file_parser = argparse.ArgumentParser(description="Takes the directory and file arguments needed for the program.")

    file_parser.add_argument("-c", "--csv_output", 
                            type=Path, required=True, 
                            help="CSV path and filename output.")
    
    file_parser.add_argument("-f", "--fasta_files", 
                            type=Path, required=True,
                            help="Directory containing FASTA files.")
    
    file_parser.add_argument("-m", "--motif_file", 
                            type=Path, required=True,
                            help="CSV file containing the motif information.")
    
    return file_parser.parse_args()

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

    # Validate FASTA directory
    fasta_dir = Path(args.fasta_files)

    if not fasta_dir.exists():
        raise ValueError("FASTA files path does not exist.")

    if not fasta_dir.is_dir():
        raise ValueError("FASTA files path must be a directory containing FASTA files.")

    # Ensure directory contains at least one FASTA file
    fasta_files = list(fasta_dir.iterdir())
    if not fasta_files:
        raise ValueError("FASTA directory is empty.")

    # Check file extensions in FASTA directory
    valid_extensions = {".fasta", ".fa", ".fna"}
    for file in fasta_files:
        if file.is_file() and file.suffix.lower() not in valid_extensions:
            raise ValueError(
                f"Invalid file detected in FASTA directory: {file.name}. "
                f"Expected FASTA extensions: {valid_extensions}"
            )

    # Validate motif file
    motif_file = Path(args.motif_file)

    if not motif_file.exists():
        raise ValueError("Motif file does not exist.")

    if not motif_file.is_file():
        raise ValueError("Motif file path must point to a file.")

    if motif_file.suffix.lower() != ".csv":
        raise ValueError("Motif file must be a CSV file (.csv).")

    # Validate output path (basic check)
    output_path = Path(args.csv_output)

    if output_path.suffix.lower() != ".csv":
        raise ValueError("CSV output file must have a .csv extension.")

    # Ensure output directory exists (optional but useful)
    if output_path.parent and not output_path.parent.exists():
        raise ValueError("Output directory does not exist.")

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