#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys
from pathlib import Path

# File handling
from src.file_reader.reader import Sequences
from src.file_reader.reader import Enzymes

# Motif locating
from src.motif_locator.locator import Motifs


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
        Directory containing FASTA files: {fasta_files}
        Motif information CSV file: {motif_file}
        Output CSV file: {csv_output}
    """.format(
        fasta_files=args.fasta_files,
        motif_file=args.motif_file,
        csv_output=args.csv_output
    ))

    # validate_args(args)

    #==================================
    # INPUT FILE HANDLING
    #==================================
    genomes = Sequences(args.fasta_files)
    genome_files = genomes.collect_fasta_files()
    # print(genome_files)

    motifs = Enzymes(args.motif_file)
    motif_info = motifs.collect_motifs()
    # print(motif_info)


    #==================================
    # MOTIF LOCATOR
    #==================================
    motif_location_results = Motifs(motif_info, genomes)
    print(motif_location_results)



    sys.stdout.write("""
    Program executed successfully.
        Output CSV File: {csv_output}
    """.format(
        csv_output=args.csv_output
    ))

    return 0

if __name__ == "__main__":
    sys.exit(main())