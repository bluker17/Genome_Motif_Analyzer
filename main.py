#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys
from pathlib import Path

# File handling
from src.file_reader.reader import Sequences
from src.file_reader.reader import Enzymes

# Motif locating
from src.motif_locator.numpy_locator import Numpy_Motif_Search
from src.motif_locator.numba_locator import Numba_Motif_Search

# Statistical analysis
from src.statistic_analysis.statistics import Statistics

# CSV output
from src.csv_output.output import CSVWriter

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
    
    # file_parser.add_argument("-s", "--strand_to_search",
    #                          type=str, required=True,
    #                          help="Distinguish which strand to investigate. Options are 'forward', 'reverse', or 'both'.")
    
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
    short_motifs, long_motifs = motifs.collect_motifs()
    # print(short_motifs)
    # print(long_motifs)


    #==================================
    # MOTIF LOCATION AND STATISTICS
    #==================================
    # Initialize the CSV output file
    csv_writer = CSVWriter(Path(args.csv_output))
    csv_writer.create_csv_file()

    # Initialize the statistics class
    stats_engine = Statistics()

    # Search for motifs in each genome. If the motif length is less than or equal to 4 bases then the sliding window method is used. Otherwise, numba is used.

    # Initialize engines once
    numpy_locator = Numpy_Motif_Search(short_motifs, genome_files)
    numba_locator = Numba_Motif_Search(long_motifs, genome_files)

    for genome in genome_files:
        print(f"\nProcessing {genome}")

        for chrom_name, sequence in numpy_locator.stream_fasta(genome):
            print(f"\tProcessing chromosome {chrom_name}")

            if short_motifs:
                chrom_stats = numpy_locator.process_chromosome(
                    chrom_name,
                    sequence
                )

                chrom_stats = stats_engine.run_proportion_test(chrom_stats)

                csv_writer.append_csv(
                    stats=chrom_stats,
                    fasta_file=genome.name,
                    chrom_name=chrom_name
                )

            if long_motifs:
                chrom_stats = numba_locator.process_chromosome(
                    chrom_name,
                    sequence
                )

                chrom_stats = stats_engine.run_proportion_test(chrom_stats)

                csv_writer.append_csv(
                    stats=chrom_stats,
                    fasta_file=genome.name,
                    chrom_name=chrom_name
                )


    sys.stdout.write("""
    Program executed successfully.
        Output CSV File: {csv_output}
    """.format(
        csv_output=args.csv_output
    ))

    return 0

if __name__ == "__main__":
    sys.exit(main())