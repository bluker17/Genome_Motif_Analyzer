#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys
from pathlib import Path

# Alphabet Initialization
from src.macromolecule_alphabet.alphabet import ALPHABETS

# File handling
from src.file_reader.reader import Sequences
from src.file_reader.reader import Enzymes

# Motif locating
from src.motif_locator.locator import AhoCorasick_Motif_Search

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
    
    file_parser.add_argument("-s", "--strand_to_search",
                             type=str, required=True,
                             choices=["forward", "reverse", "both"],
                             help="Distinguishes which strand to investigate. Options are 'forward', 'reverse', or 'both'.")
    
    file_parser.add_argument("--macromolecule",
                             type=str, required=True,
                             choices=["DNA", "RNA"],
                             help="Distinguishes which macromolecule is provided in the FASTA directory to parse. Options are 'DNA', or 'RNA'. Only one option can be present in the FASTA directory.")
    
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
    valid_extensions = [".fasta", ".fa", ".fna", ".faa", ".fas", ".ffn", ".frn"]

    if not fasta_dir.exists():
        raise ValueError("FASTA files path does not exist.")

    if not fasta_dir.is_dir():
        raise ValueError("FASTA files path must be a directory containing FASTA files.")

    if not fasta_dir.is_relative_to("data/") and not fasta_dir.is_relative_to("testing_materials/example_data/"):
        raise FileNotFoundError("FASTA directory must be located in 'data/'.")

    # Recursively find FASTA files
    fasta_files = [file for file in fasta_dir.rglob("*") if file.is_file() and file.suffix.lower() in valid_extensions]

    # Ensure at least one FASTA file exists
    if not fasta_files:
        raise ValueError(
            f"No FASTA files found in '{fasta_dir}'. "
            f"Expected extensions: {sorted(valid_extensions)}"
        )

    # Validate motif file
    motif_file = Path(args.motif_file)

    if not motif_file.exists():
        raise FileNotFoundError("Motif file does not exist.")

    if not motif_file.is_file():
        raise FileNotFoundError("Motif file path must point to a file.")

    if motif_file.suffix.lower() != ".csv":
        raise ValueError("Motif file must be a CSV file (.csv).")
    
    if not motif_file.is_relative_to("data/") and not motif_file.is_relative_to("testing_materials/example_data/"):
        raise FileNotFoundError("Motif file must be located in 'data/'.")

    # Validate CSV output file
    output_path = Path(args.csv_output)

    if output_path.suffix.lower() != ".csv":
        raise ValueError("CSV output file must have a .csv extension.")

    if not output_path.is_relative_to("output/") and not output_path.is_relative_to("testing_materials/example_outputs/"):
        raise FileNotFoundError("CSV output file must be located in 'output/'.")

def main() -> int: 
    args = parse_args()

    sys.stdout.write("""
    Arguments received:
        Directory containing FASTA file(s): {fasta_files}
        Macromolecule in FASTA file(s): {macromolecule}
        Strand to investigate: {strand_to_search}
        Motif information CSV file: {motif_file}
        Output CSV file: {csv_output}
    """.format(
        fasta_files=args.fasta_files,
        macromolecule = args.macromolecule,
        strand_to_search = args.strand_to_search,
        motif_file=args.motif_file,
        csv_output=args.csv_output
    ))

    validate_args(args)

    alphabet = ALPHABETS[args.macromolecule.upper()]

    #==================================
    # INPUT FILE HANDLING
    #==================================
    genomes = Sequences(args.fasta_files)
    genome_files = genomes.collect_fasta_files()
    # print(genome_files)

    motifs = Enzymes(args.motif_file, alphabet)
    short_motifs = motifs.collect_motifs()


    #==================================
    # MOTIF LOCATION AND STATISTICS
    #==================================
    # Initialize the CSV output file
    csv_writer = CSVWriter(Path(args.csv_output), alphabet)
    csv_writer.create_csv_file()

    # Initialize the statistics class
    stats_engine = Statistics(alphabet)

    # # Initialize engines once
    ahocorasick = AhoCorasick_Motif_Search(short_motifs, genome_files, args.strand_to_search, alphabet)
    for genome in genome_files:
        print(f"\nProcessing {genome}")

        for chrom_name, sequence in ahocorasick.stream_fasta(genome):
            # print(f"\tProcessing chromosome {chrom_name}")
            chrom_stats = ahocorasick.process_chromosome(chrom_name, sequence)
            chrom_stats = stats_engine.run_proportion_test(chrom_stats)

            csv_writer.append_csv(chrom_stats, genome.name)

    sys.stdout.write("""
    Program executed successfully.
        Output CSV File: {csv_output}
    """.format(
        csv_output=args.csv_output
    ))

    return 0

if __name__ == "__main__":
    sys.exit(main())