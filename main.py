#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, sys, traceback
from pathlib import Path

# Alphabets
from src.alphabet.macromolecule_alphabet import ALPHABETS
from src.alphabet.result_alphabet import EntryResults

# File handling
from src.file_reader.reader import Sequences
from src.file_reader.reader import Enzymes

# Motif locating
from src.motif_locator.numpy_locator import Numpy_Motif_Search
from src.motif_locator.numba_locator import Numba_Motif_Search

# Statistical analysis
from src.statistic_analysis.statistics import Statistics

# CSV stats outputs
from src.csv_output.output import CSVWriter
from src.summary_statistics.summary import SummaryStats

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
    try:

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

        macromolecule = ALPHABETS[args.macromolecule]

        #==================================
        # INPUT FILE HANDLING
        #==================================
        genomes = Sequences(args.fasta_files)
        genome_files = genomes.collect_fasta_files()
        # print(genome_files)

        motifs = Enzymes(args.motif_file, macromolecule)
        short_motifs, long_motifs = motifs.collect_motifs()
        # print(short_motifs)
        # print(long_motifs)

        #==================================
        # MOTIF LOCATION AND STATISTICS
        #==================================
        # Initialize the CSV output file
        csv_writer = CSVWriter(Path(args.csv_output), macromolecule)
        csv_writer.create_csv_file()

        # Initialize the statistics class
        stats_engine = Statistics(macromolecule)

        # Search for motifs in each genome. If the motif length is less than or equal to 4 bases then the sliding window method is used. Otherwise, numba is used.

        # Initialize engines once
        # print(short_motifs.keys())
        # print(long_motifs.keys())
        numpy_locator = Numpy_Motif_Search(short_motifs, genome_files, macromolecule)
        numba_locator = Numba_Motif_Search(long_motifs, genome_files, macromolecule)

        for genome in genome_files:
            print(f"\nProcessing {genome}")

            if short_motifs:
                for entry_name, sequence in numpy_locator.stream_fasta(genome):
                    #print(f"\tProcessing FASTA entry {entry_name}")
                    entry_data: EntryResults = numpy_locator.process_entry(entry_name, sequence)
                    entry_data = stats_engine.run_proportion_test(entry_data)

                    csv_writer.append_csv(stats=entry_data, fasta_file=genome.name, entry_name=entry_name)

            if long_motifs:
                for entry_name, sequence in numba_locator.stream_fasta(genome):
                    #print(f"\tProcessing FASTA entry {entry_name}")
                    entry_data: EntryResults = numba_locator.process_entry(entry_name, sequence)
                    entry_data = stats_engine.run_proportion_test(entry_data)

                    csv_writer.append_csv(stats=entry_data, fasta_file=genome.name,entry_name=entry_name)

        summary = SummaryStats(Path(args.csv_output))
        summary.summary_findings()

        sys.stdout.write("""
        Program executed successfully.
            Output CSV File: {csv_output}
            Summary Output CSV File: {summary_csv_output}
        """.format(
            csv_output=args.csv_output,
            summary_csv_output=f"output/{args.csv_output.stem}_summary_statistics.csv"
        ))

        return 0

    except ValueError as e:
        traceback.print_exc(file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())