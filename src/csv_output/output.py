#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path

from src.alphabet.macromolecule_alphabet import Alphabet
from src.alphabet.result_alphabet import EntryResults

class CSVWriter:
    """
    Handles initialization and incremental writing of motif analysis results to a CSV file.
    """

    def __init__(self, filename: Path, alphabet: Alphabet) -> None:
        """
        Parameters
        ----------
        filename : Path
            Path to output CSV file.
        """
        self.filename = filename
        self.alphabet = alphabet

    def create_csv_file(self) -> None:
        """
        Create and initialize CSV file with headers.
        """
        with open(self.filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "FASTA Entry",
                "FASTA File",
                "Strand",
                "Enzyme",
                "Enzyme Source Organism",
                "Motif",
                "Significance",
                "p-value",
                "z-stat",
                "Observed Matches",
                "Possible Positions",
                "Expected Matches",
                "Expected Motif Prob",
                "Genome Length",
                "GC Content",
                *self.alphabet.bases
            ])

    def append_csv(self, stats: EntryResults, fasta_file: str, entry_name: str) -> None:
        """
        Append one entry's results to the CSV file.

        Parameters
        ----------
        stats : dict
            Chromosome-level analysis output from motif + stats pipeline.
        fasta_file : str
            Source FASTA filename.
        entry_name : str
            Entry name/identifier.
        """
        with open(self.filename, "a", newline="") as f:
            writer = csv.writer(f)


            strands = [("forward", stats.forward), ("reverse", stats.reverse)]

            for strand_name, strand_stats in strands:
                probs = strand_stats.base_probs

                for motif, data in strand_stats.proportion_test.items():
                    writer.writerow([
                        entry_name,
                        fasta_file,
                        strand_name,
                        data.enzyme,
                        data.organism,
                        motif,
                        data.significance,
                        data.p_value,
                        data.z_stat,
                        data.observed,
                        data.total_positions,
                        data.expected_count,
                        data.expected_motif_prob,
                        stats.genome_length,
                        strand_stats.GC_content,
                        *[probs.get(b, 0.0) for b in self.alphabet.bases]
                    ])