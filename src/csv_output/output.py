#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path
from src.macromolecule_alphabet.alphabet import Alphabet
from src.results.result_dataclasses import ChromosomeResult

class CSVWriter:
    """
    Handles initialization and incremental writing of motif analysis
    results to a CSV file for DNA or RNA.
    """

    def __init__(self, filename: Path, alphabet: Alphabet) -> None:
        """
        Parameters
        ----------
        filename : Path
            Path to output CSV file.
        macromolecule : str
            Either 'DNA' or 'RNA'.
        """
        self.filename = filename
        self.alphabet = alphabet

    def create_csv_file(self) -> None:
        """
        Create and initialize CSV file with headers.
        """
        headers = [
            "FASTA Organism",
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
        ]

        with open(self.filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)

    # def append_csv(self, stats: dict, fasta_file: str, chrom_name: str) -> None:
    #     """
    #     Append one chromosome's results to the CSV file.
    #     """

    #     with open(self.filename, "a", newline="") as f:
    #         writer = csv.writer(f)

    #         for strand in ("forward", "reverse"):
    #             if strand not in stats:
    #                 continue

    #             probs = stats[strand]["base_probs"]
    #             for motif, data in stats[strand]["proportion_test"].items():
    #                 row = [
    #                     chrom_name,
    #                     fasta_file,
    #                     strand,
    #                     data["enzyme"],
    #                     data.get("organism"),
    #                     motif,
    #                     data.get("significance"),
    #                     data.get("p_value"),
    #                     data.get("z_stat"),
    #                     data["observed"],
    #                     data.get("total_positions"),
    #                     data.get("expected_count"),
    #                     data.get("expected_motif_prob"),
    #                     stats["genome_length"],
    #                     stats[strand]["GC_content"],
    #                 ]

    #                 # Append nucleotide probabilities dynamically
    #                 row.extend(probs[base] for base in self.alphabet.bases)

    #                 writer.writerow(row)

    def append_csv(self, chrom: ChromosomeResult, fasta_file: str) -> None:

        with open(self.filename, "a", newline="") as f:
            writer = csv.writer(f)

            for strand_name, strand in [
                ("forward", chrom.forward),
                ("reverse", chrom.reverse)
            ]:

                if strand is None:
                    continue

                probs = strand.base_probs

                for motif in strand.motifs:

                    row = [
                        chrom.chromosome,
                        fasta_file,
                        strand_name,
                        motif.enzyme,
                        motif.organism,
                        motif.motif,
                        motif.significance,
                        motif.p_value,
                        motif.z_stat,
                        motif.observed,
                        motif.total_positions,
                        motif.expected_count,
                        motif.expected_motif_prob,
                        chrom.genome_length,
                        strand.gc_content,
                    ]

                    # append nucleotide probabilities
                    row.extend(probs[base] for base in self.alphabet.bases)

                    writer.writerow(row)