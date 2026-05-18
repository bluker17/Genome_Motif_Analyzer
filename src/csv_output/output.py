#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path


class CSVWriter:
    """
    Handles initialization and incremental writing of motif analysis results to a CSV file.
    """

    def __init__(self, filename: Path) -> None:
        """
        Parameters
        ----------
        filename : Path
            Path to output CSV file.
        """
        self.filename = filename

    def create_csv_file(self) -> None:
        """
        Create and initialize CSV file with headers.
        """
        with open(self.filename, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "FASTA Organism",
                "FASTA File",
                "Strand",
                "Cas Enzyme",
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
                "A", "T", "G", "C"
            ])

    def append_csv(self, stats: dict, fasta_file: str, chrom_name: str) -> None:
        """
        Append one chromosome's results to the CSV file.

        Parameters
        ----------
        stats : dict
            Chromosome-level analysis output from motif + stats pipeline.
        fasta_file : str
            Source FASTA filename.
        chrom_name : str
            Chromosome name/identifier.
        """
        with open(self.filename, "a", newline="") as f:
            writer = csv.writer(f)

            for strand in ("forward", "reverse"):

                if strand not in stats:
                    continue
                probs = stats[strand]["base_probs"]

                for motif, data in stats[strand]["proportion_test"].items():
                    writer.writerow([
                        chrom_name,
                        fasta_file,
                        strand,
                        data["enzyme"],
                        data.get("organism"),
                        motif,
                        data.get("significance"),
                        data.get("p_value"),
                        data.get("z_stat"),
                        data["observed"],
                        data.get("total_positions"),
                        data.get("expected_count"),
                        data.get("expected_motif_prob"),
                        stats["genome_length"],
                        stats[strand]["GC_content"],
                        probs["A"], probs["T"], probs["G"], probs["C"]
                    ])