#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv


# -------------------------
# CSV INITIALIZATION
# -------------------------
def create_csv_file(filename) -> bool:
    """
    Create and initialize the output CSV file with headers.
    """
    with open(filename, "w", newline="") as f:
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
    return True


# -------------------------
# CSV APPEND
# -------------------------
def append_csv(stats, fasta_file, chrom_name, filename):
    """
    Append motif analysis results to CSV file.
    """
    with open(filename, "a", newline="") as f:
        writer = csv.writer(f)

        for strand in ["forward", "reverse"]:
            probs = stats[strand]["base_probs"]

            for motif, data in stats[strand]["proportion_test"].items():
                writer.writerow([
                    chrom_name,
                    fasta_file,
                    strand,
                    data["enzyme"],
                    data["enzyme_organism"],
                    motif,
                    data["significance"],
                    data["p_value"],
                    data["z_stat"],
                    data["observed"],
                    data["total_positions"],
                    data["expected_count"],
                    data["expected_motif_prob"],
                    stats["genome_length"],
                    stats[strand]["GC_content"],
                    probs["A"], probs["T"], probs["G"], probs["C"]
                ])