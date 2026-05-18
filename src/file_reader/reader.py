#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path
from typing import Iterator, Tuple, List
from Bio import SeqIO

class Sequences:
    def __init__(self, fasta_dir: Path) -> None:
        self.fasta_dir = fasta_dir
        self.fasta_files: List[Path] = []

    def collect_fasta_files(self) -> List[Path]:
        """
        Recursively collect FASTA files from a directory.
        """
        for ext in [".fa", ".fasta", ".fna"]:
            self.fasta_files.extend(self.fasta_dir.rglob(f"*{ext}"))
        return self.fasta_files

class Enzymes:
    def __init__(self, motif_file: str) -> None:
        self.motif_file = motif_file
        self.short_enzyme_info = {}
        self.long_enzyme_info = {}

    def collect_motifs(self) -> dict:
        with open(self.motif_file, "r", newline="") as f:
            reader = csv.DictReader(f)

            for row in reader:
                motif = row["motif_sequence"].strip()
                if len(motif) <= 4:
                    self.short_enzyme_info[motif] = {
                        "motif_sequence": row["motif_sequence"],
                        "enzyme": row["enzyme"].strip(),
                        "organism": row["organism"].strip()
                    }
                else:
                    self.long_enzyme_info[motif] = {
                        "motif_sequence": row["motif_sequence"],
                        "enzyme": row["enzyme"].strip(),
                        "organism": row["organism"].strip()
                    }
        return self.short_enzyme_info, self.long_enzyme_info



