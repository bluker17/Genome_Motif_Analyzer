#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path
from typing import List

from src.alphabet.macromolecule_alphabet import Alphabet

class Sequences:
    def __init__(self, fasta_dir: Path) -> None:
        self.fasta_dir = fasta_dir
        self.fasta_files: List[Path] = []

    def collect_fasta_files(self) -> List[Path]:
        """
        Recursively collect FASTA files from a directory.
        """
        for ext in [".fasta", ".fa", ".fna", ".faa", ".fas", ".ffn", ".frn"]:
            self.fasta_files.extend(self.fasta_dir.rglob(f"*{ext}"))
        return self.fasta_files

class Enzymes:
    def __init__(self, motif_file: str, macromolecule: Alphabet) -> None:
        self.motif_file = motif_file
        self.short_enzyme_info = {}
        self.long_enzyme_info = {}
        self.macromolecule = macromolecule

    def collect_motifs(self) -> dict:
        with open(self.motif_file, "r", newline="") as f:
            reader = csv.DictReader(f)

            for row in reader:
                motif = row["motif_sequence"].strip()
                enzyme = row["enzyme"].strip()
                composite_key = f"{enzyme}:{motif}"

                for b in motif:
                    if b not in self.macromolecule.degenerate_map:
                       raise ValueError(f"Macromolecule is set to {self.macromolecule.name}. {motif} contains '{b}'. Please update the motif to handle only {self.macromolecule.name} bases.")
                
                if len(motif) <= 4:
                    self.short_enzyme_info[composite_key] = {
                        "motif_sequence": row["motif_sequence"],
                        "enzyme": row["enzyme"].strip(),
                        "organism": row["organism"].strip()
                    }
                else:
                    self.long_enzyme_info[composite_key] = {
                        "motif_sequence": row["motif_sequence"],
                        "enzyme": row["enzyme"].strip(),
                        "organism": row["organism"].strip()
                    }
        return self.short_enzyme_info, self.long_enzyme_info