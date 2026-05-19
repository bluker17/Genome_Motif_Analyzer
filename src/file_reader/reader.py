#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
from pathlib import Path
from typing import List
from collections import defaultdict
from src.macromolecule_alphabet.alphabet import Alphabet

class Sequences:
    def __init__(self, fasta_dir: Path) -> None:
        self.fasta_dir = fasta_dir
        self.fasta_files: List[Path] = []

    def collect_fasta_files(self) -> List[Path]:
        """
        Recursively collect FASTA files from a directory.
        """
        self.fasta_files = [p for p in self.fasta_dir.rglob("*") if p.suffix in {".fasta", ".fa", ".fna", ".faa", ".fas", ".ffn", ".frn"}]
        return self.fasta_files

class Enzymes:
    def __init__(self, motif_file: Path, alphabet: Alphabet) -> None:
        self.motif_file = motif_file
        self.enzyme_info = defaultdict(list)
        self.alphabet = alphabet

    def invalid_motif_chars(self, motif: str) -> set:
        valid_bases = self.alphabet.degenerate_map.keys()
        return set(motif) - set(valid_bases)

    def collect_motifs(self) -> dict:
        with open(self.motif_file, "r", newline="") as mf:
            required_columns = {"enzyme", "motif_sequence", "organism"}
            reader = csv.DictReader(mf)

            # CSV error handling
            # CSV is not empty
            if reader.fieldnames is None:
                raise ValueError("CSV file is empty or missing headers.")

            # CSV is not missing any required fields
            headers = set(reader.fieldnames)
            missing_cols = required_columns - headers

            if missing_cols:
                raise ValueError(
                    f"CSV is missing required columns: {missing_cols}. "
                    f"Expected columns: {required_columns}"
                )

            # Checks that values provided in each field for every row
            for i, row in enumerate(reader, start=2):
                enzyme = (row.get("enzyme") or "").strip()
                motif = (row.get("motif_sequence") or "").upper().strip()
                organism = (row.get("organism") or "").strip()

                if not enzyme or not motif or not organism:
                    raise ValueError(f"Missing value in row {i}: {row}")
                
                invalid = self.invalid_motif_chars(motif)
                if invalid:
                    raise ValueError(
                        f"Invalid motif in row {i}: '{motif}'. "
                        f"Invalid characters: {invalid}"
                    )

                # Store validated data
                self.enzyme_info[motif].append({
                    "motif_sequence": motif,
                    "enzyme": enzyme,
                    "organism": organism
                })

        return self.enzyme_info



