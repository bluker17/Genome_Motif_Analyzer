#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from Bio import SeqIO


# -------------------------
# FILE DISCOVERY
# -------------------------
def collect_fasta_files(fasta_dir: Path) -> list[Path]:
    """
    Recursively collect FASTA files from a directory.

    Supported extensions: .fa, .fasta, .fna
    """
    fasta_files = []
    for ext in [".fa", ".fasta", ".fna"]:
        fasta_files.extend(fasta_dir.rglob(f"*{ext}"))
    return fasta_files


# -------------------------
# FASTA PARSING
# -------------------------
def read_chromosomes(file: Path):
    """
    Yield (chromosome_name, sequence) from a FASTA file.
    """
    for record in SeqIO.parse(file, "fasta"):
        yield record.description, str(record.seq)