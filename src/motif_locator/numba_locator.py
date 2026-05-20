#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from typing import Iterator
from Bio import SeqIO
from numba import njit

import numpy as np

degenerate_map: dict[str, list[str]] = {
    'A': ['A'],
    'T': ['T'],
    'G': ['G'],
    'C': ['C'],
    'N': ['A', 'T', 'G', 'C'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'W': ['A', 'T'],
    'S': ['G', 'C'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G']
}

RC_MAP: dict[str, str] = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'V': 'B',
    'D': 'H',
    'H': 'D',
    'N': 'N'
}

DNA_MASK = np.zeros(256, dtype=np.uint8)
DNA_MASK[ord('A')] = 1
DNA_MASK[ord('C')] = 2
DNA_MASK[ord('G')] = 4
DNA_MASK[ord('T')] = 8


DEGENERATE_MASKS: dict[str, int] = {
    'A': 1,
    'C': 2,
    'G': 4,
    'T': 8,
    'R': 1 | 4,
    'Y': 2 | 8,
    'S': 2 | 4,
    'W': 1 | 8,
    'K': 4 | 8,
    'M': 1 | 2,
    'B': 2 | 4 | 8,
    'D': 1 | 4 | 8,
    'H': 1 | 2 | 8,
    'V': 1 | 2 | 4,
    'N': 1 | 2 | 4 | 8,
}

@njit(cache=True)
def motif_search_numba(genome: np.ndarray, motif: np.ndarray) -> int:
    """
    Fast bitmask motif search using early exits and Numba JIT compilation.
    """
    n, m, count = len(genome), len(motif), 0

    for i in range(n - m + 1):
        for j in range(m):
            if (genome[i + j] & motif[j]) == 0:
                break
        else:
            count += 1

    return count

class Numba_Motif_Search:
    """
    Motif detection and genome-level sequence statistics.
    """

    def __init__(self, motif_info: dict, genomes: list[Path]) -> None:
        """
        Initialize motif search system.
        """
        self.motif_info = motif_info
        self.genomes = genomes
        self.motif_results: dict = {}

        self.motif_masks: dict[str, dict[str, np.ndarray]] = {
            enzyme: {
                "forward": self.build_motif_bitmasks(info["motif_sequence"]),
                "reverse": self.build_motif_bitmasks(self.revcomp_motif(info["motif_sequence"]))
            }
            for enzyme, info in motif_info.items()
        }

    def stream_fasta(self, file: Path) -> Iterator[tuple[str, str]]:
        """
        Stream FASTA records one chromosome at a time.
        """
        for record in SeqIO.parse(file, "fasta"):
            yield record.description, str(record.seq)

    def seq_to_array(self, sequence: str) -> np.ndarray:
        """
        Convert DNA sequence into ASCII NumPy array.
        """
        return np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    def encode_genome(self, arr: np.ndarray) -> np.ndarray:
        """
        Encode genome into bitmask representation.
        """
        return DNA_MASK[arr]

    def build_motif_bitmasks(self, motif: str) -> np.ndarray:
        """
        Convert motif into bitmask array.
        """
        return np.array([DEGENERATE_MASKS[b] for b in motif.upper()], dtype=np.uint8)

    def revcomp_motif(self, motif: str) -> str:
        """
        Generate reverse complement motif (no genome reversal needed).
        """
        return "".join(RC_MAP[b] for b in motif.upper()[::-1])

    def motif_search(self, genome: np.ndarray, motif: np.ndarray) -> int:
        """
        Search encoded genome for encoded motif.
        """
        if len(genome) < len(motif):
            return 0
        return motif_search_numba(genome, motif)

    def compute_base_probs(self, seq: str) -> dict[str, float]:
        """
        Compute nucleotide frequencies.
        """
        seq = seq.upper()
        counts = Counter(seq)
        total = sum(counts[b] for b in "ACGT")

        if total == 0:
            return {b: 0.0 for b in "ACGT"}

        return {b: counts[b] / total for b in "ACGT"}

    def reverse_base_probs(self, forward_probs: dict[str, float]) -> dict[str, float]:
        """
        Convert forward strand probabilities into reverse strand probabilities.
        """
        return {"A": forward_probs["T"], "T": forward_probs["A"], "G": forward_probs["C"], "C": forward_probs["G"]}

    def process_chromosome(self, chrom_name: str, sequence: str) -> dict:
        """
        Process chromosome for motif occurrence statistics.
        """

        sequence = sequence.upper()

        genome_encoded = self.encode_genome(self.seq_to_array(sequence))

        genome_length = len(sequence)

        fwd_base_probs = self.compute_base_probs(sequence)
        rev_base_probs = self.reverse_base_probs(fwd_base_probs)

        chrom_stats = {
            "chromosome": chrom_name,
            "genome_length": genome_length,
            "forward": {
                "base_probs": fwd_base_probs,
                "GC_content": fwd_base_probs["G"] + fwd_base_probs["C"],
                "proportion_test": {},
            },
            "reverse": {
                "base_probs": rev_base_probs,
                "GC_content": rev_base_probs["G"] + rev_base_probs["C"],
                "proportion_test": {},
            }
        }

        for enzyme, info in self.motif_info.items():

            motif = info["motif_sequence"]

            fwd_mask = self.motif_masks[enzyme]["forward"]
            rev_mask = self.motif_masks[enzyme]["reverse"]

            fwd_count = self.motif_search(genome_encoded, fwd_mask)
            rev_count = self.motif_search(genome_encoded, rev_mask)

            chrom_stats["forward"]["proportion_test"][motif] = {
                "observed": fwd_count,
                "enzyme": enzyme,
                "organism": info.get("organism")
            }

            chrom_stats["reverse"]["proportion_test"][motif] = {
                "observed": rev_count,
                "enzyme": enzyme,
                "organism": info.get("organism")
            }

        return chrom_stats