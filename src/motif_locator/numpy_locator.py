#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from typing import Iterator
from Bio import SeqIO

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

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

class Numpy_Motif_Search:
    """
    Motif detection and genome-level sequence statistics.

    Parameters
    ----------
    motif_info : dict
        Dictionary containing motif metadata and sequences.
    genomes : list[Path]
        FASTA genome files to process.
    """

    def __init__(self, motif_info: dict, genomes: list[Path]) -> None:

        self.motif_info = motif_info
        self.genomes = genomes

        self.motif_results: dict = {}

        self.motif_masks: dict[str, dict[str, list[np.ndarray]]] = {
            enzyme: {
                "forward": self.build_motif_masks(info["motif_sequence"]),
                "reverse": self.build_motif_masks(self.revcomp_motif(info["motif_sequence"]))
            }
            for enzyme, info in motif_info.items()
        }

    def stream_fasta(self, file: Path) -> Iterator[tuple[str, str]]:
        """
        Stream FASTA records one chromosome at a time.

        Parameters
        ----------
        file : Path
            FASTA file path.

        Yields
        ------
        Iterator[tuple[str, str]]
            Chromosome name and sequence.
        """

        for record in SeqIO.parse(file, "fasta"):
            yield record.description, str(record.seq)

    def seq_to_array(self, sequence: str) -> np.ndarray:
        """
        Convert DNA sequence into compact NumPy byte array.

        Parameters
        ----------
        sequence : str
            DNA sequence.

        Returns
        -------
        np.ndarray
            NumPy array of ASCII byte characters.
        """

        return np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)

    def build_motif_masks(self, motif: str) -> list[np.ndarray]:
        """
        Build ASCII lookup masks for degenerate motif matching.

        Parameters
        ----------
        motif : str
            Motif sequence.

        Returns
        -------
        list[np.ndarray]
            List of boolean lookup masks.
        """

        masks: list[np.ndarray] = []

        for char in motif.upper():

            mask = np.zeros(256, dtype=bool)

            for base in degenerate_map[char]:
                mask[ord(base)] = True

            masks.append(mask)

        return masks

    def motif_search(self, genome_arr: np.ndarray, motif_masks: list[np.ndarray]) -> tuple[np.ndarray, int]:
        """
        Locate motif occurrences in a genome array.

        Parameters
        ----------
        genome_arr : np.ndarray
            Genome sequence as NumPy byte array.
        motif_masks : list[np.ndarray]
            Precomputed motif lookup masks.

        Returns
        -------
        tuple[np.ndarray, int]
            Match positions and total count.
        """

        motif_len = len(motif_masks)

        if len(genome_arr) < motif_len:
            return np.array([], dtype=int), 0

        windows = sliding_window_view(genome_arr, motif_len)

        matches = np.ones(len(windows), dtype=bool)

        for i, mask in enumerate(motif_masks):

            matches &= mask[windows[:, i]]

            if not np.any(matches):
                return np.array([], dtype=int), 0

        positions = np.where(matches)[0]

        return positions, len(positions)

    def compute_base_probs(self, seq: str) -> dict[str, float]:
        """
        Compute nucleotide frequencies.

        Parameters
        ----------
        seq : str
            DNA sequence.

        Returns
        -------
        dict[str, float]
            Base probabilities for A/T/G/C.
        """

        seq = seq.upper()

        counts = Counter(seq)

        total = sum(counts[b] for b in "ACGT")

        if total == 0:
            return {b: 0.0 for b in "ACGT"}

        return {b: counts[b] / total for b in "ACGT"}

    def revcomp_motif(self, motif: str) -> str:
        """
        Generate reverse complement motif.

        Parameters
        ----------
        motif : str
            Motif sequence.

        Returns
        -------
        str
            Reverse complement motif.
        """

        return "".join(RC_MAP[b] for b in motif.upper()[::-1])

    def reverse_base_probs(self, forward_probs: dict[str, float]) -> dict[str, float]:
        """
        Convert forward strand probabilities into reverse strand probabilities.

        Parameters
        ----------
        forward_probs : dict[str, float]
            Forward strand probabilities.

        Returns
        -------
        dict[str, float]
            Reverse strand probabilities.
        """

        return {
            "A": forward_probs["T"],
            "T": forward_probs["A"],
            "G": forward_probs["C"],
            "C": forward_probs["G"],
        }

    def process_chromosome(self, chrom_name: str, sequence: str) -> dict:
        """
        Process a chromosome for motif occurrence statistics.

        Parameters
        ----------
        chrom_name : str
            Chromosome identifier.
        sequence : str
            Chromosome sequence.

        Returns
        -------
        dict
            Chromosome motif statistics.
        """

        sequence = sequence.upper()

        genome_arr = self.seq_to_array(sequence)

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

            fwd_masks = self.motif_masks[enzyme]["forward"]
            rev_masks = self.motif_masks[enzyme]["reverse"]

            _, fwd_count = self.motif_search(genome_arr, fwd_masks)
            _, rev_count = self.motif_search(genome_arr, rev_masks)

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