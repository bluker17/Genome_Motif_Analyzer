#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from typing import Iterator
from Bio import SeqIO
from numba import njit

from src.alphabet.macromolecule_alphabet import Alphabet
from src.alphabet.result_alphabet import EntryResults, StrandResults, MotifObservation

import numpy as np

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

    def __init__(self, motif_info: dict, genomes: list[Path], alphabet: Alphabet) -> None:
        """
        Initialize motif search system.
        """
        self.motif_info = motif_info
        self.genomes = genomes
        self.motif_results: dict = {}
        self.alphabet = alphabet

        self.motif_masks: dict[str, dict[str, np.ndarray]] = {
            enzyme: {
                "forward": self.build_motif_bitmasks(info["motif_sequence"]),
                "reverse": self.build_motif_bitmasks(
                    self.alphabet.reverse_complement(
                    info["motif_sequence"]
                    )
                )
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
        return self.alphabet.ascii_mask[arr]

    def build_motif_bitmasks(self, motif: str) -> np.ndarray:
        """
        Convert motif into bitmask array.
        """
        return np.array([self.alphabet.mask_map[b] for b in motif.upper()], dtype=np.uint8)

    def motif_search(self, genome: np.ndarray, motif: np.ndarray) -> int:
        """
        Search encoded genome for encoded motif.
        """
        if len(genome) < len(motif):
            return 0
        return motif_search_numba(genome, motif)

    def compute_base_probs(self, seq: str) -> dict[str, float]:

        seq = seq.upper()

        counts = Counter(seq)

        total = sum(counts[b] for b in self.alphabet.bases)

        if total == 0:
            return {b: 0.0 for b in self.alphabet.bases}

        return {b: counts[b] / total for b in self.alphabet.bases}

    def reverse_base_probs(self, forward_probs: dict[str, float]) -> dict[str, float]:

        return {base: forward_probs[self.alphabet.complement_map[base]] for base in self.alphabet.bases}

    def process_entry(self, entry_name: str, sequence: str) -> EntryResults:
        """
        Process chromosome for motif occurrence statistics.
        """

        sequence = sequence.upper()

        genome_encoded = self.encode_genome(self.seq_to_array(sequence))

        genome_length = len(sequence)

        fwd_base_probs = self.compute_base_probs(sequence)
        rev_base_probs = self.reverse_base_probs(fwd_base_probs)

        fwd_stats = StrandResults(
            base_probs=fwd_base_probs,
            GC_content=sum(fwd_base_probs[b] for b in self.alphabet.gc_bases)
        )

        rev_stats = StrandResults(
            base_probs=rev_base_probs,
            GC_content=sum(rev_base_probs[b] for b in self.alphabet.gc_bases)
        )

        for enzyme, info in self.motif_info.items():

            motif = info["motif_sequence"]

            fwd_mask = self.motif_masks[enzyme]["forward"]
            rev_mask = self.motif_masks[enzyme]["reverse"]

            fwd_count = self.motif_search(genome_encoded, fwd_mask)
            rev_count = self.motif_search(genome_encoded, rev_mask)

            fwd_stats.proportion_test[motif] = MotifObservation(
                observed=fwd_count,
                enzyme=enzyme,
                organism=info.get("organism")
            )

            rev_stats.proportion_test[motif] = MotifObservation(
                observed=rev_count,
                enzyme=enzyme,
                organism=info.get("organism")
            )

        return EntryResults(
            entry=entry_name,
            genome_length=genome_length,
            forward=fwd_stats,
            reverse=rev_stats
            )
