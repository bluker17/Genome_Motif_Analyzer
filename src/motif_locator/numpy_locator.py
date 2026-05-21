#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from collections import Counter
from typing import Iterator
from Bio import SeqIO

from src.alphabet.macromolecule_alphabet import Alphabet
from src.alphabet.result_alphabet import EntryResults, StrandResults, MotifObservation

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

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

    def __init__(self, motif_info: dict, genomes: list[Path], alphabet: Alphabet) -> None:

        self.motif_info = motif_info
        self.genomes = genomes
        self.alphabet = alphabet

        self.motif_results: dict = {}

        self.motif_masks: dict[str, dict[str, list[np.ndarray]]] = {
            enzyme: {
                "forward": self.build_motif_masks(info["motif_sequence"]),
                "reverse": self.build_motif_masks(
                        self.alphabet.reverse_complement(
                        info["motif_sequence"]
                        )
                    )
                }
            for enzyme, info in motif_info.items()
        }

    def stream_fasta(self, file: Path) -> Iterator[tuple[str, str]]:
        """
        Stream FASTA records one entry at a time.

        Parameters
        ----------
        file : Path
            FASTA file path.

        Yields
        ------
        Iterator[tuple[str, str]]
            Entry name and sequence.
        """

        for record in SeqIO.parse(file, "fasta"):
            yield record.description, str(record.seq)

    def seq_to_array(self, sequence: str) -> np.ndarray:
        """
        Convert genetic sequence into compact NumPy byte array.

        Parameters
        ----------
        sequence : str
            Genetic sequence.

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

            for base in self.alphabet.degenerate_map[char]: mask[ord(base)] = True

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
            Genetic sequence.

        Returns
        -------
        dict[str, float]
            Base probabilities for A/T or U/G/C.
        """

        seq = seq.upper()

        counts = Counter(seq)

        total = sum(counts[b] for b in self.alphabet.bases)

        if total == 0:
            return {b: 0.0 for b in self.alphabet.bases}

        return {b: counts[b] / total for b in self.alphabet.bases}

    # def reverse_base_probs(self, forward_probs: dict[str, float]) -> dict[str, float]:
    #     """
    #     Convert forward strand probabilities into reverse strand probabilities.

    #     Parameters
    #     ----------
    #     forward_probs : dict[str, float]
    #         Forward strand probabilities.

    #     Returns
    #     -------
    #     dict[str, float]
    #         Reverse strand probabilities.
    #     """

    #     return {base: forward_probs[self.alphabet.complement_map[base]] for base in self.alphabet.bases}
    
    def process_entry(self, entry_name: str, sequence: str) -> EntryResults:
        """
        Process a FASTA entry for motif occurrence statistics.

        Parameters
        ----------
        chrom_name : str
            Entry identifier.
        sequence : str
            Entry sequence.

        Returns
        -------
        dict
            Entry motif statistics.
        """

        sequence = sequence.upper()

        genome_arr = self.seq_to_array(sequence)

        genome_length = len(sequence)

        base_probs = self.compute_base_probs(sequence)
        # rev_base_probs = self.reverse_base_probs(fwd_base_probs)

        fwd_stats = StrandResults(
            base_probs=base_probs,
            GC_content=sum(base_probs[b] for b in self.alphabet.gc_bases)
        )

        rev_stats = StrandResults(
            base_probs=base_probs,
            GC_content=sum(base_probs[b] for b in self.alphabet.gc_bases)
        )

        for enzyme, info in self.motif_info.items():

            motif = info["motif_sequence"]

            fwd_masks = self.motif_masks[enzyme]["forward"]
            rev_masks = self.motif_masks[enzyme]["reverse"]

            _, fwd_count = self.motif_search(genome_arr, fwd_masks)
            _, rev_count = self.motif_search(genome_arr, rev_masks)

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