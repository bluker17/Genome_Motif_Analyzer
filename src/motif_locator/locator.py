#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Iterator, Dict, Tuple
from Bio import SeqIO
import ahocorasick
from collections import Counter


RC_MAP: dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N"
}


class AhoCorasick_Motif_Search:
    """
    Genome-wide motif detection using Aho–Corasick multi-pattern matching.
    """

    def __init__(self, motif_info: dict, genomes: list[Path], strand_to_search: str) -> None:
        """
        Initialize automata and metadata lookup.
        """
        self.motif_info = motif_info
        self.genomes = genomes
        self.strand_to_search = strand_to_search

        self.motif_meta: dict[str, dict] = {
            v["motif_sequence"]: {
                "enzyme": v["enzyme"],
                "organism": v.get("organism")
            }
            for v in motif_info.values()
        }

        self.automata: dict[str, ahocorasick.Automaton] = {}
        self._build_automata()

    def _build_automata(self) -> None:
        """
        Build Aho–Corasick automata for forward/reverse strands.
        """
        motifs = [v["motif_sequence"] for v in self.motif_info.values()]

        if self.strand_to_search in ("forward", "both"):
            self.automata["forward"] = self._build_automaton(motifs)

        if self.strand_to_search in ("reverse", "both"):
            rev = [self._revcomp(m) for m in motifs]
            self.automata["reverse"] = self._build_automaton(rev)

    def _build_automaton(self, motifs: list[str]) -> ahocorasick.Automaton:
        """
        Build AC automaton from motif list.
        """
        A = ahocorasick.Automaton()
        for m in motifs: A.add_word(m, m)
        A.make_automaton()
        return A

    def _revcomp(self, seq: str) -> str:
        """
        Reverse complement sequence.
        """
        return "".join(RC_MAP[b] for b in seq[::-1])

    def stream_fasta(self, file: Path) -> Iterator[Tuple[str, str]]:
        """
        Stream FASTA records.
        """
        for record in SeqIO.parse(file, "fasta"):
            yield record.description, str(record.seq)

    def compute_base_probs(self, seq: str) -> dict[str, float]:
        """
        Compute nucleotide frequencies once per chromosome.
        """
        c = Counter(seq)
        total = sum(c[b] for b in "ACGT")
        if total == 0:
            return {b: 0.0 for b in "ACGT"}
        return {b: c[b] / total for b in "ACGT"}

    def reverse_base_probs(self, p: dict[str, float]) -> dict[str, float]:
        """
        Convert forward base probabilities to reverse strand.
        """
        return {"A": p["T"], "T": p["A"], "G": p["C"], "C": p["G"]}

    def process_chromosome(self, chrom_name: str, sequence: str) -> dict:
        """
        Run AC-based motif search on one chromosome.
        """
        sequence = sequence.upper()

        base_probs = self.compute_base_probs(sequence)
        gc_content = base_probs["G"] + base_probs["C"]
        genome_length = len(sequence)

        chrom_stats: dict = {
            "chromosome": chrom_name,
            "genome_length": genome_length,
        }

        for strand, automaton in self.automata.items():

            counts: dict[str, int] = {}

            for _, motif in automaton.iter(sequence):
                counts[motif] = counts.get(motif, 0) + 1

            full_results: dict[str, dict] = {}

            for motif, meta in self.motif_meta.items():

                observed = counts.get(motif, 0)

                full_results[motif] = {
                    "observed": observed,
                    "enzyme": meta["enzyme"],
                    "organism": meta.get("organism")
                }

            chrom_stats[strand] = {
                "base_probs": base_probs if strand == "forward"
                else self.reverse_base_probs(base_probs),
                "GC_content": gc_content,
                "proportion_test": full_results,
            }

        return chrom_stats