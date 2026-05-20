# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# from pathlib import Path
# from typing import Iterator, Tuple
# from Bio import SeqIO
# from collections import Counter

# from src.macromolecule_alphabet.alphabet import Alphabet
# from src.results.result_dataclasses import MotifResult, StrandResult, ChromosomeResult

# import ahocorasick

# class AhoCorasick_Motif_Search:
#     """
#     Genome-wide motif detection using Aho–Corasick multi-pattern matching.
#     """

#     def __init__(self, motif_info: dict, genomes: list[Path], strand_to_search: str, alphabet: Alphabet) -> None:
#         """
#         Initialize automata and metadata lookup.
#         """
#         self.motif_info = motif_info
#         self.genomes = genomes
#         self.strand_to_search = strand_to_search
#         self.alphabet = alphabet

#         self.motif_meta = {}

#         for value in motif_info.values():

#             # case 1: already a dict (single record)
#             if isinstance(value, dict):
#                 value = [value]

#             # case 2: list of records
#             for v in value:
#                 self.motif_meta[v["motif_sequence"]] = {
#                     "enzyme": v["enzyme"],
#                     "organism": v.get("organism")
#                 }


#         self.automata: dict[str, ahocorasick.Automaton] = {}
#         self._build_automata()

#     def _build_automata(self) -> None:
#         """
#         Build Aho–Corasick automata for forward/reverse strands.
#         """
#         # motifs = [v["motif_sequence"] for v in self.motif_info.values()]
#         motifs = [
#             record["motif_sequence"]
#             for records in self.motif_info.values()
#             for record in records
# ]

#         if self.strand_to_search in ("forward", "both"):
#             self.automata["forward"] = self._build_automaton(motifs)

#         if self.strand_to_search in ("reverse", "both"):
#             rev = [self._revcomp(m) for m in motifs]
#             self.automata["reverse"] = self._build_automaton(rev)

#     def _build_automaton(self, motifs: list[str]) -> ahocorasick.Automaton:
#         """
#         Build AC automaton from motif list.
#         """
#         A = ahocorasick.Automaton()
#         for m in motifs: A.add_word(m, m)
#         A.make_automaton()
#         return A

#     def _revcomp(self, seq: str) -> str:
#         """
#         Reverse complement sequence.
#         """
#         return "".join(self.alphabet.complement_map[b] for b in reversed(seq))

#     def stream_fasta(self, file: Path) -> Iterator[Tuple[str, str]]:
#         """
#         Stream FASTA records.
#         """
#         for record in SeqIO.parse(file, "fasta"):
#             yield record.description, str(record.seq)

#     def compute_base_probs(self, seq: str) -> dict[str, float]:
#         """
#         Compute nucleotide frequencies once per chromosome.
#         """
#         c = Counter(seq)
#         total = sum(c[b] for b in self.alphabet.bases)
#         if total == 0:
#             return {b: 0.0 for b in self.alphabet.bases}
#         return {b: c[b] / total for b in self.alphabet.bases}

#     def reverse_base_probs(self, p: dict[str, float]) -> dict[str, float]:
#         """
#         Convert forward base probabilities to reverse strand.
#         """
#         return {base: p[self.alphabet.complement_map[base]] for base in self.alphabet.bases}

#     def process_chromosome(self, chrom_name: str, sequence: str) -> ChromosomeResult:
#         sequence = sequence.upper()

#         base_probs = self.compute_base_probs(sequence)
#         gc_content = base_probs["G"] + base_probs["C"]
#         genome_length = len(sequence)

#         strand_results = {}

#         for strand, automaton in self.automata.items():

#             counts: dict[str, int] = {}

#             for _, motif in automaton.iter(sequence):
#                 counts[motif] = counts.get(motif, 0) + 1

#             motif_objects: list[MotifResult] = []

#             for motif, meta in self.motif_meta.items():

#                 motif_objects.append(
#                     MotifResult(
#                         motif=motif,
#                         enzyme=meta["enzyme"],
#                         organism=meta.get("organism"),
#                         observed=counts.get(motif, 0)
#                     )
#                 )

#             strand_results[strand] = StrandResult(
#                 base_probs=base_probs if strand == "forward"
#                 else self.reverse_base_probs(base_probs),
#                 gc_content=gc_content,
#                 motifs=motif_objects
#             )

#             print("Example counts:", list(counts.items())[:10])
#             print("Motif meta keys:", list(self.motif_meta.keys())[:10])

#         return ChromosomeResult(
#             chromosome=chrom_name,
#             genome_length=genome_length,
#             forward=strand_results.get("forward"),
#             reverse=strand_results.get("reverse")
#         )


from pathlib import Path
from typing import Iterator, Tuple, Dict, List
from Bio import SeqIO
from collections import Counter
import ahocorasick

from alphabet.macromolecule_alphabet import Alphabet
from src.results.result_dataclasses import MotifResult, StrandResult, ChromosomeResult


class AhoCorasick_Motif_Search:
    """
    Strand-safe genome-wide motif detection using Aho–Corasick with IUPAC expansion.
    """

    # ----------------------------
    # IUPAC expansion table
    # ----------------------------
    IUPAC: Dict[str, List[str]] = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "N": ["A", "C", "G", "T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
    }

    def __init__(self, motif_info: dict, genomes: list[Path], strand_to_search: str, alphabet: Alphabet) -> None:
        self.motif_info = motif_info
        self.genomes = genomes
        self.strand_to_search = strand_to_search
        self.alphabet = alphabet

        # canonical motif metadata
        self.motif_meta: Dict[str, dict] = {}

        motifs: List[str] = []
        self.expanded_to_canonical: Dict[str, str] = {}

        for value in motif_info.values():
            if isinstance(value, dict):
                value = [value]

            for v in value:
                motif = v["motif_sequence"].upper()

                self.motif_meta[motif] = {
                    "enzyme": v["enzyme"],
                    "organism": v.get("organism")
                }

                motifs.append(motif)

        # build expanded motif universe
        expanded_motifs = self._expand_all_motifs(motifs)

        self.automata: Dict[str, ahocorasick.Automaton] = {}
        self._build_automata(expanded_motifs)

    # ----------------------------
    # IUPAC expansion
    # ----------------------------
    def _expand_iupac(self, seq: str) -> List[str]:
        pools = [self.IUPAC[b] for b in seq]
        results = [""]

        for pool in pools:
            results = [prefix + b for prefix in results for b in pool]

        return results

    def _expand_all_motifs(self, motifs: List[str]) -> List[str]:
        expanded = []

        for motif in motifs:
            variants = self._expand_iupac(motif)

            for v in variants:
                expanded.append(v)
                self.expanded_to_canonical[v] = motif

        return expanded

    # ----------------------------
    # Automaton construction
    # ----------------------------
    def _build_automata(self, expanded_motifs: List[str]) -> None:
        if self.strand_to_search in ("forward", "both"):
            self.automata["forward"] = self._build_automaton(expanded_motifs)

        if self.strand_to_search in ("reverse", "both"):
            self.automata["reverse"] = self._build_automaton(expanded_motifs)

    def _build_automaton(self, motifs: List[str]) -> ahocorasick.Automaton:
        A = ahocorasick.Automaton()

        for m in motifs:
            # store expanded motif as key, canonical motif as value
            canonical = self.expanded_to_canonical[m]
            A.add_word(m, canonical)

        A.make_automaton()
        return A

    # ----------------------------
    # Sequence utilities
    # ----------------------------
    def _revcomp(self, seq: str) -> str:
        return "".join(self.alphabet.complement_map[b] for b in reversed(seq))

    def stream_fasta(self, file: Path) -> Iterator[Tuple[str, str]]:
        for record in SeqIO.parse(file, "fasta"):
            yield record.description, str(record.seq)

    # ----------------------------
    # Stats
    # ----------------------------
    def compute_base_probs(self, seq: str) -> dict[str, float]:
        seq = seq.upper()
        c = Counter(seq)

        total = sum(c[b] for b in self.alphabet.bases)
        if total == 0:
            return {b: 0.0 for b in self.alphabet.bases}

        return {b: c[b] / total for b in self.alphabet.bases}

    def reverse_base_probs(self, p: dict[str, float]) -> dict[str, float]:
        return {
            base: p[self.alphabet.complement_map[base]]
            for base in self.alphabet.bases
        }

    # ----------------------------
    # Core processing
    # ----------------------------
    def process_chromosome(self, chrom_name: str, sequence: str) -> ChromosomeResult:
        sequence = sequence.upper()

        base_probs = self.compute_base_probs(sequence)
        gc_content = base_probs["G"] + base_probs["C"]
        genome_length = len(sequence)

        strand_results = {}

        rev_sequence = self._revcomp(sequence)

        for strand, automaton in self.automata.items():

            seq_to_scan = sequence if strand == "forward" else rev_sequence

            counts: Dict[str, int] = Counter()

            for _, canonical_motif in automaton.iter(seq_to_scan):
                counts[canonical_motif] += 1

            motif_objects: List[MotifResult] = []

            for motif, meta in self.motif_meta.items():
                motif_objects.append(
                    MotifResult(
                        motif=motif,
                        enzyme=meta["enzyme"],
                        organism=meta.get("organism"),
                        observed=counts.get(motif, 0)
                    )
                )

            strand_results[strand] = StrandResult(
                base_probs=base_probs if strand == "forward"
                else self.reverse_base_probs(base_probs),
                gc_content=gc_content,
                motifs=motif_objects
            )

        return ChromosomeResult(
            chromosome=chrom_name,
            genome_length=genome_length,
            forward=strand_results.get("forward"),
            reverse=strand_results.get("reverse")
        )