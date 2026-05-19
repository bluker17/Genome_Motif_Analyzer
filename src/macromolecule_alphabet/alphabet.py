#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass


@dataclass(frozen=True)
class Alphabet:
    """
    Defines the biological rules for a sequence type.
    """

    name: str
    bases: tuple[str, ...]
    complement_map: dict[str, str]
    degenerate_map: dict[str, list[str]]
    has_reverse_complement: bool = True

    @property
    def base_set(self) -> set[str]:
        """
        Fast lookup set for sequence validation.
        """
        return set(self.bases)
    
DNA = Alphabet(
    name="DNA",
    bases=("A", "T", "G", "C"),
    complement_map={
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",

        # IUPAC degenerate complements
        "R": "Y",
        "Y": "R",
        "W": "W",
        "S": "S",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N"
    },
    degenerate_map={
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
        'V': ['A', 'C', 'G']}
)

RNA = Alphabet(
    name="RNA",
    bases=("A", "U", "G", "C"),
    complement_map={
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",

        # IUPAC degenerate complements
        "R": "Y",
        "Y": "R",
        "W": "W",
        "S": "S",
        "K": "M",
        "M": "K",
        "B": "V",
        "V": "B",
        "D": "H",
        "H": "D",
        "N": "N"
    },


    degenerate_map={
        'A': ['A'],
        'U': ['U'],
        'G': ['G'],
        'C': ['C'],
        'N': ['A', 'U', 'G', 'C'],
        'R': ['A', 'G'],
        'Y': ['C', 'U'],
        'W': ['A', 'U'],
        'S': ['G', 'C'],
        'K': ['G', 'U'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'U'],
        'D': ['A', 'G', 'U'],
        'H': ['A', 'C', 'U'],
        'V': ['A', 'C', 'G']}
)

ALPHABETS = {
    "DNA": DNA,
    "RNA": RNA
}