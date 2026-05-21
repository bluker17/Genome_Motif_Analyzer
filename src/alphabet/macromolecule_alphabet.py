#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from functools import cached_property

@dataclass(frozen=True)
class Alphabet:
    """
    Defines biological sequence rules.
    """

    name: str

    # Canonical bases
    bases: tuple[str, ...]

    gc_bases: tuple[str, ...]

    # Complement lookup
    complement_map: dict[str, str]

    # Degenerate base expansion
    degenerate_map: dict[str, tuple[str, ...]]

    # Canonical base bit encodings
    bit_map: dict[str, int]

    has_reverse_complement: bool = True

    @cached_property
    def base_set(self) -> set[str]:
        return set(self.bases)

    @cached_property
    def all_symbols(self) -> set[str]:
        return set(self.degenerate_map)

    @cached_property
    def mask_map(self) -> dict[str, int]:
        """
        Convert all symbols into bitmasks.
        """

        return {
            symbol: sum(self.bit_map[b] for b in bases)
            for symbol, bases in self.degenerate_map.items()
        }

    @cached_property
    def ascii_mask(self):
        """
        ASCII lookup table for ultra-fast encoding.
        """

        import numpy as np

        arr = np.zeros(256, dtype=np.uint8)

        for base, mask in self.mask_map.items(): arr[ord(base)] = mask

        return arr

    def reverse_complement(self, sequence: str) -> str:
        """
        Reverse complement sequence.
        """

        return "".join(self.complement_map[b] for b in reversed(sequence.upper()))
    

DNA = Alphabet(
    name="DNA",
    bases=("A", "T", "G", "C"),
    gc_bases=("G", "C"),

    complement_map={
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
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
        "A": ('A',),
        'T': ('T',),
        'G': ('G',),
        'C': ('C',),
        'N': ('A', 'T', 'G', 'C'),
        'R': ('A', 'G'),
        'Y': ('C', 'T'),
        'W': ('A', 'T'),
        'S': ('G', 'C'),
        'K': ('G', 'T'),
        'M': ('A', 'C'),
        'B': ('C', 'G', 'T'),
        'D': ('A', 'G', 'T'),
        'H': ('A', 'C', 'T'),
        'V': ('A', 'C', 'G')
    },

    bit_map={
        "A": 1,
        "C": 2,
        "G": 4,
        "T": 8,
    }
)

RNA = Alphabet(
    name="RNA",
    bases=("A", "U", "G", "C"),
    gc_bases=("G", "C"),
    
    complement_map={
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
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
        "A": ('A',),
        'U': ('U',),
        'G': ('G',),
        'C': ('C',),
        'N': ('A', 'U', 'G', 'C'),
        'R': ('A', 'G'),
        'Y': ('C', 'U'),
        'W': ('A', 'U'),
        'S': ('G', 'C'),
        'K': ('G', 'U'),
        'M': ('A', 'C'),
        'B': ('C', 'G', 'U'),
        'D': ('A', 'G', 'U'),
        'H': ('A', 'C', 'U'),
        'V': ('A', 'C', 'G')
    },

    bit_map={
        "A": 1,
        "C": 2,
        "G": 4,
        "U": 8,
    }
)

ALPHABETS = {
    "DNA": DNA,
    "RNA": RNA
}

