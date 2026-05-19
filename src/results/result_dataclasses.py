# src/results/models.py
from dataclasses import dataclass
from typing import Optional

@dataclass(slots=True)
class MotifResult:
    motif: str
    enzyme: str
    organism: Optional[str]
    observed: int
    z_stat: float | None = None
    p_value: float | None = None
    significance: str | None = None
    expected_count: float | None = None
    total_positions: int | None = None
    expected_motif_prob: float | None = None

@dataclass(slots=True)
class StrandResult:
    base_probs: dict[str, float]
    gc_content: float
    motifs: list[MotifResult]

@dataclass(slots=True)
class ChromosomeResult:
    chromosome: str
    genome_length: int
    forward: StrandResult | None
    reverse: StrandResult | None