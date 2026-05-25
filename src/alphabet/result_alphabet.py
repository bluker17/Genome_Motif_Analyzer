from dataclasses import dataclass, field
import math

@dataclass(slots=True)
class MotifObservation:
    "Stores the search results for a specific motif."
    observed: int
    enzyme: str
    organism: str | None = None
    z_stat: float = math.nan
    p_value: float = math.nan
    significance: str | float = math.nan
    expected_count: float = math.nan
    total_positions: int = 0
    expected_motif_prob: float = 0.0


@dataclass(slots=True)
class StrandResults:
    "Contains the base probability, and proportion test information per strand."
    base_probs: dict[str, float]
    GC_content: float
    proportion_test:dict[str, MotifObservation] = field(default_factory=dict)

@dataclass(slots=True)
class EntryResults:
    "Container for entry level results"
    entry: str
    genome_length: int
    forward: StrandResults | None = None
    reverse: StrandResults | None = None