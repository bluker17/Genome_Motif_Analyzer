#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
from statsmodels.stats.proportion import proportions_ztest

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

class Statistics:
    """
    Performs one-sided proportion tests on motif occurrences per chromosome.

    This class is stateless and operates on chromosome-level dictionaries
    produced by the Motifs module.
    """

    def __init__(self) -> None:
        """
        Initializes the Statistics processor.
        """
        pass

    def expected_prob(self, base_probs: dict[str, float], motif: str) -> float:
        """
        Compute expected probability of a motif given base probabilities.

        Parameters
        ----------
        base_probs : dict[str, float]
            Base probabilities for A/T/G/C.
        motif : str
            Motif sequence (may include degenerate bases).
        degenerate_map : dict[str, list[str]]
            IUPAC degeneracy map.

        Returns
        -------
        float
            Expected probability of motif occurrence.
        """
        prob = 1.0

        for base in motif.upper():
            allowed_bases = degenerate_map[base]
            prob *= sum(base_probs[b] for b in allowed_bases)

        return prob

    def run_proportion_test(self, chrom_stats: dict) -> dict:
            """
            Apply one-sided z-test for motif enrichment/depletion.

            Parameters
            ----------
            chrom_stats : dict
                Output from Motifs.process_chromosome().
            degenerate_map : dict[str, list[str]]
                IUPAC degeneracy map.

            Returns
            -------
            dict
                Updated chrom_stats with statistical results added.
            """

            genome_length = chrom_stats["genome_length"]

            for strand in ("forward", "reverse"):

                if strand not in chrom_stats:
                    continue
                
                base_probs = chrom_stats[strand]["base_probs"]

                for motif, data in chrom_stats[strand]["proportion_test"].items():

                    observed = data["observed"]
                    possible_positions = genome_length - len(motif) + 1

                    expected_prob = self.expected_prob(base_probs=base_probs, motif=motif)

                    # Edge case protection
                    if possible_positions <= 0 or expected_prob in [0, 1]:
                        data.update({
                            "z_stat": math.nan,
                            "p_value": math.nan,
                            "significance": math.nan,
                            "expected_count": math.nan,
                            "total_positions": possible_positions,
                            "expected_motif_prob": expected_prob
                        })
                        continue

                    stat, pval = proportions_ztest(
                        count=observed,
                        nobs=possible_positions,
                        value=expected_prob
                    )

                    data.update({
                        "z_stat": float(stat),
                        "p_value": float(pval),
                        "significance": float(np.sign(stat)),
                        "expected_count": possible_positions * expected_prob,
                        "total_positions": possible_positions,
                        "expected_motif_prob": expected_prob
                    })

            return chrom_stats