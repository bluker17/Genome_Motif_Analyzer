#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
from statsmodels.stats.proportion import proportions_ztest

from src.alphabet.macromolecule_alphabet import Alphabet
from src.alphabet.result_alphabet import EntryResults, StrandResults, MotifObservation

class Statistics:
    """
    Performs one-sided proportion tests on motif occurrences per chromosome.

    This class is stateless and operates on chromosome-level dictionaries
    produced by the Motifs module.
    """

    def __init__(self, alphabet: Alphabet) -> None:
        """
        Initializes the Statistics processor.
        """
        self.alphabet = alphabet

    def expected_prob(self, base_probs: dict[str, float], motif: str) -> float:
        """
        Compute expected probability of a motif given base probabilities.

        Parameters
        ----------
        base_probs : dict[str, float]
            Base probabilities for A/T or U/G/C.
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
            allowed_bases = self.alphabet.degenerate_map[base]
            prob *= sum(base_probs[b] for b in allowed_bases)

        return prob

    def run_proportion_test(self, entry_stats: EntryResults) -> EntryResults:
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

            genome_length = entry_stats.genome_length

            for strand in [entry_stats.forward, entry_stats.reverse]:
                base_probs = strand.base_probs

                # print(strand.proportion_test.items())

                for motif, data in strand.proportion_test.items():
                    # print(f"Motif - {motif}")
                    # print(f"Data - {data}")
                    observed = data.observed
                    possible_positions =  genome_length - len(motif) + 1

                    expected_prob = self.expected_prob(base_probs=base_probs, motif=motif)

                    # Edge case protection
                    if possible_positions <= 0 or expected_prob in [0, 1]:
                        data.z_stat = math.nan
                        data.p_value = math.nan
                        data.significance = math.nan
                        data.expected_count = math.nan
                        data.total_positions = possible_positions
                        data.expected_motif_prob = expected_prob
                        continue

                    stat, pval = proportions_ztest(
                        count=observed,
                        nobs=possible_positions,
                        value=expected_prob
                    )
                     
                    data.z_stat = float(stat)
                    data.p_value = float(pval)
                    data.significance = "+" if float(np.sign(stat)) > 0 else "-" if float(np.sign(stat)) < 0 else 0
                    data.expected_count = possible_positions * expected_prob
                    data.total_positions = possible_positions
                    data.expected_motif_prob = expected_prob

            return entry_stats