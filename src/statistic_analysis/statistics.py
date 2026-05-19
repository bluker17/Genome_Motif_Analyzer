#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import numpy as np
from src.macromolecule_alphabet.alphabet import Alphabet

from statsmodels.stats.proportion import proportions_ztest
from src.results.result_dataclasses import MotifResult, ChromosomeResult


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
            allowed_bases = self.alphabet.degenerate_map[base]
            prob *= sum(base_probs[b] for b in allowed_bases)

        return prob

    def run_proportion_test(self, chrom_stats: ChromosomeResult) -> ChromosomeResult:
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

            genome_length = chrom_stats.genome_length

            for strand in (chrom_stats.forward, chrom_stats.reverse):
                if strand is None:
                    continue

                base_probs = strand.base_probs

                for motif_obj in strand.motifs:

                    observed = motif_obj.observed
                    motif = motif_obj.motif

                    possible_positions = genome_length - len(motif) + 1
                    expected_prob = self.expected_prob(base_probs, motif)

                    # Edge case protection
                    if possible_positions <= 0 or expected_prob in (0, 1):
                        motif_obj.z_stat = math.nan
                        motif_obj.p_value = math.nan
                        motif_obj.significance = ("+" if motif_obj.z_stat > 0 else "-" if motif_obj.z_stat < 0 else None)
                        motif_obj.expected_count = math.nan
                        motif_obj.total_positions = possible_positions
                        motif_obj.expected_motif_prob = expected_prob
                        continue

                    stat, pval = proportions_ztest(
                        count=observed,
                        nobs=possible_positions,
                        value=expected_prob
                    )

                    motif_obj.z_stat = float(stat)
                    motif_obj.p_value = float(pval)
                    motif_obj.significance = ("+" if motif_obj.z_stat > 0 else "-" if motif_obj.z_stat < 0 else None)
                    motif_obj.expected_count = possible_positions * expected_prob
                    motif_obj.total_positions = possible_positions
                    motif_obj.expected_motif_prob = expected_prob

            return chrom_stats