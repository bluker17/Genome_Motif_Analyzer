#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path

class SummaryStats:
    """
    Initializes the SummaryStats class with input and output file paths.
    """
    def __init__(self, input_file: Path):
        self.input_file = input_file
        self.output_file = f"{input_file.stem}_summary_statistics.csv"
        self.findings = None
        return
    
    def summary_findings(self):
        """
        Calculates and exports summary statistics for the dataset.

        Groups the input data by 'FASTA file', 'Strand', and 'Enzyme'. For each group, it calculates the entry count, as well as the mean and median for both 'p-value' and 'z-stat' columns. The resulting DataFrame is flattened and saved to a new CSV file.
        """
        df = pd.read_csv(self.input_file)
        groups = df.groupby(['FASTA File', 'Strand', 'Enzyme'])
        self.findings = groups.agg({
                'FASTA File': 'count',
                'p-value':['mean','median'],
                'z-stat':['mean','median']
        }).reset_index()

        self.findings.columns = [
            'FASTA File', 'Strand', 'Enzyme', 
            'FASTA Entry Count', 
            'p-value_mean', 'p-value_median', 
            'z-stat_mean', 'z-stat_median'
        ]

        self.findings.to_csv(self.output_file, index=False)
        return