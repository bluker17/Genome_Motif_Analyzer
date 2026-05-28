#!/usr/bin/env nextflow

// Process that calculates summary statistics based upon the previous statistical analysis.
process Summarize_Stats {
    tag "Generating summary statistics."
    publishDir "${file(params.csv_output).getParent()}", mode: 'copy'

    input:
    path csv_output

    output:
    path "*_summary_statistics.csv", emit: summary_csv_output

    script:
    """
    run_summary.py --initial_csv ${csv_output}
    """
}