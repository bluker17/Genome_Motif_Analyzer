#!/usr/bin/env nextflow

// Base process that searches for motifs and performs statistical analysis on the findings.
process Run_Analyzer {
    tag "Running motif analysis."
    publishDir "${file(params.csv_output).getParent()}", mode: 'copy'

    input:
    val csv_output
    path fasta_files
    path motif_file
    val strand
    val macromolecule

    output:
    path "${file(params.csv_output).getName()}", emit: initial_csv_output

    script:
    def output_filename = file(params.csv_output).getName()

    """
    main.py \\
        --csv_output ${output_filename} \\
        --fasta_files ${params.fasta_files} \\
        --motif_file ${params.motif_file} \\
        --strand ${params.strand} \\
        --macromolecule ${params.macromolecule}
    """
}