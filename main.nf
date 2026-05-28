#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {Run_Analyzer} from './modules/run_analyzer.nf'
include {Summarize_Stats} from './modules/summarize_stats.nf'

workflow {
    main:
    // Collect the files from the parameters
    fasta_dir_ch = file(params.fasta_files)
    motif_file_ch = file(params.motif_file)

    // Run the pipeline
    analysis_outputs = Run_Analyzer(params.csv_output, fasta_dir_ch, motif_file_ch, params.strand, params.macromolecule)
    Summarize_Stats(analysis_outputs.initial_csv_output)
}

