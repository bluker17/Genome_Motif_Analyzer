#!/usr/bin/env bash

echo "This is a test run of the program that searches the forward, reverse, and both strands of DNA bacteriophage genomes for the motifs found in 'test_RNA_motifs.csv'." 

strands=(forward reverse both)

for s in "${strands[@]}"; do
    /usr/bin/time ./main.py \
        -f ./testing_materials/example_data/DNA_test_phages \
        -m ./testing_materials/example_data/test_DNA_motifs.csv \
        -c ./testing_materials/example_outputs/DNA_test/DNA_phage_${s}_strand_test.csv \
        --macromolecule DNA \
        -s both
done