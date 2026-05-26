#!/usr/bin/env bash

echo "This is a test run of the program that searches the forward, reverse, and both strands of an RNA bacteriophage genomes for the motifs found in 'test_RNA_motifs.csv'." 

strands=(forward reverse both)

for s in "${strands[@]}"; do
    /usr/bin/time ./main.py \
        -f ./testing_materials/example_data/RNA_test_phage \
        -m ./testing_materials/example_data/test_RNA_motifs.csv \
        -c ./testing_materials/example_outputs/RNA_test/RNA_phage_${s}_strand_test.csv \
        --macromolecule RNA \
        -s both
done