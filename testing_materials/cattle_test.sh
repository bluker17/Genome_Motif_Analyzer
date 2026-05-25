#!/usr/bin/env bash

echo "This is a test run of the program that searches the forward and reverse strands of Bos taurus genome for the motifs found in 'test_cattle_motifs.csv'."

/usr/bin/time ./main.py \
    -f ./testing_materials/example_data/cattle_test_genome \
    -m ./testing_materials/example_data/test_cattle_motifs.csv \
    -c ./testing_materials/example_outputs/cattle_test/cattle_both_strands_test.csv \
    --macromolecule DNA \
    -s both