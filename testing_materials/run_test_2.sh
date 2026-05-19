#!/usr/bin/env bash

# Test 2
echo "This is a test run of the program that searches the forward and reverse strands of cattle genomes for the NNGRRT motif of SaCas9 from Staphylococcus aureus." 

/usr/bin/time ./main.py \
    -f ./testing_materials/example_data/cattle \
    -m ./testing_materials/example_data/sa_motif.csv \
    -c ./testing_materials/example_outputs/cattle_sa_test.csv \
    -s both