#!/usr/bin/env bash

# Test 1
echo "This is a test run of the program that searches the forward and reverse strands of cattle genomes for the NGG motif of SpCas9 from Streptococcus pyogenes." 

/usr/bin/time ./main.py \
    -f ./testing_materials/example_data/cattle \
    -m ./testing_materials/example_data/sp_motif.csv \
    -c ./testing_materials/example_outputs/cattle_sp_test.csv \
    -s both