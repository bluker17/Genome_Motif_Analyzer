#!/usr/bin/env bash

# Test 3
echo "This is a test run of the program that searches the forward and reverse strands of cattle genomes for the NNGRRN motif of SpCas9 from Streptococcus pyogenes and NGG from Staphylococcus aureus." 

/usr/bin/time ./main.py \
    -f ./testing_materials/example_data/cattle \
    -m ./testing_materials/example_data/combined_motif.csv \
    -c ./testing_materials/example_outputs/cattle_combined_test.csv \
    -s both