#!/usr/bin/env bash

# Test 1
echo "This is a test run of the program that searches cattle genomes for the NGG motif of SpCas9 from Streptococcus pyogenes." 

python3 main.py \
    -f testing_materials/example_data/cattle \
    -m testing_materials/example_data/motif.csv \
    -c testing_materials/example_outputs/cattle_test.csv

# Test 2
# echo "This is a test run of th eprogram that searches a prokaryotic genome for the _________ motif of ________ from __________."

# python3 main.py \
#     -f testing_materials/example_data/cattle \
#     -m testing_materials/example_data/motif.csv \
#     -c testing_materials/example_outputs/cattle_test.csv