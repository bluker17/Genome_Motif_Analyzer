#!/usr/bin/env bash

# Test 4
echo "This is a test run of the program that searches forward, reverse, and both strands of bacteriophage genomes for the NNGRRN motif of SpCas9 from Streptococcus pyogenes and NGG from Staphylococcus aureus." 

strands=(forward reverse both)

for s in strands; do
    /usr/bin/time ./main.py \
        -f ./testing_materials/example_data/cattle \
        -m ./testing_materials/example_data/sp_motif.csv \
        -c ./testing_materials/example_outputs/bacteriophage_${s}_test.csv \
        -s ${s}
done