#!/bin/bash

python3 ../../main.py \
    --biom_file ../laurie_additions/bangladesh_casey_otu.biom \
    --meta_file ../laurie_additions/metadata_full_w_times_modified.tsv \
    --tree_file ../laurie_additions/dada2.tre \
    --subject M01 \
    --filter_num 1000 \
    --rared 7500 \
    --precision 1000000 \
    --output permutation_results.csv
