#!/bin/bash
python3 workflow/main.py \
    --biom_file results/251006_wrangling/M01_interp_biom_clr_prev20_abd0p001x5_mad20.csv \
    --meta_file results/251006_wrangling/M01_avg_meta_prev20_abd0p001x5_mad20.csv \
    --precision 10000 \
    --output_csv results/lagged_regression/251006_M01_test_10000_prev20_abd0p001x5_mad20.csv