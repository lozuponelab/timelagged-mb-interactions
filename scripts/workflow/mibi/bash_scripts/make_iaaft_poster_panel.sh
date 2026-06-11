#!/bin/bash
#run from bangladesh_tlc_organized directory
python3 workflow/mibi/make_iaaft_poster_panel_for_pair.py \
  --biom_file "results/251006_wrangling/M01_interp_biom_clr_prev20_abd0p001x5_mad20.csv" \
  --meta_file "results/251006_wrangling/M01_avg_meta_prev20_abd0p001x5_mad20.csv" \
  --otu_j 4426438	\
  --otu_i 4476527 \
  --n_surr 10000 \
  --n_iter 20 \
  --outdir "plots/panels_pair_90_4426438_4476527_10k"