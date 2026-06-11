#!/bin/bash
python3 workflow/mibi/time_series_data_wrangling_251006_update.py \
  -i data/dada2_otu_table_w_tax_no_pynast_failures_rare7500.biom \
  -m data/metadata_full_w_times.csv \
  -n results/251006_wrangling \
  -o results/251006_wrangling \
  -s M01 \
  --prevalence 0.20 \
  --min_rel 0.001 \
  --min_samples 5 \
  --mad_quantile 0.2