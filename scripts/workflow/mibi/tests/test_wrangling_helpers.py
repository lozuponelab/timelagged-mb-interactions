# tests/test_wrangling_helpers.py
import numpy as np
import pandas as pd
import pytest

# Import the functions from your script
from bangladesh_tlc_organized.workflow.mibi.time_series_data_wrangling_251006_update import (
    avg_sample_by_day_clr,
    interp_missing_day_clr,
)

def test_avg_sample_by_day_clr_means_correctly():
    # Make a CLR table with two samples on the same day
    # Rows=OTUs, Cols=Samples (CLR)
    clr = pd.DataFrame(
        {"S1": [ 0.6, -0.3, -0.3],   # day=1
         "S2": [ 0.4, -0.2, -0.2],   # day=1 (replicate)
         "S3": [-0.1,  0.05, 0.05]}, # day=2
        index=["OTU1","OTU2","OTU3"]
    )
    meta = pd.DataFrame({
        "#SampleID": ["S1","S2","S3"],
        "ANONYMIZED_NAME": ["M01","M01","M01"],
        "epoch_time": [1,1,2]
    })
    avg_meta, avg = avg_sample_by_day_clr(meta, clr)
    # Expect two day-columns: 1 and 2
    assert list(avg.columns) == [1,2]
    # Day 1 column should be the mean of S1 and S2
    expected_day1 = (clr["S1"] + clr["S2"]) / 2.0
    assert np.allclose(avg[1].values, expected_day1.values, atol=1e-12)


def test_interp_missing_day_clr_fills_internal_gap(daily_clr_small):
    # Convert day indices (1,2,3) → epoch seconds at daily spacing (0, 86400, 172800)
    df = daily_clr_small.copy()
    df.columns = (df.columns - 1) * 86400  # day 1 -> 0s, day 2 -> 86400s, day 3 -> 172800s

    # Remove the middle day to create a gap between day0 and day2
    df_gap = df[[0, 172800]]
    out = interp_missing_day_clr(df_gap)

    # Should now include 0, 86400, 172800
    assert list(out.columns) == [0, 86400, 172800]

    # Interpolated middle should be between neighbors for each OTU
    mid = out[86400]
    d0  = out[0]
    d2  = out[172800]
    for v0, vm, v2 in zip(d0.values, mid.values, d2.values):
        assert (min(v0, v2) - 1e-9) <= vm <= (max(v0, v2) + 1e-9)