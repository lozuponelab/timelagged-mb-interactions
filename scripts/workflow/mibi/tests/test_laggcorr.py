# tests/test_lagcorr.py
import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from mibi.timeseries_0625_update import MiBiTimeSeries

@pytest.fixture
def tiny_tlc_data():
    """
    Build a small-but-not-too-small CLR matrix (T=24) with smooth variation.
    Longer T avoids FFT/IAAFT edge cases that produce NaNs on very short series.
    """
    T = 24
    times = np.arange(1, T + 1)

    # driver j: smooth sinusoid with trend -> non-constant, non-pathological spectrum
    j = 0.3 * np.sin(2 * np.pi * times / 12.0) + 0.02 * (times - T/2)

    # responder i: make Δi(t+1) somewhat anti-correlated with j(t)
    # construct i cumulatively from a base trend minus a scaled j plus small noise
    rng = np.random.default_rng(123)
    di = 0.03 * np.ones(T) - 0.1 * j + 0.01 * rng.normal(size=T)  # increments
    i = np.cumsum(di) - np.mean(np.cumsum(di))  # centered cumulative sum

    biom = pd.DataFrame(
        [i, j],
        index=["OTUi", "OTUj"],
        columns=times
    )  # rows=OTUs, cols=time (this matches MiBiTimeSeries expectation)
    meta = pd.DataFrame({
        "#SampleID": [str(t) for t in times],
        "ANONYMIZED_NAME": ["M01"] * T,
        "epoch_time": times
    })
    return meta, biom

def manual_tlc(i_series, j_series):
    # Δi(t+1)
    delta_i = i_series.shift(-1) - i_series
    j_aligned = j_series.loc[delta_i.index]
    mask = delta_i.notna() & j_aligned.notna()
    rho, _ = spearmanr(j_aligned[mask].values, delta_i[mask].values)
    return rho

def test_lag_regress_matches_manual(tiny_tlc_data):
    meta, biom = tiny_tlc_data
    ts = MiBiTimeSeries(meta, biom)
    # ts.biom has rows=OTUs, cols=sorted times
    i_series = ts.biom.loc["OTUi"]
    j_series = ts.biom.loc["OTUj"]

    rho_manual = manual_tlc(i_series, j_series)
    out = ts.lag_regress("OTUi", "OTUj", ts.biom)
    assert np.isfinite(out["rho"])
    assert abs(out["rho"] - rho_manual) < 1e-12

@pytest.mark.skipif(pytest.importorskip("pyunicorn", reason="pyunicorn not installed") is None, reason="no pyunicorn")
def test_iaaft_generator_and_one_surrogate(tiny_tlc_data):
    """
    Smoke test: ensure we can make a generator and compute at least one surrogate rho.
    Keep it FAST; we do not test significance here.
    """
    meta, biom = tiny_tlc_data
    ts = MiBiTimeSeries(meta, biom)
    ts.make_aaft_generator("OTUj", silence_level=2)
    ts.gen_null_distribution_iaaft("OTUi", "OTUj", precision=5, n_iterations=5)
    # Should have an array of 5 surrogate rhos
    arr = ts._null_distribution_iaaft
    assert arr.shape == (5,)
    assert np.isfinite(arr).all()