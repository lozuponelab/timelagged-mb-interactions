# tests/conftest.py
import numpy as np
import pandas as pd
import pytest

@pytest.fixture
def tiny_counts():
    """
    OTUs x samples count table. We pretend it's rarefied to 10 reads/sample
    to keep numbers small and human-checkable.
    """
    # 4 OTUs, 6 samples (two subjects, M01 & M02)
    data = pd.DataFrame(
        [[5, 6, 7, 2, 1, 2],   # OTU1
         [3, 2, 1, 6, 7, 6],   # OTU2
         [2, 1, 1, 2, 1, 1],   # OTU3 (rare)
         [0, 1, 1, 0, 1, 1]],  # OTU4 (sparse)
        index=["OTU1","OTU2","OTU3","OTU4"],
        columns=["S1","S2","S3","S4","S5","S6"]
    )
    # make column sums equal (simulate rarefaction): here sums=10 already
    return data

@pytest.fixture
def tiny_metadata():
    """
    Metadata including #SampleID, ANONYMIZED_NAME, epoch_time.
    We make S1..S3 belong to M01, S4..S6 to M02. Epoch_time strictly increases.
    """
    df = pd.DataFrame({
        "#SampleID": ["S1","S2","S3","S4","S5","S6"],
        "ANONYMIZED_NAME": ["M01","M01","M01","M02","M02","M02"],
        "epoch_time": [1, 2, 3, 1, 2, 3]  # simple ints for clarity
    })
    return df

@pytest.fixture
def biom_table(tiny_counts):
    """Build a BIOM Table from the tiny counts so biom2df can be tested."""
    biom = pytest.importorskip("biom")  # skip test if biom isn’t installed
    from biom.table import Table
    mat = tiny_counts.values
    obs_ids = list(tiny_counts.index)      # OTU IDs
    samp_ids = list(tiny_counts.columns)   # Sample IDs
    return Table(mat, obs_ids, samp_ids)

@pytest.fixture
def daily_clr_small():
    """
    A tiny daily CLR table for wrangling helper tests.
    Rows=OTUs, Cols=epoch days (ints). We pretend CLR is already done.
    """
    # Three days (1,2,3), three OTUs
    # Column means are ~0 in CLR (we enforce exact sums to 0 here)
    df = pd.DataFrame(
        [[ 0.4,  0.1, -0.2],
         [-0.2, -0.1,  0.3],
         [-0.2,  0.0, -0.1]],
        index=["OTU1","OTU2","OTU3"],
        columns=[1,2,3]
    )
    # Make column means exactly zero
    df = df.sub(df.mean(axis=0), axis=1)
    return df