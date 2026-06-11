# tests/test_biom_utils.py
import numpy as np
import pandas as pd
import pytest

from mibi import biom_utils_updated as U

def test_biom2df_roundtrip(biom_table, tiny_counts):
    df = U.biom2df(biom_table)
    assert list(df.index) == list(tiny_counts.index)
    assert list(df.columns) == list(tiny_counts.columns)
    assert (df.values == tiny_counts.values).all()

def test_clr_transform_centering():
    # Simple proportions (2 OTUs, 3 samples)
    props = pd.DataFrame(
        [[0.5, 0.2, 0.9],
         [0.5, 0.8, 0.1]],
        index=["A","B"],
        columns=["S1","S2","S3"]
    )
    clr = U.clr_transform(props, pseudocount=1e-12)
    # CLR columns should have mean ~ 0
    assert np.allclose(clr.mean(axis=0).values, 0.0, atol=1e-10)
    # If a column has equal parts (0.5,0.5), clr should be [0,0]
    assert np.allclose(clr["S1"].values, [0.0, 0.0], atol=1e-10)

def test_meta_biom_filter_preserves_subject_and_order(tiny_metadata, tiny_counts):
    md, biom = U.meta_biom_filter(tiny_metadata, tiny_counts, "M01")
    # Should only have S1,S2,S3
    assert set(md["#SampleID"]) == {"S1","S2","S3"}
    assert set(biom.columns) == {"S1","S2","S3"}
    # Order in biom columns should match metadata order as returned
    assert list(biom.columns) == list(md["#SampleID"])

def test_filter_by_prevalence_after_clr_behavior():
    # Make a tiny CLR table where OTU3 is often at column-min (treated as absent)
    clr = pd.DataFrame(
        [[ 0.3, 0.1, 0.2,  0.2],
         [ 0.0, 0.0, 0.0,  0.0],
         [-0.3,-0.1,-0.2, -0.2]],
        index=["OTU1","OTU2","OTU3"],
        columns=["S1","S2","S3","S4"]
    )
    # In this heuristic, OTU3 equals the min in all columns → prevalence ~ 0
    filtered = U.filter_by_prevalence(clr, threshold=0.75)
    assert "OTU3" not in filtered.index
    assert "OTU1" in filtered.index