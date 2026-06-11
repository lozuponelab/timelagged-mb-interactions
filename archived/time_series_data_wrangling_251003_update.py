#!/usr/bin/env python3
"""
time_series_data_wrangling_clr_filtered_per_subject.py

Purpose
-------
Prepare a single subject's 16S OTU table and metadata for
time-lagged correlation (TLC) and plotting, using principled
compositional preprocessing:

1) Filter to the subject of interest (keeping sample order from metadata).
2) Convert counts → proportions (closure).
3) Apply a prevalence filter on *presence in proportions* (pre-CLR).
4) Replace zeros multiplicatively, re-close, and apply CLR (skbio).
5) Average replicate samples within the same day *in CLR space* (== compositional mean).
6) (Optional) Interpolate missing days *in CLR space* for smooth plots (avoid for TLC stats).
7) Save daily CLR matrices + daily metadata.

Why this order?
---------------
- Prevalence belongs on presence (proportions > 0), not on CLR values.
- Multiplicative zero-replacement + closure is more faithful than adding a flat pseudocount.
- Averaging in CLR is equivalent to the compositional (geometric) mean on proportions.

Outputs
-------
- <intermediate_dir>/<subject>_avg_biom_clr_<thr>.csv : rows=OTUs, cols=epoch_day (daily, CLR)
- <intermediate_dir>/<subject>_avg_meta_<thr>.csv     : one row per day, includes epoch_time
- <final_dir>/<subject>_interp_biom_clr_<thr>.csv     : (optional) daily CLR with gaps filled

Author: Laurie Lyon (commentary expanded) via ChatGPT-5 (vibe coding)
Last updated: 2025-10-03
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from biom import load_table
from scipy.interpolate import PchipInterpolator
from skbio.stats.composition import multiplicative_replacement, closure, clr
#multiplicative_replacement is depricated in latest skbio, can replace with multi_replace from skbio.stats.composition

# Utilities you already have; we do NOT modify these.
from .biom_utils_updated import biom2df, meta_biom_filter

# ---------- Assertions & Sanity Checks ----------

def _assert_nonempty(df: pd.DataFrame, name: str):
    if df.shape[0] == 0 or df.shape[1] == 0:
        raise AssertionError(f"{name} is empty (shape={df.shape}). Upstream filtering may be too strict.")

def _assert_prevalence_range(thresh: float):
    if not (0.0 < thresh <= 1.0):
        raise AssertionError(f"prevalence_threshold must be in (0,1], got {thresh}")

def _assert_columns_align(df: pd.DataFrame, cols: list, name: str):
    if list(df.columns) != list(cols):
        raise AssertionError(f"{name} columns do not align with expected order from metadata.\n"
                             f"First few expected={list(cols)[:5]}\nFirst few actual={list(df.columns)[:5]}")

def _assert_increasing(times: pd.Series, name: str):
    if not (times.values[:-1] < times.values[1:]).all():
        raise AssertionError(f"{name} must be strictly increasing. Found non-increasing sequence.")

def _assert_proportions(df_prop: pd.DataFrame, tol: float = 1e-6):
    col_sums = df_prop.sum(axis=0)
    if not np.allclose(col_sums.values, 1.0, atol=tol):
        bad = np.where(~np.isclose(col_sums.values, 1.0, atol=tol))[0][:5]
        raise AssertionError(f"Proportions columns must sum to 1 (±{tol}). Offenders (first 5 idx): {bad}, "
                             f"sums={col_sums.iloc[bad].tolist()}")

def _assert_no_nans(df: pd.DataFrame, name: str):
    if df.isna().any().any():
        where = np.argwhere(df.isna().values)
        i, j = where[0]
        raise AssertionError(f"{name} contains NaNs. Example at (row={df.index[i]}, col={df.columns[j]}).")

def _assert_clr_centered(df_clr: pd.DataFrame, tol: float = 1e-8):
    # CLR per sample should have mean 0; check across columns (samples)
    col_means = df_clr.mean(axis=0)
    if not np.allclose(col_means.values, 0.0, atol=tol):
        bad = np.where(~np.isclose(col_means.values, 0.0, atol=tol))[0][:5]
        raise AssertionError(f"CLR columns must be centered (~0 mean). Offenders (first 5 idx): {bad}, "
                             f"means={col_means.iloc[bad].tolist()}")

def _assert_daily_alignment(avg_meta: pd.DataFrame, avg_biom_clr: pd.DataFrame, time_col: str = "epoch_time"):
    # avg_biom_clr has columns = epoch days; avg_meta has a row per day with epoch_time
    times = avg_meta[time_col].astype(int).tolist()
    cols  = [int(c) for c in avg_biom_clr.columns]
    if times != cols:
        raise AssertionError("Daily CLR columns do not match avg_meta epoch_time order.\n"
                             f"First few meta={times[:5]}\nFirst few clr cols={cols[:5]}")

def _report_interpolation(before_cols: list, after_cols: list):
    before_set = set(before_cols)
    after_set  = set(after_cols)
    added = sorted(list(after_set - before_set))
    return len(added), added[:5]

# ---------- Core Functions ----------
def avg_sample_by_day_clr(subject_metadata: pd.DataFrame,
                          subj_clr: pd.DataFrame,
                          sample_col: str = "#SampleID",
                          time_col: str = "epoch_time"):
    """
    Average replicate samples collected on the same day in *CLR space*.

    Why in CLR?
    -----------
    Arithmetic mean in CLR equals the compositional (geometric) mean
    in proportions, then closed back to the simplex. So this is the
    compositionally correct way to collapse replicates.

    Parameters
    ----------
    subject_metadata : DataFrame
        Metadata filtered to the subject, one row per *sample*.
        Must include sample_col and integer epoch_time (seconds).
    subj_clr : DataFrame
        CLR-transformed OTU table, rows=OTUs, cols=samples.
    sample_col : str
        Column name for sample IDs.
    time_col : str
        Column name for epoch timestamps.

    Returns
    -------
    avg_meta : DataFrame
        One metadata row per day (first sample's metadata kept).
    avg_biom_clr : DataFrame
        Daily CLR matrix (rows=OTUs, cols=epoch_day).
    """
    meta = subject_metadata.copy()
    # Ensure epoch_time is int (seconds since epoch) for consistent sorting/indexing
    meta[time_col] = meta[time_col].astype(int)

    # Collect unique days we actually have samples for
    days = sorted(meta[time_col].unique())

    # We'll create a new matrix with one column per day
    out_cols = []
    out_vals = []

    for d in days:
        # All samples collected on day d (possibly multiple samples/day)
        cols = meta.loc[meta[time_col] == d, sample_col].tolist()

        if len(cols) == 1:
            # Single sample that day → use its CLR vector
            out_vals.append(subj_clr[cols[0]].values)
        else:
            # Multiple samples that day → arithmetic mean in CLR space
            # (== compositional mean in proportions)
            out_vals.append(subj_clr[cols].mean(axis=1).values)

        out_cols.append(d)

    # Build OTUs x days DataFrame from the collected columns
    avg_biom_clr = pd.DataFrame(
        np.column_stack(out_vals),
        index=subj_clr.index,  # OTUs
        columns=out_cols       # epoch days (int)
    )

    # One metadata row per day. We keep the first entry's info for that day.
    avg_meta = meta.groupby(time_col).first().reset_index()

    return avg_meta, avg_biom_clr


def interp_missing_day_clr(ts_clr_df: pd.DataFrame) -> pd.DataFrame:
    """
    Interpolate missing days in *CLR space* using PCHIP (shape-preserving).

    Use-case
    --------
    - Great for *figures* (smooth daily curves).
    - Avoid using interpolated values for *statistical inference* (e.g., TLC),
      or at least mark/omit imputed deltas in the correlation.

    Parameters
    ----------
    ts_clr_df : DataFrame
        rows=OTUs, cols=epoch days (int seconds). Daily but with gaps.

    Returns
    -------
    out : DataFrame
        rows=OTUs, cols=epoch days (int seconds), with daily coverage from min..max.
    """
    df = ts_clr_df.copy()

    # Ensure day columns are integers (epoch seconds)
    if not np.issubdtype(df.columns.dtype, np.integer):
        df.columns = df.columns.astype(int)

    # Convert integer epoch seconds to pandas datetime for range construction
    dt_cols = pd.to_datetime(df.columns.astype(np.int64), unit='s')

    # Daily frequency between min and max observed day
    full = pd.date_range(start=dt_cols.min(), end=dt_cols.max(), freq='D')

    # Prepare output table of same rows (OTUs) and full daily columns
    out = pd.DataFrame(index=df.index, columns=full, dtype=float)

    for otu in df.index:
        # Build a time series (datetime index) for a single OTU
        s = pd.Series(df.loc[otu].values, index=dt_cols).sort_index()

        # Align to full date range (introduces NaNs for missing days)
        s_full = s.reindex(full)
        valid = s_full.notna()

        if valid.sum() > 1:
            # PCHIP requires at least 2 points
            xk = s_full.index[valid].view(np.int64).astype(float)  # known x (ns since epoch)
            yk = s_full[valid].values                               # known y (CLR)
            xf = s_full.index.view(np.int64).astype(float)          # all x in range

            # Shape-preserving interpolation
            p = PchipInterpolator(xk, yk)
            out.loc[otu] = p(xf)
        else:
            # Not enough points → just copy what we have (NaNs remain otherwise)
            out.loc[otu] = s_full.values

    # Convert back to epoch seconds as integer columns (consistent with rest of pipeline)
    out.columns = (out.columns.view(np.int64) // 10**9).astype(int)

    return out


def main(input_biom: str,
         metadata_csv: str,
         intermediate_dir: str,
         final_dir: str,
         subject_id: str,
         prevalence_threshold: float = 0.90,
         do_interpolate: bool = False):
    """
    Orchestrates the full preprocessing for one subject.

    Parameters
    ----------
    input_biom : str
        Path to BIOM file (counts).
    metadata_csv : str
        Path to metadata CSV containing #SampleID, ANONYMIZED_NAME, epoch_time, etc.
    intermediate_dir : str
        Directory to save daily CLR and daily metadata.
    final_dir : str
        Directory to save (optional) interpolated daily CLR.
    subject_id : str
        The subject's anonymized name to filter on (e.g., 'M01').
    prevalence_threshold : float
        Proportion of days in which an OTU must be present (>0 proportion) to keep it.
        Typical: 0.75–0.90. Higher → fewer OTUs, more stable TLC.
    do_interpolate : bool
        If True, fill missing days (for plots). Avoid for TLC inference.
    """
    # Ensure output directories exist
    Path(intermediate_dir).mkdir(parents=True, exist_ok=True)
    Path(final_dir).mkdir(parents=True, exist_ok=True)

    # --- Load inputs ---
    biom_obj = load_table(input_biom)  # parses BIOM format into a BIOM Table
    meta_df  = pd.read_csv(metadata_csv)  # sample-wise metadata
    otu_df   = biom2df(biom_obj)          # convert BIOM to DataFrame (OTUs x samples, counts)

    # --- Filter to subject (keeps only samples for subject_id) ---
    # Note: original utility used set() which can scramble order; we will
    # enforce chronological order immediately after this call.
    subj_meta, subj_otu = meta_biom_filter(meta_df, otu_df, subject_id)

    # --- Enforce chronological order (critical for time series) ---
    # epoch_time must be integer seconds for sorting and for consistent columns downstream
    subj_meta = subj_meta.copy()
    subj_meta["epoch_time"] = subj_meta["epoch_time"].astype(int)

    # Sort rows (samples) by time, then reorder columns of the OTU table to match
    subj_meta = subj_meta.sort_values("epoch_time")
    subj_otu  = subj_otu.loc[:, subj_meta["#SampleID"]]  # align OTU columns to the sorted sample IDs

    # --- Counts → proportions per sample (closure) ---
    # We move from raw or rarefied counts to compositions (sum 1 per sample).
    subj_prop = subj_otu.div(subj_otu.sum(axis=0), axis=1)

    # --- Prevalence on presence BEFORE CLR ---
    # Presence = proportion > 0 (post-rarefaction zeros remain zeros).
    # We keep OTUs that are present in at least `prevalence_threshold` fraction of samples.
    presence = (subj_prop > 0).sum(axis=1) / subj_prop.shape[1]
    keep_mask = presence >= prevalence_threshold
    subj_prop = subj_prop.loc[keep_mask]

    # --- Multiplicative zero-replacement + closure + CLR (skbio) ---
    # Why? CLR requires strictly positive compositions. Multiplicative replacement
    # preserves ratios better than adding a flat pseudocount; closure re-normalizes.
    # skbio expects "samples x features" → transpose, apply, then transpose back.
    replaced = multiplicative_replacement(subj_prop.T.values)  # shape: samples x OTUs
    closed   = closure(replaced)                                # safety: ensure rows sum to 1
    clr_vals = clr(closed)                                      # CLR per sample (still samples x OTUs)

    # Build a DataFrame in your expected orientation: rows=OTUs, cols=samples
    subj_clr = pd.DataFrame(
        clr_vals,
        index=subj_prop.columns,   # sample IDs
        columns=subj_prop.index    # OTU IDs
    ).T

    # --- Average replicate samples per day in CLR space ---
    # This collapses multiple samples collected on the same epoch day to one column.
    avg_meta, avg_biom_clr = avg_sample_by_day_clr(subj_meta, subj_clr)

    # --- Optional: interpolate missing days (CLR space) for PLOTTING ONLY ---
    # TLC inference should either (a) use only observed days, or (b) mask/omit
    # imputed deltas to avoid inflating correlation.
    if do_interpolate:
        interp_biom_clr = interp_missing_day_clr(avg_biom_clr)
    else:
        interp_biom_clr = avg_biom_clr.copy()

    # --- Save outputs ---
    # Suffix indicates prevalence threshold (e.g., 90 => 0.90)
    suffix = f"{int(prevalence_threshold * 100):02d}"

    # Daily CLR (no interpolation) → good input to TLC (columns = epoch days)
    avg_biom_clr.to_csv(f"{intermediate_dir}/{subject_id}_avg_biom_clr_{suffix}.csv", index=True)

    # One row per day with the original metadata fields (epoch_time retained)
    avg_meta.to_csv(f"{intermediate_dir}/{subject_id}_avg_meta_{suffix}.csv", index=False)

    # Interpolated daily CLR (for figures); rows=OTUs, cols=epoch days
    interp_biom_clr.index.name = "#OTU ID"
    interp_biom_clr.to_csv(f"{final_dir}/{subject_id}_interp_biom_clr_{suffix}.csv", index=True)

    # --- Console summary ---
    print(f"[OK] Saved: {subject_id}_avg_biom_clr_{suffix}.csv | {subject_id}_avg_meta_{suffix}.csv "
          f"| {subject_id}_interp_biom_clr_{suffix}.csv")
    print(f"OTUs kept after prevalence {prevalence_threshold:.2f}: "
          f"{keep_mask.sum()} of {len(keep_mask)} ({keep_mask.mean()*100:.1f}%)")
    print("NOTE: be aware when using avg_biom_clr for TLC; use interp_biom_clr for plots only.")

def _build_argparser():
    """CLI argument parser with clear help strings."""
    p = argparse.ArgumentParser(
        description="Compositional wrangling (per subject) for TLC: prevalence→zero-replacement→closure→CLR→daily aggregation."
    )
    p.add_argument("-i", "--input_biom",       required=True,
                   help="Path to BIOM table (counts).")
    p.add_argument("-m", "--metadata_csv",     required=True,
                   help="CSV with metadata, including #SampleID, ANONYMIZED_NAME, epoch_time (seconds).")
    p.add_argument("-n", "--intermediate_dir", required=True,
                   help="Directory to write daily CLR + daily metadata CSVs.")
    p.add_argument("-o", "--final_dir",        required=True,
                   help="Directory to write interpolated daily CLR (for plots).")
    p.add_argument("-s", "--subject_id",       required=True,
                   help="Subject ID to filter on (e.g., 'M01').")
    p.add_argument("--prevalence", type=float, default=0.90,
                   help="Prevalence threshold on presence in proportions before CLR (default 0.90).")
    p.add_argument("--interpolate", action="store_true",
                   help="If set, fill missing days in CLR space (for plots only; avoid for TLC stats).")
    return p

if __name__ == "__main__":
    parser = _build_argparser()
    args = parser.parse_args()

    main(input_biom=args.input_biom,
         metadata_csv=args.metadata_csv,
         intermediate_dir=args.intermediate_dir,
         final_dir=args.final_dir,
         subject_id=args.subject_id,
         prevalence_threshold=args.prevalence,
         do_interpolate=args.interpolate)