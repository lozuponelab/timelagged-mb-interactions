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
3) Apply feature filtering on *presence in proportions* (pre-CLR):
   3a) Prevalence filter (presence > 0 in ≥ X% of samples).
   3b) Abundance gate (reaches ≥ min_rel abundance in ≥ min_samples).
   3c) Variance gate (retain by MAD quantile or top-N by MAD).
4) Replace zeros multiplicatively, re-close, and apply CLR (skbio).
5) Average replicate samples within the same day *in CLR space* (== compositional mean).
6) (Optional) Interpolate missing days *in CLR space* for smooth plots (avoid for TLC stats).
7) Save daily CLR matrices + daily metadata.

Why this order?
---------------
- Filtering belongs on presence/abundance/variance *before* CLR so the CLR
  denominator (geometric mean) is computed on stable, meaningful features.
- Multiplicative zero-replacement + closure is more faithful than adding a flat pseudocount.
- Averaging in CLR is equivalent to the compositional (geometric) mean on proportions.

Outputs
-------
- <intermediate_dir>/<subject>_avg_biom_clr_<suffix>.csv : rows=OTUs, cols=epoch_day (daily, CLR)
- <intermediate_dir>/<subject>_avg_meta_<suffix>.csv     : one row per day, includes epoch_time
- <final_dir>/<subject>_interp_biom_clr_<suffix>.csv     : (optional) daily CLR with gaps filled

Suffix encodes key filter params (see _suffix_from_params).

Author: Laurie Lyon (commentary expanded) via ChatGPT-5 (vibe coding with review by me)
Last updated: 2025-10-06
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from biom import load_table
from scipy.interpolate import PchipInterpolator
from scipy.stats import median_abs_deviation
from skbio.stats.composition import multiplicative_replacement, closure, clr
#multiplicative_replacement is depricated in latest skbio, can replace with multi_replace from skbio.stats.composition

# Utilities you already have; we do NOT modify these.
from biom_utils_updated import biom2df, meta_biom_filter

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

# ---------- Feature Filters (Pre-CLR) ----------
# added 10-5-25 to deal with the issue of setting too high of a prevalenced threshold
# this was to reduce compute time but would only allow the most prevalent taxa to be kept for clr

def prevalence_filter(subj_prop: pd.DataFrame, min_prop_present: float) -> pd.DataFrame:
    """
    Keep OTUs present (>0 proportion) in ≥ min_prop_present fraction of samples.
    """
    presence = (subj_prop > 0).sum(axis=1) / subj_prop.shape[1]
    keep = presence >= min_prop_present
    return subj_prop.loc[keep]

def abundance_filter(subj_prop: pd.DataFrame, min_rel: float, min_samples: int) -> pd.DataFrame:
    """
    Keep OTUs that reach at least min_rel proportion in ≥ min_samples samples.
    This removes 'present but vanishing' taxa that destabilize CLR and TLC.
    """
    hits = (subj_prop >= min_rel).sum(axis=1)
    keep = hits >= int(min_samples)
    return subj_prop.loc[keep]

def variance_filter_mad(subj_prop: pd.DataFrame,
                        mad_quantile: float | None = 0.50,
                        top_n: int | None = None) -> pd.DataFrame:
    """
    Retain dynamic taxa by Median Absolute Deviation (MAD) across time.
    Use either a quantile threshold OR keep top-N by MAD.
    """
    # Compute MAD across samples (time) for each OTU in *proportion space*
    mads = subj_prop.apply(lambda s: median_abs_deviation(s, scale='normal'), axis=1)
    if top_n is not None:
        keep_idx = mads.sort_values(ascending=False).head(int(top_n)).index
        return subj_prop.loc[keep_idx]
    if mad_quantile is not None:
        thr = mads.quantile(mad_quantile)
        keep = mads >= thr
        return subj_prop.loc[keep]
    return subj_prop

#dry run function to see how many taxa would be kept under different filter parameters
def dry_run_filter_counts(subj_prop, prev_vals=(0.3, 0.4, 0.5),
                          min_rel_vals=(0.001, 0.002),
                          min_samples_vals=(5, 7, 10),
                          mad_quantiles=(0.5, 0.6),
                          top_ns=(120, 150)):
    import pandas as pd
    from scipy.stats import median_abs_deviation

    def prevalence_filter(df, thr):
        presence = (df > 0).sum(axis=1) / df.shape[1]
        return df.loc[presence >= thr]

    def abundance_filter(df, min_rel, min_samples):
        hits = (df >= min_rel).sum(axis=1)
        return df.loc[hits >= min_samples]

    def count_after_all(df, prev, min_rel, min_s, mad_q=None, top_n=None):
        A = prevalence_filter(df, prev)
        B = abundance_filter(A, min_rel, min_s)
        mads = B.apply(lambda s: median_abs_deviation(s, scale='normal'), axis=1)
        if top_n is not None:
            keep = mads.sort_values(ascending=False).head(top_n).index
            C = B.loc[keep]
        else:
            thr = mads.quantile(mad_q)
            C = B.loc[mads >= thr]
        return len(C)

    rows = []
    for p in prev_vals:
        for ar in min_rel_vals:
            for ms in min_samples_vals:
                # MAD quantile route
                for q in mad_quantiles:
                    rows.append({"mode":"mad_q", "prev":p, "min_rel":ar, "min_samples":ms,
                                 "mad_q":q, "top_n":None,
                                 "n_taxa":count_after_all(subj_prop, p, ar, ms, mad_q=q)})
                # Top-N route
                for N in top_ns:
                    rows.append({"mode":"topN", "prev":p, "min_rel":ar, "min_samples":ms,
                                 "mad_q":None, "top_n":N,
                                 "n_taxa":count_after_all(subj_prop, p, ar, ms, top_n=N)})
    return pd.DataFrame(rows).sort_values("n_taxa", ascending=False)

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

#added 10-6-25
def _suffix_from_params(prevalence_threshold: float,
                        min_rel: float,
                        min_samples: int,
                        mad_quantile: float | None,
                        top_n: int | None) -> str:
    """
    Build a concise suffix encoding filter parameters for file names.
    Examples:
      prev40_abd0p1x5_mad50         (40% prev, 0.1% in ≥5, MAD ≥ median)
      prev50_abd0p1x5_top150        (50% prev, 0.1% in ≥5, top 150 by MAD)
    """
    prev = f"prev{int(prevalence_threshold*100)}"
    abd  = f"abd{str(min_rel).replace('0.', '0p')}x{min_samples}"
    if top_n is not None:
        var = f"top{int(top_n)}"
    else:
        var = f"mad{int(mad_quantile*100)}"
    return f"{prev}_{abd}_{var}"


# changed 10-6-25 to add extra filtering steps to fine tune which samples go into TLC
def main(input_biom: str,
         metadata_csv: str,
         intermediate_dir: str,
         final_dir: str,
         subject_id: str,
         prevalence_threshold: float = 0.40,
         min_rel: float = 0.001,
         min_samples: int = 5,
         mad_quantile: float | None = 0.50,
         top_n_by_mad: int | None = None,
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
        Fraction of samples requiring >0 presence to keep an OTU (default 0.40).
    min_rel : float
        Minimum relative abundance gate (e.g., 0.001 = 0.1%) in at least min_samples.
    min_samples : int
        Number of samples in which min_rel must be reached.
    mad_quantile : float | None
        Retain OTUs with MAD ≥ this quantile (e.g., 0.50 = median). Ignored if top_n_by_mad is set.
    top_n_by_mad : int | None
        If provided, override mad_quantile and keep the top-N OTUs by MAD.
    do_interpolate : bool
        If True, fill missing days (for plots). Avoid for TLC inference.
    """
    # Ensure output directories exist
    Path(intermediate_dir).mkdir(parents=True, exist_ok=True)
    Path(final_dir).mkdir(parents=True, exist_ok=True)

    # --- Load inputs ---
    biom_obj = load_table(input_biom)
    meta_df  = pd.read_csv(metadata_csv)
    otu_df   = biom2df(biom_obj)  # OTUs x samples (counts)

    # --- Filter to subject; enforce chronological order ---
    subj_meta, subj_otu = meta_biom_filter(meta_df, otu_df, subject_id)
    subj_meta = subj_meta.copy()
    subj_meta["epoch_time"] = subj_meta["epoch_time"].astype(int)
    subj_meta = subj_meta.sort_values("epoch_time")
    subj_otu  = subj_otu.loc[:, subj_meta["#SampleID"]]  # reorder columns to sample order

    # --- Counts → proportions per sample (closure) ---
    # Work in proportion space for presence/abundance/variance filters.
    col_sums = subj_otu.sum(axis=0).replace(0, np.nan)
    subj_prop = subj_otu.div(col_sums, axis=1).fillna(0.0)
    _assert_nonempty(subj_prop, "subj_prop (proportions)")
    _assert_proportions(subj_prop)

    # --- Pre-CLR feature filtering (three gates) ---
    n0 = subj_prop.shape[0]

    # 3a) Prevalence
    prop_prev = prevalence_threshold
    subj_prop = prevalence_filter(subj_prop, min_prop_present=prop_prev)
    n_prev = subj_prop.shape[0]

    # 3b) Abundance gate
    subj_prop = abundance_filter(subj_prop, min_rel=min_rel, min_samples=min_samples)
    n_abd = subj_prop.shape[0]

    # 3c) Variance gate (MAD)
    subj_prop = variance_filter_mad(subj_prop,
                                    mad_quantile=None if top_n_by_mad is not None else mad_quantile,
                                    top_n=top_n_by_mad)
    n_var = subj_prop.shape[0]

    _assert_nonempty(subj_prop, "Filtered subj_prop (after prevalence/abundance/variance)")

    # --- Multiplicative zero-replacement + closure + CLR (skbio expects samples x features) ---
    replaced = multiplicative_replacement(subj_prop.T.values)  # samples x OTUs
    closed   = closure(replaced)
    clr_vals = clr(closed)                                     # samples x OTUs

    # Rebuild DataFrame as rows=OTUs, cols=samples (match your downstream expectations)
    subj_clr = pd.DataFrame(
        clr_vals,
        index=subj_prop.columns,   # sample IDs
        columns=subj_prop.index    # OTU IDs
    ).T

    _assert_no_nans(subj_clr, "subj_clr")
    _assert_clr_centered(subj_clr)

    # --- Average replicate samples per day in CLR space ---
    avg_meta, avg_biom_clr = avg_sample_by_day_clr(subj_meta, subj_clr)
    _assert_daily_alignment(avg_meta, avg_biom_clr)

    # --- Optional: interpolate missing days (CLR space) for PLOTTING ONLY ---
    if do_interpolate:
        before_cols = avg_biom_clr.columns.tolist()
        interp_biom_clr = interp_missing_day_clr(avg_biom_clr)
        added_n, added_preview = _report_interpolation(before_cols, interp_biom_clr.columns.tolist())
    else:
        interp_biom_clr = avg_biom_clr.copy()
        added_n, added_preview = 0, []

    # --- Save outputs ---
    suffix = _suffix_from_params(prevalence_threshold, min_rel, min_samples, mad_quantile, top_n_by_mad)

    avg_biom_clr.to_csv(f"{intermediate_dir}/{subject_id}_avg_biom_clr_{suffix}.csv", index=True)
    avg_meta.to_csv(f"{intermediate_dir}/{subject_id}_avg_meta_{suffix}.csv", index=False)

    interp_biom_clr.index.name = "#OTU ID"
    interp_biom_clr.to_csv(f"{final_dir}/{subject_id}_interp_biom_clr_{suffix}.csv", index=True)

    # --- Console summary ---
    print(f"[OK] Saved: {subject_id}_avg_biom_clr_{suffix}.csv | {subject_id}_avg_meta_{suffix}.csv "
          f"| {subject_id}_interp_biom_clr_{suffix}.csv")
    print(f"Start OTUs: {n0}")
    print(f"After prevalence (≥{int(prop_prev*100)}% present): {n_prev}")
    print(f"After abundance (≥{min_rel:.4f} in ≥{min_samples} samples): {n_abd}")
    if top_n_by_mad is not None:
        print(f"After variance (top {top_n_by_mad} by MAD): {n_var}")
    else:
        print(f"After variance (MAD ≥ {int(mad_quantile*100)}th pct): {n_var}")
    if do_interpolate:
        print(f"Interpolated days added: {added_n} (e.g., {added_preview})")
    #print("NOTE: Use avg_biom_clr for TLC; interp_biom_clr is for plots only.")

def _build_argparser():
    """CLI argument parser with clear help strings."""
    p = argparse.ArgumentParser(
        description="Compositional wrangling (per subject) for TLC: subject filter → closure → prevalence/abundance/MAD filters → zero-replace → closure → CLR → daily aggregation (optional interpolate for plots)."
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

    # --- Filtering params (tuned defaults for TLC robustness) ---
    p.add_argument("--prevalence", type=float, default=0.40,
                   help="Prevalence threshold on presence in proportions before CLR (default 0.40 = 40%%).")
    p.add_argument("--min_rel", type=float, default=0.001,
                   help="Abundance gate: minimum relative abundance per sample (default 0.001 = 0.1%%).")
    p.add_argument("--min_samples", type=int, default=5,
                   help="Samples required to meet --min_rel (default 5).")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--mad_quantile", type=float, default=0.50,
                       help="Keep OTUs with MAD ≥ this quantile (default 0.50 = median). Ignored if --top_n_by_mad is set.")
    group.add_argument("--top_n_by_mad", type=int, default=None,
                       help="Override quantile: keep top-N taxa by MAD (e.g., 150).")

    # --- Interpolation flag ---
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
         min_rel=args.min_rel,
         min_samples=args.min_samples,
         mad_quantile=None if args.top_n_by_mad is not None else args.mad_quantile,
         top_n_by_mad=args.top_n_by_mad,
         do_interpolate=args.interpolate)