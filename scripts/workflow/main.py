import argparse
import pandas as pd
from functools import partial
from multiprocessing import Pool, Queue
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

# import your updated class
from mibi.timeseries_0625_update import MiBiTimeSeries

def _process_pair(pair, ts, precision):
    i, j = pair
    # 4) get shuffled‐null summary (includes spear_rho_real, mean & pseudo_pval)
    #stats_shuf = ts.get_summary_stats_shuffled(i, j, precision=precision)
    # 5) get iaaft‐null summary (includes mean & pseudo_pval)
    stats_iaaft = ts.get_summary_stats_iaaft(i, j, precision=precision)

    # 6) assemble the row
    return {
        "j":                      j,
        "i":                      i,
        "spear_rho_real":         stats_iaaft["spear_rho_real"],
        #"mean_rho_shuffled":      stats_shuf["mean_rho_shuffled"],
        #"pseudo_pval_shuffled":   stats_shuf["pseudo_pval_shuffled"],
        "mean_rho_iaaft":         stats_iaaft["mean_rho_iaaft"],
        "pseudo_pval_iaaft":      stats_iaaft["pseudo_pval_iaaft"],
    }

def main():
    p = argparse.ArgumentParser(
        description="Compute observed and null‐model stats for every OTU pair"
    )
    p.add_argument("--biom_file",   required=True,
                   help="CSV of OTU table (rows=OTUs, cols=epoch_time)")
    p.add_argument("--meta_file",   required=True,
                   help="CSV of metadata with epoch_time column")
    p.add_argument("--precision",   type=int, default=1000,
                   help="Number of permutations/surrogates per pair")
    p.add_argument("--output_csv",  required=True,
                   help="Where to write the combined results CSV")
    args = p.parse_args()

    # 1) load inputs
    biom_df = pd.read_csv(args.biom_file, index_col=0)
    meta_df = pd.read_csv(args.meta_file)


    # init your timeseries object 
    ts = MiBiTimeSeries(meta_df, biom_df)

    # loop over all OTU‐pairs
    otus = list(ts.biom.index)

    pairs = [(i, j) for i in otus for j in otus if i != j]
    process_pair = partial(_process_pair, ts=ts, precision=args.precision)

    # 2) use multiprocessing to process pairs
    with Pool(10) as pool:
        rows = list(tqdm(pool.imap(process_pair, pairs), desc="Processing OTU pairs", total=len(pairs)))


    # for i, j in tqdm(pairs, desc="Processing OTU pairs", total=len(pairs)):

    #     # 4) get shuffled‐null summary (includes spear_rho_real, mean & pseudo_pval)
    #     stats_shuf  = ts.get_summary_stats_shuffled(
    #         i, j, precision=args.precision
    #         )
    #     # 5) get iaaft‐null summary (includes mean & pseudo_pval)
    #     stats_iaaft = ts.get_summary_stats_iaaft(
    #         i, j, precision=args.precision
    #      )

    #     # 6) assemble the row
    #     rows.append({
    #         "j":                      j,
    #         "i":                      i,
    #         "spear_rho_real":         stats_shuf["spear_rho_real"],
    #         "mean_rho_shuffled":      stats_shuf["mean_rho_shuffled"],
    #         "pseudo_pval_shuffled":   stats_shuf["pseudo_pval_shuffled"],
    #         "mean_rho_iaaft":         stats_iaaft["mean_rho_iaaft"],
    #         "pseudo_pval_iaaft":      stats_iaaft["pseudo_pval_iaaft"],
    #      })

    # 7) write out
    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(out_df)} pairwise results to {args.output_csv}")

    #pcols = ["pseudo_pval_shuffled", "pseudo_pval_iaaft"]
    #commending out shuffled null model for now
    pcols = ["pseudo_pval_iaaft"]
    alpha = 0.05

    for col in pcols:
        # multipletests returns: reject_array, pvals_corrected, _, _
        reject, pvals_corrected, _, _ = multipletests(
            out_df[col].values,
            alpha=alpha,
            method="fdr_bh"
        )
        out_df[f"{col}_fdr"]    = pvals_corrected
        out_df[f"{col}_signif"] = reject

    # 10) write out
    out_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(out_df)} rows (with FDR columns) to {args.output_csv}")

if __name__ == "__main__":
    main()