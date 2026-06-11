
import pandas as pd
import numpy as np
from biom import load_table
from scipy.stats import gmean

## UTILITIES FOR TIME SERIES DATA WRANGLING

def biom2df(biom_table):
    """Convert a BIOM table into a pandas DataFrame (OTUs as rows, samples as columns)."""
    data_matrix = biom_table.matrix_data.toarray()
    df = pd.DataFrame(data_matrix,
                      index=biom_table.ids(axis="observation"),
                      columns=biom_table.ids(axis="sample"))
    return df

def clr_transform(df, pseudocount=1e-6):
    """Apply centered log-ratio (CLR) transformation."""
    log_df = np.log(df + pseudocount)
    # compute the geometric mean across rows (OTUs)
    geometric_mean = log_df.mean(axis=0)
    # subtract the geometric mean from each log-transformed value to center the data within each sample
    clr_df = log_df.subtract(geometric_mean, axis=1)
    return clr_df

def meta_biom_filter(metadata, biom_df, anon_name, sample_col="#SampleID", subject_col="ANONYMIZED_NAME"):
    """Filter metadata and BIOM table for a specific subject with ID alignment."""
    metadata[sample_col] = metadata[sample_col].astype(str).str.strip()
    biom_df.columns = biom_df.columns.astype(str).str.strip()

    metadata.index = metadata[sample_col]
    # filter metadata for the subject of interest (anon_name)
    subject_metadata = metadata[metadata[subject_col] == anon_name]

    # list of the set intersection of sample_col in metadata and columns in biom_df
    shared_indices = list(set(subject_metadata[sample_col]) & set(biom_df.columns))
    # subject-specific biom table and metadata with only the shared indices
    subject_biom = biom_df[shared_indices]
    subject_metadata = subject_metadata.loc[shared_indices]

    return subject_metadata, subject_biom

def filter_by_prevalence(df, threshold=0.9):
    """
    Filter OTUs that were originally non‐zero in at least `threshold` proportion
    of samples, *after* CLR transform.  We assume any OTU that was zero
    before CLR becomes the *minimum* value in that column, so we count
    an OTU as 'present' if its CLR value > column‐min.

    Parameters
    ----------
    df : pandas.DataFrame
        CLR‐transformed OTU table (OTUs × samples).
    threshold : float
        Proportion of samples in which an OTU must be present to keep it.

    Returns
    -------
    filtered_df : pandas.DataFrame
        Subset of rows (OTUs) passing the prevalence threshold.
    """
    #identify the min CLR value in each sample (column)
    col_mins = df.min(axis=0)
    #building a boolean mask for each OTU: True if any CLR value is greater than the column min
    #axis=1 tells pandas to broadcast the column min across all rows
    presence_mask = df.gt(col_mins, axis=1)
    #calculate the prevalence as the proportion of samples where the OTU is present
    prevalence = presence_mask.sum(axis=1) / df.shape[1]
    #return the filtered df of all rows that meet or exceed the threshold
    filtered_df = df.loc[prevalence >= threshold]
    return filtered_df

# not currently using these that Casey had, but leaving them in for now just in case

def norm_biom_table(biom_filepath, sep="\t", skiprows=1, index_col=0):
    """Normalize TSV-formatted biom table: rows=samples, cols=OTUs."""
    biom_df = pd.read_csv(biom_filepath, sep=sep, skiprows=skiprows, index_col=index_col)
    biom_df_wide = biom_df.T
    biom_df_wide.index.rename("#SampleID", inplace=True)
    return biom_df_wide

def wide_to_long(biom_df_wide):
    """Convert wide-format BIOM table to long format (for plotting)."""
    biom_df_long = biom_df_wide.reset_index().melt(id_vars="#SampleID", 
                                                   var_name="#OTU ID", 
                                                   value_name="abundance")
    return biom_df_long

def merge_biom_metadata(biom_df_long, meta_df, biom_col="#SampleID", meta_col="#SampleID"):
    """Merge long-format biom data with metadata."""
    biom_meta_df = pd.merge(biom_df_long, meta_df, how="left", left_on=biom_col, right_on=meta_col)
    return biom_meta_df
