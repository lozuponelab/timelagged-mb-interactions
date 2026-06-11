from skbio import TreeNode
import pandas as pd
import os

#PSEUDOCODE
# Load in all interp bioms from each subject 
# make combined list of OTUs that are present in any subject's biom
# load that as the df to go into the pruning 

def find_all_otus(biom_path, subject_ids):
    """
    Finds shared indices (OTU IDs) across multiple subject-specific BIOM files
    
    Args:
        biom_path (str): Path to the directory containing subject-specific BIOM files.
            should be "interpolated_bioms" if using the output of the time series data wrangling script.
        subject_ids (list): List of subject IDs to look for in the BIOM files.
        
    Returns:
        biom_df (pd.DataFrame): DataFrame containing the shared OTUs across all subjects.
        """
    """"""

    all_otus = []

    for subject in subject_ids:
        subject_file = os.path.join(biom_path, f"{subject}_interp_biom_50.csv") #changed from 90 to 50 to match the new prevalence filter
        if os.path.exists(subject_file):
            biom_df = pd.read_csv(subject_file, index_col=0)
            otus = [str(x).strip().strip('"').strip("'") for x in biom_df.columns.tolist()]
            all_otus.extend(otus)
        else:
            print(f"Warning: {subject_file} does not exist.")
    
    return list(set(all_otus))

# call find_shared_otus to find shared OTUs across subjects
#all_subj_otus  = find_all_otus("results/interpolated_bioms", ["F01", "M01", "M02"])
# 09-24-25 changing this line to pull the interpolated bioms with the 50% prevalence filter to make a larger tree
all_subj_otus  = find_all_otus("results/250924_prev_wrangling", ["F01", "M01", "M02"])
print(all_subj_otus)

tree = TreeNode.read("data/qiime/sepp-tree/tree.nwk")
all_tips = [tip.name for tip in tree.tips()]
all_tips = [t.strip().strip('"').strip("'") for t in all_tips]  # normalize

matched = set(all_subj_otus).intersection(all_tips)
print(matched)

# Prune the tree to only those matched tips (TreeNode shear command)
pruned_tree = tree.shear(list(matched))

# (C) Optionally, write out the pruned tree to Newick
with open("results/allsubj_50pctprev_pruned_tree.nwk", "w") as out:
    pruned_tree.write(out, format="newick")

print("Original tip count:", len(all_tips))
print("Pruned tip count:  ", len([t.name for t in pruned_tree.tips()]))