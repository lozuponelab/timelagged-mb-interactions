"""
SILVA Taxonomy Mapping Script
-----------------------------

Purpose:
This script processes subject-specific OTU tables and maps OTUs to taxonomy using QIIME2 and the SILVA database.
The workflow includes converting OTU tables to QIIME2-compatible formats, performing taxonomy classification,
and exporting results for downstream analysis.

Workflow:
1. Convert OTU tables (TSV) to BIOM format.
2. Import BIOM files into QIIME2 as FeatureTable artifacts.
3. Classify OTUs using the SILVA Naive Bayes taxonomy classifier.
4. Export taxonomy classifications to TSV format.

Inputs:
- OTU tables in TSV format, with sample IDs as headers and OTU IDs as rows.
- Representative sequences file (QIIME2 `.qza` artifact).
- SILVA taxonomy classifier (`.qza` artifact).

Outputs:
For each OTU table, the script generates:
- `<subject_id>_combined_otu.biom`: BIOM-formatted OTU table.
- `<subject_id>_feature_table.qza`: QIIME2 FeatureTable artifact.
- `<subject_id>_taxonomy.qza`: QIIME2 Taxonomy artifact.
- `<subject_id>_taxonomy.tsv`: Exported taxonomy in TSV format.

Usage:
    python silva_taxonomy_mapping.py \
        -i <otu_table_dir> \
        -r <rep_seq_path> \
        -c <classifier_path> \
        -o <output_dir>

Arguments:
- `-i, --otu_table_dir`: Directory containing subject-specific OTU tables (e.g., *_combined_otu_qiime2.tsv).
- `-r, --rep_seq_path`: Path to the QIIME2 artifact for representative sequences (e.g., .qza file).
- `-c, --classifier_path`: Path to the SILVA classifier artifact (e.g., silva-138-99-nb-classifier.qza).
- `-o, --output_dir`: Directory to store QIIME2 outputs (feature tables, taxonomy, TSV exports).

Example:
    python silva_taxonomy_mapping.py \
        -i ../data/combined_meta_and_otu_outputs/ \
        -r ../data/uclust_casey_rep_set.qza \
        -c ../data/silva-138-99-nb-classifier.qza \
        -o ../data/qiime_outputs/

Notes:
- Ensure the QIIME2 environment is activated before running the script.
- Input files must align with the requirements for QIIME2 taxonomy classification.

Author: Laurie Lyon
Date: 01/21/2025

"""

# Import required libraries
from tqdm import tqdm  # For creating progress bars
import os  # For file and directory operations
import subprocess  # For running shell commands
from pathlib import Path  # For handling directory paths
import argparse  # For command-line argument parsing

def main(final_output_dir, rep_seq_path, classifier_path, qiime_output_dir):
    """
    Process OTU tables and map OTUs to taxonomy using QIIME2 and the SILVA database.

    Parameters:
        - final_output_dir (str): Directory containing subject-specific OTU tables.
        - rep_seq_path (str): Path to the QIIME2 artifact for representative sequences (e.g., .qza file).
        - classifier_path (str): Path to the SILVA classifier (e.g., silva-138-99-nb-classifier.qza).
        - qiime_output_dir (str): Directory to store output files (e.g., feature tables, taxonomy, and TSV exports).
    """
    # Ensure the output directory exists
    Path(qiime_output_dir).mkdir(parents=True, exist_ok=True)

    #Check if the rep seq file is .fna and convert to .qza if needed 
    if rep_seq_path.endswith(".fna"):
        print("Converting representative sequences to QIIME2 artifact...")
        req_seq_qza = rep_seq_path.replace(".fna", ".qza")
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "FeatureData[Sequence]",
            "--input-path", rep_seq_path,
            "--output-path", req_seq_qza
        ], check=True)
        rep_seq_path = req_seq_qza
    
    # Retrieve a list of OTU files to process
    otu_files = [f for f in os.listdir(final_output_dir) if f.endswith("_combined_otu_qiime2.tsv")]

    # Process each OTU table using a progress bar
    for otu_file in tqdm(otu_files, desc="Processing OTU Tables", unit="file"):
        # Extract subject ID (e.g., "F01" from "F01_combined_otu_qiime2.tsv")
        subject_id = otu_file.split("_")[0]
        print(f"\nProcessing SILVA taxonomy mapping for subject {subject_id}...")

        # Define paths for intermediate and output files
        input_tsv = os.path.join(final_output_dir, otu_file)  # Input OTU table (TSV format)
        biom_file = os.path.join(qiime_output_dir, f"{subject_id}_combined_otu.biom")  # Intermediate BIOM file
        feature_table = os.path.join(qiime_output_dir, f"{subject_id}_feature_table.qza")  # QIIME2 FeatureTable artifact
        taxonomy_output = os.path.join(qiime_output_dir, f"{subject_id}_taxonomy.qza")  # QIIME2 taxonomy artifact
        taxonomy_export = os.path.join(qiime_output_dir, f"{subject_id}_taxonomy.tsv")  # Exported taxonomy TSV file

        # Step 1: Convert TSV to BIOM format
        print(f"Converting {input_tsv} to BIOM format...")
        subprocess.run([
            "biom", "convert",
            "-i", input_tsv,               # Input TSV file
            "-o", biom_file,               # Output BIOM file
            "--table-type", "OTU table",   # Specify table type
            "--to-hdf5"                    # Use HDF5 format for BIOM
        ], check=True)

        # Step 2: Import BIOM file as a QIIME2 FeatureTable artifact
        print(f"Importing BIOM file {biom_file} into QIIME2...")
        subprocess.run([
            "qiime", "tools", "import",
            "--type", "FeatureTable[Frequency]",  # Specify QIIME2 artifact type
            "--input-path", biom_file,            # Path to the BIOM file
            "--output-path", feature_table        # Output QIIME2 FeatureTable artifact
        ], check=True)

        # Step 3: Classify OTUs using the SILVA taxonomy classifier
        print(f"Classifying OTUs for subject {subject_id} using SILVA...")
        subprocess.run([
            "qiime", "feature-classifier", "classify-sklearn",
            "--i-classifier", classifier_path,    # SILVA classifier artifact
            "--i-reads", rep_seq_path,            # Representative sequences artifact
            "--o-classification", taxonomy_output  # Output taxonomy classification artifact
        ], check=True)

        # Step 4: Export the taxonomy classification to a TSV file
        print(f"Exporting taxonomy for subject {subject_id}...")
        subprocess.run([
            "qiime", "tools", "export",
            "--input-path", taxonomy_output,  # QIIME2 taxonomy artifact
            "--output-path", qiime_output_dir # Export directory
        ], check=True)

        print(f"Taxonomy mapping for subject {subject_id} saved to {taxonomy_export}.")

# Command-line interface for the script
if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Map OTUs to SILVA taxonomy using QIIME2.")
    
    # Input directory containing subject-specific OTU tables
    parser.add_argument(
        "-i", "--final_output_dir", required=True, 
        help="Directory containing subject-specific OTU tables (e.g., *_combined_otu_qiime2.tsv)."
    )
    
    # Path to the representative sequences QIIME2 artifact
    parser.add_argument(
        "-r", "--rep_seq_path", required=True, 
        help="Path to the QIIME2 artifact for representative sequences (e.g., .qza file)."
    )
    
    # Path to the SILVA taxonomy classifier
    parser.add_argument(
        "-c", "--classifier_path", required=True, 
        help="Path to the SILVA classifier artifact (e.g., silva-138-99-nb-classifier.qza)."
    )
    
    # Output directory for QIIME2 artifacts and exports
    parser.add_argument(
        "-o", "--qiime_output_dir", required=True, 
        help="Directory to store QIIME2 outputs (e.g., feature tables, taxonomy, and TSV exports)."
    )

    # parser.add_argument(
    #    "-n", "--n_jobs", required=False, default=4, type=int,
    #    help="Number of jobs to run in parallel for taxonomy classification."
    # )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(
        final_output_dir=args.final_output_dir,
        rep_seq_path=args.rep_seq_path,
        classifier_path=args.classifier_path,
        qiime_output_dir=args.qiime_output_dir,
        # n_jobs=args.n_jobs
    )