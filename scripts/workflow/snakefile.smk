## Global variables
# environments 
QIIME_CONDA = "qiime2-amplicon-2024.5"
MIBI_CONDA = "mibi_madi"

# wildcards
SUBJECTS = ["F01", "M01", "M02"]

## Workflow
rule all:
    input:
        expand("results/avg_otu_meta_pre_interp/{subject}_avg_biom.csv", subject=SUBJECTS),
        expand("results/avg_otu_meta_pre_interp/{subject}_avg_meta.csv", subject=SUBJECTS),
        expand("results/interpolated_bioms/{subject}_interp_biom.csv", subject=SUBJECTS),
        "data/qiime/sepp-tree/tree.qza",
        "data/qiime/sepp-tree/placements.qza",
        "data/qiime/sepp-tree/tree.nwk",
        "data/qiime/taxonomy.qza", 
        "results/filtered-table.qza",
        "results/removed-table.qza"
        #expand("results/subject_biom/{subject}_biom.biom", subject=SUBJECTS)
        #expand("results/lagged_regression/{subject}_permutation_results.csv", subject=SUBJECTS)



#step 1: data wrangling
#input biom and metadata files to split by subject, average relative abundances
#for >1 sample/day, identify days with missing samples, and interpolate relative abundances
#also updates metadata with entries for interpolated timepoints 
#!!NEED TO UPDATE TO FIX SUBJECT FOR LOOP VS. SNAKEMAKE WILDCARDS!!
#Currently, running the UPDATED wrangling script works on its own but not when the wrapper is called via snakemake 
rule wrangle_by_subject:
    input:
        biom = "workflow/mibi/laurie_additions/dada2_otu_table_w_tax_no_pynast_failures_rare7500.biom",
        metadata = "workflow/mibi/laurie_additions/metadata_full_w_times.csv"
    output:
        avg_biom = "results/avg_otu_meta_pre_interp/{subject}_avg_biom.csv",
        avg_meta = "results/avg_otu_meta_pre_interp/{subject}_avg_meta.csv",
        interp_biom = "results/interpolated_bioms/{subject}_interp_biom.csv"
    conda:
        MIBI_CONDA
    params:
        n_threads = 8,
        intermediate_dir = "results/avg_otu_meta_pre_interp",
        final_dir = "results/interpolated_bioms"
    shell:
        """
        python3 workflow/mibi/laurie_additions/time_series_data_wrangling_clr_filtered_per_subject.py \
            --input_biom {input.biom} \
            --metadata_csv {input.metadata} \
            --intermediate_dir {params.intermediate_dir} \
            --final_dir {params.final_dir} \
            --subject_id {wildcards.subject}
        """

#step 2: generate tree file from rep_seqs
#needed for later lagged regression rule
rule sepp_phylo_tree:
    input:
        rep_seq = "data/cleaned_rep_seqs.qza",
        silva_ref = "data/sepp-refs-silva-128.qza"
    output:
        tree = "data/qiime/sepp-tree/tree.qza",
        placements = "data/qiime/sepp-tree/placements.qza"
    conda:
        QIIME_CONDA
    shell:
        """
        qiime fragment-insertion sepp \
            --i-representative-sequences {input.rep_seq} \
            --i-reference-database {input.silva_ref} \
            --o-tree {output.tree} \
            --o-placements {output.placements}
        """


#step 2.5 exporting tree from .qza to .nwk for downstream use
rule export_tree: 
    input: 
        tree_file = "data/qiime/sepp-tree/tree.qza"
    output: 
        tree_nwk = "data/qiime/sepp-tree/tree.nwk"
    conda:
        QIIME_CONDA
    params: 
        out_path = "data/qiime/sepp-tree"
    shell: 
        """
        qiime tools export \
            --input-path {input.tree_file} \
            --output-path {params.out_path}
        """




# step 3: SILVA classifier 
rule taxonomy_classifier: 
    input: 
        classifier = "data/silva-138-99-nb-classifier.qza",
        seqs = "data/rep_seqs.qza"
    output: 
        taxonomy = "data/qiime/taxonomy.qza"
    conda: 
        QIIME_CONDA
    shell: 
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.seqs} \
            --o-classification {output.taxonomy}
        """

#step 2.75 prune tree (trim_tree.py) once for all subjects

# run lagged regression with MiBiTimeSeries for each subject
# rule lagged_regression:
#     input:
#         biom = "results/subject_biom/{subject}_biom.biom",
#         meta = "results/combined_meta_and_otu_outputs/{subject}_updated_metadata.csv",
#         tree = "results/{subject}_pruned_tree.nwk"
#     output:
#         result_csv = "results/lagged_regression/{subject}_permutation_results.csv"
#     params: 
#         n_threads = 10
#     conda:
#         MIBI_CONDA
#     shell:
#         """
#         mkdir -p results/lagged_regression
#         python3 workflow/main.py \
#             --biom_file {input.biom} \
#             --meta_file {input.meta} \
#             --tree_file {input.tree} \
#             --subject {wildcards.subject} \
#             --precision 10 \
#             --output {output.result_csv}
#         """
#change precision to a higher number after testing