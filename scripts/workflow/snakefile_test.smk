## Global variables
# environments 
QIIME_CONDA = "qiime2-amplicon-2024.5"
MIBI_CONDA = "mibi"

# wildcards
SUBJECTS = ["F01", "M01", "M02"]

## Workflow
rule all:
    input:
        #telling snakemake the order of outputs that you want
        "results/mafft-tree/aligned_rep_seqs.qza",
        "results/mafft-tree/masked_aligned_rep_seqs.qza",
        "results/mafft-tree/unrooted_tree.qza",
        "results/mafft-tree/rooted_tree.qza",
        "results/sepp-tree/tree.qza",
        "results/sepp-tree/placements.qza", 
        "results/tree/tree.nwk"
        "results/qiime/taxonomy.qza", 
        expand("results/lagged_regression/{subject}.done", subject=SUBJECTS)

#step 1: data wrangling
rule wrangle_by_subject: 
    input: 
        biom = "data/dada2_otu_rare7500_lml.biom",
        metadata = "data/metadata_full_w_times.csv",
    output: 
        tsv = "results/processed_otus/{subject}_combined_otu.tsv"
    conda:
        MIBI_CONDA
    params: 
        n_threads = 8
    shell: 
    """
    python3 workflow/mibi/laurie_additions/time_series_data_wrangling_w_parser.py \
    # -i {input.biom} \
    # -m {input.metadata} \
    # -n results/interpolated_bioms \
    # -o results/combined_meta_and_otu_outputs \
    # -l {wildcards.subject}
    """


#step whatever
rule generate_tree:
    input:
        rep_seq = "data/cleaned_rep_seqs.qza" #referencing from directory where you run snakefile (run from mibi_longitude in this case)
    output:
        alignment = "data/qiime/fasttree/aligned_rep_seqs.qza",
        masked_alignment = "data/qiime/fasttree/masked_aligned_rep_seqs.qza",
        tree = "data/qiime/fasttree/unrooted_tree.qza",
        rooted_tree = "data/qiime/fasttree/rooted_tree.qza"
    conda:
        QIIME_CONDA
    params:
        n_threads = 8
    shell:
        """
        qiime phylogeny align-to-tree-mafft-fasttree \
            --i-sequences {input.rep_seq} \
            --o-alignment {output.alignment} \
            --o-masked-alignment {output.masked_alignment} \
            --o-tree {output.tree} \
            --o-rooted-tree {output.rooted_tree} \
            --p-n-threads {params.n_threads}
        """

#trying to generate tree by another method
rule sepp_phylo_tree:
    input:
        rep_seq = "data/cleaned_rep_seqs.qza",
        silva_ref = "data/sepp-refs-silva-128.qza"
    output:
        tree = "data/sepp-tree/tree.qza",
        placements = "data/sepp-tree/placements.qza"
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

rule export_tree: 
    input: 
        tree = "data/tree/rooted_tree.qza"
    output: 
        tree_nwk = "data/tree/tree.nwk"
    shell: 
        """
        qiime tools export \
            --input-path {input.tree} \
            --output-path {ouput.tree_nwk}
        """

# SILVA classifier 
rule taxonomy_classifier: 
    input: 
        classifier = "data/silva-138-99-nb-classifier.qza"
        seqs = "data/rep_seqs.qza"
    output: 
        taxonomy = "data/qiime/taxonomy.qza"
    conda: 
        QIIME_CONDA
    shell: 
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier}
            --i-reads {input.seqs}
            --o-classification {output.taxomony}
        """

# convert merged tsv back to biom for each subject
rule tsv_to_biom:

# # convert merged tsv back to biom for each subject
# rule tsv_to_biom:
#     input:
#         tsv = "results/combined_meta_and_otu_outputs/{subject}_combined_otu_clr_90pct.tsv"
#     output:
#         biom= "results/subject_biom/{subject}_biom.biom"
#     conda:
#         MIBI_CONDA
#     shell: 
#         """
#         mkdir -p results/subject_biom
#         python3 workflow/mibi/laurie_additions/convert_tsv_to_biom.py \
#             --input_tsv {input.tsv} \
#             --output_biom {output.biom}
#         """

# run lagged regression with MiBiTimeSeries for each subject
rule lagged_regression: 