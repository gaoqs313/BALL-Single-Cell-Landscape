library(infercnv)

## Load rds
count_matrix <- readRDS("raw_count.matrix.rds")

## Create the InferCNV Object
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=count_matrix,
  annotations_file="annotations.tsv",
  delim="\t",
  gene_order_file="genes_ordered.tsv",
  chr_exclude = c("chrY", "chrM"),
  ref_group_names=c("B_lineage", "DC", "erythroid", "MK", "monocyte", "T_NK"))

## run inferCNV
infercnv_obj = infercnv::run(
    infercnv_obj,
    cutoff=0.1,
    out_dir='./subclusters/',
    plot_steps=F,
    HMM=TRUE,
    analysis_mode = "subclusters",
    num_threads = 16,
    k_obs_groups = 8,
    cluster_by_groups = F,
    denoise=TRUE,
    output_format="pdf"
)

