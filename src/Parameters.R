# title: "Parameters.R"
# author: "Pedro Sodré"
# date: "2025-06-04"

# --- Input/Output Paths ----
# Path to the input RDS file containing Seurat objects (e.g., a list of Seurat objects).
input_rds_path <- "/data04/projects04/PatriciaPossik/singlecell_am/bin/objects/prj_list.rds"

# Directory to save the processed Seurat object and intermediate files.
output_directory <- "/data04/projects04/PatriciaPossik/singlecell_am/bin/objects/"

# Path to the R script containing custom functions (e.g., functions.R).
functions_script_path <- "/data04/projects04/PatriciaPossik/singlecell_am/bin/scripts/Functions.R"

# Custom R library path to load packages from.
custom_lib_paths <- "/home/pedrosrb/singlecell_am/lib/Rpackages/" 

# Key in the input RDS list specifying the Seurat object to process.
list_key_to_process <- "All_Primary"

# --- Quality Control (QC) Parameters ----
# Maximum mitochondrial percentage allowed for cells during filtering.
mt_threshold <- 5

# Random seed for doublet detection (scDblFinder) to ensure reproducibility.
doublet_detection_seed <- 28

# --- Integration Parameters ----
# Metadata variable to group by for batch correction during Harmony integration
group_by_variable_for_integration <- "orig.ident"

# Comma-separated metadata variables to regress out during SCTransform.
variables_to_regress <- c("percent.mt")

# Number of principal components (PCs) to compute and use for PCA and Harmony integration.
num_principal_components <- 50

# --- Clustering Parameters ----
# Dimensions from the integrated reduction (e.g., Harmony) to use for clustering.
clustering_dimensions <- "1:30"

# Numeric vector of resolutions for Seurat's FindClusters function.
clustering_resolutions_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# --- UMAP Parameters ----
# Dimensions from the integrated reduction (e.g., Harmony) to use for UMAP embedding.
umap_dimensions <- "1:30"

# Reduction to use for UMAP (e.g., "harmony", "pca").
umap_reduction_method <- "harmony"

# --- Annotation Parameters ----
# Reference dataset for SingleR annotation.
singler_reference_dataset <- "BlueprintEncodeData"

# --- Subsetting Parameters ----
# Cell type to isolate from the annotated object. Set to NULL or "" to skip this step.
cell_type_to_subset <- "Melanocytes"

# Suffix to add to the output filename for the subsetted object.
subset_output_suffix <- "_melanocytes_only"

