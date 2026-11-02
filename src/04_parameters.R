# Parameters of the single-cell pipeline: paths, QC, integration, panels

# ---- Load here library to define paths ----
library("here")

# ---- Input/Output Paths ----

output_directory <- here("results")

cellranger_directory <- here("data", "count10x_outputs")

functions_script_path <- here("src", "03_functions.R")

custom_lib_paths <- NULL

# ---- Quality Control (QC) Parameters ----
mt_threshold <- 5

doublet_detection_seed <- 28

# ---- Integration Parameters ----
group_by_variable_for_integration <- "orig.ident"

variables_to_regress <- c("percent.mt")

num_principal_components <- 50

# ---- Clustering Parameters ----
clustering_dimensions <- "1:30"

clustering_resolutions_to_test <- c(0.1, 0.2, 0.3, 0.4, 0.5)

# ---- UMAP Parameters ----
umap_dimensions <- "1:30"

umap_reduction_method <- "harmony"

# ---- Annotation Parameters ----
singler_reference_dataset <- "BlueprintEncodeData"

# ---- Subsetting Parameters ----
cell_type_to_subset <- "Melanocytes"

# ---- Samples ----
tumour_labels <- c(SRR20791206 = "AM1", SRR20791205 = "AM2", SRR20791203 = "AM3",
                   SRR20791201 = "AM4", SRR20791200 = "AM5", SRR20791199 = "AM6",
                   SRR20791198 = "CM1", SRR20791197 = "CM2", SRR20791204 = "CM3")

# ---- Plotting ----
genes_to_highlight <- c("CRKL", "DPAGT1", "SF3B4", "THOC2",
                        "SOX10", "MITF", "TYRP1")

cell_type_labels <- c("B-cells" = "B cells", "CD4+ T-cells" = "CD4+ T cells",
                      "CD8+ T-cells" = "CD8+ T cells", "Endothelial cells" = "Endothelial cells",
                      "Fibroblasts" = "Fibroblasts", "Melanocytes" = "Melanocytes",
                      "Monocytes" = "Monocytes", "NK cells" = "NK cells")

custom_colors <- NULL

# ---- Differential expression and internal AM/CM score ----
de_test <- "MAST"
de_min_pct <- 0.05
de_logfc_threshold <- 0.25
de_log2fc_cutoff <- 1
de_padj_cutoff <- 0.05
score_top_genes <- 100

# ---- Panels e, f ----
focus_gene <- "SF3B4"
subtype_colors <- c(AM = "#F8766D", CM = "#00BFC4")
volcano_colors <- c("n.s" = "grey", "down" = "lightblue", "up" = "firebrick2")
dotplot_colors <- c(low = "blue", mid = "grey", high = "firebrick2")
ci_z <- 1.96

# ---- Figure sizes (inches): width, height ----
figure_sizes <- list(a = c(6, 5), b = c(6, 5), c = c(6, 5), d = c(6, 5), e = c(3.2, 3.25), f = c(3.45, 3.25))
figure_dpi <- 300

# ---- Object creation ----
min_cells <- 3
min_features <- 200
mt_pattern <- "^MT-"
