# title: "Processing and Integration Script"
# author: "Pedro Sodré"
# date: "2025-06-04"

# --- Setup Project Root ----
library("here")
here::i_am("src/02_main.R")

# --- Activate renv Environment ----
tryCatch({
  source(here("renv", "activate.R"))
}, error = function(e) {
  stop("renv environment not found. Using default libraries.")
})

# --- Load Parameters ----
tryCatch({
  source(here("src", "04_parameters.R")) 
}, error = function(e) {
  stop("Error: Could not find or source 'src/04_Parameters.R'. Make sure the file exists.")
})

# --- Validate and Process Parameters ---
if (is.null(functions_script_path) || !file.exists(functions_script_path)) {
  stop(paste0("Error: Functions script not found. Check the 'functions_script_path' in your parameters file. Path provided: '", functions_script_path, "'"))
}

if (length(variables_to_regress) == 0) {
  variables_to_regress <- NULL
}

clustering_dimensions_numeric <- eval(parse(text = clustering_dimensions))
umap_dimensions_numeric <- eval(parse(text = umap_dimensions))

# --- Setup Environment ----
message(paste0("--- Starting Single-Cell Processing Script (", Sys.time(), ") ---"))

# Set custom R library paths if provided in parameters.R
if (!is.null(custom_lib_paths) && all(dir.exists(custom_lib_paths))) {
  .libPaths(c(custom_lib_paths, .libPaths()))
  message(paste0("Custom R library paths set to: ", paste(custom_lib_paths, collapse = ", ")))
}

# Create output directory
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
message(paste0("Output directory set to: '", output_directory, "'"))

# Source custom functions (path taken from parameters.R)
message(paste0("Sourcing custom functions from: '", functions_script_path, "'..."))
tryCatch({
  source(functions_script_path)
}, error = function(e) {
  stop(paste0("Error sourcing functions script: ", e$message))
})
message("Custom functions sourced successfully.")

# Load required packages (using function from functions.R)
message("Loading required R packages...")
tryCatch({
  load_pkgs() 
}, error = function(e) {
  stop(paste0("Error loading packages via load_pkgs(): ", e$message))
})
message("Required packages loaded.")

# --- Main Processing Workflow ----

# Step 1: Load Seurat Objects List
message("Step 1: Creating Seurat objects list from Cell Ranger output'")
seurat_objs_list <- tryCatch({
  create_seurat_list(
    sample_dirs = list.dirs(cellranger_directory, full.names = TRUE, recursive = FALSE), 
    output_csv = file.path(output_directory, "initial_qc_metadata.csv") 
  )
}, error = function(e) {
  stop(paste0("Error creating seurat objects list: ", e$message))
})
message("Seurat objects list loaded.")

# Step 2: Filtering (Mitochondrial percentage and Doublets)
message(paste0("Step 2: Filtering Seurat object (MT% < ", mt_threshold, ", doublet detection with seed ", doublet_detection_seed, ")..."))

aux_obj_filtered_list <- tryCatch({
  filter_seurat_list(seurat_objs_list, mt_threshold = mt_threshold, seed = doublet_detection_seed)
}, error = function(e) {
  stop(paste0("Error during filtering: ", e$message))
})
message("Filtering complete.")

# Step 3: Integration (SCTransform, PCA, Harmony)
message(paste0("Step 3: Integrating Seurat object using Harmony (grouping by '", group_by_variable_for_integration, "', regressing '",
               paste(variables_to_regress, collapse = ", "), "', using ", num_principal_components, " PCs)..."))

aux_obj_integrated <- tryCatch({
  integrate_seurat_list(aux_obj_filtered_list,
                        gp = group_by_variable_for_integration,
                        vars.to.regress = variables_to_regress,
                        npcs = num_principal_components)
}, error = function(e) {
  stop(paste0("Error during integration: ", e$message))
})
aux_obj_integrated$Type <- ifelse(seurat_obj$orig.ident %in% c("SRR20791198", "SRR20791197", "SRR20791204"), "CM", "AM")
message("Integration complete.")

# Step 4: Clustering
message(paste0("Step 4: Performing clustering (using dims ", paste(clustering_dimensions_numeric, collapse = ","),
               ", resolutions: ", paste(clustering_resolutions_to_test, collapse = ","), ")..."))
aux_obj_clustered <- tryCatch({
  cluster_seurat_obj(aux_obj_integrated,
                     dims = clustering_dimensions_numeric,
                     resolutions = clustering_resolutions_to_test,
                     reduction = umap_reduction_method)
}, error = function(e) {
  stop(paste0("Error during clustering: ", e$message))
})
message("Clustering complete.")

# Step 5: UMAP 
message(paste0("Step 5: Running UMAP embedding (using reduction '", umap_reduction_method, "', dims ",
               paste(umap_dimensions_numeric, collapse = ","), ")..."))
aux_obj_umap <- tryCatch({
  RunUMAP(aux_obj_clustered,
          reduction = umap_reduction_method,
          dims = umap_dimensions_numeric)
}, error = function(e) {
  stop(paste0("Error during UMAP embedding: ", e$message))
})
message("UMAP embedding complete.")

# Step 6: Cell Type Annotation (SingleR)
message(paste0("Step 6: Annotating cell types using SingleR with reference: '", singler_reference_dataset, "'..."))
aux_obj_annotated <- tryCatch({
  annotate_singler(aux_obj_umap, ref_dataset_name = singler_reference_dataset)
}, error = function(e) {
  stop(paste0("Error during SingleR annotation: ", e$message))
})
message("Cell type annotation complete.")

# Step 7: Subset Specific Cell Type 
subset_obj <- NULL 
if (!is.null(cell_type_to_subset) && nchar(cell_type_to_subset) > 0) {
  message(paste0("--- Step 7: Subsetting for '", cell_type_to_subset, "' cells ---"))
  
  if (!"SingleR.labels" %in% colnames(aux_obj_annotated@meta.data)) {
    stop("Error: 'SingleR.labels' column not found in metadata. Cannot subset.")
  }
  
  cells_before_subset <- ncol(aux_obj_annotated)
  subset_obj <- tryCatch({
    subset(aux_obj_annotated, subset = SingleR.labels == cell_type_to_subset)
  }, error = function(e) {
    stop(paste0("Error during subsetting: ", e$message))
  })
  cells_after_subset <- ncol(subset_obj)
  
  message(paste0("Subsetting complete. Kept ", cells_after_subset, " '", cell_type_to_subset, "' cells out of ", cells_before_subset, "."))
  
} else {
  message("--- Step 7: No cell type specified for subsetting. Skipping. ---")
}

# --- Plotting and Saving Results ---
# Step 8: Plotting - Figure 7A (UMAP of all cell types)
fig7a <- DimPlot(aux_obj_annotated,
                 group.by = "SingleR.labels", 
                 label = TRUE, 
                 repel = TRUE, 
                 raster = FALSE,
                 cols = custom_colors) +
  labs(
    title = "UMAP - Primary AM and CM", 
    subtitle = "All Cell Types") +
  theme_classic2()

ggsave(filename = file.path(output_directory, "Fig7A_UMAP_All_Cell_Types.pdf"),
       plot = fig7a, width = 12, height = 8, dpi = 300)
message("Saved Figure 7A.")

# Step 9: Plotting - Figure 7B (DotPlot of marker genes)
fig7b <- DotPlot(aux_obj_annotated,
                 features = genes_to_highlight,
                 assay = "SCT",
                 dot.scale = 8,
                 group.by = "SingleR.labels") +
  scale_color_gradient2(low = "blue", mid = "grey", high = "firebrick2") + 
  labs(
    title = "DotPlot - Marker Gene Expression",
    subtitle = "All Cell Types",
    x = "Features",
    y = "Cell Types") +
  theme_classic2() +
  RotatedAxis()

ggsave(filename = file.path(output_directory, "Fig7B_DotPlot_Markers.pdf"),
       plot = fig7b, width = 10, height = 8, dpi = 300)
message("Saved Figure 7B.")

# Step 10: Plotting - Figure 7C (UMAP of subset)
if (!is.null(subset_obj)) {
  fig7c <- DimPlot(subset_obj,
                   group.by = "Type", 
                   label = FALSE, 
                   repel = TRUE, 
                   raster = FALSE) +
    labs(
      title = "UMAP - Primary AM and CM", 
      subtitle = paste(cell_type_to_subset, "Subset")) +
    theme_classic2()
  
  ggsave(filename = file.path(output_directory, "Fig7C_UMAP_Subset.pdf"),
         plot = fig7c, width = 12, height = 8, dpi = 300)
  message("Saved Figure 7C.")
}

# Step 11: Plotting - Figure 7D (Volcano Plot)
if (!is.null(subset_obj)) {
  de_results <- FindMarkers(
    object = subset_obj,
    group.by = "Type", 
    ident.1 = "AM",     
    ident.2 = "CM",     
    assay = "SCT",
    test.use = "MAST",
    min.pct = 0.05,
    logfc.threshold = 0.25
  )
  
  de_results$threshold <- "Not significant"
  de_results$threshold[de_results$avg_log2FC > 1 & de_results$p_val_adj < 0.05] <- "Upregulated"
  de_results$threshold[de_results$avg_log2FC < -1 & de_results$p_val_adj < 0.05] <- "Downregulated"
  
  highlight_genes_data <- de_results[rownames(de_results) %in% genes_to_highlight, ]
  
  fig7d <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = threshold), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red")) +
    geom_point(data = highlight_genes_data, aes(x = avg_log2FC, y = -log10(p_val_adj)), color = "black", size = 2, shape = 1) +
    ggrepel::geom_text_repel(data = highlight_genes_data, aes(label = rownames(highlight_genes_data)), box.padding = 0.5) +
    labs(title = "Volcano Plot: DEGs AM vs CM", x = "Log2 Fold Change", y = "-Log10(Adjusted p-value)") +
    theme_classic2()
  
  ggsave(filename = file.path(output_directory, "Fig7D_Volcano_Plot.pdf"),
         plot = fig7d, width = 12, height = 8, dpi = 300)
  message("Saved Figure 7D.")
}

# # Step 12: Plotting - Figure 7E (Violin Plot)
am_like_genes <- rownames(head(de_results[de_results$threshold == "Upregulated", ][order(de_results[de_results$threshold == "Upregulated", ]$p_val_adj), ], 100))
  
cm_like_genes <- rownames(head(de_results[de_results$threshold == "Downregulated", ][order(de_results[de_results$threshold == "Downregulated", ]$p_val_adj), ], 100))

subset_obj <- AddModuleScore(subset_obj, features = list(am_like_genes), name = "am_like_score", assay = "SCT")
subset_obj <- AddModuleScore(subset_obj, features = list(cm_like_genes), name = "cm_like_score", assay = "SCT")

subset_obj$am_cm_like_ratio <- subset_obj$am_like_score1 - subset_obj$cm_like_score1

max_CM <- max(subset_obj@meta.data$am_cm_like_ratio[subset_obj@meta.data$Type == "CM"], na.rm = TRUE)
min_AM <- min(subset_obj@meta.data$am_cm_like_ratio[subset_obj@meta.data$Type == "AM"], na.rm = TRUE)

fig7e <- VlnPlot(subset_obj, 
                 features = "am_cm_like_ratio", 
                 group.by = "Type", 
                 assay = "SCT", 
                 pt.size = 0.0) +
  labs(title = "Top DEGs Ratio Score") +
  geom_hline(yintercept = max_CM, color = "black", linetype = "dashed") +
  geom_hline(yintercept = min_AM, color = "black", linetype = "dashed") +
  theme_classic2()

ggsave(filename = file.path(output_directory, "Fig7E_VlnPlot_DEGs_Score.pdf"),
       plot = fig7e, width = 8, height = 6, dpi = 300)
message("Saved Figure 7E.")
 
# # Step 13: Plotting - Figure 7F (Linear Model)
model_data <- data.frame(
  am_cm_like_ratio = subset_obj@meta.data$am_cm_like_ratio,
  SF3B4_Expression = GetAssayData(subset_obj, assay = "SCT", layer = "data")["SF3B4", ], 
  Type = factor(subset_obj@meta.data$Type) 
)

lm_model <- lm(am_cm_like_ratio ~ SF3B4_Expression * Type, data = model_data)

ggsave(filename = file.path(output_directory, "Linear_Model_Validation.pdf"),
       plot = summary(lm_fit), width = 8, height = 6, dpi = 300)


fig7f <- ggplot(model_data, aes(x = SF3B4_Expression, y = am_cm_like_ratio, color = Type, fill = Type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  labs(
    title = "Top DEGs Score Ratio vs SF3B4 Expression",
    x = "SF3B4 Expression (Normalized)",
    y = "AM/CM-like ratio"
  ) +
  theme_classic2()

ggsave(filename = file.path(output_directory, "Fig7F_LM_Plot.pdf"),
       plot = fig7f, width = 10, height = 6, dpi = 300)
message("Saved Figure 7F.")

message(paste0("--- Single-Cell Processing Script Finished (", Sys.time(), ") ---"))
