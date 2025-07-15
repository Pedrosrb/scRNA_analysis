# title: "Processing and Integration Script"
# author: "Pedro Sodré"
# date: "2025-06-04"

# --- Load Parameters ----
tryCatch({
  source("scripts/Parameters.R") 
}, error = function(e) {
  stop(paste0("Error: Could not find or source 'scripts/Parameters.R'. ..."))
})

# --- Validate and Process Parameters ---
if (is.null(input_rds_path) || !file.exists(input_rds_path)) {
  stop(paste0("Error: Input RDS file not found or not specified at '", input_rds_path, "'"))
}
if (is.null(functions_script_path) || !file.exists(functions_script_path)) {
  stop(paste0("Error: Functions script not found or not specified at '", functions_script_path, "'"))
}

if (length(variables_to_regress) == 0) {
  variables_to_regress <- NULL
}

clustering_dimensions_numeric <- eval(parse(text = clustering_dimensions))
umap_dimensions_numeric <- eval(parse(text = umap_dimensions))


# --- Setup Environment ----
message(paste0("--- Starting Single-Cell Processing Script (", Sys.time(), ") ---"))

# Set custom R library paths if provided in parameters.R
if (!is.null(custom_lib_paths)) {
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
  stop(paste0("Error loading packages via inst_load_pkgs(): ", e$message))
})
message("Required packages loaded.")

# --- Main Processing Workflow ----

# Step 1: Load Seurat Objects List
message("Step 1: Creating Seurat objects list from Cell Ranger output'")
seurat_objs_list <- tryCatch({
  create_seurat_list(
    sample_dirs = list.dirs(cellrange_directory,
                            full.names = TRUE, recursive = FALSE), output_csv = paste0(output_directory, "zhang_metadata.csv"))
}, error = function(e) {
  stop(paste0("Error creating seurat objects list: ", e$message))
})
message("Seurat objects list loaded.")

# Extract the Seurat object to be processed
aux_obj <- seurat_objs_list
message("Processing Seurat objects list'")


# Step 2: Filtering (Mitochondrial percentage and Doublets)
message(paste0("Step 2: Filtering Seurat object (MT% < ", mt_threshold, ", doublet detection with seed ", doublet_detection_seed, ")..."))
aux_obj_filtered_list <- tryCatch({
  filter_seurat_list(list(aux_obj), mt_threshold = mt_threshold, seed = doublet_detection_seed)
}, error = function(e) {
  stop(paste0("Error during filtering: ", e$message))
})
aux_obj_filtered <- aux_obj_filtered_list[[1]]
message("Filtering complete. Number of cells after filtering: ", ncol(aux_obj_filtered))


# Step 3: Integration (SCTransform, PCA, Harmony)
message(paste0("Step 3: Integrating Seurat object using Harmony (grouping by '", group_by_variable_for_integration, "', regressing '",
               paste(variables_to_regress, collapse = ", "), "', using ", num_principal_components, " PCs)..."))
aux_obj_integrated_list <- tryCatch({
  integrate_seurat_list(list(aux_obj_filtered),
                        gp = group_by_variable_for_integration,
                        vars.to.regress = variables_to_regress,
                        npcs = num_principal_components)
}, error = function(e) {
  stop(paste0("Error during integration: ", e$message))
})
aux_obj_integrated <- aux_obj_integrated_list[[1]]
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
if (!is.null(cell_type_to_subset) && nchar(cell_type_to_subset) > 0) {
  message(paste0("--- Step 7: Subsetting for '", cell_type_to_subset, "' cells ---"))
  
  # Check if the annotation column exists
  if (!"SingleR.labels" %in% colnames(aux_obj_annotated@meta.data)) {
    stop("Error: 'SingleR.labels' column not found in metadata. Cannot subset.")
  }
  
  # Perform the subsetting
  cells_before_subset <- ncol(aux_obj_annotated)
  subset_obj <- tryCatch({
    subset(aux_obj_annotated, subset = SingleR.labels == cell_type_to_subset)
  }, error = function(e) {
    stop(paste0("Error during subsetting: ", e$message))
  })
  cells_after_subset <- ncol(subset_obj)
  
  message(paste0("Subsetting complete. Kept ", cells_after_subset, " '", cell_type_to_subset, "' cells out of ", cells_before_subset, "."))
  
} else {
  message("--- No cell type specified for subsetting. ---")
}

# Step 8: Plotting - Figure 7A
fig7a <- DimPlot(aux_obj_annotated,
                 group.by = "SingleR.labels", 
                 label = T, 
                 repel = T, 
                 raster = F,
                 cols = custom_colors) +
  labs(
    title = "UMAP - Primary AM and CM", 
    subtitle = "All Cell Types") +
  theme_classic2() 

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "UMAP - Cell Types.pdf"),
       plot = fig7a, width = 12, height = 8)

# Step 9: Plotting - Figure 7B
fig7b <- DotPlot(seurat_obj,
                 features = genes,
                 assay = "SCT",
                 dot.scale = 20,
                 group.by = "SingleR.labels") +
  
  scale_color_gradient2(low = "blue", mid = "grey", high = "firebrick2") + 
  
  labs(
    title = "DotPlot - Primary Dataset",
    subtitle = "All Cell Types",
    x = "Features",
    y = "Cell Types") +
  theme_classic2() +
  RotatedAxis()

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "Dot Plot - Features by Cell Types .pdf"),
       plot = fig7b, width = 10, height = 8)

# Step 10: Plotting - Figure 7C
fig7c <- DimPlot(subset_obj,
                group.by = "Type", 
                label = F, 
                repel = T, 
                raster = F) +
  labs(
    title = "UMAP - Primary AM and CM", 
    subtitle = "Melanocytes Subset") +
  theme_classic2() 

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "UMAP - Melanocytes Subset.pdf"),
       plot = fig7c, width = 12, height = 8)

# Step 11: Plotting - Figure 7D
de_results <- FindMarkers(
  object = subset_obj,
  group.by = "Type",
  ident.1 = "AM",
  ident.2 = "CM",
  assay = "SCT",
  test.use = "MAST",
  min.pct = 0.05,
  logfc.threshold = 0.25,
  only.pos = FALSE
)

highlight_genes_data <- de_results[genes_to_highlight, ]
current_logfc_threshold <- 1 
current_p_val_cutoff <- 0.05

fig7d <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.3) + # All genes as points
  scale_color_manual(
    values = c("Downregulated" = "lightblue", "Not significant" = "grey", "Upregulated" = "firebrick2")
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = paste0("Volcano Plot: DEGs AM vs CM"),
    subtitle = paste0("Highlighted Genes: ", paste(rownames(highlight_genes_data), collapse = ", ")),
    x = paste0("Log2 Fold Change (AM vs CM)"),
    y = "-Log10(Adjusted p-value)",
    color = "Significance"
  ) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic")
  ) +
  # Add dashed lines for logFC and p-value cutoffs
  geom_vline(xintercept = c(-current_logfc_threshold, current_logfc_threshold),
             linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
  geom_hline(yintercept = -log10(current_p_val_cutoff),
             linetype = "dashed", color = "darkgrey", linewidth = 0.5) +
  
  # --- Step 4: Highlight specified genes ---
  # Add highlighted genes as distinct points
  geom_point(data = highlight_genes_data, aes(x = avg_log2FC, y = -log10(p_val_adj)),
             color = "black", size = 0.8, shape = 21, fill = "black", stroke = 1) + # Use a distinct color and shape
  
  # Add labels for highlighted genes
  geom_text_repel(data = highlight_genes_data, aes(label = rownames(highlight_genes_data)),
                  box.padding = 0.5, point.padding = 0.5,
                  min.segment.length = 0, # Draw segment line even if label is close
                  segment.color = 'black',
                  color = "black", fontface = "bold", max.overlaps = Inf) + theme_classic2()

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "Volcano Plot - DEGs AM vs CM.pdf"),
       plot = fig7d, width = 12, height = 8)

# Step 12: Plotting - Figure 7E
fig7e <- VlnPlot(subset_obj, 
                features = "top_degs_score_ratio", 
                group.by = "Type", 
                assay = "SCT", 
                pt.size = 0.0, 
                raster = FALSE) +
  ggtitle("Top DEGs Ratio Score") +
  geom_hline(yintercept = max_CM, color = "black", linetype = "dashed") +
  geom_hline(yintercept = min_AM, color = "black", linetype = "dashed")

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "VlnPlot - Top DEGs Ratio Score.pdf"),
       plot = fig7e, width = 8, height = 6)

# Step 13: Plotting - Figure 7F
lm_model <- lm(top_degs_score_ratio ~ SF3B4_Expression * Type, data = model_data)

fig7f <- ggplot(model_data, aes(x = SF3B4_Expression, y = top_degs_score_ratio, color = Type, fill = Type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  labs(
    title = "Top DEGs Score Ratio vs SF3B4 Expression by Melanoma Type",
    x = "SF3B4 Expression (Normalized)",
    y = "Top DEGs Ratio Score",
    color = "Melanoma Type",
    fill = "Melanoma Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave(paste0("/data04/projects04/PatriciaPossik/singlecell_am/results/", "Linear Regression - Top DEGs Score Ratio vs SF3B4 Expression.pdf"),
       plot = fig7f, width = 10, height = 6)

