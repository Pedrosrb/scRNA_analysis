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
  inst_load_pkgs() 
}, error = function(e) {
  stop(paste0("Error loading packages via inst_load_pkgs(): ", e$message))
})
message("Required packages loaded.")

# --- Main Processing Workflow ----

# Step 1: Load Seurat Objects List
message(paste0("Step 1: Loading Seurat objects list from '", input_rds_path, "'..."))
seurat_objs_list <- tryCatch({
  readRDS(input_rds_path)
}, error = function(e) {
  stop(paste0("Error loading RDS file: ", e$message))
})
message("Seurat objects list loaded.")

# Check if the specified key exists in the list
if (!list_key_to_process %in% names(seurat_objs_list)) {
  stop(paste0("Error: Key '", list_key_to_process, "' not found in the input RDS list. Available keys: ",
              paste(names(seurat_objs_list), collapse = ", ")))
}

# Extract the Seurat object to be processed
aux_obj <- seurat_objs_list[[list_key_to_process]]
message(paste0("Processing Seurat object corresponding to key: '", list_key_to_process, "'"))


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
  message("--- Step 7: No cell type specified for subsetting. ---")
}

# --- Save Final Processed Object ----
message(paste0("--- Saving final processed Seurat object to '", file.path(output_directory, paste0(list_key_to_process, "_processed.rds")), "' ---"))

seurat_objs_list[[list_key_to_process]] <- aux_obj_annotated 

tryCatch({
  saveRDS(seurat_objs_list, file = file.path(output_directory, paste0(list_key_to_process, "_processed.rds")))
  message("Final processed Seurat object saved successfully.")
}, error = function(e) {
  stop(paste0("Error saving final RDS file: ", e$message))
})

# Save the subsetted object to a separate file
message(paste0("--- Saving subsetted object to '", subset_filename, "' ---"))

tryCatch({
  saveRDS(subset_obj, file = file.path(output_directory, paste0(list_key_to_process, subset_output_suffix, ".rds")))
  message("Subsetted Seurat object saved successfully.")
}, error = function(e) {
  stop(paste0("Error saving subsetted RDS file: ", e$message))
})

message(paste0("--- Script finished successfully (", Sys.time(), ") ---"))

