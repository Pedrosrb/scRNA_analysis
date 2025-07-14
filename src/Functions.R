# title: "AP1 - Proj. Analysis"
# author: "Pedro Sodré"
# date: "2025-06-04"

# Load Required Packages
load_pkgs <- function(){
  required_packages <- c("BiocManager", "Seurat", "dplyr", "tidyr", "purrr", "ggplot2", 
                         "SingleR", "celldex", "glmGamPoi", "SingleCellExperiment", "scDblFinder",
                         "clustree", "gridExtra", "harmony", "patchwork",
                         "parallel", "MAST", "CellChat", "optparse")
  
  lapply(required_packages, library, character.only = TRUE)
}
    
# Create Seurat Objects from CellRanger outputs
create_seurat_list <- function(sample_dirs, output_csv = "qc_metrics.csv") {

  dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)
  
  message("Starting processing of ", length(sample_dirs), " samples for Seurat object creation...")
  
  seurat_list <- lapply(sample_dirs, function(dir) {
    sample_name <- basename(dir) 
    message("Processing: ", sample_name)
    
    seurat_obj <- CreateSeuratObject(
      counts = Read10X(file.path(dir, "outs/filtered_feature_bc_matrix")),
      min.cells = 3, 
      min.features = 200, 
      project = sample_name 
    )
    
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, "^MT-")
    return(seurat_obj)
  })
  
  qc_metrics <- do.call(rbind, lapply(seurat_list, function(obj) {
    with(obj[[]], data.frame(
      Sample = obj@project.name,
      nCells = ncol(obj),
      Mean_nFeatures = mean(nFeature_RNA),
      Median_nFeatures = median(nFeature_RNA),
      Mean_nCounts = mean(nCount_RNA),
      Median_nCounts = median(nCount_RNA),
      Mean_percent_mt = mean(percent.mt),
      Median_percent_mt = median(percent.mt)
    ))
  }))
  
  write.csv(qc_metrics, output_csv, row.names = FALSE)
  
  message("\nSuccessfully processed all ", length(seurat_list), " samples and saved QC metrics.")
  
  return(seurat_list)
}

# Process Doublets and filter based on mitochondrial content
filter_seurat_list <- function(seurat_objs_list, seed = 28, mt_threshold = 5) {
  message("Processing ", length(seurat_objs_list), " samples for filtering...")
  
  filtered_list <- lapply(seurat_objs_list, function(seurat_obj) {
    set.seed(seed) # Set seed for reproducibility of scDblFinder
    sample_name <- seurat_obj@project.name
    message("\nProcessing filtering for sample: ", sample_name)
    
    # Add mitochondrial percentage if missing (redundant if create_seurat_list is used, but good safeguard)
    if (!"percent.mt" %in% colnames(seurat_obj[[]])) {
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    }
    
    # Convert to SingleCellExperiment (SCE) and find doublets using scDblFinder
    # scDblFinder automatically adds doublet_score and doublet_class to the SCE object
    sce <- scDblFinder(as.SingleCellExperiment(seurat_obj))
    
    # Add doublet information from SCE back to Seurat object metadata
    seurat_obj <- AddMetaData(
      seurat_obj,
      metadata = list(
        doublet_score = sce$scDblFinder.score,
        doublet_label = sce$scDblFinder.class
      )
    )
    
    # Filter cells based on doublet label and mitochondrial percentage
    initial_cells <- ncol(seurat_obj)
    seurat_obj_filtered <- subset(seurat_obj,
                                  doublet_label == "singlet" & percent.mt < mt_threshold)
    filtered_cells <- ncol(seurat_obj_filtered)
    message(paste0("  Sample '", sample_name, "': Filtered ", initial_cells - filtered_cells,
                   " cells (", round((initial_cells - filtered_cells)/initial_cells * 100, 2), "%) out of ", initial_cells,
                   " total. Remaining cells: ", filtered_cells))
    
    return(seurat_obj_filtered)
  })
  
  message("\nAll samples filtered successfully.")
  return(filtered_list)
}

# Integrate Seurat Objects (SCTransform, PCA, Harmony)
integrate_seurat_list <- function(seurat_objs_list, gp = "orig.ident", vars.to.regress = "percent.mt", npcs = 50, verbose = FALSE) {
  
  if (length(seurat_objs_list) == 0) {
    stop("Input 'seurat_objs_list' is empty. Cannot integrate.")
  }
  
  if (length(seurat_objs_list) == 1) {
    message("Only one Seurat object provided. Skipping merge step, performing SCTransform, PCA, and Harmony directly.")
    merged_obj <- seurat_objs_list[[1]]
  } else {
    message("Merging ", length(seurat_objs_list), " Seurat objects...")

    merged_obj <- merge(x = seurat_objs_list[[1]],
                        y = seurat_objs_list[2:length(seurat_objs_list)])
    message("Objects merged.")
  }
  
  message("Performing SCTransform...")
  merged_obj <- SCTransform(merged_obj,
                            vars.to.regress = vars.to.regress,
                            verbose = verbose,
                            method = "glmGamPoi") 
  message("SCTransform complete.")
  
  message("Running PCA...")
  merged_obj <- RunPCA(merged_obj,
                       npcs = npcs,
                       verbose = verbose)
  message("PCA complete.")
  
  message("Running Harmony integration...")
  if (!gp %in% colnames(merged_obj@meta.data)) {
    stop(paste("Grouping variable '", gp, "' not found in metadata. Cannot perform Harmony integration."))
  }
  
  integrated_obj <- RunHarmony(merged_obj,
                               group.by.vars = gp,
                               assay.use = "SCT", 
                               verbose = verbose)
  message("Harmony integration complete.")
  
  return(integrated_obj)
}

# Clustering of Seurat Object
cluster_seurat_obj <- function(seurat_obj, dims = 1:30, resolutions = seq(0.1, 0.5, by = 0.1), reduction = "harmony") {
  
  message("Performing clustering...")
  message(paste0("  Finding neighbors using reduction: '", reduction, "' and dims: ", paste(range(dims), collapse = "-"), "..."))
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, reduction = reduction)
  message("  Neighbors found.")
  
  for (res in resolutions) {
    message(paste0("  Finding clusters for resolution: ", res, "..."))
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  }
  message("Clustering complete for all specified resolutions.")
  return(seurat_obj)
}

# Annotate Cells using SingleR
annotate_singler <- function(seurat_obj, ref_dataset_name = "BlueprintEncodeData") {
  message(paste0("Annotating cells using SingleR with reference: '", ref_dataset_name, "'..."))

  ref <- tryCatch({
    eval(parse(text = paste0("celldex::", ref_dataset_name, "()")))
  }, error = function(e) {
    stop(paste0("Error loading SingleR reference '", ref_dataset_name, "'. Ensure it's a valid celldex dataset name and celldex is installed: ", e$message))
  })

  if ("SCT" %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- "SCT"
    message("Using 'SCT' assay for SingleR annotation.")
  } else {
    message("SCT assay not found. Using default assay for SingleR annotation.")
  }
  
  annotations <- SingleR(test = GetAssayData(seurat_obj, slot = "data"),
                         ref = ref,
                         labels = ref$label.main)

  seurat_obj@meta.data$SingleR.labels <- annotations$labels
  message("SingleR annotation complete.")
  
  return(seurat_obj)
}

# Calculate annotation statistics (cell counts and percentages)
annotation_stats <- function(seurat_obj, group_cols) {
  message("Calculating annotation statistics...")
  valid_cols <- group_cols[group_cols %in% colnames(seurat_obj@meta.data)]
  
  if (length(valid_cols) == 0) {
    stop("None of the provided columns for 'group_cols' are present in the Seurat object metadata.")
  }
  
  annotation_stats <- seurat_obj@meta.data %>%
    group_by(across(all_of(valid_cols))) %>% 
    summarise(
      Cell_Count = n(), 
      Percentage = (Cell_Count / nrow(seurat_obj@meta.data)) * 100, 
      .groups = "drop" 
    ) %>%
    arrange(desc(Cell_Count)) 
  
  message("Annotation statistics generated.")
  return(annotation_stats)
}

