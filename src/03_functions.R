# title: "Single-Cell Analysis Functions"
# author: "Pedro Sodré"
# date: "2025-06-04"

# Load Required Packages
load_pkgs <- function(){
  required_packages <- c("BiocManager", "Seurat", "dplyr", "tidyr", "purrr", "ggplot2", 
                         "SingleR", "celldex", "glmGamPoi", "SingleCellExperiment", "scDblFinder",
                         "clustree", "gridExtra", "harmony", "patchwork", "ggrepel",
                         "parallel", "MAST", "CellChat", "optparse", "here", "ggpubr", "renv")
  
  suppressPackageStartupMessages({
    lapply(required_packages, library, character.only = TRUE)
  })
}

# Create a list of Seurat Objects from CellRanger output directories
create_seurat_list <- function(sample_dirs, output_csv) {
  
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
    
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
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
  
  message("\nSuccessfully processed all ", length(seurat_list), " samples and saved QC metrics to '", output_csv, "'.")
  
  return(seurat_list)
}

# Filter a list of Seurat objects for doublets and high mitochondrial content
filter_seurat_list <- function(seurat_objs_list, seed = 28, mt_threshold = 5) {
  message("Processing ", length(seurat_objs_list), " Seurat object(s) for filtering...")
  
  filtered_list <- lapply(seurat_objs_list, function(seurat_obj) {
    set.seed(seed)
    sample_name <- seurat_obj@project.name
    message("\nProcessing filtering for sample: ", sample_name)

    if (!"percent.mt" %in% colnames(seurat_obj[[]])) {
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    }
    
    sce <- scDblFinder(as.SingleCellExperiment(seurat_obj))
    
    seurat_obj$doublet_score <- sce$scDblFinder.score
    seurat_obj$doublet_label <- sce$scDblFinder.class

    initial_cells <- ncol(seurat_obj)
    seurat_obj_filtered <- subset(seurat_obj,
                                  subset = doublet_label == "singlet" & percent.mt < mt_threshold)
    filtered_cells <- ncol(seurat_obj_filtered)
    message(paste0("  Sample '", sample_name, "': Filtered ", initial_cells - filtered_cells,
                   " cells (", round((initial_cells - filtered_cells)/initial_cells * 100, 2), "%) out of ", initial_cells,
                   " total. Remaining cells: ", filtered_cells))
    
    return(seurat_obj_filtered)
  })
  
  message("\nAll samples filtered successfully.")
  return(filtered_list)
}

# Integrate a list of Seurat Objects using SCTransform, PCA, and Harmony
integrate_seurat_list <- function(seurat_objs_list, gp = "orig.ident", vars.to.regress = "percent.mt", npcs = 50, verbose = FALSE) {
  
  if (length(seurat_objs_list) == 0) {
    stop("Input 'seurat_objs_list' is empty. Cannot integrate.")
  }
  
  if (length(seurat_objs_list) == 1) {
    message("Only one Seurat object provided. Skipping merge step, performing SCTransform, PCA, and Harmony directly.")
    merged_obj <- seurat_objs_list[[1]]
  } else {
    message("Merging ", length(seurat_objs_list), " Seurat objects...")
    merged_obj <- merge(x = seurat_objs_list[[1]], y = seurat_objs_list[2:length(seurat_objs_list)])
    message("Objects merged.")
  }
  
  message("Performing SCTransform...")
  merged_obj <- SCTransform(merged_obj,
                            vars.to.regress = vars.to.regress,
                            verbose = verbose,
                            method = "glmGamPoi") 
  message("SCTransform complete.")
  
  message("Running PCA...")
  merged_obj <- RunPCA(merged_obj, npcs = npcs, verbose = verbose)
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

# Perform clustering on a Seurat Object
cluster_seurat_obj <- function(seurat_obj, dims = 1:30, resolutions = seq(0.1, 0.5, by = 0.1), reduction = "harmony") {
  
  message("Performing clustering...")
  message(paste0("  Finding neighbors using reduction: '", reduction, "' and dims: ", paste(range(dims), collapse = "-"), "..."))
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims, reduction = reduction)
  message("  Neighbors found.")
  
  # Find clusters for each specified resolution
  for (res in resolutions) {
    message(paste0("  Finding clusters for resolution: ", res, "..."))
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  }
  message("Clustering complete for all specified resolutions.")
  return(seurat_obj)
}

# Annotate cells using SingleR
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
  
  seurat_obj$SingleR.labels <- annotations$labels
  message("SingleR annotation complete.")
  
  return(seurat_obj)
}
