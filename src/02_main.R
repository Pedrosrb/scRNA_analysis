# =============================================================================
# Single-cell pipeline: Extended Data Figure 18, panels a-f
# Run:      source("src/02_main.R")
# Outputs:  output_directory (04_parameters.R)
# Inputs:   src/04_parameters.R
#   cellranger_directory  <SRR>/outs/filtered_feature_bc_matrix, one folder per sample
#   tumour_labels         sample folder -> tumour label (AM1-AM6, CM1-CM3)
# =============================================================================

# ---- Setup Project Root ----
library("here")
here::i_am("src/02_main.R")

# ---- Activate renv Environment ----
tryCatch({
  source(here("renv", "activate.R"))
}, error = function(e) {
  stop("renv environment not found. Using default libraries.")
})

# ---- Load Parameters ----
tryCatch({
  source(here("src", "04_parameters.R")) 
}, error = function(e) {
  stop("Error: Could not find or source 'src/04_Parameters.R'. Make sure the file exists.")
})

# ---- Validate and Process Parameters ----
if (is.null(functions_script_path) || !file.exists(functions_script_path)) {
  stop(paste0("Error: Functions script not found. Check the 'functions_script_path' in your parameters file. Path provided: '", functions_script_path, "'"))
}

if (length(variables_to_regress) == 0) {
  variables_to_regress <- NULL
}

clustering_dimensions_numeric <- eval(parse(text = clustering_dimensions))
umap_dimensions_numeric <- eval(parse(text = umap_dimensions))

# ---- Setup Environment ----
message(paste0("--- Starting Single-Cell Processing Script (", Sys.time(), ") ---"))

if (!is.null(custom_lib_paths) && all(dir.exists(custom_lib_paths))) {
  .libPaths(c(custom_lib_paths, .libPaths()))
  message(paste0("Custom R library paths set to: ", paste(custom_lib_paths, collapse = ", ")))
}

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
message(paste0("Output directory set to: '", output_directory, "'"))

message(paste0("Sourcing custom functions from: '", functions_script_path, "'..."))
tryCatch({
  source(functions_script_path)
}, error = function(e) {
  stop(paste0("Error sourcing functions script: ", e$message))
})
message("Custom functions sourced successfully.")

message("Loading required R packages...")
tryCatch({
  load_pkgs() 
}, error = function(e) {
  stop(paste0("Error loading packages via load_pkgs(): ", e$message))
})
message("Required packages loaded.")

# ---- Main Processing Workflow ----

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

message(paste0("Step 2: Filtering Seurat object (MT% < ", mt_threshold, ", doublet detection with seed ", doublet_detection_seed, ")..."))

aux_obj_filtered_list <- tryCatch({
  filter_seurat_list(seurat_objs_list, mt_threshold = mt_threshold, seed = doublet_detection_seed)
}, error = function(e) {
  stop(paste0("Error during filtering: ", e$message))
})
message("Filtering complete.")

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
unknown_samples <- setdiff(unique(aux_obj_integrated$orig.ident), names(tumour_labels))
if (length(unknown_samples) > 0) {
  stop(paste0("Samples missing from 'tumour_labels' in 04_parameters.R: ", paste(unknown_samples, collapse = ", ")))
}
aux_obj_integrated$tumour <- unname(tumour_labels[aux_obj_integrated$orig.ident])
aux_obj_integrated$Type <- substr(aux_obj_integrated$tumour, 1, 2)
message("Integration complete.")

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

message(paste0("Step 6: Annotating cell types using SingleR with reference: '", singler_reference_dataset, "'..."))
aux_obj_annotated <- tryCatch({
  annotate_singler(aux_obj_umap, ref_dataset_name = singler_reference_dataset)
}, error = function(e) {
  stop(paste0("Error during SingleR annotation: ", e$message))
})
message("Cell type annotation complete.")

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

if (!is.null(subset_obj)) {
  subset_obj <- RunUMAP(subset_obj,
                        reduction = umap_reduction_method,
                        dims = umap_dimensions_numeric)
}

# ---- Plotting and Saving Results: Extended Data Figure 18 ----
save_panel <- function(name, plot, panel) {
  ggsave(filename = file.path(output_directory, name), plot = plot,
         width = figure_sizes[[panel]][1], height = figure_sizes[[panel]][2], dpi = figure_dpi)
  message("Saved Figure 18", panel, ".")
}
fig18_theme <- theme_classic2() +
  theme(legend.position = "top", legend.justification = "left", legend.title = element_blank())

cells_in_figure <- colnames(aux_obj_annotated)[aux_obj_annotated$SingleR.labels %in% names(cell_type_labels)]
aux_obj_figure <- subset(aux_obj_annotated, cells = cells_in_figure)
aux_obj_figure$cell_type_label <- factor(unname(cell_type_labels[aux_obj_figure$SingleR.labels]),
                                         levels = unname(cell_type_labels))

fig18a <- DimPlot(aux_obj_figure,
                  group.by = "cell_type_label",
                  label = FALSE,
                  raster = FALSE,
                  cols = custom_colors) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
  fig18_theme

save_panel("Fig18a_UMAP_All_Cell_Types.pdf", fig18a, "a")

fig18b <- DotPlot(aux_obj_figure,
                  features = genes_to_highlight,
                  assay = "SCT",
                  dot.scale = 8,
                  group.by = "cell_type_label") +
  scale_color_gradient2(low = dotplot_colors[["low"]], mid = dotplot_colors[["mid"]], high = dotplot_colors[["high"]]) +
  guides(color = guide_colorbar(title = "expression\nlevel"), size = guide_legend(title = "% of cells")) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_classic2() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

save_panel("Fig18b_DotPlot_Markers.pdf", fig18b, "b")

if (is.null(subset_obj)) {
  stop("Figures 18c-f need the melanocyte subset: set 'cell_type_to_subset' in 04_parameters.R.")
}

fig18c <- DimPlot(subset_obj,
                  group.by = "Type",
                  label = FALSE,
                  raster = FALSE,
                  cols = subtype_colors) +
  labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
  fig18_theme

save_panel("Fig18c_UMAP_Subset.pdf", fig18c, "c")

subset_obj <- PrepSCTFindMarkers(subset_obj, verbose = FALSE)
de_results <- FindMarkers(
  object = subset_obj,
  group.by = "Type",
  ident.1 = "AM",
  ident.2 = "CM",
  assay = "SCT",
  test.use = de_test,
  min.pct = de_min_pct,
  logfc.threshold = de_logfc_threshold,
  recorrect_umi = FALSE
)
de_results$gene <- rownames(de_results)

de_results$threshold <- "n.s"
de_results$threshold[de_results$avg_log2FC > de_log2fc_cutoff & de_results$p_val_adj < de_padj_cutoff] <- "up"
de_results$threshold[de_results$avg_log2FC < -de_log2fc_cutoff & de_results$p_val_adj < de_padj_cutoff] <- "down"
de_results$threshold <- factor(de_results$threshold, levels = c("n.s", "down", "up"))
write.table(de_results, file.path(output_directory, "Fig18d_markers_MAST_AM_vs_CM.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

highlight_genes_data <- de_results[rownames(de_results) %in% genes_to_highlight, ]

fig18d <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = threshold), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-de_log2fc_cutoff, de_log2fc_cutoff), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(de_padj_cutoff), linetype = "dashed", color = "grey60") +
  geom_point(data = highlight_genes_data, color = "black", size = 1.5) +
  ggrepel::geom_text_repel(data = highlight_genes_data, aes(label = gene), fontface = "italic",
                           box.padding = 0.5, min.segment.length = 0, seed = 1) +
  labs(title = NULL, x = expression(Log[2]~"Fold Change (AM vs CM)"), y = expression(-Log[10]~"(adj. p value)")) +
  fig18_theme

save_panel("Fig18d_Volcano_Plot.pdf", fig18d, "d")

up_genes <- de_results[de_results$threshold == "up", ]
down_genes <- de_results[de_results$threshold == "down", ]
am_like_genes <- head(up_genes$gene[order(up_genes$p_val_adj)], score_top_genes)
cm_like_genes <- head(down_genes$gene[order(down_genes$p_val_adj)], score_top_genes)
write.table(rbind(data.frame(direction = "AM-like", gene = am_like_genes),
                  data.frame(direction = "CM-like", gene = cm_like_genes)),
            file.path(output_directory, "internal_score_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

subset_obj <- AddModuleScore(subset_obj, features = list(am_like_genes), name = "am_like_score", assay = "SCT")
subset_obj <- AddModuleScore(subset_obj, features = list(cm_like_genes), name = "cm_like_score", assay = "SCT")
subset_obj$internal_score <- subset_obj$am_like_score1 - subset_obj$cm_like_score1

rna_assay <- JoinLayers(subset_obj[["RNA"]])
cell_data <- data.frame(
  tumour = subset_obj$tumour,
  Type = factor(subset_obj$Type, levels = c("AM", "CM")),
  gene_umi = as.numeric(LayerData(rna_assay, layer = "counts")[focus_gene, colnames(subset_obj)]),
  gene_sct = as.numeric(GetAssayData(subset_obj, assay = "SCT", layer = "data")[focus_gene, ]),
  internal_score = subset_obj$internal_score,
  log10_umi = log10(subset_obj$nCount_RNA)
)
write.table(cell_data, file.path(output_directory, "cell_level_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
expressing <- cell_data[cell_data$gene_umi > 0, ]

panel_e <- expressing %>%
  group_by(tumour, Type) %>%
  summarise(mean = mean(gene_umi), sd = sd(gene_umi), n_expressing = n(), .groups = "drop") %>%
  mutate(lo = mean - ci_z * sd / sqrt(n_expressing), hi = mean + ci_z * sd / sqrt(n_expressing)) %>%
  left_join(cell_data %>% group_by(tumour) %>% summarise(melanocytes = n(), .groups = "drop"), by = "tumour") %>%
  mutate(detected_pct = 100 * n_expressing / melanocytes,
         tumour = factor(tumour, levels = rev(unname(tumour_labels))))
write.table(panel_e, file.path(output_directory, "Fig18e_gene_per_tumour.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fig18e <- ggplot(panel_e, aes(mean, tumour, color = Type)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.3) +
  geom_point(size = 2.2) +
  scale_color_manual(values = subtype_colors) +
  labs(x = bquote(italic(.(focus_gene))~"UMIs per expressing cell"), y = NULL) +
  fig18_theme

save_panel("Fig18e_UMIs_per_tumour.pdf", fig18e, "e")

fig18f <- ggplot(expressing, aes(gene_sct, internal_score, color = Type, fill = Type)) +
  geom_point(size = 0.15, alpha = 0.25, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.8, alpha = 0.2, show.legend = FALSE) +
  scale_color_manual(values = subtype_colors) +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  labs(x = bquote(italic(.(focus_gene))~"expression (SCT)"), y = "Internal AM/CM score") +
  fig18_theme

save_panel("Fig18f_LM_Plot.pdf", fig18f, "f")

fit_lm <- function(data, label) {
  fit <- summary(lm(internal_score ~ gene_sct * Type, data = data))
  am <- summary(lm(internal_score ~ gene_sct, data = data[data$Type == "AM", ]))$coefficients
  cm <- summary(lm(internal_score ~ gene_sct, data = data[data$Type == "CM", ]))$coefficients
  data.frame(cells = label, n = nrow(data),
             AM_slope = am["gene_sct", 1], AM_p = am["gene_sct", 4],
             CM_slope = cm["gene_sct", 1], CM_p = cm["gene_sct", 4],
             interaction = fit$coefficients["gene_sct:TypeCM", 1],
             interaction_p = fit$coefficients["gene_sct:TypeCM", 4],
             adj_R2 = fit$adj.r.squared, F = fit$fstatistic[[1]],
             df1 = fit$fstatistic[[2]], df2 = fit$fstatistic[[3]])
}
write.table(rbind(fit_lm(expressing, paste(focus_gene, "> 0")), fit_lm(cell_data, "all cells")),
            file.path(output_directory, "Fig18f_regression.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

per_tumour <- do.call(rbind, lapply(split(expressing, expressing$tumour), function(d) {
  co <- summary(lm(internal_score ~ gene_sct + log10_umi, data = d))$coefficients
  data.frame(tumour = d$tumour[1], Type = d$Type[1], n = nrow(d),
             slope = co["gene_sct", 1], p = co["gene_sct", 4])
}))
write.table(per_tumour, file.path(output_directory, "Fig18f_per_tumour_slopes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message(paste0("--- Single-Cell Processing Script Finished (", Sys.time(), ") ---"))
