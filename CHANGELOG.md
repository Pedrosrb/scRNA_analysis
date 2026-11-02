# Changelog

## 2.0 (November 2026), manuscript revision

### Single-cell pipeline (`src/`)
- Plots now match Extended Data Figure 18, panels a–f. Output names: `Fig18a` to `Fig18f`.
- Panel e is new: *SF3B4* UMIs per expressing melanocyte, per tumour (mean, 95% CI).
- Panel f now uses melanocytes with detectable *SF3B4*, labelled "Internal AM/CM score" and SCT background.
- Removed: the violin plot of the score (old panel E).
- `SOX2` and `SOX9` are no longer highlighted; `THOC2` is.
- New tables: markers, score genes, per-cell scores, panel e values, panel f regression (correct degrees of freedom) and per-tumour slopes.

### Fixes (`src/`)
- `seurat_obj$orig.ident` was undefined in Step 3.
- The `ggsave()` of an undefined `lm_fit` was removed; the regression is written to `Fig18f_regression.tsv`.
- `merge()` receives `add.cell.ids`, `RunHarmony()` receives `reduction.use = "pca"` and `PrepSCTFindMarkers()` runs before `FindMarkers()`.

### New
- `revision/`: analyses of the response to the referees (see `revision/README.md`), using DepMap Public 26Q1.

### Code style
- Analysis values moved to `src/04_parameters.R` and `revision/config/params.R`.


## 1.0 (July 2025), first submission
