# Comparative Meta-analysis of Acral and Cutaneous Melanoma at Single-Cell Resolution

This repository is part of the [Acral Melanoma PDX Models from Latin America](https://github.com/team113sanger/Acral_Melanoma_PDX_models_from_Latin_America) project. It contains the code and documentation for a meta-analysis of single-cell RNA sequencing (scRNA-seq) data from primary Acral Melanoma (AM) and Cutaneous Melanoma (CM) samples.

---

## Overview

This workflow processes raw scRNA-seq data, integrates individual samples into a unified dataset, and performs comparative analyses to identify transcriptional differences between AM and CM. The analysis includes cell type annotation, clustering, and the generation of publication-ready figures.

---

## Analysis Workflow

The analysis is structured as a modular pipeline, with each stage implemented in dedicated scripts located in the `src/` directory.

### 1. SRA Data Processing (`01_process_sra_data.sh`)

- Downloads raw sequencing data using a list of SRA accession numbers (`sra_list.txt`).
- Converts `.sra` files to gzipped FASTQ format using `fastq-dump`.
- Runs `cellranger count` to perform alignment, barcode assignment, and UMI counting, generating feature-barcode matrices for each sample.

### 2. Seurat Object Preparation (`02_main.R`, `03_functions.R`)

- **Loading**: `create_seurat_list()` reads Cell Ranger outputs and creates a list of Seurat objects, computing initial QC metrics.
- **Filtering**: `filter_seurat_list()` removes low-quality cells (e.g., cells with high mitochondrial content, `percent.mt > 5%`) and identifies doublets using `scDblFinder`.

### 3. Integration and Clustering (`02_main.R`, `03_functions.R`)

- **Integration**: `integrate_seurat_list()` merges filtered Seurat objects and uses `SCTransform` and `RunHarmony` to normalize and correct for batch effects across samples.
- **Clustering**: `cluster_seurat_obj()` performs cell clustering over multiple resolutions (0.1–0.5).

### 4. Annotation and Subsetting (`02_main.R`, `03_functions.R`)

- **Dimensionality Reduction**: Executes `RunUMAP` on Harmony-corrected dimensions for 2D visualization.
- **Annotation**: `annotate_singler()` labels cell types using the `BlueprintEncodeData` reference dataset via `SingleR`.
- **Subsetting**: Isolates the "Melanocyte" population for downstream differential analysis.

### 5. Figure Generation and Analysis (`02_main.R`)

Generates and saves plots in the `results/single_cell_analysis/` directory, corresponding to manuscript figures:

- **Fig 7A**: UMAP of all annotated cell types.
- **Fig 7B**: Dot plot of key marker gene expression.
- **Fig 7C**: UMAP of melanocytes, colored by condition (AM vs. CM).
- **Fig 7D**: Volcano plot of differentially expressed genes (DEGs) between AM and CM melanocytes using `MAST`.
- **Fig 7E**: Violin plot showing AM/CM signature scores per cell.
- **Fig 7F**: Linear regression plot correlating *SF3B4* expression with AM/CM signature scores.

---

## Software and Dependencies

- R version: `4.4.3`
- CellRanger: `8.0.1`
- Dependencies managed via [`renv`](https://rstudio.github.io/renv/)
- Key R packages:
  - `Seurat`
  - `harmony`
  - `SingleR`
  - `scDblFinder`
  - `here`
  - `renv`

To reproduce the analysis environment, use the `renv.lock` file in this repository and refer to the [renv documentation](https://rstudio.github.io/renv/).

---

## How to Run

### 1. Configuration

Edit `src/04_parameters.R` to set all file paths and analysis parameters. Ensure the `sra_list.txt` file contains valid SRA accession numbers.

---

### 2. SRA Data Processing

Run the shell script to download and process raw scRNA-seq data.

```bash
bash 01_process_sra_data.sh \
    -s path/to/your/sra_list.txt \
    -t path/to/your/transcriptome \
    -o path/to/your/output_directory
```

This step will generate Cell Ranger outputs in the specified directory.

---

### 3. Run the Core Analysis

Once the count matrices are generated, run the main R script:

```R
# In an R session
source("src/02_main.R")
```

---

### 4. Outputs

All intermediate and final outputs (QC reports, plots, differential analysis results) will be saved in the directory specified by `output_directory` in `04_parameters.R`.

---

## Contact

If you have any questions, suggestions, or feedback about this repository, please contact:

- **Pedro Sodré**: [pedrosodrerb@gmail.com](mailto\:pedrosodrerb@gmail.com)

