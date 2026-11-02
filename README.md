# Comparative Meta-analysis of Acral and Cutaneous Melanoma at Single-Cell Resolution

This repository is part of the [Acral Melanoma PDX Models from Latin America](https://github.com/team113sanger/Acral_Melanoma_PDX_models_from_Latin_America) project. It contains the code and documentation for a meta-analysis of single-cell RNA sequencing (scRNA-seq) data from primary Acral Melanoma (AM) and Cutaneous Melanoma (CM) samples.

---

## Data

Public scRNA-seq data of 6 AM and 3 CM primary tumours from [GSE215121](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215121) (Zhang et al., *Nat Commun* 2022). The runs are listed in `data/sra_list.txt`; their tumour labels (AM1–AM6, CM1–CM3) are in `tumour_labels`, in `src/04_parameters.R`.

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

Generates and saves the panels of **Extended Data Figure 18** in the directory set by `output_directory`:

- **Fig 18a**: UMAP of all annotated cell types.
- **Fig 18b**: Dot plot of *CRKL*, *DPAGT1*, *SF3B4*, *THOC2*, *SOX10*, *MITF* and *TYRP1* across cell types.
- **Fig 18c**: UMAP of melanocytes, colored by subtype (AM vs. CM).
- **Fig 18d**: Volcano plot of differentially expressed genes between AM and CM melanocytes using `MAST`.
- **Fig 18e**: *SF3B4* UMIs per expressing melanocyte in each tumour (mean, 95% CI).
- **Fig 18f**: *SF3B4* expression (SCT) vs. internal AM/CM score in melanocytes with detectable *SF3B4*, with one linear fit per subtype.

`results/Extended_Data_Figure_18.pdf` is the assembled figure.

### 6. Revision analyses (`revision/`)

Analyses of the response to the referees, as an independent modular pipeline with its own runner and parameters. See [`revision/README.md`](revision/README.md).

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
  - `MAST`
  - `here`
  - `renv`

To reproduce the analysis environment, run `Rscript src/00_setup_env.R` from the repository root. It installs the versions in `renv.lock`; see the [renv documentation](https://rstudio.github.io/renv/).

---

## How to Run

### 1. Configuration

All paths are relative to the repository root. `src/04_parameters.R` holds the paths and analysis parameters. Ensure the `sra_list.txt` file contains valid SRA accession numbers.

---

### 2. SRA Data Processing

Run the shell script to download and process raw scRNA-seq data.

```bash
bash 01_process_sra_data.sh \
    -s path/to/your/sra_list.txt \
    -t path/to/your/transcriptome \
    -o path/to/your/output_directory
```

This step will generate Cell Ranger outputs in the specified directory. `02_main.R` expects one folder per sample, named by its SRA accession, directly under `cellranger_directory` (`<cellranger_directory>/<SRR>/outs/filtered_feature_bc_matrix`). The accessions must match `tumour_labels` in `04_parameters.R`.

---

### 3. Run the Core Analysis

Once the count matrices are generated, run the main R script:

```R
# In an R session
source("src/02_main.R")
```

---

### 4. Outputs

All intermediate and final outputs (QC reports, plots, tables) will be saved in the directory specified by `output_directory` in `04_parameters.R`.

---

## Versions

See [`CHANGELOG.md`](CHANGELOG.md). Version 1.0 (tag `1.0.0`) is the code of the first submission; version 2.0 (tag `2.0.0`) is the code of the revision.

---

## Contact

If you have any questions, suggestions, or feedback about this repository, please contact:

- **Pedro Sodré**: [pedro.barros@ensino.inca.gov.br](mailto\:pedro.barros@ensino.inca.gov.br)

