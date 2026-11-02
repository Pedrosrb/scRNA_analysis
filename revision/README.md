# Revision analyses

Analyses by Pedro Sodré in the response to the referees. Each script lists its inputs, required columns and outputs at the top; all parameters are in `config/params.R`.

## Letter items

| Item | Module | Output |
|---|---|---|
| Tables R2, R4 | `ancestry` | `tables/` |
| Figures R3, R8 | `crispr_cn` | `figures/harmonised/CN_call_per_cell_line.pdf` |
| Figure R9 | `crispr_cn` | `figures/correlation/hit_correlation_chronos.pdf` |
| Figure R10 | `crispr_cn` | `figures/hit_waterfalls/hit_genes_waterfall_chronos.pdf` |
| Figure R16 | `crispr_cn` | `figures/harmonised/CN_call_CRKL.pdf` |
| Figure R17 | `crispr_cn` | `figures/correlation/CRKL_cn_vs_chronos.pdf` |
| Figure R18 | `crispr_cn` | `figures/hit_waterfalls/CRKL_waterfall_chronos.pdf` |
| Figure R6 | `cnapp_mex` | `figures/mex_hailstorm_frequency_by_arm.pdf` |
| Figure R13 | `cnapp_mex` | `figures/mex_hailstorm_braf_v600_tiles.pdf` |
| Figure R14, Table R15 | `cnapp_mex` | `figures/mex_cnapp_by_braf.pdf` |
| Figure R15, Table R15 | `cnapp_bra` | `figures/bra_cnapp_by_braf.pdf` |
| Table R9, MEX rows | `atm_chek1_mex` | `tables/atm_chek1_3arms_mex_results.tsv` |
| Figure R12 | `cn_profile` | `figures/M040416_genome.pdf` |

Script names keep the tags used during the analysis, so they may differ from the numbering in the final letter.

## Run

From this folder:

```bash
Rscript R/setup_env.R           # install the versions in renv.lock
Rscript run.R                   # list modules
Rscript run.R --preflight all   # check inputs and packages
Rscript run.R --jobs=3 all      # run everything
Rscript run.R crispr_cn         # run one module
```

Results go to `results/<module>/`, with the parameters and input checksums of each run in `provenance/`.

## Input data

```
data/
  depmap_26Q1/   DepMap Public 26Q1: Model.csv, CRISPRGeneEffect.csv, OmicsCNGeneWGS.csv,
                 OmicsCNSegmentsWGS.csv, OmicsGlobalSignatures.csv
  bra/           Supplementary Tables S1, S2, S13; genotyping; GISTIC all_lesions.conf;
                 detected hailstorms; CNApp scores; BRAF sequencing;
                 bayes_combined_AM_depmap.csv
  mex/           Supplementary_Table_1.csv, All_samples_hailstorms_mexican_cohort.tsv,
                 detection_run_log.tsv, sequenza/<Sample>/*segments.txt
  cnvkit/        PD53336l and PD53344l .cns files (AM007, AM016)
```
