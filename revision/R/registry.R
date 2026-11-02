# ---- Module registry: key -> script, reviewer questions, output folder ----

MODULES <- list(
  ancestry = list(
    script    = "R1-6_ancestry_cna_hailstorm_bra.R",
    questions = "R1-6",
    title     = "Ancestry vs CNA and hailstorm events, Firth logistic (Table R4)",
    packages  = c("tidyverse", "logistf", "readxl")),
  cnapp_bra = list(
    script    = "R3-1_R3-13_cnapp_hailstorm_braf_bra.R",
    questions = c("R3-1", "R3-13"),
    title     = "CNApp scores x BRAF, Brazilian cohort (Figure R15)",
    packages  = c("dplyr", "ggplot2", "readr", "tidyr", "purrr", "stringr",
                  "tibble", "patchwork", "readxl")),
  atm_chek1_mex = list(
    script    = "R3-1_atm_chek1_hailstorm_mex.R",
    questions = "R3-1",
    title     = "ATM/CHEK1 loss vs hailstorms, Mexican cohort (Table R9)",
    packages  = c("readr", "dplyr", "stringr")),
  cnapp_mex = list(
    script    = "R3-1_R3-13_cnapp_hailstorm_braf_mex.R",
    questions = c("R3-1", "R3-13"),
    title     = "Hailstorm frequency, BRAF tiles, CNApp x BRAF, Mexican cohort (Figures R6, R13, R14)",
    packages  = c("tidyverse", "patchwork")),
  crispr_cn = list(
    script    = "R3-12_R1-13_crispr_dependency_vs_cn_bra_depmap.R",
    questions = c("R3-12", "R1-13"),
    title     = "CRISPR dependency vs copy number (Figures R3, R8, R9, R10, R16-R18)",
    packages  = c("tidyverse", "data.table", "patchwork", "ggrepel", "scales")),
  cn_profile = list(
    script    = "R3-13_cn_profile_drivers_hailstorm_depmap.R",
    questions = "R3-13",
    title     = "DepMap line CN profile at AM drivers (Figure R12)",
    packages  = c("tidyverse", "data.table", "ggrepel", "scales",
                  "GenomicFeatures", "GenomicRanges",
                  "TxDb.Hsapiens.UCSC.hg38.knownGene"))
)

module_id <- function(key) sub("\\.R$", "", MODULES[[key]]$script)

all_packages <- function() sort(unique(unlist(lapply(MODULES, `[[`, "packages"))))
