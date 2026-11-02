# All input paths and analysis choices. Each module reads PARAMS$modules$<key>.

P <- list()

# ---- ROOTS ----
# Relative to revision/ (see "Inputs" in README.md)
P$roots <- list(
  data    = file.path(Sys.getenv("AMREV_ROOT"), "data"),
  results = file.path(Sys.getenv("AMREV_ROOT"), "results")
)

# ---- SHARED INPUTS ----
P$depmap <- list(
  release = "DepMap Public 26Q1",
  dir     = file.path(P$roots$data, "depmap_26Q1")
)
P$depmap$files <- list(
  model       = file.path(P$depmap$dir, "Model.csv"),
  crispr      = file.path(P$depmap$dir, "CRISPRGeneEffect.csv"),
  cn_gene     = file.path(P$depmap$dir, "OmicsCNGeneWGS.csv"),
  cn_segments = file.path(P$depmap$dir, "OmicsCNSegmentsWGS.csv"),
  signatures  = file.path(P$depmap$dir, "OmicsGlobalSignatures.csv")
)

P$bayes_combined <- file.path(P$roots$data, "bra", "bayes_combined_AM_depmap.csv")

P$cnvkit_dir <- file.path(P$roots$data, "cnvkit")
P$cns <- c(
  AM007 = file.path(P$cnvkit_dir, "PD53336l.sample.dupmarked.call.median_centred.weight_ge_30.0.cns"),
  AM016 = file.path(P$cnvkit_dir, "PD53344l.sample.dupmarked.call.median_centred.weight_ge_30.0.cns")
)
cns_inputs <- structure(as.list(P$cns), names = paste0("cns_", names(P$cns)))

P$bra <- list(dir = file.path(P$roots$data, "bra"))
find_in <- function(dir, pattern) {
  hits <- if (dir.exists(dir))
    list.files(dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE) else character(0)
  hits <- hits[!grepl("^~\\$", basename(hits))]
  if (length(hits) > 1) {
    xl <- hits[grepl("\\.xlsx$", hits, ignore.case = TRUE)]
    if (length(xl)) hits <- xl
  }
  if (length(hits) == 1) return(hits)
  file.path(dir, sprintf("<%s: %s>", if (length(hits)) "AMBIGUOUS" else "NOT FOUND",
                         if (length(hits)) paste(basename(hits), collapse = " | ") else pattern))
}
P$bra$files <- list(
  samples    = find_in(P$bra$dir, "Supplementary[ _-]*Table[ _-]*S13"),
  clinical   = find_in(P$bra$dir, "Supplementary[ _-]*Table[ _-]*S1([^0-9]|$)"),
  histology  = find_in(P$bra$dir, "Supplementary[ _-]*Table[ _-]*S2([^0-9]|$)"),
  genotyping = find_in(P$bra$dir, "GENOTYPING"),
  braf       = find_in(P$bra$dir, "BRAF[ _-]*SEQUENCING"),
  cnapp      = find_in(P$bra$dir, "Cnapp[ _-]*Scores"),
  hailstorms = find_in(P$bra$dir, "Detected[ _-]*hailstorms"),
  gistic     = find_in(P$bra$dir, "all_lesions\\.conf")
)

P$bra$min_cellularity <- 0.35
P$bra$hailstorm_pass  <- c("Pass", "Manual_Pass")

P$mex <- list(dir = file.path(P$roots$data, "mex"))
P$mex$files <- list(
  supp = file.path(P$mex$dir, "Supplementary_Table_1.csv"),
  hailstorms = file.path(P$mex$dir, "All_samples_hailstorms_mexican_cohort.tsv"),
  detection_log = file.path(P$mex$dir, "detection_run_log.tsv")
)

# ---- Copy-number boundaries ----
P$cn <- list(
  harmonised   = list(deep_del = -1.10, loss = -0.40, gain = 0.585, amp = 1.00),
  floor_linear = 0.05
)

P$braf_definition <- "v600_only"

P$cell_lines <- list(
  am_bra    = c(AM007 = "PD53336", AM016 = "PD53344"),
  am_depmap = c(WM3211 = "ACH-002508", M040416 = "ACH-002510",
                MM160113 = "ACH-002512")
)

P$run <- list(
  hash_inputs = TRUE
)


# ---- MODULE PARAMETERS ----
P$modules <- list()

P$modules$ancestry <- list(
  seed = 1L,
  inputs = P$bra$files[c("genotyping", "clinical", "histology", "samples",
                         "gistic", "hailstorms")],
  min_cellularity = P$bra$min_cellularity,
  hailstorm_pass  = P$bra$hailstorm_pass,
  ancestry_map = c(eur = "K1", amr = "K2", sas = "K3", afr = "K4", eas = "K5"),
  min_events   = 3,
  firth_maxit  = 1000,
  covariates   = c("sex", "age", "stage_group"),
  early_stages = c("Ib", "IIa", "IIb", "IIc"),
  late_stages  = c("III", "IV"),
  stage_from_metastasis = TRUE,
  gistic_q_max     = 0.1,
  gistic_min_level = 1
)

P$modules$cnapp_bra <- list(
  seed = 1L,
  inputs = P$bra$files[c("samples", "cnapp", "hailstorms", "braf")],
  min_cellularity  = P$bra$min_cellularity,
  hailstorm_pass   = P$bra$hailstorm_pass,
  analysis_unit    = "patient",
  collapse_fun     = "median",
  min_group_n      = 3,
  min_arm_cases    = 3,
  braf_source      = "consensus",
  braf_definition  = P$braf_definition,
  braf_v600_pattern = "^V600",
  freq_denominator = "profiled",
  arms_to_exclude  = character(0),
  score_names      = c("BCS", "FCS", "GCS")
)

P$modules$cnapp_mex <- list(
  seed = 1L,
  inputs = P$mex$files[c("supp", "hailstorms", "detection_log")],
  braf_definition        = P$braf_definition,
  braf_v600_pattern      = "V600",
  hs_pass_values         = c("Yes"),
  analysis_unit          = "patient",
  collapse_fun           = "median",
  min_group_n            = 3,
  freq_denominator       = "detected",
  exclude_no_breakpoints = TRUE,
  min_arm_cases          = 3,
  arms_to_exclude        = c("7p")
)

P$modules$atm_chek1_mex <- list(
  seed = 1L,
  inputs = P$mex$files[c("supp", "hailstorms", "detection_log")],
  sequenza_dir = file.path(P$mex$dir, "sequenza"),
  hs_pass_values = c("Yes"),
  genes = list(
    ATM   = list(chr = "11", start = 108222484, end = 108369102),
    CHEK1 = list(chr = "11", start = 125625135, end = 125676255)),
  arms = c("5p", "11q", "22q"),
  autosomes = as.character(1:22)
)

P$modules$crispr_cn <- list(
  seed = 28L,
  inputs = c(P$depmap$files[c("model", "crispr", "cn_gene")],
             list(bayes_combined = P$bayes_combined), cns_inputs),
  pd_map          = c(PD53336 = "AM007", PD53344 = "AM016"),
  bf_cols_inhouse = c(AM007 = "AM007_D21", AM016 = "AM016_D21"),
  am_bra          = names(P$cell_lines$am_bra),
  am_depmap       = P$cell_lines$am_depmap,
  hit_genes = c("MITF", "SOX10", "CDK4", "MDM2", "CCND1", "CCNB1", "ARF4",
                "SF3B4", "CRKL", "MAPK1", "THOC2", "DPAGT1", "FERMT2",
                "SBNO1", "USP7"),
  focus_genes   = c("CRKL", "MITF", "SOX10"),
  cn_ratio      = P$cn$harmonised,
  cn_harmonised = P$cn$harmonised,
  cn_scheme     = "harmonised",
  min_models_per_metric = 10,
  corr_scopes   = c("pancancer", "melanoma"),
  melanoma_groups = c("AM_DepMap", "AM_DepMap_other", "Cutaneous_Melanoma",
                      "Uveal_Melanoma", "Melanoma_NOS"),
  detail_gene   = "CRKL",
  common_essential_chronos = -1,
  exact_spearman_max_n = 7,
  depmap_id_prefix = "ACH",
  am_bra_annotation = list(subtype = "Acral Melanoma", disease = "Melanoma", lineage = "Skin"),
  cns_non_genes = c("", "-", ".", "BACKGROUND", "ANTITARGET"),
  regression_check = list(
    gene = "SF3B4",
    expected = c(AM007 = 0.90, AM016 = 0.37, WM3211 = 0.18,
                 M040416 = 0.11, MM160113 = 0.85),
    tolerance = 0.05),
  complete_cases_max_models = 10,
  min_models_per_gene       = 20,
  min_frac_models           = 0.50,
  fdr_method     = "BH",
  fdr_alpha      = 0.05,
  fig_dpi = 300
)

P$modules$cn_profile <- list(
  seed = 1L,
  inputs = P$depmap$files,
  model = "M040416",
  gene_coords = "txdb",
  drivers = c("TERT", "BRAF", "CCND1", "GAB2", "CDK4", "MDM2",
              "CRKL", "EP300", "SOX10", "MITF"),
  cn = P$cn,
  label_only_gained   = TRUE,
  default_ploidy = 2,
  chromosomes     = paste0("chr", c(1:22, "X")),
  chromosomes_all = paste0("chr", c(1:22, "X", "Y"))
)

PARAMS <- P
rm(P, cns_inputs, find_in)
