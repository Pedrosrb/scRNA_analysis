# =============================================================================
# Ancestry vs CNA and hailstorms, Firth logistic regression
# Run:      Rscript run.R ancestry
# Outputs:
#   tables/cna_hailstorm_ancestry.tsv    Table R4
#   tables/patient_features.tsv, tables/ancestry_scaling.tsv
# Inputs:   config/params.R -> PARAMS$modules$ancestry$inputs
#   genotyping    SUPPLEMENTARY_TABLE_GENOTYPING: Id lab, Sample id, Sample type, K1-K5
#   clinical      Supplementary Table S1: Sample_ID, Sex, Age
#   histology     Supplementary Table S2: Sample_ID, Tumor_type, Clinical_Stage
#   samples       Supplementary Table S13
#   gistic        GISTIC all_lesions.conf
#   hailstorms    Detected hailstorms: Sample_ID, Sanger_DNA_ID, Arm, Filter
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("ancestry")

suppressPackageStartupMessages({
  library(tidyverse)
  library(logistf)
})

config <- c(M$params[setdiff(names(M$params), "inputs")],
            M$inputs,
            list(out_dir = M$tables))

load_ancestry <- function(cfg) {
  a <- bra_ancestry(cfg$genotyping)
  as_tibble(a) |> transmute(patient_id = patient, !!!syms(cfg$ancestry_map))
}

load_covariates <- function(cfg) {
  cl <- bra_clinical(cfg$clinical, cfg$histology, cfg$early_stages, cfg$late_stages,
                     cfg$stage_from_metastasis)
  as_tibble(cl) |> transmute(patient_id = patient, sex, age, stage_group, stage_source)
}

load_features <- function(cfg) {
  smp <- bra_samples(cfg$samples, cfg$min_cellularity)
  an  <- smp[smp$analysable, , drop = FALSE]
  flow_add("analysable patient tumours (S13)", nrow(an), "samples",
           paste0("cellularity >= ", cfg$min_cellularity))
  flow_add("analysable patients (S13)", n_distinct(an$patient))

  g <- bra_gistic_long(cfg$gistic, cfg$gistic_q_max, cfg$gistic_min_level)
  out <- setdiff(unique(g$sample_id), an$sample_id)
  if (length(out))
    message("[gistic] tumour(s) not analysable in S13, dropped: ", paste(out, collapse = ", "))
  g <- g[g$sample_id %in% an$sample_id, , drop = FALSE]
  g$patient_id <- an$patient[match(g$sample_id, an$sample_id)]
  if (any(tapply(g$sample_id, g$patient_id, n_distinct) > 1))
    stop("GISTIC has more than one tumour for a patient", call. = FALSE)
  flow_add("GISTIC patients", n_distinct(g$patient_id), note = paste0(
    n_distinct(g$feature), " peaks with q < ", cfg$gistic_q_max))
  cna <- as_tibble(g) |>
    group_by(patient_id, feature) |>
    summarise(present = as.integer(any(present == 1)), .groups = "drop") |>
    mutate(type = "CNA")

  calls <- bra_hailstorm_calls(cfg$hailstorms, cfg$hailstorm_pass)
  calls <- calls[calls$am_id == am_lesion(calls$am_id), , drop = FALSE]
  arms  <- sort(unique(calls$arm[calls$sample_id %in% an$sample_id]))
  st <- hailstorm_status(calls, an$sample_id, arms)
  st$patient_id <- an$patient[match(st$sample_id, an$sample_id)]
  hs <- as_tibble(st) |>
    group_by(patient_id, feature = arm) |>
    summarise(present = as.integer(any(present == 1)), .groups = "drop") |>
    mutate(type = "Hailstorm")
  flow_add("hailstorm-evaluated patients", n_distinct(hs$patient_id), note = paste0(
    n_distinct(hs$patient_id[hs$present == 1]), " with >= 1 hailstorm"))

  bind_rows(cna, hs)
}

build_model_data <- function(df, ancestry, cfg, scaling) {
  sc <- scaling[scaling$ancestry == ancestry, ]
  df |>
    mutate(ancestry_scaled = (.data[[ancestry]] - sc$mean) / sc$sd) |>
    select(present, ancestry_scaled, all_of(cfg$covariates)) |>
    drop_na()
}

fit_firth <- function(model_df, covariates) {
  active <- covariates[map_lgl(covariates, \(cv) n_distinct(model_df[[cv]]) > 1)]
  logistf(
    reformulate(c("ancestry_scaled", active), response = "present"),
    data      = model_df,
    pl        = TRUE,
    control   = logistf.control(maxit = M$params$firth_maxit),
    plcontrol = logistpl.control(maxit = M$params$firth_maxit)
  )
}

test_feature <- function(df, ancestry, cfg, scaling) {
  model_df <- build_model_data(df, ancestry, cfg, scaling)
  n_pos <- sum(model_df$present)
  if (n_pos < cfg$min_events || nrow(model_df) - n_pos < cfg$min_events) return(NULL)
  fit <- fit_firth(model_df, cfg$covariates)
  tibble(
    n_pos      = n_pos,
    n          = nrow(model_df),
    odds_ratio = exp(coef(fit)[["ancestry_scaled"]]),
    ci_lower   = exp(fit$ci.lower[["ancestry_scaled"]]),
    ci_upper   = exp(fit$ci.upper[["ancestry_scaled"]]),
    p_value    = fit$prob[["ancestry_scaled"]]
  )
}

ancestry_scaling <- function(anc_df, cfg) {
  comps <- names(cfg$ancestry_map)
  tibble(ancestry = comps,
         mean = map_dbl(comps, \(a) mean(anc_df[[a]])),
         sd   = map_dbl(comps, \(a) sd(anc_df[[a]])),
         n    = nrow(anc_df))
}

run_association <- function(features, ancestry, covariates, cfg, scaling) {
  complete <- covariates |> drop_na(all_of(cfg$covariates)) |> pull(patient_id)
  in_model <- intersect(intersect(unique(features$patient_id), ancestry$patient_id), complete)
  flow_add("with features, ancestry and complete covariates", length(in_model))

  features |>
    inner_join(ancestry, by = "patient_id") |>
    left_join(covariates, by = "patient_id") |>
    nest(data = -c(type, feature)) |>
    expand_grid(ancestry = names(cfg$ancestry_map)) |>
    mutate(result = map2(data, ancestry, \(d, a) test_feature(d, a, cfg, scaling))) |>
    select(-data) |>
    unnest(result) |>
    group_by(type) |>
    mutate(p_fdr = p.adjust(p_value, method = "BH")) |>
    ungroup()
}

summarise_ancestry <- function(ancestry, cfg) {
  ancestry |>
    pivot_longer(-patient_id, names_to = "ancestry", values_to = "proportion") |>
    group_by(ancestry) |>
    summarise(
      N      = n(),
      Mean   = mean(proportion),
      SD     = sd(proportion),
      Median = median(proportion),
      Q1     = quantile(proportion, 0.25),
      Q3     = quantile(proportion, 0.75),
      Min    = min(proportion),
      Max    = max(proportion),
      .groups = "drop"
    ) |>
    arrange(desc(Mean))
}

format_proportions <- function(summary, cfg) {
  summary |>
    transmute(
      Ancestry  = toupper(ancestry),
      Component = cfg$ancestry_map[ancestry],
      N,
      across(Mean:Max, \(x) round(x, 3))
    )
}

write_tables <- function(tables, cfg) {
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  iwalk(tables, \(tbl, name) write_tsv(tbl, file.path(cfg$out_dir, paste0(name, ".tsv"))))
  tables
}

format_associations <- function(results) {
  results |>
    arrange(type, feature, ancestry) |>
    transmute(
      Type             = type,
      Feature          = feature,
      Ancestry         = toupper(ancestry),
      `Positive/N`     = paste0(n_pos, "/", n),
      `Odds Ratio`     = signif(odds_ratio, 3),
      `Lower 95% CI`   = signif(ci_lower, 3),
      `Upper 95% CI`   = signif(ci_upper, 3),
      `P value`        = signif(p_value, 2),
      `P value FDR`    = signif(p_fdr, 2)
    )
}

main <- function(cfg) {
  ancestry <- load_ancestry(cfg)
  flow_add("genotyped (blood, ADMIXTURE K5, QC pass)", nrow(ancestry))
  covariates <- load_covariates(cfg)
  flow_add("in clinical tables (S1/S2)", nrow(covariates), note = paste0(
    sum(covariates$stage_source %in% "nodal_or_in_transit_metastasis"),
    " staged from nodal/in-transit metastasis"))
  features <- load_features(cfg)

  summary <- summarise_ancestry(ancestry, cfg)
  scaling <- ancestry_scaling(ancestry, cfg)
  associations <- run_association(features, ancestry, covariates, cfg, scaling) |>
    format_associations()
  write_tables(
    list(
      ancestry_proportions   = format_proportions(summary, cfg),
      ancestry_scaling       = scaling,
      patient_features       = features |>
        pivot_wider(id_cols = patient_id, names_from = c(type, feature),
                    values_from = present) |>
        left_join(covariates, by = "patient_id") |>
        left_join(ancestry, by = "patient_id"),
      cna_hailstorm_ancestry = associations
    ),
    cfg
  )
}

main(config)
module_finish(M)
