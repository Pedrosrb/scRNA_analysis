# =============================================================================
# CNApp scores x BRAF, Brazilian cohort
# Run:      Rscript run.R cnapp_bra
# Outputs:
#   figures/bra_cnapp_by_braf.pdf        Figure R15, Table R15
#   tables/bra_braf_hailstorm_arms.tsv    Table R15, hailstorm row
# Inputs:   config/params.R -> PARAMS$modules$cnapp_bra$inputs
#   samples       Supplementary Table S13: Sample, Lab_id, Sample_type, Cellularity, Ploidy, Will we use it?
#   cnapp         Cnapp_Scores: AM_id, Sanger_id, BCS, FCS, GCS
#   hailstorms    Detected hailstorms: Sample_ID, Sanger_DNA_ID, Arm, Filter
#   braf          BRAF sequencing: Samples, Sanger call, WES call
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("cnapp_bra")

if (requireNamespace("tidyverse", quietly = TRUE)) {
  suppressPackageStartupMessages(library(tidyverse))
} else {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(readr); library(tidyr)
    library(purrr); library(stringr); library(tibble)
  })
}

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

if (!requireNamespace("patchwork", quietly = TRUE))
  stop("patchwork is required to assemble the panel grids: install.packages('patchwork')")

p <- M$params
samples_path <- M$inputs$samples
hs_path      <- M$inputs$hailstorms
braf_path    <- M$inputs$braf
cnapp_path   <- M$inputs$cnapp

analysis_unit    <- p$analysis_unit
collapse_fun     <- match.fun(p$collapse_fun)
min_group_n      <- p$min_group_n
min_arm_cases    <- p$min_arm_cases
braf_source      <- p$braf_source
braf_definition  <- p$braf_definition
freq_denominator <- p$freq_denominator
arms_to_exclude  <- p$arms_to_exclude
score_names      <- p$score_names

bar_fill <- "#E8452F"
braf_pal <- c("Non-mutant" = "grey75", "Mutant" = bar_fill)
hs_pal   <- c("Hailstorm-" = "grey75", "Hailstorm+" = bar_fill)

panel_w_arm    <- 3.2
panel_h_arm    <- 2.6
panel_w_global <- 3.6
panel_h_global <- 3.8

out_dir <- M$dir
out_tab <- M$tables
out_fig <- M$figures

# ---- HELPERS  (shared with the MX module) ----

section <- function(txt) cat("\n", strrep("=", 72), "\n", txt, "\n",
                             strrep("=", 72), "\n", sep = "")

fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.4f", p)

save_fig <- function(plot, stem, width, height) {
  ggsave(file.path(out_fig, paste0(stem, ".pdf")), plot, width = width, height = height)
  cat(sprintf("  figure -> %s.pdf\n", stem))
}

norm_braf <- function(x) {
  x <- str_squish(str_replace_all(as.character(x), "[\r\n]", ""))
  case_when(
    is.na(x) | x == ""                                   ~ NA_character_,
    str_detect(x, regex("^(wt|wild[- _]?type)$", TRUE))  ~ "WT",
    str_detect(x, regex("^not[ _]sequenced$", TRUE))     ~ NA_character_,
    str_detect(x, "^[A-Z]\\d+[A-Z*]?(fs)?$")             ~ str_to_upper(x),
    TRUE                                                 ~ NA_character_)
}

compare_groups <- function(values, group, feature, group_levels, min_n = min_group_n) {
  r <- wilcox_two_group(values, group, feature, group_levels, min_n)
  if (is.null(r)) NULL else as_tibble(r)
}

make_arm_panel <- function(d, sc, a, lim, p_lab) {
  ggplot(d, aes(x = arm_hs, y = value, fill = arm_hs)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, colour = "black",
                 linewidth = 0.4, alpha = 0.9) +
    geom_jitter_fixed(width = 0.16, size = 1.4, shape = 21,
                fill = "white", colour = "black", stroke = 0.35, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.1,
                 fill = "white", colour = "black", stroke = 0.4) +
    scale_fill_manual(values = c(No = "grey75", Yes = bar_fill), guide = "none") +
    coord_cartesian(ylim = c(lim$ymin, lim$ymax)) +
    labs(title = if (length(p_lab) == 1 && !is.na(p_lab))
      paste0("p.adj = ", signif(p_lab, 4)) else "",
      x = paste("Hailstorm", a), y = sc) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 9.5, hjust = 0.5, colour = "black"),
          axis.title.x = element_text(face = "bold", size = 10.5, colour = "black"),
          axis.title.y = element_text(face = "bold", size = 10.5, colour = "black"),
          axis.text = element_text(colour = "black", size = 9),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          plot.margin = margin(4, 6, 4, 4))
}

make_box <- function(df, value_col, group_col, title_str, y_lab,
                     p_value = NA_real_, palette) {
  ggplot(df, aes(x = .data[[group_col]], y = .data[[value_col]],
                 fill = .data[[group_col]])) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = "black",
                 linewidth = 0.4, alpha = 0.85) +
    geom_jitter_fixed(width = 0.14, size = 1.8, shape = 21,
                fill = "white", colour = "black", stroke = 0.4) +
    scale_fill_manual(values = palette, guide = "none") +
    labs(title = if (is.na(p_value)) NULL else paste0("p.adj = ", signif(p_value, 4)),
         x = NULL, y = y_lab) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 9.5, hjust = 0.5, colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(face = "bold"))
}

# ---- LOAD ----
section("LOADING")

# ---- analysable patient tumours (S13) ----
smp <- bra_samples(samples_path, p$min_cellularity)
an  <- smp[smp$analysable, , drop = FALSE]
flow_add("analysable patient tumours (S13)", nrow(an), "samples",
         paste0("Patient, cellularity >= ", p$min_cellularity))
flow_add("analysable patients (S13)", n_distinct(an$patient))
write_tsv(as_tibble(smp), file.path(out_tab, "bra_sample_inclusion.tsv"))

# ---- CNApp scores ----
cn_raw <- bra_cnapp(cnapp_path, score_names)
not_an <- setdiff(cn_raw$sample_id[cn_raw$am_id == am_lesion(cn_raw$am_id)], an$sample_id)
if (length(not_an))
  cat(sprintf("CNApp patient tumours not analysable in S13 (excluded): %s\n",
              paste(not_an, collapse = ", ")))
cnapp <- as_tibble(cn_raw) |>
  filter(sample_id %in% an$sample_id) |>
  transmute(Sample = sample_id, Patient = an$patient[match(sample_id, an$sample_id)],
            across(all_of(score_names)))
cat(sprintf("CNApp: %d samples | %d patients\n",
            nrow(cnapp), n_distinct(cnapp$Patient)))
flow_add("with CNApp scores", nrow(cnapp), "samples")
flow_add("with CNApp scores", n_distinct(cnapp$Patient), "patients")
calls <- bra_hailstorm_calls(hs_path, p$hailstorm_pass)
calls <- calls[calls$am_id == am_lesion(calls$am_id), , drop = FALSE]
hs_arms <- sort(unique(calls$arm[calls$sample_id %in% an$sample_id]))
hs_long <- as_tibble(hailstorm_status(calls, an$sample_id, hs_arms)) |>
  transmute(Sample = sample_id, Arm = arm,
            hs = if_else(present == 1L, "Yes", "No"),
            Patient = an$patient[match(sample_id, an$sample_id)])

hs_calls <- hs_long |> filter(hs == "Yes") |> distinct(Sample, Patient, Arm)
flow_add("hailstorm-evaluated", n_distinct(hs_long$Patient), "patients",
         paste0(n_distinct(hs_calls$Patient), " with >= 1 hailstorm"))

cat(sprintf("Hailstorms: %d samples | %d patients | %d arms | %d positive calls\n",
            n_distinct(hs_long$Sample), n_distinct(hs_long$Patient),
            n_distinct(hs_long$Arm), nrow(hs_calls)))

id_col <- if (analysis_unit == "patient") "Patient" else "Sample"

unit_tbl <- cnapp |>
  mutate(unit_id = .data[[id_col]]) |>
  group_by(unit_id) |>
  summarise(n_cnapp_samples = n(),
            across(all_of(score_names), ~ collapse_fun(.x, na.rm = TRUE)),
            .groups = "drop")

cat(sprintf("Analysis unit = %s | n = %d (%d with >1 sample)\n",
            analysis_unit, nrow(unit_tbl), sum(unit_tbl$n_cnapp_samples > 1)))
unit_lab_flow <- paste0(analysis_unit, "s")
flow_add("CNApp analysis units", nrow(unit_tbl), unit_lab_flow,
         paste0("multi-sample units collapsed by ", p$collapse_fun))

# ---- global hailstorm status per unit (any arm) ----
hs_unit <- hs_long |>
  mutate(unit_id = .data[[id_col]]) |>
  group_by(unit_id, Arm) |>
  summarise(hs = any(hs == "Yes"), .groups = "drop") |>
  group_by(unit_id) |>
  summarise(n_hs_arms = sum(hs), .groups = "drop") |>
  mutate(hs_status = factor(if_else(n_hs_arms > 0, "Hailstorm+", "Hailstorm-"),
                            levels = c("Hailstorm-", "Hailstorm+")))

# ---- BRAF by patient (any mutant lesion wins) ----
braf_unit <- NULL
if (file.exists(braf_path)) {
  braf <- as_tibble(bra_braf(braf_path)) |>
    mutate(BRAF_Sanger = norm_braf(sanger), BRAF_WES = norm_braf(wes),
           BRAF_call = case_when(
             braf_source == "wes"    ~ BRAF_WES,
             braf_source == "sanger" ~ BRAF_Sanger,
             TRUE                    ~ coalesce(BRAF_WES, BRAF_Sanger)),
           BRAF = case_when(
             is.na(BRAF_call)                                 ~ NA_character_,
             braf_definition == "v600_only" &
               str_detect(BRAF_call, p$braf_v600_pattern)                 ~ "Mutant",
             braf_definition == "v600_only"                   ~ "Non-mutant",
             BRAF_call == "WT"                                ~ "Non-mutant",
             TRUE                                             ~ "Mutant"))
  braf_by_patient <- braf |>
    group_by(patient) |>
    summarise(BRAF_call = paste(sort(unique(na.omit(BRAF_call))), collapse = ";"),
              BRAF = if (all(is.na(BRAF))) NA_character_
              else if (any(BRAF == "Mutant", na.rm = TRUE)) "Mutant"
              else "Non-mutant",
              .groups = "drop")

  unit_patient <- unit_tbl |>
    mutate(patient = if (analysis_unit == "patient") unit_id
           else an$patient[match(unit_id, an$sample_id)])
  braf_unit <- unit_patient |>
    select(unit_id, patient) |>
    left_join(braf_by_patient, by = "patient") |>
    mutate(BRAF = factor(BRAF, levels = c("Non-mutant", "Mutant")))

  unmatched <- braf_unit |> filter(is.na(BRAF)) |> pull(patient) |> unique()
  if (length(unmatched) > 0)
    cat(sprintf("No BRAF call for %d patient(s), excluded from figure R15: %s\n",
                length(unmatched), paste(unmatched, collapse = ", ")))
} else {
  message("BRAF sheet not found — figure R15 will be skipped: ", braf_path)
}

# ---- FIGURE R15 — scores x BRAF ----
section("FIGURE R15 — scores x BRAF")

global_one <- function(d, group_col, group_levels, palette, fig_title, tag) {
  n_by_grp <- table(factor(d[[group_col]], levels = group_levels))
  if (any(n_by_grp < min_group_n)) {
    message(sprintf("  %s: groups are %s — below min_group_n, figure skipped", tag,
                    paste(sprintf("%s = %d", names(n_by_grp), n_by_grp),
                          collapse = ", ")))
    return(NULL)
  }

  res <- map_dfr(score_names, function(sc)
    compare_groups(d[[sc]], d[[group_col]], sc, group_levels))
  if (nrow(res) == 0) return(NULL)

  res <- res |> mutate(p_fdr = p.adjust(p_value, "BH")) |> mutate(contrast = tag, .before = 1)

  cat(sprintf("\n-- %s --\n", tag))
  print(res |> select(contrast, feature, n_ref, n_comp, median_ref, median_comp,
                      rank_biserial, p_value, p_fdr), width = Inf)
  write_tsv(res, file.path(out_tab, sprintf("bra_cnapp_by_%s.tsv", tag)))

  panels <- map(score_names, function(sc) {
    p  <- res |> filter(feature == sc) |> pull(p_fdr)
    dd <- d |> filter(!is.na(.data[[sc]]), !is.na(.data[[group_col]]))
    make_box(dd, sc, group_col, sc,
             paste0(sc, " score"), if (length(p) == 1) p else NA_real_, palette)
  })

  fig <- patchwork::wrap_plots(panels, nrow = 1) +
    patchwork::plot_annotation(
      title = fig_title,
      caption = "Wilcoxon rank-sum; p.adj: Benjamini-Hochberg across the three scores. Points: patients.",
      theme = theme(plot.caption = element_text(size = 8.5, hjust = 0),
                    plot.title = element_text(face = "bold", size = 12.5,
                                              colour = "black", hjust = 0,
                                              margin = margin(b = 8))))

  save_fig(fig, sprintf("bra_cnapp_by_%s", tag),
           width  = panel_w_global * length(score_names),
           height = panel_h_global + 0.3)
  cat(sprintf("  (n = %d)\n", sum(!is.na(d[[group_col]]))))

  res
}

if (is.null(braf_unit)) {
  res_braf <- NULL
} else {
  braf_df <- unit_tbl |> left_join(select(braf_unit, unit_id, BRAF), by = "unit_id")
  cat(sprintf("BRAF funnel: %d %ss -> %d genotyped -> %d Mutant / %d Non-mutant\n",
              nrow(braf_df), analysis_unit, sum(!is.na(braf_df$BRAF)),
              sum(braf_df$BRAF == "Mutant", na.rm = TRUE),
              sum(braf_df$BRAF == "Non-mutant", na.rm = TRUE)))
  flow_add("CNApp + BRAF genotype (figure R15)", sum(!is.na(braf_df$BRAF)),
           unit_lab_flow, paste0("BRAF_source = ", braf_source))
  res_braf <- global_one(filter(braf_df, !is.na(BRAF)), "BRAF",
                         c("Non-mutant", "Mutant"), braf_pal,
                         "CNApp copy-number scores by BRAF status (Brazilian cohort)",
                         "braf")
  write_tsv(braf_df |> filter(BRAF == "Mutant") |> select(unit_id) |>
              left_join(hs_calls |> mutate(unit_id = .data[[id_col]]) |>
                          distinct(unit_id, Arm), by = "unit_id"),
            file.path(out_tab, "bra_braf_hailstorm_arms.tsv"))
}

# ---- SUMMARY ----
section("SUMMARY")

cat(sprintf("CNApp source:     %s\n", basename(cnapp_path)))
cat(sprintf("Unit:             %s (n = %d)\n", analysis_unit, nrow(unit_tbl)))
cat(sprintf("BRAF definition:  %s / source %s\n", braf_definition, braf_source))
cat(sprintf("\nTables:  %s\nFigures: %s (PDF only)\nDone.\n", out_tab, out_fig))

module_finish(M)
