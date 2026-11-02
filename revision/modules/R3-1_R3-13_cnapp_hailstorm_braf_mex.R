# =============================================================================
# Hailstorm frequency and CNApp x BRAF, Mexican cohort (Basurto-Lozada et al., 2026)
# Run:      Rscript run.R cnapp_mex
# Outputs:
#   figures/mex_hailstorm_frequency_by_arm.pdf   Figure R6
#   figures/mex_hailstorm_braf_v600_tiles.pdf    Figure R13
#   figures/mex_cnapp_by_braf.pdf                Figure R14, Table R15
# Inputs:   config/params.R -> PARAMS$modules$cnapp_mex$inputs
#   supp          Supplementary Table 1: Sample, Patient, BCS, FCS, GCS, BRAF_status, Ploidy, Purity, WGD
#   hailstorms    hailstorm calls: Sample, Arm, pass
#   detection_log detector run log: Sample, status
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("cnapp_mex")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter

p <- M$params
supp_path          <- M$inputs$supp
hailstorm_path     <- M$inputs$hailstorms
detection_log_path <- M$inputs$detection_log

out_dir <- M$dir
out_tab <- M$tables
out_fig <- M$figures

braf_definition <- p$braf_definition

hs_pass_values <- p$hs_pass_values

analysis_unit <- p$analysis_unit
collapse_fun  <- match.fun(p$collapse_fun)
min_group_n   <- p$min_group_n

freq_denominator <- p$freq_denominator

exclude_no_breakpoints <- p$exclude_no_breakpoints

min_arm_cases   <- p$min_arm_cases
arms_to_exclude <- p$arms_to_exclude

bar_fill <- "#E8452F"
pal_hs   <- c(No = "grey75", Yes = "#E8452F")
pal_braf <- c("grey75", bar_fill)

# ---- HELPERS ----

fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.4f", p)

section <- function(txt) cat("\n", strrep("=", 72), "\n", txt, "\n",
                             strrep("=", 72), "\n", sep = "")

save_fig <- function(plot, stem, width, height) {
  ggsave(file.path(out_fig, paste0(stem, ".pdf")), plot, width = width, height = height)
  cat(sprintf("  figure -> %s.pdf\n", stem))
}

compare_groups <- function(values, group, feature, group_levels, min_n = min_group_n) {
  r <- wilcox_two_group(values, group, feature, group_levels, min_n)
  if (is.null(r)) NULL else as_tibble(r)
}

make_box <- function(df, value_col, group_col, title_str, y_lab,
                     p_value = NA_real_, palette = pal_braf) {
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

fisher_report <- function(tab, label) {
  if (!all(dim(tab) == c(2, 2))) return(NULL)
  ft <- fisher.test(tab)
  cat(sprintf("%s: OR = %.3f (95%% CI %.3f-%s), P = %s\n",
              label, unname(ft$estimate), ft$conf.int[1],
              ifelse(is.infinite(ft$conf.int[2]), "Inf", sprintf("%.3f", ft$conf.int[2])),
              fmt_p(ft$p.value)))
  tibble(comparison = label, n_total = sum(tab),
         odds_ratio = unname(ft$estimate),
         ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
         p_value = ft$p.value)
}

# ---- LOAD ----
section("LOADING")

supp <- read_csv(supp_path, show_col_types = FALSE) |>
  mutate(across(c(Sample, Patient, BRAF_status), ~ str_trim(as.character(.x))))

supp <- supp |>
  mutate(BRAF = case_when(
    is.na(BRAF_status)                       ~ NA_character_,
    braf_definition == "v600_only" &
      str_detect(BRAF_status, p$braf_v600_pattern)        ~ "Mutant",
    BRAF_status == "WT"                      ~ "Non-mutant",
    braf_definition == "v600_only"           ~ NA_character_,
    TRUE                                     ~ "Mutant"))

cat(sprintf("supp:         %d samples, %d patients | BRAF-mutant %d samples\n",
            nrow(supp), n_distinct(supp$Patient),
            sum(supp$BRAF == "Mutant", na.rm = TRUE)))
flow_add("supplementary table", n_distinct(supp$Patient), "patients")

# ---- CNApp-evaluable subset ----
cn <- supp |> filter(!is.na(GCS) & !is.na(BCS) & !is.na(FCS))
cat(sprintf("cn_evaluable: %d samples, %d patients\n",
            nrow(cn), n_distinct(cn$Patient)))
flow_add("CNApp-evaluable (non-NA BCS/FCS/GCS)", n_distinct(cn$Patient), "patients")

# ---- detection run log ----
det_log <- read_tsv(detection_log_path, show_col_types = FALSE) |>
  mutate(Sample = str_trim(Sample), status = str_trim(status)) |>
  left_join(select(supp, Sample, Patient), by = "Sample")

if (any(is.na(det_log$Patient)))
  warning("Run-log samples absent from the supplementary table: ",
          paste(unique(det_log$Sample[is.na(det_log$Patient)]), collapse = ", "),
          call. = FALSE)

cat("run log status: ",
    paste(names(table(det_log$status)), table(det_log$status),
          sep = "=", collapse = "  "), "\n", sep = "")

det_eval <- if (exclude_no_breakpoints)
  filter(det_log, status != "no_breakpoints") else det_log

cat(sprintf("detected:     %d samples, %d patients (no_breakpoints %s)\n",
            nrow(det_eval), n_distinct(det_eval$Patient),
            if (exclude_no_breakpoints) "excluded" else "kept"))
flow_add("hailstorm detection evaluable", n_distinct(det_eval$Patient), "patients",
         paste0("no_breakpoints ", if (exclude_no_breakpoints) "excluded" else "kept"))

# ---- hailstorm calls ----
hs_raw <- read_tsv(hailstorm_path, show_col_types = FALSE) |>
  mutate(Sample = str_trim(Sample),
         Arm    = str_remove(str_trim(Arm), "^chr"))

if ("pass" %in% names(hs_raw)) {
  cat("pass values: ",
      paste(names(table(hs_raw$pass)), table(hs_raw$pass), sep = "=", collapse = "  "),
      "\n", sep = "")
  hs_raw <- hs_raw |> filter(pass %in% hs_pass_values)
}

hs_calls <- hs_raw |>
  distinct(Sample, Arm) |>
  left_join(select(supp, Sample, Patient), by = "Sample")

if (any(is.na(hs_calls$Patient)))
  warning("Hailstorm samples absent from the supplementary table: ",
          paste(unique(hs_calls$Sample[is.na(hs_calls$Patient)]), collapse = ", "),
          call. = FALSE)

cat(sprintf("hailstorms:   %d calls | %d samples | %d patients | %d arms\n",
            nrow(hs_calls), n_distinct(hs_calls$Sample),
            n_distinct(hs_calls$Patient), n_distinct(hs_calls$Arm)))

log_n <- det_log |>
  filter(str_detect(status, "^hailstorm:")) |>
  transmute(Sample, n_log = as.integer(str_remove(status, "^hailstorm:")))
obs_n <- hs_raw |> count(Sample, name = "n_file")
disc  <- full_join(log_n, obs_n, by = "Sample") |>
  replace_na(list(n_log = 0L, n_file = 0L)) |>
  filter(n_log != n_file)
if (nrow(disc) > 0) {
  cat("\nNOTE: run log and calls file disagree (expected if pass filtering is on):\n")
  print(disc, n = Inf)
}
hs_burden <- hs_calls |> count(Sample, name = "n_hs_arms")

analysis_df <- cn |>
  select(Sample, Patient, BCS, FCS, GCS, BRAF, Ploidy, Purity, WGD,
         any_of(c("Sample_type", "Tumour_stage", "TMB"))) |>
  left_join(hs_burden, by = "Sample") |>
  replace_na(list(n_hs_arms = 0L)) |>
  mutate(hs_status = factor(if_else(n_hs_arms > 0, "Hailstorm+", "Hailstorm-"),
                            levels = c("Hailstorm-", "Hailstorm+")))

write_tsv(analysis_df, file.path(out_tab, "mex_analysis_table_sample_level.tsv"))

if (analysis_unit == "patient") {
  n_before <- nrow(analysis_df)
  analysis_df <- analysis_df |>
    group_by(Patient) |>
    summarise(
      n_samples = n(),
      Sample    = paste(sort(Sample), collapse = ";"),
      across(any_of(c("BCS", "FCS", "GCS", "Ploidy")), ~ collapse_fun(.x, na.rm = TRUE)),
      Purity = max(Purity, na.rm = TRUE),
      WGD    = max(WGD, na.rm = TRUE),
      across(any_of("TMB"), ~ collapse_fun(.x, na.rm = TRUE)),
      n_hs_arms = max(n_hs_arms, na.rm = TRUE),
      hs_status = if (any(hs_status == "Hailstorm+")) "Hailstorm+" else "Hailstorm-",
      BRAF = if (all(is.na(BRAF))) NA_character_
      else if (any(BRAF == "Mutant", na.rm = TRUE)) "Mutant" else "Non-mutant",
      .groups = "drop") |>
    mutate(hs_status = factor(hs_status, levels = c("Hailstorm-", "Hailstorm+")),
           BRAF      = factor(BRAF, levels = c("Non-mutant", "Mutant")))
  cat(sprintf("\nCollapsed %d samples -> %d patients (%d multi-region)\n",
              n_before, nrow(analysis_df), sum(analysis_df$n_samples > 1)))
} else {
  analysis_df <- analysis_df |>
    mutate(BRAF = factor(BRAF, levels = c("Non-mutant", "Mutant")))
}

write_tsv(analysis_df, file.path(out_tab, "mex_analysis_table.tsv"))
cat(sprintf("Analysis unit = %s | n = %d | Hailstorm+ = %d | BRAF-mutant = %d\n",
            analysis_unit, nrow(analysis_df),
            sum(analysis_df$hs_status == "Hailstorm+"),
            sum(analysis_df$BRAF == "Mutant", na.rm = TRUE)))
flow_add("score analyses (CNApp-evaluable)", nrow(analysis_df), paste0(analysis_unit, "s"))

score_names <- c("BCS", "FCS", "GCS")
section("FIGURE 1 — hailstorm frequency per arm")

universe <- switch(
  freq_denominator,
  "detected"     = det_eval,
  "cn_evaluable" = cn,
  "all_samples"  = supp,
  stop("unknown freq_denominator: ", freq_denominator))

if (analysis_unit == "patient") {
  all_units <- unique(universe$Patient); unit_col <- "Patient"; unit_lab <- "patients"
} else {
  all_units <- unique(universe$Sample);  unit_col <- "Sample";  unit_lab <- "samples"
}
all_units <- all_units[!is.na(all_units)]
n_units   <- length(all_units)

freq_df <- hs_calls |>
  filter(.data[[unit_col]] %in% all_units) |>
  distinct(.data[[unit_col]], Arm) |>
  count(Arm, name = "n_cases") |>
  mutate(pct = 100 * n_cases / n_units) |>
  filter(n_cases > 0) |>
  mutate(chr_num = suppressWarnings(as.integer(str_extract(Arm, "^\\d+"))),
         chr_num = if_else(is.na(chr_num), 99L, chr_num),
         arm_let = str_extract(Arm, "[pq]$")) |>
  arrange(pct, desc(chr_num), arm_let) |>
  mutate(Arm = factor(Arm, levels = Arm))

cat(sprintf("Denominator: %d %s (%s) | %d arms\n",
            n_units, unit_lab, freq_denominator, nrow(freq_df)))

freq_tab <- freq_df |>
  transmute(Arm = as.character(Arm), n_cases, n_total = n_units,
            pct = round(pct, 1),
            label = sprintf("%d/%d (%.1f%%)", n_cases, n_total, pct)) |>
  arrange(desc(n_cases), Arm)
print(freq_tab, n = Inf)
write_tsv(freq_tab, file.path(out_tab, "mex_hailstorm_frequency_by_arm.tsv"))

x_max <- max(freq_df$pct)
x_lim <- if (x_max <= 5) ceiling(x_max * 2) / 2 else ceiling(x_max / 10) * 10
x_brk <- pretty(c(0, x_lim), n = 5); x_brk <- x_brk[x_brk <= x_lim]

fig1 <- ggplot(freq_df, aes(x = pct, y = Arm)) +
  geom_col(fill = bar_fill, colour = "black", linewidth = 0.45, width = 0.72) +
  scale_x_continuous(limits = c(0, x_lim), breaks = x_brk,
                     labels = format(x_brk, trim = TRUE),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Cases with hailstorm (%)", y = "Chromosome arms") +
  theme_classic(base_size = 13) +
  theme(axis.text  = element_text(colour = "black", face = "bold", size = 11),
        axis.title = element_text(colour = "black", size = 13),
        axis.line  = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.ticks.length = unit(0.15, "cm"),
        plot.margin = margin(8, 12, 8, 8))

save_fig(fig1, "mex_hailstorm_frequency_by_arm",
         width = 3.6, height = max(4, 0.26 * nrow(freq_df) + 1.2))

# ---- missingness diagnostic ----
cat("\n-- CNApp availability vs hailstorm status (detected universe) --\n")
miss <- tibble(Patient = unique(det_eval$Patient)) |>
  filter(!is.na(Patient)) |>
  mutate(hailstorm = if_else(Patient %in% hs_calls$Patient, "Hailstorm+", "Hailstorm-"),
         cnapp     = if_else(Patient %in% cn$Patient, "Scored", "No scores"))
tab_miss <- table(miss$hailstorm, miss$cnapp)
print(tab_miss)
miss_out <- fisher_report(tab_miss, "hailstorm status x CNApp availability")
if (!is.null(miss_out))
  write_tsv(miss_out, file.path(out_tab, "mex_cnapp_missingness_vs_hailstorm.tsv"))
write_tsv(miss, file.path(out_tab, "mex_cnapp_missingness_by_patient.tsv"))
section("FIGURE 4 — hailstorm distribution across BRAF-mutant cases")

tile_universe <- "tested"
tile_unit     <- "patient"

braf_mut <- supp |>
  filter(!is.na(BRAF_status), BRAF_status != "WT") |>
  mutate(is_v600 = str_detect(BRAF_status, p$braf_v600_pattern)) |>
  select(Sample, Patient, BRAF_status, is_v600)

n_non_v600 <- sum(!braf_mut$is_v600)
if (n_non_v600 > 0)
  cat(sprintf("Non-V600 BRAF mutations shown and labelled (%d): %s\n", n_non_v600,
              paste(sprintf("%s %s", braf_mut$Sample[!braf_mut$is_v600],
                            braf_mut$BRAF_status[!braf_mut$is_v600]), collapse = ", ")))

tested_samples <- unique(det_eval$Sample)
untested <- setdiff(braf_mut$Sample, tested_samples)

if (length(untested) > 0)
  cat(sprintf("BRAF-mutant samples never run through the detector (%d): %s\n",
              length(untested), paste(sort(untested), collapse = ", ")))

tile_samples <- if (tile_universe == "tested")
  filter(braf_mut, Sample %in% tested_samples) else braf_mut

unit_of <- tile_samples |>
  mutate(unit = if (tile_unit == "patient") Patient else Sample) |>
  select(Sample, unit)

tile_df_base <- tile_samples |>
  mutate(unit = if (tile_unit == "patient") Patient else Sample) |>
  group_by(unit) |>
  summarise(Patient   = first(Patient),
            n_samples = n(),
            is_v600   = all(is_v600),
            mut_lab   = paste(sort(unique(BRAF_status[!is_v600])), collapse = "/"),
            .groups   = "drop") |>
  mutate(unit_lab = if_else(is_v600, unit, paste0(unit, " (", mut_lab, ")")))

n_collapsed <- nrow(tile_samples) - nrow(tile_df_base)
cat(sprintf("Universe = %s | unit = %s | %d samples -> %d rows (%d replicate(s) collapsed)\n",
            tile_universe, tile_unit, nrow(tile_samples), nrow(tile_df_base), n_collapsed))
if (n_collapsed > 0)
  cat(sprintf("  multi-sample: %s\n",
              paste(sprintf("%s (n=%d)",
                            tile_df_base$unit[tile_df_base$n_samples > 1],
                            tile_df_base$n_samples[tile_df_base$n_samples > 1]),
                    collapse = ", ")))

draw_tiles <- function(dat, stem, label, title_str) {

  if (nrow(dat) == 0) {
    cat(sprintf("%-24s no cases — skipped\n", label)); return(invisible(NULL))
  }

  hits <- hs_calls |>
    inner_join(unit_of, by = "Sample") |>
    filter(unit %in% dat$unit) |>
    distinct(unit, Arm)

  arms_shown <- hits |>
    distinct(Arm) |>
    mutate(chr_num = suppressWarnings(as.integer(str_extract(Arm, "^\\d+"))),
           chr_num = if_else(is.na(chr_num), 99L, chr_num),
           arm_let = str_extract(Arm, "[pq]$")) |>
    arrange(chr_num, arm_let) |>
    pull(Arm)

  if (length(arms_shown) == 0) {
    cat(sprintf("%-24s no hailstorm in any case — skipped\n", label))
    return(invisible(NULL))
  }

  hit_u <- unique(hits$unit)

  cat(sprintf("%-24s %2d rows (%s) | arms: %s\n",
              label, nrow(dat), tile_unit, paste(arms_shown, collapse = ", ")))
  cat(sprintf("%-24s hailstorm+ %d/%d (%.1f%%)\n",
              "", length(hit_u), nrow(dat), 100 * length(hit_u) / nrow(dat)))

  lab_levels <- dat |> arrange(unit) |> pull(unit_lab)

  tile_df <- tidyr::expand_grid(unit = sort(dat$unit), Arm = arms_shown) |>
    left_join(mutate(hits, hit = TRUE), by = c("unit", "Arm")) |>
    left_join(select(dat, unit, unit_lab, is_v600, n_samples), by = "unit") |>
    mutate(Status = factor(if_else(is.na(hit), "No Hit", "Hailstorm Hit"),
                           levels = c("No Hit", "Hailstorm Hit")),
           Arm      = factor(Arm, levels = arms_shown),
           unit_lab = factor(unit_lab, levels = rev(lab_levels)))

  write_tsv(tile_df |> select(unit, unit_lab, is_v600, n_samples, Arm, Status),
            file.path(out_tab, paste0(stem, ".tsv")))

  is_flagged <- str_detect(levels(tile_df$unit_lab), "\\(")
  lab_face   <- if_else(is_flagged, "italic", "plain")
  lab_colour <- if_else(is_flagged, "#8B1A1A", "black")

  fig <- ggplot(tile_df, aes(x = Arm, y = unit_lab, fill = Status)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    scale_fill_manual(values = c("No Hit" = "grey88", "Hailstorm Hit" = bar_fill)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = title_str, x = "Arm",
         y = if (tile_unit == "patient") "Patient" else "Sample", fill = "Status") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 11.5, colour = "black",
                                    hjust = 0, margin = margin(b = 8)),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = lab_face, colour = lab_colour),
          axis.title = element_text(colour = "black", size = 12),
          legend.position = "right")

  save_fig(fig, stem,
           width  = max(5.5, 0.55 * length(arms_shown) + 4.4),
           height = max(3.5, 0.28 * nrow(dat) + 2.0))
}

draw_tiles(filter(tile_df_base, is_v600), "mex_hailstorm_braf_v600_tiles", "V600E only:",
           "Hailstorm distribution across BRAF V600E cases (Mexican cohort)")
section("FIGURE 2 — BRAF x CNApp")

braf_df <- analysis_df |> filter(!is.na(BRAF))
flow_add("BRAF comparison (figure 2)", nrow(braf_df), paste0(analysis_unit, "s"),
         paste0("braf_definition = ", braf_definition))
cat(sprintf("n = %d (%d mutant, %d non-mutant)\n", nrow(braf_df),
            sum(braf_df$BRAF == "Mutant"), sum(braf_df$BRAF == "Non-mutant")))

braf_res <- map_dfr(score_names, function(sc)
  compare_groups(braf_df[[sc]], braf_df$BRAF, sc, c("Non-mutant", "Mutant")))

if (nrow(braf_res) > 0) {
  braf_res <- braf_res |> mutate(p_fdr = p.adjust(p_value, "BH")) |> arrange(p_value)
  print(braf_res, width = Inf)
  write_tsv(braf_res, file.path(out_tab, "mex_cnapp_by_braf.tsv"))

  panels <- map(score_names, function(sc) {
    p <- braf_res |> filter(feature == sc) |> pull(p_fdr)
    make_box(braf_df, sc, "BRAF", sc, paste0(sc, " score"),
             if (length(p) == 1) p else NA_real_, palette = pal_braf)
  })
  fig2 <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = "CNApp copy-number scores by BRAF status (Mexican cohort)",
      theme = theme(plot.title = element_text(face = "bold", size = 12.5,
                                              colour = "black", hjust = 0,
                                              margin = margin(b = 8))))
  save_fig(fig2, "mex_cnapp_by_braf", width = 9.5, height = 4.1)
} else {
  message("No BRAF comparison reached min_group_n — fig2 not produced.")
}

# ---- SUMMARY ----
section("SUMMARY")
cat(sprintf("Universes:          supp %d | detected %d | cn_evaluable %d (patients)\n",
            n_distinct(supp$Patient), n_distinct(det_eval$Patient),
            n_distinct(cn$Patient)))
cat(sprintf("Frequency figure:   %s, n = %d %s\n", freq_denominator, n_units, unit_lab))
cat(sprintf("Most frequent arm:  %s (%d/%d, %.1f%%)\n",
            freq_tab$Arm[1], freq_tab$n_cases[1], freq_tab$n_total[1], freq_tab$pct[1]))
cat(sprintf("Score analyses:     %s (n = %d) | Hailstorm+ %d | BRAF-mutant %d\n",
            analysis_unit, nrow(analysis_df),
            sum(analysis_df$hs_status == "Hailstorm+"),
            sum(analysis_df$BRAF == "Mutant", na.rm = TRUE)))
cat(sprintf("\nTables:  %s\nFigures: %s\n", out_tab, out_fig))
cat("Figures written: mex_hailstorm_frequency_by_arm, mex_hailstorm_braf_v600_tiles,\n",
    "                 mex_cnapp_by_braf\nDone.\n", sep = "")

module_finish(M)
