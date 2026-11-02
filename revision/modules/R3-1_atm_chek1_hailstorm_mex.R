# =============================================================================
# ATM/CHEK1 loss vs hailstorms, Mexican cohort (Fisher, BH within gene)
# Run:      Rscript run.R atm_chek1_mex
# Outputs:
#   tables/atm_chek1_3arms_mex_results.tsv   Table R9, MEX rows
#   tables/atm_chek1_3arms_mex_{samples,patients,gain_excluded}.tsv
# Inputs:   config/params.R -> PARAMS$modules$atm_chek1_mex$inputs
#   supp          Supplementary Table 1: Sample, Patient, Ploidy
#   hailstorms    hailstorm calls: Sample, Arm, pass
#   detection_log detector run log: Sample, status
#   sequenza_dir  (parameter) <dir>/<Sample>/*segments.txt: chromosome, start.pos, end.pos, CNt
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("atm_chek1_mex")

suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr) })

p        <- M$params
SUPP     <- M$inputs$supp
CALLS    <- M$inputs$hailstorms
RUN_LOG  <- M$inputs$detection_log
SEQ_DIR  <- p$sequenza_dir
HS_PASS  <- p$hs_pass_values
GENES    <- p$genes
ARMS     <- p$arms
OUT_DIR  <- M$tables

# ---- inputs ----
supp <- read_csv(SUPP, show_col_types = FALSE) |>
  transmute(Sample = str_trim(Sample), Patient = str_trim(Patient),
            supp_ploidy = suppressWarnings(as.numeric(Ploidy)))

run <- read_tsv(RUN_LOG, show_col_types = FALSE) |>
  mutate(Sample = str_trim(Sample), status = str_trim(status)) |>
  filter(status != "no_breakpoints") |>
  left_join(supp, by = "Sample")
cat(sprintf("evaluated samples: %d | patients: %d\n",
            nrow(run), n_distinct(run$Patient)))

calls <- read_tsv(CALLS, show_col_types = FALSE) |>
  mutate(Sample = str_trim(Sample), Arm = str_remove(str_trim(Arm), "^chr"))
if ("pass" %in% names(calls)) calls <- filter(calls, pass %in% HS_PASS)

# ---- gene status from Sequenza, per sample ----
find_one <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(f)) f[1] else NA_character_
}

status_one <- function(s) {
  out <- tibble(Sample = s, ploidy = NA_real_)
  for (g in names(GENES)) out[[paste0(g, "_cn")]] <- NA_real_
  for (g in names(GENES)) out[[g]] <- NA_character_
  seg_f <- find_one(file.path(SEQ_DIR, s), "segments\\.txt$")
  if (is.na(seg_f)) return(out)
  seg <- read.table(seg_f, header = TRUE, stringsAsFactors = FALSE)
  seg$chromosome <- sub("^chr", "", seg$chromosome)
  auto <- seg[seg$chromosome %in% p$autosomes & is.finite(seg$CNt), ]
  w <- auto$end.pos - auto$start.pos + 1
  out$ploidy <- sum(auto$CNt * w) / sum(w)
  ref <- round(out$ploidy)
  for (g in names(GENES)) {
    G <- GENES[[g]]
    hit <- seg[seg$chromosome == G$chr & seg$start.pos <= G$end &
                 seg$end.pos >= G$start & is.finite(seg$CNt), ]
    if (!nrow(hit)) next
    cn <- min(hit$CNt)
    out[[paste0(g, "_cn")]] <- cn
    out[[g]] <- if (cn < ref) "Loss" else if (cn > ref) "Gain" else "Neutral"
  }
  out
}

st <- bind_rows(lapply(run$Sample, status_one))

# ---- sample table ----
samples <- run |>
  select(Patient, Sample, supp_ploidy) |>
  left_join(st, by = "Sample") |>
  arrange(Patient, Sample)
for (a in ARMS)
  samples[[paste0("hs_", a)]] <- samples$Sample %in% calls$Sample[calls$Arm == a]
write_tsv(samples, file.path(OUT_DIR, "atm_chek1_3arms_mex_samples.tsv"))

chk <- samples |> filter(is.finite(supp_ploidy), is.finite(ploidy))
cat(sprintf(paste("ploidy check vs Supplementary Table 1: %d samples,",
                  "median |difference| = %.2f\n"),
            nrow(chk), median(abs(chk$ploidy - chk$supp_ploidy))))
for (g in names(GENES)) {
  if (any(is.na(samples[[g]])))
    cat("WARNING:", g, "status missing for:",
        paste(samples$Sample[is.na(samples[[g]])], collapse = ", "), "\n")
  cat(g, "sample-level status:",
      paste(names(table(samples[[g]])), table(samples[[g]]), sep = "=", collapse = "  "), "\n")
}

# ---- patient table ----
collapse_status <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) NA_character_ else if (any(x == "Loss")) "Loss" else
    if (any(x == "Neutral")) "Neutral" else "Gain"
}
patients <- samples |>
  group_by(patient = Patient) |>
  summarise(samples = paste(Sample, collapse = ";"),
            ATM = collapse_status(ATM), CHEK1 = collapse_status(CHEK1),
            across(all_of(paste0("hs_", ARMS)), ~ if (any(.x)) "Yes" else "No"),
            .groups = "drop")
write_tsv(patients, file.path(OUT_DIR, "atm_chek1_3arms_mex_patients.tsv"))

# ---- Fisher: each gene x 3 arms, BH within gene ----
fisher_one <- function(gene, arm) {
  d  <- patients[patients[[gene]] %in% c("Loss", "Neutral"), ]
  hs <- d[[paste0("hs_", arm)]] == "Yes"
  ls <- d[[gene]] == "Loss"
  ft <- fisher.test(table(factor(hs, c(TRUE, FALSE)), factor(ls, c(TRUE, FALSE))))
  tibble(gene = gene, arm = arm, n = nrow(d),
         n_gain_excluded = sum(patients[[gene]] %in% "Gain"),
         loss_with_hs = sum(hs & ls), n_with_hs = sum(hs),
         loss_without_hs = sum(!hs & ls), n_without_hs = sum(!hs),
         pct_with_hs = 100 * sum(hs & ls) / sum(hs),
         pct_without_hs = 100 * sum(!hs & ls) / sum(!hs),
         odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1],
         ci_high = ft$conf.int[2], p_value = ft$p.value)
}
res <- bind_rows(lapply(names(GENES), function(g)
  bind_rows(lapply(ARMS, function(a) fisher_one(g, a))))) |>
  group_by(gene) |>
  mutate(padj = p.adjust(p_value, "BH")) |>
  ungroup()
write_tsv(res, file.path(OUT_DIR, "atm_chek1_3arms_mex_results.tsv"))

excl <- bind_rows(lapply(names(GENES), function(g) {
  pts <- patients$patient[patients[[g]] %in% "Gain"]
  if (!length(pts)) return(NULL)
  samples |> filter(Patient %in% pts) |>
    transmute(gene = g, Patient, Sample, copy_number = .data[[paste0(g, "_cn")]],
              ploidy = round(ploidy, 2), status = .data[[g]])
}))
write_tsv(excl, file.path(OUT_DIR, "atm_chek1_3arms_mex_gain_excluded.tsv"))

cat(sprintf("\npatients: %d\n", nrow(patients)))
res |>
  transmute(gene, arm, n,
            with_hs = sprintf("%d/%d (%.0f%%)", loss_with_hs, n_with_hs, pct_with_hs),
            without_hs = sprintf("%d/%d (%.0f%%)", loss_without_hs, n_without_hs, pct_without_hs),
            OR = sprintf("%.2f (%.2f-%.2f)", odds_ratio, ci_low, ci_high),
            P = signif(p_value, 3), padj = signif(padj, 3)) |>
  as.data.frame() |> print()
if (nrow(excl)) { cat("\nExcluded (gain/amplification only):\n"); print(as.data.frame(excl)) }

flow_add("Patients in the Fisher tests", nrow(patients), "patients")
module_finish(M)
