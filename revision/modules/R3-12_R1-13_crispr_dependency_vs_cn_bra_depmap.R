# =============================================================================
# CRISPR dependency vs copy number, AM cell lines + DepMap
# Run:      Rscript run.R crispr_cn
# Outputs:
#   figures/harmonised/CN_call_per_cell_line.pdf      Figures R3, R8
#   figures/harmonised/CN_call_<detail_gene>.pdf      Figure R16
#   figures/correlation/hit_correlation_chronos.pdf   Figure R9
#   figures/correlation/<detail_gene>_cn_vs_chronos.pdf  Figure R17
#   figures/hit_waterfalls/hit_genes_waterfall_chronos.pdf  Figure R10
#   figures/hit_waterfalls/<detail_gene>_waterfall_chronos.pdf  Figure R18
#   tables/  rho, FDR, background FDR, percentiles, log2 ratios
# Inputs:   config/params.R -> PARAMS$modules$crispr_cn$inputs
#   model         DepMap Model.csv
#   crispr        DepMap CRISPRGeneEffect.csv (models x genes, Chronos)
#   cn_gene       DepMap OmicsCNGeneWGS.csv (models x genes, relative ratio)
#   bayes_combinedBAGEL2 Bayes Factors: Gene + one column per line
#   cns_AM007     CNVkit .cns: chromosome, start, end, gene, log2
#   cns_AM016     CNVkit .cns: chromosome, start, end, gene, log2
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("crispr_cn")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
  library(ggrepel)
})

options(stringsAsFactors = FALSE, datatable.fread.datatable = FALSE)

# ---- CONFIGURATION ----

p <- M$params
CFG <- c(p[setdiff(names(p), c("inputs", "optional", "seed"))],
         list(out_root = M$dir,
              files = list(bayes_combined  = M$inputs$bayes_combined,
                           model           = M$inputs$model,
                           cn_depmap       = M$inputs$cn_gene,
                           crispr_chronos  = M$inputs$crispr,
                           cns             = M$shared$cns)))

# ---- UTILITIES ----

DIRS <- list(root   = M$reports,
             tables = M$tables,
             figs   = M$figures,
             qc     = M$provenance)
invisible(lapply(DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))


.LOG <- character(0)
say <- function(...) {
  m <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  .LOG <<- c(.LOG, m); message(m)
}

need <- function(p, what) {
  if (!file.exists(p)) stop(sprintf("[missing %s] %s", what, p), call. = FALSE)
  invisible(TRUE)
}

out_path <- function(base, name, sub = NULL) {
  d <- if (is.null(sub)) base else file.path(base, sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, name)
}

wtsv <- function(x, name, sub = NULL) {
  p <- out_path(DIRS$tables, name, sub)
  readr::write_tsv(x, p, na = "")
  say("[table] ", if (is.null(sub)) name else file.path(sub, name),
      " (", nrow(x), " rows)")
  invisible(p)
}

sfig <- function(pl, name, w = 10, h = 7, sub = NULL) {
  p <- out_path(DIRS$figs, name, sub)
  ggsave(p, pl, width = w, height = h, dpi = CFG$fig_dpi, limitsize = FALSE)
  say("[figure] ", if (is.null(sub)) name else file.path(sub, name))
  invisible(p)
}

strip_entrez <- function(x) trimws(sub("\\s*\\(\\d+\\)\\s*$", "", x))


# ---- readr parsing diagnostics ----
read_csv_checked <- function(path, label) {
  need(path, label)
  x <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  pb <- readr::problems(x)
  if (nrow(pb)) {
    say("[parse] ", label, ": ", nrow(pb), " parsing issue(s) in column(s) ",
        paste(unique(pb$col), collapse = ", "),
        " | expected ", paste(unique(pb$expected), collapse = "; "))
    wtsv(head(pb, 500), paste0("parsing_issues_", label, ".tsv"))
  }
  x
}
norm_chr     <- function(v) sub("^chr", "", as.character(v))

# ---- CN scale handling ----

force_log2 <- function(v, floor_val = 0.05, label = "") {
  sc <- detect_cn_scale(v)
  say("[cn-scale] ", label, ": detected ", sc,
      " | min=", sprintf("%.4f", min(v, na.rm = TRUE)),
      " median=", sprintf("%.4f", median(v, na.rm = TRUE)),
      " max=", sprintf("%.4f", max(v, na.rm = TRUE)),
      " | frac<0 = ", sprintf("%.4f", mean(v < 0, na.rm = TRUE)))
  if (identical(sc, "log2")) return(v)
  n_floor <- sum(v < floor_val, na.rm = TRUE)
  v[v < floor_val] <- floor_val
  say("[cn-scale] ", label, ": converting linear -> log2, floored ",
      n_floor, " values at ", floor_val)
  log2(v)
}

# ---- CN classification ----


classify_cn_ratio <- function(log2r, th = CFG$cn_ratio) cn_call_log2(log2r, th)

is_gain_amp <- function(x) as.character(x) %in% c("Gain", "Amp")


# ---- Dependency harmonisation ----
dep_percentile <- function(bf) {
  n <- sum(!is.na(bf))
  if (n < 2) return(rep(NA_real_, length(bf)))
  100 * (rank(bf, na.last = "keep", ties.method = "average") - 1) / (n - 1)
}

# ---- Exact Spearman for tiny n ----
perms_n <- function(n) {
  if (n == 1) return(matrix(1L))
  sub <- perms_n(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i)
    cbind(i, matrix(seq_len(n)[-i][sub], nrow = nrow(sub)))))
}
.PC <- new.env(parent = emptyenv())
get_perms <- function(n) {
  k <- as.character(n)
  if (is.null(.PC[[k]])) .PC[[k]] <- perms_n(n)
  .PC[[k]]
}
spearman_exact <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; n <- length(x)
  if (n < 3) return(tibble(n = n, r = NA_real_, p = NA_real_))
  rx <- rank(x); ry <- rank(y)
  r <- suppressWarnings(cor(rx, ry))
  if (!is.finite(r)) return(tibble(n = n, r = NA_real_, p = NA_real_))
  p <- if (n <= 7) {
    P <- get_perms(n)
    null_r <- apply(P, 1, function(idx) cor(rx, ry[idx]))
    mean(abs(null_r) >= abs(r) - 1e-12)
  } else {
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  }
  tibble(n = n, r = r, p = p)
}

row_spearman_pairwise <- function(A_mat, B_mat, min_n) {
  stopifnot(identical(dim(A_mat), dim(B_mat)))
  n <- nrow(A_mat)
  r <- rep(NA_real_, n); k <- integer(n)
  for (i in seq_len(n)) {
    a <- A_mat[i, ]; b <- B_mat[i, ]
    m <- is.finite(a) & is.finite(b)
    k[i] <- sum(m)
    if (k[i] < min_n) next
    ra <- rank(a[m]); rb <- rank(b[m])
    if (sd(ra) == 0 || sd(rb) == 0) next
    r[i] <- cor(ra, rb)
  }
  list(r = r, n = k)
}

read_cns <- function(path, line) {
  need(path, paste("CNVkit .cns for", line))
  x <- fread(path)
  nm <- tolower(names(x))
  gc <- function(p) { i <- which(nm == p); if (length(i)) names(x)[i[1]] else NA }
  c_chr <- gc("chromosome"); c_s <- gc("start"); c_e <- gc("end")
  c_l2  <- gc("log2");       c_g <- gc("gene")
  c_w   <- gc("weight");     c_p  <- gc("probes")
  if (any(is.na(c(c_chr, c_s, c_e, c_l2))))
    stop("unexpected .cns layout (need chromosome/start/end/log2): ", path)

  seg <- tibble(
    cell_line = line,
    chrom = norm_chr(x[[c_chr]]),
    start = as.numeric(x[[c_s]]),
    end   = as.numeric(x[[c_e]]),
    log2r = as.numeric(x[[c_l2]]),
    gene_field = if (!is.na(c_g)) as.character(x[[c_g]]) else NA_character_,
    weight = if (!is.na(c_w)) as.numeric(x[[c_w]]) else NA_real_,
    probes = if (!is.na(c_p)) as.numeric(x[[c_p]]) else NA_real_
  )
  say("[cnvkit] ", line, ": ", nrow(seg), " segments | log2 range ",
      sprintf("%.2f", min(seg$log2r, na.rm = TRUE)), " to ",
      sprintf("%.2f", max(seg$log2r, na.rm = TRUE)))
  seg
}

genes_from_gene_field <- function(seg) {
  seg %>%
    filter(!is.na(gene_field), gene_field != "") %>%
    separate_rows(gene_field, sep = ",") %>%
    mutate(gene = toupper(trimws(gene_field)),
           width = end - start) %>%
    filter(!gene %in% CFG$cns_non_genes) %>%
    group_by(cell_line, gene) %>%
    slice_max(width, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(cell_line, gene, log2r, chrom, start, end)
}

load_am_bra_cn <- function() {
  segs <- imap_dfr(CFG$files$cns, ~ read_cns(.x, .y))
  wtsv(segs, "cnvkit_segments_raw.tsv")

  gl <- genes_from_gene_field(segs)
  say("[cnvkit] gene-level CN from the .cns gene column: ", nrow(gl),
      " gene x line entries, ", n_distinct(gl$gene), " distinct symbols")

  cover <- expand_grid(cell_line = CFG$am_bra, gene = CFG$hit_genes) %>%
    left_join(gl %>% select(cell_line, gene, log2r), by = c("cell_line", "gene"))
  absent <- cover %>% filter(!is.finite(log2r))
  if (nrow(absent))
    stop("panel genes absent from the .cns gene column: ",
         paste(sprintf("%s/%s", absent$cell_line, absent$gene), collapse = ", "),
         call. = FALSE)
  say("[cnvkit] all ", length(CFG$hit_genes), " panel genes annotated in both lines")

  sc <- detect_cn_scale(gl$log2r)
  if (!identical(sc, "log2"))
    warning("CNVkit values do not look like log2 ratios (detected: ", sc, ")")

  .SEG$gene_coord <- gl %>%
    filter(is.finite(start), is.finite(end)) %>%
    group_by(gene) %>% slice_head(n = 1) %>% ungroup() %>%
    select(gene, chrom, start, end)

  gl %>% mutate(cohort = "AM_BRA")
}
regression_test <- function(long) {
  rc <- CFG$regression_check
  obs <- long %>%
    filter(gene == toupper(rc$gene)) %>%
    select(cell_line, log2r) %>%
    distinct()
  if (!nrow(obs)) {
    say("[check] WARNING: ", rc$gene, " absent; regression test skipped")
    return(invisible(NULL))
  }
  cmp <- tibble(cell_line = names(rc$expected),
                expected_log2 = unname(rc$expected)) %>%
    left_join(obs, by = "cell_line") %>%
    mutate(observed_log2 = round(log2r, 3),
           delta = round(log2r - expected_log2, 3),
           pass = is.finite(log2r) & abs(delta) <= rc$tolerance) %>%
    select(-log2r)
  wtsv(cmp, "regression_check.tsv")
  nfail <- sum(!cmp$pass, na.rm = TRUE) + sum(is.na(cmp$pass))
  if (nfail == 0) {
    say("[check] ", rc$gene, " matches all ", nrow(cmp),
        " manually verified log2 values")
  } else {
    say("[check] WARNING: ", rc$gene, " differs from the verified values in ",
        nfail, "/", nrow(cmp), " lines")
  }
  invisible(cmp)
}

# ---- DEPMAP COPY NUMBER ----

.CN_CACHE <- new.env(parent = emptyenv())

report_missing_cn <- function(missing, requested) {
  if (!length(missing)) return(invisible(NULL))
  say("[depmap] no copy number for ", length(missing), " of ", requested,
      " requested models (identifiers in missing_cn_models.tsv)")
  prev <- if (file.exists(file.path(DIRS$tables, "missing_cn_models.tsv")))
    readr::read_tsv(file.path(DIRS$tables, "missing_cn_models.tsv"),
                    show_col_types = FALSE)$model_id else character(0)
  wtsv(tibble(model_id = sort(unique(c(prev, missing)))),
       "missing_cn_models.tsv")
}

load_depmap_cn <- function(model_ids) {
  if (!is.null(.CN_CACHE$mat)) {
    keep <- intersect(model_ids, colnames(.CN_CACHE$mat))
    report_missing_cn(setdiff(model_ids, colnames(.CN_CACHE$mat)),
                      length(model_ids))
    return(.CN_CACHE$mat[, keep, drop = FALSE])
  }

  p <- CFG$files$cn_depmap; need(p, "DepMap gene CN")
  say("[depmap] reading ", basename(p), " ...")
  cn <- fread(p)

  idcol <- intersect(c("ModelID", "V1", "...1"), names(cn))[1]
  if (is.na(idcol)) idcol <- names(cn)[1]

  meta_cols <- intersect(
    c("...1", "V1", "SequencingID", "ModelConditionID", "ProfileID",
      "IsDefaultEntryForMC", "IsDefaultEntryForModel"), names(cn))

  if ("IsDefaultEntryForModel" %in% names(cn)) {
    before <- nrow(cn)
    cn <- cn[cn$IsDefaultEntryForModel %in% c("Yes", TRUE, "TRUE"), , drop = FALSE]
    say("[depmap] IsDefaultEntryForModel filter: ", before, " -> ", nrow(cn), " rows")
  }

  ids <- cn[[idcol]]
  mat <- as.matrix(cn[, setdiff(names(cn), c(idcol, meta_cols)), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  colnames(mat) <- strip_entrez(colnames(mat))

  # ---- the critical step ----
  samp <- as.numeric(mat[sample(seq_len(nrow(mat)), min(50, nrow(mat))), ])
  qc <- tibble(
    file = basename(p),
    detected_scale = detect_cn_scale(samp),
    frac_negative = mean(samp < 0, na.rm = TRUE),
    min = min(samp, na.rm = TRUE),
    q25 = quantile(samp, .25, na.rm = TRUE),
    median = median(samp, na.rm = TRUE),
    q75 = quantile(samp, .75, na.rm = TRUE),
    max = max(samp, na.rm = TRUE))
  wtsv(qc, "depmap_cn_scale_QC.tsv")

  mat_log2 <- matrix(force_log2(as.numeric(mat), label = basename(p)),
                     nrow = nrow(mat), dimnames = dimnames(mat))

  .CN_CACHE$mat <- t(mat_log2)
  say("[depmap] CN matrix cached: ", nrow(.CN_CACHE$mat), " genes x ",
      ncol(.CN_CACHE$mat), " models")

  keep <- intersect(model_ids, colnames(.CN_CACHE$mat))
  report_missing_cn(setdiff(model_ids, colnames(.CN_CACHE$mat)),
                    length(model_ids))
  .CN_CACHE$mat[, keep, drop = FALSE]
}

.CHRONOS <- new.env(parent = emptyenv())

load_chronos <- function() {
  if (!is.null(.CHRONOS$mat)) return(.CHRONOS$mat)
  p <- CFG$files$crispr_chronos
  say("[chronos] reading ", basename(p), " ...")
  x <- fread(p)
  m <- omics_matrix(x, "chronos")
  if (!is.matrix(m)) { .CHRONOS$mat <- NA; return(NA) }
  say("[chronos] ", nrow(m), " genes x ", ncol(m), " models | median = ",
      sprintf("%.3f", median(m, na.rm = TRUE)),
      " (more negative = more dependent)")
  .CHRONOS$mat <- m
  m
}

.SEG <- new.env(parent = emptyenv())

omics_matrix <- function(x, tag) {
  nm <- names(x)
  meta_known <- c("ProfileID", "ModelID", "ModelConditionID", "SequencingID",
                  "IsDefaultEntryForModel", "IsDefaultEntryForMC",
                  "ProfileType", "V1", "...1")
  meta <- intersect(meta_known, nm)
  lead <- character(0)
  for (cnm in head(nm, 20)) {
    if (is.numeric(x[[cnm]]) && !cnm %in% meta_known) break
    lead <- c(lead, cnm)
  }
  meta <- unique(c(meta, lead))
  say("[", tag, "] metadata columns detected: ",
      paste(head(meta, 8), collapse = ", "),
      if (length(meta) > 8) " ..." else "")

  idcol <- if ("ModelID" %in% nm) "ModelID" else {
    hit <- meta[vapply(meta, function(cc)
      mean(grepl("^ACH-", as.character(x[[cc]]))) > 0.5, logical(1))]
    if (length(hit)) hit[1] else NA_character_
  }
  if (is.na(idcol)) {
    say("[", tag, "] WARNING: no ModelID column found. First columns are ",
        paste(head(nm, 5), collapse = ", "),
        ". A profile-to-model mapping (OmicsProfiles.csv) would be needed.")
    return(NA)
  }
  say("[", tag, "] keyed on ", idcol)

  if ("IsDefaultEntryForModel" %in% nm) {
    before <- nrow(x)
    x <- x[x$IsDefaultEntryForModel %in% c("Yes", "TRUE", TRUE, 1), ,
           drop = FALSE]
    say("[", tag, "] IsDefaultEntryForModel filter: ", before, " -> ",
        nrow(x), " rows")
  }
  ids <- as.character(x[[idcol]])
  dup <- duplicated(ids)
  if (any(dup)) {
    say("[", tag, "] dropping ", sum(dup), " duplicate model rows")
    x <- x[!dup, , drop = FALSE]; ids <- ids[!dup]
  }

  gene_cols <- setdiff(nm, meta)
  m <- as.matrix(x[, gene_cols, drop = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  colnames(m) <- strip_entrez(colnames(m))
  t(m)
}

.EXPR <- new.env(parent = emptyenv())

load_bayes <- function() {
  p <- CFG$files$bayes_combined
  bf <- read_csv_checked(p, "bayes_combined")
  say("[bf] combined matrix: ", nrow(bf), " genes x ", ncol(bf) - 1, " cell lines")

  expected <- c(unname(CFG$bf_cols_inhouse), unname(CFG$am_depmap))
  absent <- setdiff(expected, colnames(bf))
  if (length(absent))
    stop("combined Bayes Factor matrix is missing required columns: ",
         paste(absent, collapse = ", "), call. = FALSE)

  bf %>% mutate(Gene = toupper(Gene))
}

load_annotation <- function(bf) {
  mi <- read_csv_checked(CFG$files$model, "Model")
  ach <- keep(colnames(bf), ~ str_starts(.x, CFG$depmap_id_prefix))

  ann <- mi %>%
    filter(ModelID %in% ach) %>%
    select(ModelID, CellLineName, OncotreeSubtype,
           OncotreePrimaryDisease, OncotreeLineage) %>%
    mutate(mel_subtype = melanoma_subtype(OncotreeSubtype, OncotreePrimaryDisease),
           group = case_when(
      ModelID %in% unname(CFG$am_depmap)   ~ "AM_DepMap",
      mel_subtype %in% "Acral"             ~ "AM_DepMap_other",
      mel_subtype %in% "Cutaneous"         ~ "Cutaneous_Melanoma",
      mel_subtype %in% "Uveal"             ~ "Uveal_Melanoma",
      mel_subtype %in% "Melanoma NOS"      ~ "Melanoma_NOS",
      TRUE ~ NA_character_)) %>%
    filter(!is.na(group)) %>%
    select(-mel_subtype)

  ours <- tibble(
    ModelID = unname(CFG$bf_cols_inhouse),
    CellLineName = names(CFG$bf_cols_inhouse),
    OncotreeSubtype = CFG$am_bra_annotation$subtype,
    OncotreePrimaryDisease = CFG$am_bra_annotation$disease,
    OncotreeLineage = CFG$am_bra_annotation$lineage,
    group = "AM_BRA")

  out <- bind_rows(ours, ann)
  for (g in sort(unique(out$group)))
    flow_add(paste0("cell lines: ", g), sum(out$group == g), "models")
  say("[ann] ", paste(sprintf("%s n=%d", names(table(out$group)),
                              as.integer(table(out$group))), collapse = " | "))
  wtsv(out, "cell_line_annotation.tsv")

  all_models <- mi %>%
    filter(ModelID %in% ach) %>%
    transmute(model_id = ModelID, CellLineName, OncotreeSubtype,
              OncotreePrimaryDisease, OncotreeLineage) %>%
    left_join(out %>% select(model_id = ModelID, group), by = "model_id") %>%
    mutate(group = ifelse(is.na(group), "Other", group))
  say("[ann] pan-cancer scope: ", nrow(all_models), " models across ",
      n_distinct(all_models$OncotreeLineage), " lineages")
  wtsv(all_models, "all_models_annotation.tsv")

  attr(out, "all_models") <- all_models

  all_models_full <- mi %>%
    transmute(model_id = ModelID, CellLineName, OncotreeSubtype,
              OncotreePrimaryDisease, OncotreeLineage,
              mel_subtype = melanoma_subtype(OncotreeSubtype, OncotreePrimaryDisease)) %>%
    mutate(group = case_when(
      model_id %in% unname(CFG$am_depmap) ~ "AM_DepMap",
      mel_subtype %in% "Acral"            ~ "AM_DepMap_other",
      mel_subtype %in% "Cutaneous"        ~ "Cutaneous_Melanoma",
      mel_subtype %in% "Uveal"            ~ "Uveal_Melanoma",
      mel_subtype %in% "Melanoma NOS"     ~ "Melanoma_NOS",
      TRUE                                ~ "Other")) %>%
    select(-mel_subtype)
  wtsv(all_models_full, "all_models_annotation_full.tsv")
  attr(out, "all_models_full") <- all_models_full
  out
}

build_long <- function(bf, cn_bra, cn_dm) {
  bf_long <- bf %>%
    pivot_longer(-Gene, names_to = "bf_col", values_to = "BF") %>%
    group_by(bf_col) %>%
    mutate(dep_pct = dep_percentile(BF)) %>%
    ungroup() %>%
    rename(gene = Gene)
  lab <- c(setNames(names(CFG$bf_cols_inhouse), unname(CFG$bf_cols_inhouse)),
           setNames(names(CFG$am_depmap), unname(CFG$am_depmap)))
  bf_long <- bf_long %>%
    mutate(cell_line = ifelse(bf_col %in% names(lab), lab[bf_col], bf_col))
  cn_dm_long <- as.data.frame(cn_dm) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "model_id", values_to = "log2r") %>%
    mutate(gene = toupper(gene),
           cell_line = setNames(names(CFG$am_depmap), unname(CFG$am_depmap))[model_id],
           cohort = "AM_DepMap")
  cn_bra_long <- cn_bra %>%
    mutate(model_id = unname(CFG$bf_cols_inhouse)[cell_line]) %>%
    select(gene, model_id, log2r, cell_line, cohort)
  bf_long %>%
    inner_join(bind_rows(select(cn_dm_long, gene, model_id, log2r, cell_line, cohort),
                         cn_bra_long), by = c("gene", "cell_line")) %>%
    mutate(cn_call = cn_call_log2(log2r, CFG$cn_harmonised))
}

call_tables <- function(long) {
  lv <- c(CFG$focus_genes, setdiff(CFG$hit_genes, CFG$focus_genes))
  tab <- long %>%
    filter(gene %in% CFG$hit_genes) %>%
    mutate(gene = factor(gene, levels = lv),
           cell_line = factor(cell_line, levels = c(CFG$am_bra, names(CFG$am_depmap)))) %>%
    arrange(gene, cell_line) %>%
    transmute(gene, cell_line, cohort, model_id,
              log2_ratio = round(log2r, 3), relative_ratio = round(2^log2r, 3),
              cn_call = as.character(cn_call), bayes_factor = BF,
              dependency_percentile = round(dep_pct, 1))
  wtsv(tab, "per_cellline_CN_and_dependency.tsv", sub = CFG$cn_scheme)
  invisible(tab)
}

# ---- correlation: copy number x Chronos ----
# ---- assemble the matrices for one scope ----
build_scope <- function(scope, ann_all) {
  metric <- "chronos"
  ids <- switch(
    scope,
    melanoma  = ann_all %>% filter(group %in% CFG$melanoma_groups) %>% pull(model_id),
    pancancer = ann_all$model_id)
  ch <- load_chronos()
  if (!is.matrix(ch)) stop("Chronos matrix unavailable")
  ids <- intersect(ids, colnames(ch))
  cn <- load_depmap_cn(ids)
  ids <- intersect(ids, colnames(cn))
  if (length(ids) < CFG$min_models_per_metric) {
    say("[scope:", scope, "] skipped: only ", length(ids),
        " models with both dependency and copy number")
    return(NULL)
  }
  shared <- intersect(rownames(ch), rownames(cn))
  A <- ch[shared, ids, drop = FALSE]
  B <- cn[shared, ids, drop = FALSE]
  models <- ann_all %>% filter(model_id %in% ids) %>%
    transmute(model_id, cell_line = CellLineName, group,
              lineage = OncotreeLineage)

  n_models <- ncol(A)
  if (n_models <= CFG$complete_cases_max_models) {
    ok <- complete.cases(A) & complete.cases(B)
    rule <- sprintf("complete data in all %d cell lines", n_models)
  } else {
    pair_n <- rowSums(is.finite(A) & is.finite(B))
    min_n <- max(CFG$min_models_per_gene,
                 ceiling(CFG$min_frac_models * n_models))
    ok <- pair_n >= min_n
    rule <- sprintf("data in at least %d of %d cell lines", min_n, n_models)
    say("[scope:", scope, "] gene coverage: median ",
        median(pair_n), " of ", n_models, " models; ",
        sum(pair_n == n_models), " genes complete in all")
  }
  say("[scope:", scope, "] ", n_models, " cell lines | ",
      sum(ok), " genes retained (", rule, ")")
  if (sum(ok) < 1000)
    say("[scope:", scope, "] WARNING: fewer than 1000 genes retained; ",
        "the background null may be a biased gene set")

  keep <- intersect(CFG$hit_genes, rownames(A))

  list(scope = scope, metric = metric,
       dep = A[ok, , drop = FALSE], cn = B[ok, , drop = FALSE],
       dep_keep = A[keep, , drop = FALSE], cn_keep = B[keep, , drop = FALSE],
       models = models, n_models = n_models, rule = rule)
}

# ---- run the correlation for one scope ----
correlation_analysis <- function(S) {
  sub <- file.path("correlation", S$metric, S$scope)
  n <- S$n_models
  small <- n <= CFG$exact_spearman_max_n

  wtsv(S$models, "models.tsv", sub = sub)

  genes <- intersect(CFG$hit_genes, rownames(S$dep))
  res <- map_dfr(genes, function(g) {
    x <- S$cn[g, ]; y <- S$dep[g, ]
    st <- spearman_exact(x, y)
    tibble(gene = g, n_models = st$n, spearman_r = st$r, p = st$p)
  }) %>%
    mutate(p_adj = p.adjust(p, CFG$fdr_method))

  min_n <- if (n <= CFG$complete_cases_max_models) n else
    max(CFG$min_models_per_gene, ceiling(CFG$min_frac_models * n))
  rs <- row_spearman_pairwise(S$dep, S$cn, min_n = min_n)
  bg <- tibble(gene = rownames(S$dep), spearman_r = rs$r,
               n_models = rs$n) %>%
    filter(is.finite(spearman_r))
  wtsv(bg, "background_null_r.tsv", sub = sub)
  say("[corr:", S$scope, "/", S$metric, "] ", "genome-wide",
      " null built from ", nrow(bg), " genes | ",
      "models per gene: median ", median(bg$n_models),
      ", range ", min(bg$n_models), "-", max(bg$n_models))

  res <- res %>%
    mutate(background_percentile = 100 * ecdf(bg$spearman_r)(spearman_r),
           empirical_p = map_dbl(spearman_r,
                                 ~ mean(abs(bg$spearman_r) >= abs(.x),
                                        na.rm = TRUE)),
           background_n = nrow(bg),
           empirical_p_adj = p.adjust(empirical_p, CFG$fdr_method),
           interpretation = case_when(
             is.na(spearman_r) ~ "not evaluable",
             small ~ "descriptive only (n too small for inference)",
             p_adj >= 0.05 ~ "no relationship (FDR n.s.)",
             (S$metric == "chronos") == (spearman_r < 0) ~
               "dependency increases with copy number (FDR significant)",
             TRUE ~ "dependency decreases with copy number (FDR significant)")) %>%
    arrange(spearman_r)
  wtsv(res, "hit_correlation.tsv", sub = sub)

  say("[corr:", S$scope, "/", S$metric, "] median ", "genome-wide",
      " r = ",
      sprintf("%.3f", median(bg$spearman_r)),
      " | IQR ", sprintf("%.3f", quantile(bg$spearman_r, .25)), " to ",
      sprintf("%.3f", quantile(bg$spearman_r, .75)))
  if (small)
    say("[corr:", S$scope, "/", S$metric, "] n = ", n,
        "; smallest attainable two-sided exact p = ",
        sprintf("%.4f", 2 / factorial(n)),
        " - results are descriptive, not inferential")
  else
    say("[corr:", S$scope, "/", S$metric, "] ",
        sum(res$p_adj < CFG$fdr_alpha, na.rm = TRUE), "/",
        nrow(res), " hits significant at FDR < CFG$fdr_alpha")

  list(scope = S$scope, metric = S$metric, res = res, bg = bg,
       n_models = n, small = small, sub = sub, rule = S$rule)
}

# ---- FIGURES ----

PAL <- c(Del = "#2166ac", Loss = "#92c5de", Neutral = "gray60",
         Gain = "#f4a582", Amp = "#d6604d")

figures_cn_calls <- function(long) {
  scheme <- CFG$cn_scheme
  lv <- rev(c(CFG$focus_genes, setdiff(CFG$hit_genes, CFG$focus_genes)))
  d1 <- long %>%
    filter(gene %in% CFG$hit_genes, !is.na(cn_call)) %>%
    mutate(gene = factor(gene, levels = lv),
           cell_line = factor(cell_line,
                              levels = c(CFG$am_bra, names(CFG$am_depmap))))
  tile <- function(d, title) ggplot(d, aes(cell_line, gene, fill = cn_call)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", log2r)), size = 2.9) +
    scale_fill_manual(values = PAL, name = "CN call", drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  sfig(tile(d1, "Copy-number status of CRISPR hits"),
       "CN_call_per_cell_line.pdf", 8, 7, sub = scheme)
  sfig(tile(droplevels(filter(d1, gene == CFG$detail_gene)),
            sprintf("Copy-number status of %s", CFG$detail_gene)),
       sprintf("CN_call_%s.pdf", CFG$detail_gene), 6, 2.2, sub = scheme)
}

FDR_LABEL <- sprintf("FDR < %s", CFG$fdr_alpha)

figure_hit_correlation <- function(corr) {
  panel <- function(R, tag) {
    d <- R$res %>% filter(is.finite(spearman_r)) %>%
      mutate(gene = reorder(gene, spearman_r),
             sig = factor(ifelse(!is.na(p_adj) & p_adj < CFG$fdr_alpha, "FDR < CFG$fdr_alpha", "n.s."),
                          levels = c("n.s.", "FDR < CFG$fdr_alpha")))
    ggplot(d, aes(spearman_r, gene, colour = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "gray55") +
      geom_point(size = 3) +
      scale_colour_manual(values = c(`n.s.` = "gray55", `FDR < CFG$fdr_alpha` = "#d6604d"),
                          name = NULL, drop = FALSE) +
      labs(title = sprintf("%s (n = %d cell lines)", tag, R$n_models),
           x = expression(paste("Spearman ", rho)), y = NULL) +
      theme_classic(base_size = 11)
  }
  p <- panel(corr$pancancer, "A  Pan-cancer") + panel(corr$melanoma, "B  Melanoma") +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  sfig(p, "hit_correlation_chronos.pdf", 10, 5.5, sub = "correlation")
}

figure_gene_scatter <- function(scopes, corr, gene = CFG$detail_gene) {
  labels <- setNames(names(CFG$am_depmap), unname(CFG$am_depmap))
  panel <- function(sp, tag) {
    S <- scopes[[sp]]
    d <- tibble(model_id = colnames(S$dep_keep), dep = S$dep_keep[gene, ],
                log2r = S$cn_keep[gene, ]) %>%
      filter(is.finite(dep), is.finite(log2r)) %>%
      mutate(cn_call = factor(as.character(cn_call_log2(log2r, CFG$cn_harmonised)),
                              levels = setdiff(names(HIT_CN_PAL), "No CN data")),
             label = unname(labels[model_id]))
    r <- corr[[sp]]$res %>% filter(gene == !!gene)
    lab <- filter(d, !is.na(label))
    ggplot(d, aes(log2r, dep)) +
      geom_hline(yintercept = CFG$common_essential_chronos, linetype = "dashed", colour = "gray40", linewidth = 0.3) +
      geom_point(aes(colour = cn_call), size = 1.2, alpha = 0.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
      geom_point(data = lab, shape = 21, size = 3, colour = "black", fill = NA) +
      geom_text_repel(data = lab, aes(label = label), size = 2.6, seed = 1,
                      min.segment.length = 0, max.overlaps = Inf) +
      scale_colour_manual(values = HIT_CN_PAL, name = "Copy number", drop = FALSE) +
      labs(title = sprintf("%s (n = %d cell lines)", tag, nrow(d)),
           caption = sprintf("rho = %.2f, FDR = %.2g", r$spearman_r, r$p_adj),
           x = sprintf("log2 CN ratio (%s)", gene),
           y = sprintf("%s Chronos gene effect", gene)) +
      theme_classic(base_size = 11)
  }
  p <- panel("pancancer", "A  Pan-cancer") + panel("melanoma", "B  Melanoma") +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  sfig(p, sprintf("%s_cn_vs_chronos.pdf", gene), 10, 5, sub = "correlation")
  wtsv(bind_rows(lapply(names(scopes), function(sp) corr[[sp]]$res %>%
                          filter(gene == !!gene) %>% mutate(scope = sp))),
       sprintf("%s_correlation.tsv", gene), sub = "correlation")
}

HIT_CN_PAL <- c(Del = "#2166ac", Loss = "#92c5de", Neutral = "gray80",
                Gain = "#f4a582", Amp = "#b2182b", `No CN data` = "gray45")

hit_waterfall_data <- function(dep, cn_long, labels, higher_is_dependent, th) {
  d <- as.data.frame(dep) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "model_id", values_to = "dep") %>%
    filter(is.finite(dep)) %>%
    left_join(cn_long, by = c("gene", "model_id")) %>%
    filter(is.finite(log2r)) %>%
    mutate(cn_call = as.character(cn_call_log2(log2r, th)),
           cn_call = factor(replace_na(cn_call, "No CN data"), levels = names(HIT_CN_PAL)),
           label = unname(labels[model_id]),
           score = if (higher_is_dependent) dep else -dep) %>%
    group_by(gene) %>%
    mutate(rank = rank(score, ties.method = "first"),
           n_models = n(),
           pct_of_lines_less_dependent = 100 * (rank(score, ties.method = "average") - 1) /
             pmax(n() - 1, 1)) %>%
    ungroup()
  d$gene <- factor(d$gene, levels = intersect(CFG$hit_genes, unique(d$gene)))
  d
}

hit_waterfall_plot <- function(d, y_lab, title, ref_line = NULL) {
  lab <- d %>% filter(!is.na(label))
  p <- ggplot(d, aes(rank, dep)) +
    geom_hline(yintercept = 0, colour = "gray60", linewidth = 0.3)
  if (!is.null(ref_line))
    p <- p + geom_hline(yintercept = ref_line, colour = "gray40",
                        linewidth = 0.3, linetype = "dashed")
  p +
    geom_segment(aes(xend = rank, y = 0, yend = dep, colour = cn_call),
                 linewidth = 0.25) +
    geom_point(data = lab, shape = 21, size = 1.6, fill = "black",
               colour = "white", stroke = 0.3) +
    geom_text_repel(data = lab, aes(label = label), size = 2.4, seed = 1,
                    min.segment.length = 0, box.padding = 0.3,
                    max.overlaps = Inf) +
    scale_colour_manual(values = HIT_CN_PAL, name = "Copy number",
                        drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.01)) +
    facet_wrap(~ gene, ncol = 5, scales = "free") +
    labs(title = title, x = "Cell lines, ranked (most dependent on the right)",
         y = y_lab) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(size = 11, face = "bold"),
          plot.caption = element_text(size = 8, hjust = 0),
          legend.position = "bottom") +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)))
}

hit_waterfalls <- function() {
  sub <- "hit_waterfalls"
  genes <- CFG$hit_genes
  labels <- setNames(names(CFG$am_depmap), unname(CFG$am_depmap))
  ch <- load_chronos()
  if (!is.matrix(ch)) stop("Chronos matrix unavailable")
  dep <- ch[intersect(genes, rownames(ch)), , drop = FALSE]
  cn <- load_depmap_cn(colnames(dep))
  cn_l <- as.data.frame(cn[intersect(genes, rownames(cn)), , drop = FALSE]) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "model_id", values_to = "log2r")
  d <- hit_waterfall_data(dep, cn_l, labels, higher_is_dependent = FALSE,
                          th = CFG$cn_harmonised)
  caption <- paste0(
    "Each bar is one cell line. Colour: copy-number call for that gene ",
    "(Del <= ", CFG$cn_harmonised$deep_del, ", Loss < ", CFG$cn_harmonised$loss,
    ", Gain >= ", CFG$cn_harmonised$gain, ", Amp >= ", CFG$cn_harmonised$amp,
    " log2). Labelled: acral lines. Dashed line: -1, median of common essential genes.")
  p <- hit_waterfall_plot(
    d, "Chronos gene effect",
    sprintf("Dependency on the AM CRISPR hits across %s DepMap cell lines",
            format(max(d$n_models), big.mark = ",")), ref_line = CFG$common_essential_chronos) +
    labs(caption = caption)
  sfig(p, "hit_genes_waterfall_chronos.pdf", 16, 10, sub = sub)
  d1 <- d %>% filter(gene == CFG$detail_gene) %>% mutate(gene = droplevels(gene))
  p1 <- hit_waterfall_plot(
    d1, sprintf("%s Chronos gene effect", CFG$detail_gene),
    sprintf("Dependency on %s across %s DepMap cell lines", CFG$detail_gene,
            format(max(d1$n_models), big.mark = ",")), ref_line = CFG$common_essential_chronos) +
    facet_wrap(~ gene, ncol = 1, scales = "free") + labs(caption = caption)
  sfig(p1, sprintf("%s_waterfall_chronos.pdf", CFG$detail_gene), 7, 4, sub = sub)

  acral <- d %>% filter(!is.na(label)) %>%
    transmute(gene = as.character(gene), cell_line = label, model_id,
              dependency = dep, rank, n_models,
              pct_of_lines_less_dependent = round(pct_of_lines_less_dependent, 1),
              cn_log2 = log2r, cn_call = as.character(cn_call)) %>%
    arrange(gene, cell_line)
  wtsv(acral, "acral_lines_rank_chronos.tsv", sub = sub)
  wtsv(d %>% count(gene, cn_call, name = "n_models") %>%
         mutate(gene = as.character(gene), cn_call = as.character(cn_call)),
       "cn_call_counts_chronos.tsv", sub = sub)
  invisible(d)
}

# ---- MAIN ----

main <- function() {
  say("=== CRISPR dependency vs copy number ===")
  bf  <- load_bayes()
  ann <- load_annotation(bf)
  cn_bra <- load_am_bra_cn()
  cn_dm <- load_depmap_cn(unname(CFG$am_depmap))
  long <- build_long(bf, cn_bra, cn_dm)
  wtsv(long %>% filter(gene %in% CFG$hit_genes), "long_matrix_hits.tsv")
  regression_test(long)

  ann_full <- attr(ann, "all_models_full")
  scopes <- list(); corr <- list()
  for (sp in CFG$corr_scopes) {
    scopes[[sp]] <- build_scope(sp, ann_full)
    if (is.null(scopes[[sp]])) stop("scope ", sp, " could not be built")
    corr[[sp]] <- correlation_analysis(scopes[[sp]])
  }
  figure_hit_correlation(corr)
  figure_gene_scatter(scopes, corr)
  hit_waterfalls()

  call_tables(long)
  figures_cn_calls(long)

  dir.create(DIRS$qc, recursive = TRUE, showWarnings = FALSE)
  writeLines(.LOG, file.path(DIRS$qc, "run.log"))
  say("=== done -> ", normalizePath(DIRS$root), " ===")
}

if (!interactive()) {
  main()
  module_finish(M)
}
