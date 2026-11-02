# =============================================================================
# Copy-number view of one DepMap cell line at AM driver loci
# Run:      Rscript run.R cn_profile
# Outputs:
#   figures/<model>_genome.pdf       Figure R12
#   tables/<model>_driver_genes.tsv
# Inputs:   config/params.R -> PARAMS$modules$cn_profile$inputs
#   model         DepMap Model.csv
#   crispr        DepMap CRISPRGeneEffect.csv
#   cn_gene       DepMap OmicsCNGeneWGS.csv
#   cn_segments   DepMap OmicsCNSegmentsWGS.csv
#   signatures    DepMap OmicsGlobalSignatures.csv (ModelID, Ploidy)
# =============================================================================

# ---- bootstrap ----
if (!nzchar(Sys.getenv("AMREV_ROOT"))) local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  Sys.setenv(AMREV_ROOT = if (length(f)) dirname(dirname(normalizePath(f[1]))) else getwd())
})
source(file.path(Sys.getenv("AMREV_ROOT"), "R", "setup.R"))
M <- module_init("cn_profile")

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggrepel)
})

p  <- M$params
hb <- p$cn$harmonised
CFG <- list(
  model       = p$model,
  gene_coords = p$gene_coords,
  out_dir     = M$dir,
  drivers     = p$drivers,
  gain_log2 = hb$gain, amp_log2 = hb$amp, loss_log2 = hb$loss, del_log2 = hb$deep_del,
  floor_linear = p$cn$floor_linear,
  label_only_gained   = p$label_only_gained,
  default_ploidy = p$default_ploidy,
  chromosomes = p$chromosomes,
  chromosomes_all = p$chromosomes_all
)

# ---- Helpers ----

say <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", ...)
norm_chr <- function(x) paste0("chr", sub("^chr", "", as.character(x)))
strip_entrez <- function(x) trimws(sub("\\s*\\(\\d+\\)\\s*$", "", x))
to_log2  <- function(v) linear_to_log2(v, CFG$floor_linear)
classify <- function(l) as.character(cn_call_log2(l, p$cn$harmonised))

dir.create(file.path(CFG$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CFG$out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
wtsv <- function(x, name) {
  readr::write_tsv(x, file.path(CFG$out_dir, "tables", name), na = "")
  say("[table] ", name)
}
sfig <- function(p, name, w, h) {
  ggsave(file.path(CFG$out_dir, "figures", name), p, width = w, height = h)
  say("[figure] ", name)
}

# ---- Inputs ----

resolve_model <- function(name) {
  m <- as.data.frame(fread(M$inputs$model))
  key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
  cols <- intersect(c("ModelID", "CellLineName", "StrippedCellLineName"), names(m))
  hit <- m[Reduce(`|`, lapply(cols, function(cc) key(m[[cc]]) == key(name))), ,
           drop = FALSE]
  if (nrow(hit) != 1)
    stop("model '", name, "' matched ", nrow(hit), " rows in Model.csv")
  say("[model] ", name, " -> ", hit$ModelID,
      if ("OncotreeSubtype" %in% names(hit)) paste0(" (", hit$OncotreeSubtype, ")"))
  hit
}

read_model_row <- function(file, model_id) {
  x <- as.data.frame(fread(file))
  idcol <- if ("ModelID" %in% names(x)) "ModelID" else names(x)[1]
  if ("IsDefaultEntryForModel" %in% names(x))
    x <- x[x$IsDefaultEntryForModel %in% c("Yes", "TRUE", TRUE), , drop = FALSE]
  gene_cols <- grep("\\(\\d+\\)$", names(x), value = TRUE)
  row <- x[x[[idcol]] == model_id, , drop = FALSE]
  if (!nrow(row)) { say("[", file, "] no row for ", model_id); return(NULL) }
  list(this = tibble(gene = strip_entrez(gene_cols),
                     entrez = sub(".*\\((\\d+)\\)$", "\\1", gene_cols),
                     value = as.numeric(unlist(row[1, gene_cols]))),
       all = x, gene_cols = gene_cols)
}

load_gene_coords <- function(src, ids) {
  std <- CFG$chromosomes_all
  if (identical(src, "txdb")) {
    if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
      stop("install with BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')")
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    g <- suppressMessages(GenomicFeatures::genes(txdb,
                                                 single.strand.genes.only = FALSE))
    g <- unlist(g)
    gr <- tibble(entrez = names(g),
                 chrom = as.character(GenomicRanges::seqnames(g)),
                 start = GenomicRanges::start(g),
                 end = GenomicRanges::end(g)) %>%
      filter(chrom %in% std) %>%
      group_by(entrez, chrom) %>%
      summarise(start = min(start), end = max(end), .groups = "drop") %>%
      group_by(entrez) %>% slice_max(end - start, n = 1, with_ties = FALSE) %>%
      ungroup()
    out <- ids %>% inner_join(gr, by = "entrez") %>% select(gene, chrom, start, end)
    say("[coords] TxDb hg38: ", nrow(out), " of ", nrow(ids),
        " DepMap genes placed by Entrez ID")
    return(out)
  }
  if (!file.exists(src)) stop("gene coordinate file not found: ", src)
  if (grepl("\\.cns$", src))
    stop(".cns gene coordinates are segment-wide; use 'txdb' or a GTF")
  g <- fread(src, sep = "\t", header = FALSE, skip = "chr")
  g <- g[V3 == "gene"]
  tibble(chrom = norm_chr(g$V1), start = g$V4, end = g$V5,
         gene = sub('.*gene_name "([^"]+)".*', "\\1", g$V9)) %>%
    filter(chrom %in% std) %>%
    group_by(gene) %>% slice_max(end - start, n = 1, with_ties = FALSE) %>%
    ungroup()
}

load_segments <- function(model_id) {
  s <- as.data.frame(fread(M$inputs$cn_segments))
  idcol <- intersect(c("ModelID", "DepMap_ID"), names(s))[1]
  if (is.na(idcol))
    stop("segments are not keyed by ModelID; columns: ",
         paste(head(names(s), 8), collapse = ", "))
  s <- s[s[[idcol]] == model_id, , drop = FALSE]
  if ("IsDefaultEntryForModel" %in% names(s))
    s <- s[s$IsDefaultEntryForModel %in% c("Yes", "TRUE", TRUE), , drop = FALSE]
  if ("ProfileID" %in% names(s) && n_distinct(s$ProfileID) > 1) {
    keep <- names(which.max(table(s$ProfileID)))
    say("[segments] ", n_distinct(s$ProfileID), " profiles; keeping ", keep)
    s <- s[s$ProfileID == keep, , drop = FALSE]
  }
  vcol <- intersect(c("RelCopyRatio", "SegmentMean", "CN"), names(s))[1]
  v <- as.numeric(s[[vcol]])
  log2v <- if (identical(detect_cn_scale(v), "log2")) v else to_log2(v)
  out <- tibble(chrom = norm_chr(s$Chromosome), start = as.numeric(s$Start),
                end = as.numeric(s$End), log2 = log2v) %>%
    filter(is.finite(log2))
  say("[segments] ", nrow(out), " segments from ", vcol,
      if (identical(log2v, v)) " (already log2)" else " (linear, converted)")
  out
}

# ---- Plots (layout follows plot_chr_from_cnvkit.R) ----

get_y_scale <- function(v, ploidy, step = 1, pad = 0.3, clamp = c(-2, 6)) {
  v <- v[is.finite(v)]
  lo <- max(floor((min(v) - pad) / step) * step, clamp[1])
  hi <- min(ceiling((max(v) + pad) / step) * step, clamp[2])
  br <- seq(lo, hi, by = step)
  list(limits = c(lo, hi), breaks = br, labels = round(ploidy * 2^br, 1))
}

plot_genome <- function(genes, segs, drivers_tbl, ploidy, title) {
  ord <- CFG$chromosomes
  lens <- genes %>% filter(chrom %in% ord) %>%
    group_by(chrom) %>% summarise(len = max(end), .groups = "drop") %>%
    arrange(factor(chrom, ord)) %>%
    mutate(offset = lag(cumsum(as.numeric(len)), default = 0),
           mid = offset + len / 2, shade = row_number() %% 2 == 0)
  place <- function(d, cs = "start") d %>% inner_join(lens, by = "chrom") %>%
    mutate(gx = (.data[[cs]] + offset) / 1e6)

  g  <- place(genes %>% mutate(start = (start + end) / 2))
  s  <- segs %>% inner_join(lens, by = "chrom") %>%
    mutate(x0 = (start + offset) / 1e6, x1 = (end + offset) / 1e6)
  hl <- place(drivers_tbl %>% filter(is.finite(log2)) %>%
                mutate(start = (start + end) / 2,
                       gained = is_gain_or_amp(call)))
  if (CFG$label_only_gained) hl <- hl %>% filter(gained)
  ys <- get_y_scale(c(g$log2, s$log2), ploidy)

  ggplot() +
    geom_rect(data = lens %>% filter(shade),
              aes(xmin = offset / 1e6, xmax = (offset + len) / 1e6,
                  ymin = -Inf, ymax = Inf), fill = "gray95") +
    geom_hline(yintercept = c(CFG$gain_log2, CFG$amp_log2),
               linetype = "dashed", colour = "gray80", linewidth = 0.3) +
    geom_point(data = g, aes(gx, log2), colour = "gray", size = 0.2) +
    geom_segment(data = s, aes(x = x0, xend = x1, y = log2, yend = log2),
                 colour = "gold", linewidth = 0.7) +
    geom_point(data = hl, aes(gx, log2, colour = gained), size = 1) +
    geom_text_repel(data = hl, aes(gx, log2, label = gene, colour = gained),
                    size = 2.8, min.segment.length = 0, max.overlaps = Inf) +
    scale_colour_manual(values = c(`TRUE` = "darkred", `FALSE` = "gray40"),
                        guide = "none") +
    scale_x_continuous(breaks = lens$mid / 1e6, labels = sub("chr", "", lens$chrom),
                       expand = c(0, 0), name = "Chromosome") +
    scale_y_continuous(limits = ys$limits, breaks = ys$breaks,
                       oob = scales::squish, name = "log2 ratio",
                       sec.axis = sec_axis(~ ., breaks = ys$breaks,
                                           labels = ys$labels,
                                           name = sprintf("Copy number (ploidy %.2f)", ploidy))) +
    ggtitle(title) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
}

# ---- Main ----

model   <- resolve_model(CFG$model)
mid     <- model$ModelID
label   <- if ("CellLineName" %in% names(model)) model$CellLineName else CFG$model

sig <- as.data.frame(fread(M$inputs$signatures))
ploidy <- as.numeric(sig$Ploidy[sig$ModelID == mid][1])
if (!is.finite(ploidy)) { ploidy <- CFG$default_ploidy; say("[ploidy] not found; assuming ", ploidy) }
say("[ploidy] ", label, " = ", sprintf("%.2f", ploidy))

cn <- read_model_row(M$inputs$cn_gene, mid)
if (is.null(cn)) stop("no copy number for ", mid)

coords <- load_gene_coords(CFG$gene_coords, cn$this %>% distinct(gene, entrez))
miss <- setdiff(CFG$drivers, coords$gene)
if (length(miss)) say("[coords] WARNING: no coordinates for ",
                      paste(miss, collapse = ", "))

genes <- cn$this %>% mutate(log2 = to_log2(value)) %>% select(-entrez) %>%
  inner_join(coords, by = "gene") %>% filter(is.finite(log2))

segs <- load_segments(mid)

drivers_tbl <- tibble(gene = CFG$drivers) %>%
  left_join(coords, by = "gene") %>%
  left_join(cn$this %>% transmute(gene, log2 = to_log2(value)), by = "gene") %>%
  mutate(abs_cn = ploidy * 2^log2, call = classify(log2))

drivers_tbl$segment_log2 <- map_dbl(seq_len(nrow(drivers_tbl)), function(i) {
  r <- drivers_tbl[i, ]
  hit <- segs %>% filter(chrom == r$chrom, start <= r$end, end >= r$start)
  if (!nrow(hit)) NA_real_ else hit$log2[which.max(pmin(hit$end, r$end) -
                                                     pmax(hit$start, r$start))]
})

dep <- read_model_row(M$inputs$crispr, mid)
if (!is.null(dep))
  drivers_tbl <- drivers_tbl %>%
  left_join(dep$this %>% select(gene, chronos = value), by = "gene")

drivers_tbl <- drivers_tbl %>% mutate(model = label, ploidy = ploidy, .before = 1)
wtsv(drivers_tbl, paste0(label, "_driver_genes.tsv"))
print(drivers_tbl %>% select(gene, chrom, log2, abs_cn, call, any_of(
  c("segment_log2", "chronos"))),
  n = Inf)

sfig(plot_genome(genes, segs, drivers_tbl, ploidy,
                 sprintf("%s - genome", label)),
     sprintf("%s_genome.pdf", label), 14, 4)

say("done -> ", normalizePath(CFG$out_dir))
module_finish(M)
