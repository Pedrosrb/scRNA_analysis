
# ---- generic readers ----

read_sheet_raw <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop("file not found: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path)))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("readxl is required to read ", basename(path), call. = FALSE)
    x <- readxl::read_excel(path, sheet = sheet %||% 1, col_names = FALSE,
                            col_types = "text", .name_repair = "minimal")
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- utils::read.delim(path, sep = if (ext == "csv") "," else "\t",
                           header = FALSE, colClasses = "character",
                           quote = if (ext == "csv") "\"" else "",
                           na.strings = character(0), fill = TRUE,
                           comment.char = "", check.names = FALSE,
                           blank.lines.skip = TRUE)
  }
  names(x) <- paste0("V", seq_along(x))
  x[] <- lapply(x, function(v) trimws(as.character(v)))
  x
}

promote_header <- function(raw, first_col) {
  hit <- which(tolower(raw[[1]]) == tolower(first_col))
  if (!length(hit))
    stop("header row starting with '", first_col, "' not found", call. = FALSE)
  h <- hit[1]
  nm <- as.character(unlist(raw[h, ]))
  body <- raw[-seq_len(h), , drop = FALSE]
  keep_col <- !is.na(nm) & nzchar(nm)
  body <- body[, keep_col, drop = FALSE]
  names(body) <- make.unique(nm[keep_col])
  empty <- apply(body, 1, function(r) all(is.na(r) | !nzchar(r)))
  rownames(body) <- NULL
  body[!empty, , drop = FALSE]
}

parse_num <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x %in% c("", "NA", "NS", "na", "-")] <- NA
  comma <- !is.na(x) & grepl(",", x) & !grepl("\\.", x)
  x[comma] <- sub(",", ".", x[comma], fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

am_patient <- function(x) {
  out <- sub("^\\s*(AM\\d+).*$", "\\1", x)
  out[!grepl("^\\s*AM\\d+", x)] <- NA_character_
  out
}

am_lesion <- function(x) {
  out <- sub("^\\s*(AM\\d+[a-z]?).*$", "\\1", x)
  out[!grepl("^\\s*AM\\d+", x)] <- NA_character_
  out
}

strip_hum <- function(x) sub("_hum$", "", trimws(x))

# ---- S13: cellularity, ploidy, sample type ----

bra_samples <- function(path, min_cellularity = 0.35) {
  d <- promote_header(read_sheet_raw(path), "Sample")
  need <- c("Sample", "Lab_id", "Sample_type", "Cellularity", "Ploidy")
  if (length(miss <- setdiff(need, names(d))))
    stop(basename(path), ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  use_col <- grep("^will we use", names(d), ignore.case = TRUE, value = TRUE)[1]
  out <- data.frame(
    sample_id   = strip_hum(d$Sample),
    lab_id      = d$Lab_id,
    patient     = am_patient(d$Lab_id),
    lesion      = am_lesion(d$Lab_id),
    sample_type = d$Sample_type,
    histology   = d$Tissue_histology %||% NA_character_,
    cellularity = parse_num(d$Cellularity),
    ploidy      = parse_num(d$Ploidy),
    marked_use  = if (is.na(use_col)) NA else tolower(d[[use_col]]) == "yes",
    stringsAsFactors = FALSE)
  out$is_patient <- out$sample_type == "Patient"

  bad <- !is.na(out$cellularity) & (out$cellularity < 0 | out$cellularity > 1)
  if (any(bad & out$is_patient))
    message("[S13] patient cellularity outside [0, 1], set to NA and excluded: ",
            paste(out$lab_id[bad & out$is_patient], collapse = ", "))
  out$cellularity[bad] <- NA
  out$analysable <- out$is_patient & !is.na(out$cellularity) &
    out$cellularity >= min_cellularity & (is.na(out$marked_use) | out$marked_use)
  disagree <- out$is_patient & !is.na(out$marked_use) & !is.na(out$cellularity) &
    out$marked_use != (out$cellularity >= min_cellularity)
  if (any(disagree))
    message("[S13] 'Will we use it?' disagrees with cellularity >= ", min_cellularity,
            " for: ", paste(out$sample_id[disagree], collapse = ", "),
            " (excluded unless both pass)")
  if (anyDuplicated(out$sample_id))
    stop("duplicated sample IDs in ", basename(path), call. = FALSE)
  out
}

# ---- S1 + S2: sex, age, clinical stage ----

bra_clinical <- function(s1_path, s2_path, early_stages, late_stages,
                         stage_from_metastasis = TRUE) {
  s1 <- promote_header(read_sheet_raw(s1_path), "Sample_ID")
  s2 <- promote_header(read_sheet_raw(s2_path), "Sample_ID")
  s2 <- s2[grepl("^AM\\d+", s2$Sample_ID), , drop = FALSE]

  sex <- toupper(substr(s1$Sex, 1, 1))
  pat <- data.frame(
    patient = am_patient(s1$Sample_ID),
    sex = factor(ifelse(sex == "F", "Female", ifelse(sex == "M", "Male", NA)),
                 levels = c("Female", "Male")),
    age = parse_num(s1$Age),
    stringsAsFactors = FALSE)
  if (anyDuplicated(pat$patient)) stop("duplicated patients in S1", call. = FALSE)

  st <- data.frame(patient = am_patient(s2$Sample_ID), lesion = s2$Sample_ID,
                   tumor_type = s2$Tumor_type,
                   stage = ifelse(s2$Clinical_Stage %in% c("", "NA"), NA, s2$Clinical_Stage),
                   stringsAsFactors = FALSE)
  unknown <- setdiff(na.omit(st$stage), c(early_stages, late_stages))
  if (length(unknown))
    message("[S2] Clinical_Stage values not in early/late sets (treated as NA): ",
            paste(unknown, collapse = ", "))
  grp <- ifelse(st$stage %in% late_stages, "III-IV",
                ifelse(st$stage %in% early_stages, "I-II", NA))
  by_pat <- tapply(grp, st$patient, function(g)
    if (any(g == "III-IV", na.rm = TRUE)) "III-IV"
    else if (any(g == "I-II", na.rm = TRUE)) "I-II" else NA_character_)
  has_met <- tapply(grepl("metastasis", st$tumor_type, ignore.case = TRUE),
                    st$patient, any)

  pat$stage_group  <- as.vector(by_pat[pat$patient])
  pat$stage_source <- ifelse(is.na(pat$stage_group), NA, "primary")
  if (stage_from_metastasis) {
    fill <- is.na(pat$stage_group) & pat$patient %in% names(has_met)[has_met]
    pat$stage_group[fill]  <- "III-IV"
    pat$stage_source[fill] <- "nodal_or_in_transit_metastasis"
  }
  pat$stage_group <- factor(pat$stage_group, levels = c("I-II", "III-IV"))
  pat
}

# ---- genotyping (ADMIXTURE K5) ----

bra_ancestry <- function(path) {
  d <- promote_header(read_sheet_raw(path), "Id lab")
  ks <- paste0("K", 1:5)
  if (length(miss <- setdiff(c("Sample id", "Sample type", ks), names(d))))
    stop(basename(path), ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  d <- d[grepl("^AM\\d+", d[["Id lab"]]), , drop = FALSE]
  failed <- grepl("\\*", d[["Sample id"]])
  if (any(failed))
    message("[genotyping] failed QC, dropped: ", paste(d[["Id lab"]][failed], collapse = ", "))
  d <- d[!failed & tolower(d[["Sample type"]]) == "blood", , drop = FALSE]
  out <- data.frame(patient = am_patient(d[["Id lab"]]),
                    sample_id = d[["Sample id"]], stringsAsFactors = FALSE)
  for (k in ks) out[[k]] <- parse_num(d[[k]])
  sums <- rowSums(out[ks])
  if (any(abs(sums - 1) > 0.01))
    stop("ancestry proportions do not sum to 1 for: ",
         paste(out$patient[abs(sums - 1) > 0.01], collapse = ", "), call. = FALSE)
  if (anyDuplicated(out$patient)) stop("duplicated patients in genotyping", call. = FALSE)
  out
}

# ---- BRAF ----

bra_braf <- function(path) {
  d <- promote_header(read_sheet_raw(path), "Samples")
  d <- d[grepl("^AM\\d+", d$Samples), , drop = FALSE]
  norm <- function(x) {
    x <- toupper(gsub("[[:space:]]", "", x))
    ifelse(x %in% c("", "NA", "NS"), NA, ifelse(x %in% c("WT", "WILDTYPE"), "WT", x))
  }
  data.frame(lesion = d$Samples, patient = am_patient(d$Samples),
             sanger = norm(d[[2]]), wes = norm(d[[3]]), stringsAsFactors = FALSE)
}

# ---- CNApp scores ----

bra_cnapp <- function(path, scores = c("BCS", "FCS", "GCS")) {
  d <- promote_header(read_sheet_raw(path), "AM_id")
  if (length(miss <- setdiff(c("Sanger_id", scores), names(d))))
    stop(basename(path), ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  out <- data.frame(sample_id = strip_hum(d$Sanger_id), am_id = d$AM_id,
                    stringsAsFactors = FALSE)
  for (s in scores) out[[s]] <- parse_num(d[[s]])
  out
}

# ---- hailstorm calls ----

bra_hailstorm_calls <- function(path, pass_values = c("Pass", "Manual_Pass")) {
  d <- promote_header(read_sheet_raw(path), "Sample_ID")
  if (length(miss <- setdiff(c("Sanger_DNA_ID", "Arm", "Filter"), names(d))))
    stop(basename(path), ": missing column(s) ", paste(miss, collapse = ", "), call. = FALSE)
  unknown <- setdiff(unique(d$Filter), pass_values)
  if (length(unknown))
    message("[hailstorm] Filter values not counted as positive: ",
            paste(unknown, collapse = ", "))
  d <- d[d$Filter %in% pass_values, , drop = FALSE]
  unique(data.frame(sample_id = strip_hum(d$Sanger_DNA_ID), am_id = d$Sample_ID,
                    arm = sub("^chr", "", d$Arm), stringsAsFactors = FALSE))
}

hailstorm_status <- function(calls, universe, arms = sort(unique(calls$arm))) {
  outside <- setdiff(unique(calls$sample_id), universe)
  if (length(outside))
    message("[hailstorm] calls on ", length(outside), " sample(s) outside the ",
            "evaluated universe ignored: ", paste(outside, collapse = ", "))
  g <- expand.grid(sample_id = universe, arm = arms, stringsAsFactors = FALSE)
  key <- paste(calls$sample_id, calls$arm)
  g$present <- as.integer(paste(g$sample_id, g$arm) %in% key)
  g
}

# ---- GISTIC ----

bra_gistic_long <- function(path, q_max = 0.1, min_level = 1) {
  raw <- utils::read.delim(path, header = TRUE, check.names = FALSE,
                           colClasses = "character", quote = "",
                           comment.char = "", strip.white = TRUE)
  raw <- raw[, nzchar(names(raw)) & !grepl("^V\\d+$|^\\.\\.\\.", names(raw)), drop = FALSE]
  name <- trimws(raw[["Unique Name"]])
  thr <- raw[nzchar(name) & !grepl("CN values", name), , drop = FALSE]
  fixed <- c("Unique Name", "Descriptor", "Wide Peak Limits", "Peak Limits",
             "Region Limits", "q values",
             "Residual q values after removing segments shared with higher peaks",
             "Broad or Focal", "Amplitude Threshold")
  samples <- setdiff(names(thr), fixed)
  samples <- samples[grepl("^PD\\d+", samples)]
  if (!length(samples)) stop("no PD sample columns in ", basename(path), call. = FALSE)

  type <- ifelse(grepl("^Amplification", thr[["Unique Name"]]), "Amp", "Del")
  feature <- make.unique(paste0(type, "_", trimws(thr$Descriptor)))
  q <- parse_num(thr[["q values"]])
  keep <- !is.na(q) & q < q_max
  thr <- thr[keep, , drop = FALSE]; feature <- feature[keep]; q <- q[keep]

  long <- do.call(rbind, lapply(seq_len(nrow(thr)), function(i) {
    lv <- parse_num(unlist(thr[i, samples]))
    data.frame(sample_id = samples, feature = feature[i], q_value = q[i],
               level = lv, present = as.integer(lv >= min_level),
               stringsAsFactors = FALSE)
  }))
  rownames(long) <- NULL
  long
}
