# ---- DepMap helpers (base R; data.table used when installed) ----

read_depmap <- function(path) {
  if (!file.exists(path)) stop("DepMap file not found: ", path, call. = FALSE)
  if (requireNamespace("data.table", quietly = TRUE))
    data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  else
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

melanoma_subtype <- function(subtype, primary_disease) {
  s <- tolower(trimws(as.character(subtype)))
  d <- tolower(trimws(as.character(primary_disease)))
  is_mel <- grepl("melanoma", d) | grepl("melanoma", s)
  out <- ifelse(grepl("acral", s), "Acral",
         ifelse(grepl("uveal", s), "Uveal",
         ifelse(grepl("mucosal", s), "Mucosal",
         ifelse(grepl("cutaneous", s), "Cutaneous",
         ifelse(s %in% "melanoma", "Melanoma NOS", "Other")))))
  out[!is_mel | is.na(s)] <- NA_character_
  out
}
