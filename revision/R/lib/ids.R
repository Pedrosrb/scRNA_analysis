# ---- Sample / patient identifiers (base R) ----

strip_region_suffix <- function(x) sub("[a-z]$", "", trimws(x))

lacks_region_suffix <- function(x) !grepl("[a-z]$", trimws(x))

pd_patient <- function(x) sub("^(PD\\d+).*$", "\\1", x)
