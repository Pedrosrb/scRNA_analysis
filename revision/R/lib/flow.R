# ---- Sample flow: how many samples/patients enter and leave each step. ----

.FLOW <- new.env(parent = emptyenv())
.FLOW$rows <- list()

flow_reset <- function() .FLOW$rows <- list()

flow_add <- function(step, n, unit = "patients", note = "") {
  .FLOW$rows[[length(.FLOW$rows) + 1]] <- data.frame(
    order = length(.FLOW$rows) + 1, step = step, unit = unit,
    n = as.integer(n), note = note, stringsAsFactors = FALSE)
  message(sprintf("[flow] %-40s %6d %s%s", step, as.integer(n), unit,
                  if (nzchar(note)) paste0("  (", note, ")") else ""))
  invisible(NULL)
}

flow_table <- function()
  if (length(.FLOW$rows)) do.call(rbind, .FLOW$rows) else NULL
