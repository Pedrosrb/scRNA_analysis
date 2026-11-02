# ---- Plot helpers (ggplot2 is resolved only when these are called) ----

geom_jitter_fixed <- function(..., width = 0.1, height = 0, seed = 1) {
  ggplot2::geom_point(..., position = ggplot2::position_jitter(
    width = width, height = height, seed = seed))
}
