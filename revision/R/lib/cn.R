# ---- Copy-number helpers shared by every module (base R) ----

CN_LEVELS <- c("Del", "Loss", "Neutral", "Gain", "Amp")

cn_call_log2 <- function(log2r, th) {
  stopifnot(th$deep_del < th$loss, th$loss < 0, 0 < th$gain, th$gain < th$amp)
  out <- rep(NA_character_, length(log2r))
  ok  <- !is.na(log2r)
  out[ok & log2r <= th$deep_del]                   <- "Del"
  out[ok & log2r >  th$deep_del & log2r < th$loss] <- "Loss"
  out[ok & log2r >= th$loss     & log2r < th$gain] <- "Neutral"
  out[ok & log2r >= th$gain     & log2r < th$amp]  <- "Gain"
  out[ok & log2r >= th$amp]                        <- "Amp"
  factor(out, levels = CN_LEVELS)
}

cn_call_integer <- function(cn_int, th) {
  out <- rep(NA_character_, length(cn_int))
  ok  <- !is.na(cn_int)
  out[ok & cn_int >= th$amp]     <- "Amp"
  out[ok & cn_int == th$gain]    <- "Gain"
  out[ok & cn_int == th$neutral] <- "Neutral"
  out[ok & cn_int == th$loss]    <- "Loss"
  out[ok & cn_int <= th$del]     <- "Del"
  factor(out, levels = CN_LEVELS)
}

is_gain_or_amp <- function(call) as.character(call) %in% c("Gain", "Amp")

linear_to_log2 <- function(v, floor_linear) log2(pmax(v, floor_linear))

detect_cn_scale <- function(v, min_frac_neg = 0.01) {
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_character_)
  if (mean(v < 0) > min_frac_neg) "log2" else "linear_ratio"
}
