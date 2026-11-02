# ---- Two-group comparison shared by the CNApp modules (base R) ----

wilcox_two_group <- function(values, group, feature, group_levels, min_n = 3) {
  keep <- !is.na(values) & !is.na(group) & as.character(group) %in% group_levels
  values <- values[keep]
  group  <- as.character(group[keep])
  v1 <- values[group == group_levels[1]]
  v2 <- values[group == group_levels[2]]
  if (length(v1) < min_n || length(v2) < min_n) {
    message(sprintf("  Skipped %-26s (n = %d vs %d < min_n)", feature,
                    length(v1), length(v2)))
    return(NULL)
  }
  wt <- tryCatch(
    suppressWarnings(wilcox.test(v2, v1, exact = FALSE, correct = TRUE,
                                 conf.int = TRUE)),
    error = function(e)
      suppressWarnings(wilcox.test(v2, v1, exact = FALSE, correct = TRUE)))
  auc <- mean(outer(v2, v1, ">") + 0.5 * outer(v2, v1, "=="))
  q <- function(v, p) unname(quantile(v, p))
  data.frame(
    feature = feature,
    group_ref = group_levels[1], group_comp = group_levels[2],
    n_ref = length(v1), n_comp = length(v2),
    median_ref = median(v1), median_comp = median(v2),
    iqr_ref  = sprintf("%.2f-%.2f", q(v1, .25), q(v1, .75)),
    iqr_comp = sprintf("%.2f-%.2f", q(v2, .25), q(v2, .75)),
    delta_median  = median(v2) - median(v1),
    rank_biserial = 2 * auc - 1,
    hl_shift   = if (is.null(wt$estimate)) NA_real_ else unname(wt$estimate),
    hl_ci_low  = if (is.null(wt$conf.int)) NA_real_ else wt$conf.int[1],
    hl_ci_high = if (is.null(wt$conf.int)) NA_real_ else wt$conf.int[2],
    p_value = unname(wt$p.value),
    stringsAsFactors = FALSE)
}

