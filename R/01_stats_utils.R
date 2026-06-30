# Statistical helpers for the minimal public workflow.

bh_adjust <- function(p) {
  p.adjust(p, method = "BH")
}

spearman_one <- function(data, x, y) {
  keep <- stats::complete.cases(data[, c(x, y)])
  n <- sum(keep)
  if (n < 3) {
    return(data.frame(x = x, y = y, n = n, rho = NA_real_, p_value = NA_real_))
  }

  test <- suppressWarnings(stats::cor.test(data[[x]][keep], data[[y]][keep], method = "spearman"))
  data.frame(
    x = x,
    y = y,
    n = n,
    rho = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

spearman_grid <- function(data, x_vars, y_vars) {
  results <- do.call(
    rbind,
    lapply(x_vars, function(x) {
      do.call(rbind, lapply(y_vars, function(y) spearman_one(data, x, y)))
    })
  )
  results$p_adj <- bh_adjust(results$p_value)
  results
}

fit_linear_model <- function(data, formula) {
  model <- stats::lm(formula, data = data)
  stats::coef(summary(model))
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x - mean(x, na.rm = TRUE))
  as.numeric(scale(x))
}

baseline_delta <- function(data,
                           participant_col = "participant_id",
                           time_col = "time",
                           baseline_time = "W0",
                           value_cols) {
  split_data <- split(data, data[[participant_col]])
  out <- lapply(split_data, function(df) {
    baseline <- df[df[[time_col]] == baseline_time, , drop = FALSE]
    if (nrow(baseline) != 1) return(NULL)
    for (col in value_cols) {
      df[[paste0("delta_", col)]] <- df[[col]] - baseline[[col]][1]
    }
    df
  })
  do.call(rbind, out)
}
