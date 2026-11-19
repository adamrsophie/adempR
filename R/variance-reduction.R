# Variance reduction ------------------------------------------------------

# Independent sampling units of a set of replications: with an antithetic
# generator, replications 2k - 1 and 2k form one unit (pair k).
rep_units <- function(study, replication) {
  if (isTRUE(study$generator$antithetic)) (as.integer(replication) + 1L) %/% 2L else
    as.integer(replication)
}

# Bootstrap resample of positions that keeps units together. With
# independent replications it is the ordinary resample.
boot_index <- function(units) {
  n <- length(units)
  if (!anyDuplicated(units)) return(sample.int(n, n, replace = TRUE))
  u <- unique(units)
  pick <- u[sample.int(length(u), length(u), replace = TRUE)]
  pos <- split(seq_len(n), factor(units, levels = u))
  unlist(pos[match(pick, u)], use.names = FALSE)
}

# MCSE of a mean of contributions that are dependent within units (cluster
# robust); with one contribution per unit it is sd(x) / sqrt(n).
cluster_mcse <- function(x, units) {
  n <- length(x)
  G <- length(unique(units))
  if (G < 2L) return(NA_real_)
  s <- tapply(x - mean(x), units, sum)
  sqrt(sum(s^2) * G / (G - 1)) / n
}

# Delete-one-pair jackknife MCSE of a measure for antithetic replications:
# valid for linear and nonlinear measures alike, and equal to
# sd(pair means) / sqrt(pairs) for a mean over complete pairs.
pair_jackknife <- function(m, rows, study, estimand, r) {
  if (!r$status %in% c("ok", "boundary", "unstable") || !is.finite(r$value)) return(r)
  units <- rep_units(study, rows$replication)
  u <- unique(units)
  np <- length(u)
  note <- function(txt) paste(c(if (!is.na(r$note)) r$note, txt), collapse = "; ")
  if (np < 2L) {
    r$mcse <- r$lower <- r$upper <- NA_real_
    r$note <- note("antithetic pairs: fewer than two pairs, no MCSE")
    return(r)
  }
  th <- vapply(u, function(k) {
    v <- measure_cell(m, rows[units != k, , drop = FALSE], study, estimand)$value
    if (is.null(v)) NA_real_ else v
  }, numeric(1))
  if (anyNA(th)) {
    r$mcse <- r$lower <- r$upper <- NA_real_
    r$note <- note("antithetic pairs: the measure is undefined without some pair, no MCSE")
    return(r)
  }
  r$mcse <- sqrt((np - 1) / np * sum((th - mean(th))^2))
  z <- stats::qnorm(1 - (1 - study$level) / 2)
  r$lower <- r$value - z * r$mcse
  r$upper <- r$value + z * r$mcse
  if (identical(m, "coverage") || (is.character(m) && m %in% PROPORTION_MEASURES)) {
    r$lower <- max(0, r$lower)
    r$upper <- min(1, r$upper)
  }
  r$note <- note(sprintf("MCSE from %d antithetic pairs (delete-one-pair jackknife)", np))
  r
}

vr_contributions <- function(res, measure) {
  e <- res$estimates
  if (measure == "coverage") {
    e <- e[e$interval_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(e$truth), , drop = FALSE]
    e$x <- as.numeric(e$lower <= e$truth & e$truth <= e$upper)
  } else {
    e <- e[e$valid & is.finite(e$estimate) & !is.na(e$truth), , drop = FALSE]
    e$x <- switch(measure, bias = e$estimate - e$truth, mse = (e$estimate - e$truth)^2,
                  mean = e$estimate)
  }
  e
}

#' Summaries for antithetic pairs
#'
#' With `sim_generator(antithetic = TRUE)`, replications 2k - 1 and 2k use
#' antithetic draws and are dependent. The pair average of each
#' contribution is the independent unit, so the value is the mean over
#' complete pairs and its MCSE is `sd(pair means) / sqrt(pairs)`. The
#' result also gives the MCSE that independent replications with the same
#' variance would have had, and the ratio of variances (above 1 when the
#' antithetic design helped).
#'
#' @param run A run or results.
#' @param measure `"bias"`, `"mse"`, `"coverage"` or `"mean"` (the mean
#'   estimate).
#' @return A data frame per scenario, method and estimand.
#' @export
sim_summarise_antithetic <- function(run, measure = c("bias", "mse", "coverage", "mean")) {
  measure <- match.arg(measure)
  res <- sim_results(run)
  study <- res$manifest$study
  if (!isTRUE(study$generator$antithetic)) stopf("The generator does not use antithetic pairs.")
  e <- vr_contributions(res, measure)
  out <- list()
  for (grp in split(e, list(e$scenario_key, e$method, e$estimand), drop = TRUE)) {
    pair <- (grp$replication + 1L) %/% 2L
    full <- names(which(table(pair) == 2L))
    g <- grp[as.character(pair) %in% full, , drop = FALSE]
    pm <- tapply(g$x, (g$replication + 1L) %/% 2L, mean)
    np <- length(pm)
    if (np < 2L) next
    mcse <- stats::sd(pm) / sqrt(np)
    naive <- stats::sd(g$x) / sqrt(2 * np)
    out[[length(out) + 1L]] <- data.frame(
      scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
      method = grp$method[1], estimand = grp$estimand[1], measure = measure,
      value = mean(pm), mcse = mcse, pairs = np, mcse_independent = naive,
      variance_ratio = (naive / mcse)^2, stringsAsFactors = FALSE)
  }
  out <- bind_rows(out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

#' Control variate adjustment of a performance measure
#'
#' A quantity recorded for every data set whose mean is known exactly (for
#' example the sample mean of a generated covariate, recorded with
#' `sim_generator(record = )`) can reduce the Monte Carlo error of a
#' measure correlated with it. With contributions `x` and control values
#' `c` of known mean `mu`, the adjusted estimate is
#' `mean(x) - b (mean(c) - mu)` with `b = cov(x, c) / var(c)`, and its MCSE
#' is the standard deviation of `x - b c` over `sqrt(n)` (the usual
#' large-sample approximation that ignores the estimation of `b`).
#'
#' @param run A run or results.
#' @param control Name of a recorded quantity (see the `records` of
#'   [sim_results()]).
#' @param control_mean Its known expectation: a number, or a function of the
#'   scenario.
#' @param measure `"bias"`, `"mse"`, `"coverage"` or `"mean"`.
#' @return A data frame per scenario, method and estimand with the
#'   unadjusted and adjusted values and MCSEs, the coefficient and the
#'   ratio of variances.
#' @export
sim_control_variate <- function(run, control, control_mean,
                                measure = c("bias", "mse", "coverage", "mean")) {
  measure <- match.arg(measure)
  if (!is_scalar_character(control)) stopf("`control` must name a recorded quantity.")
  res <- sim_results(run)
  study <- res$manifest$study
  rec <- res$records
  if (is.null(rec) || !control %in% rec$name) {
    stopf("No recorded quantity `%s`; record it with sim_generator(record = ).", control)
  }
  rec <- rec[rec$name == control, , drop = FALSE]
  e <- vr_contributions(res, measure)
  e$c <- rec$value[match(paste(e$scenario_key, e$replication), paste(rec$scenario_key, rec$replication))]
  e <- e[is.finite(e$c), , drop = FALSE]
  out <- list()
  for (grp in split(e, list(e$scenario_key, e$method, e$estimand), drop = TRUE)) {
    n <- nrow(grp)
    if (n < 3L) next
    sc <- study$scenarios$values[[scenario_index(study$scenarios, grp$scenario_key[1])]]
    mu <- if (is.function(control_mean)) call_flexible(control_mean, list(scenario = sc)) else control_mean
    if (!is_scalar_number(mu)) stopf("`control_mean` must give one number per scenario.")
    vc <- stats::var(grp$c)
    b <- if (vc > 0) stats::cov(grp$x, grp$c) / vc else 0
    adj <- grp$x - b * (grp$c - mu)
    out[[length(out) + 1L]] <- data.frame(
      scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
      method = grp$method[1], estimand = grp$estimand[1], measure = measure, n = n,
      value = mean(grp$x), mcse = stats::sd(grp$x) / sqrt(n),
      value_adjusted = mean(adj), mcse_adjusted = stats::sd(adj) / sqrt(n), coefficient = b,
      variance_ratio = stats::var(grp$x) / stats::var(adj), stringsAsFactors = FALSE)
  }
  out <- bind_rows(out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}
