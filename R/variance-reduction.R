# Variance reduction ------------------------------------------------------

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
