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
  # a mean of one contribution per replication: closed form from the pair
  # sums (the same jackknife, without refitting the measure per pair)
  lc <- linear_contributions(m, rows, study, estimand)
  if (!is.null(lc) && length(lc$x) && all(is.finite(lc$x)) &&
      isTRUE(all.equal(mean(lc$x), r$value, tolerance = 1e-10))) {
    return(linear_unit_result(r, lc$x, rep_units(study, lc$replication), study$level,
                              is_proportion_measure(m)))
  }
  units <- rep_units(study, rows$replication)
  u <- unique(units)
  th <- if (length(u) < 2L) numeric(length(u)) else vapply(u, function(k) {
    v <- measure_cell(m, rows[units != k, , drop = FALSE], study, estimand)$value
    if (is.null(v)) NA_real_ else v
  }, numeric(1))
  unit_jackknife(r, th, study$level, is_proportion_measure(m))
}

vr_contributions <- function(res, measure) {
  e <- res$estimates
  if (measure == "coverage") {
    e <- e[e$interval_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(e$truth), , drop = FALSE]
    e$x <- as.numeric(e$lower <= e$truth & e$truth <= e$upper)
  } else {
    # the mean estimate needs no truth
    e <- e[e$valid & is.finite(e$estimate) & (measure == "mean" | !is.na(e$truth)), , drop = FALSE]
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
#' complete pairs and its MCSE is `sd(pair means) / sqrt(pairs)`. Pairs
#' with only one eligible member are left out and counted in
#' `incomplete_pairs`; when there are any, the value is the mean over the
#' pairs where both members succeeded, not over all pairs. The
#' result also gives the MCSE that independent replications with the same
#' variance would have had, and the ratio of variances (above 1 when the
#' antithetic design helped).
#'
#' @param run A run or results.
#' @param measure `"bias"`, `"mse"`, `"coverage"` or `"mean"` (the mean
#'   estimate).
#' @return A data frame per scenario, method and estimand: `value`,
#'   `mcse`, `pairs`, `incomplete_pairs`, `mcse_independent` and
#'   `variance_ratio`.
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
    incomplete <- length(unique(pair)) - length(full)
    g <- grp[as.character(pair) %in% full, , drop = FALSE]
    pm <- tapply(g$x, (g$replication + 1L) %/% 2L, mean)
    np <- length(pm)
    if (np < 2L) next
    mcse <- stats::sd(pm) / sqrt(np)
    naive <- stats::sd(g$x) / sqrt(2 * np)
    out[[length(out) + 1L]] <- data.frame(
      scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
      method = grp$method[1], estimand = grp$estimand[1], measure = measure,
      value = mean(pm), mcse = mcse, pairs = np, incomplete_pairs = incomplete,
      mcse_independent = naive, variance_ratio = (naive / mcse)^2, stringsAsFactors = FALSE)
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
#' large-sample approximation that ignores the estimation of `b`). With
#' antithetic pairs the independent unit is the pair: `b` is estimated from
#' pair means and both MCSEs come from the pairs (delete-one-pair
#' jackknife), as in [sim_summarise()].
#'
#' The known mean is the control's expectation over all generated data
#' sets. The adjustment is therefore made only when every generated data
#' set contributes (no failed fit, no missing control); otherwise the
#' selected data sets can have another control mean, the adjustment would
#' change the target, and the row gets status `selected` with no adjusted
#' value. The mean estimate (`measure = "mean"`) needs no truth.
#'
#' @param run A run or results.
#' @param control Name of a recorded quantity (see the `records` of
#'   [sim_results()]).
#' @param control_mean Its known expectation: a number, or a function of the
#'   scenario.
#' @param measure `"bias"`, `"mse"`, `"coverage"` or `"mean"`.
#' @return A data frame per scenario, method and estimand with the number
#'   of contributions, of independent units (`n_units`) and of generated
#'   data sets, the unadjusted and
#'   adjusted values and MCSEs, the coefficient, the ratio of variances, a
#'   status and a note.
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
  reps <- res$replications
  out <- list()
  for (grp in split(e, list(e$scenario_key, e$method, e$estimand), drop = TRUE)) {
    n <- nrow(grp)
    if (n < 3L) next
    n_gen <- sum(reps$scenario_key == grp$scenario_key[1] & reps$generation_status == "ok")
    sc <- study$scenarios$values[[scenario_index(study$scenarios, grp$scenario_key[1])]]
    mu <- if (is.function(control_mean)) call_flexible(control_mean, list(scenario = sc)) else control_mean
    if (!is_scalar_number(mu)) stopf("`control_mean` must give one number per scenario.")
    # The known expectation is that of the control over all generated data
    # sets. When contributions were selected (failed fits, missing
    # controls), the control's mean on the selected data sets is unknown and
    # the adjustment would move the estimate towards another target.
    complete <- n == n_gen
    # the coefficient and both MCSEs come from the independent units
    # (antithetic pairs or single replications)
    units <- rep_units(study, grp$replication)
    xb <- tapply(grp$x, units, mean)
    cb <- tapply(grp$c, units, mean)
    vc <- if (length(cb) > 1L) stats::var(cb) else 0
    b <- if (vc > 0) stats::cov(xb, cb) / vc else 0
    adj <- grp$x - b * (grp$c - mu)
    unit_mcse <- function(v) linear_unit_result(kernel_result(value = mean(v), n = length(v)),
                                                v, units, study$level)$mcse
    mc_x <- unit_mcse(grp$x)
    mc_adj <- unit_mcse(adj)
    out[[length(out) + 1L]] <- data.frame(
      scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
      method = grp$method[1], estimand = grp$estimand[1], measure = measure, n = n,
      n_generated = n_gen,
      n_units = length(unique(units)),
      value = mean(grp$x), mcse = mc_x,
      value_adjusted = if (complete) mean(adj) else NA_real_,
      mcse_adjusted = if (complete) mc_adj else NA_real_,
      coefficient = if (complete) b else NA_real_,
      variance_ratio = if (complete && isTRUE(mc_adj > 0)) (mc_x / mc_adj)^2 else NA_real_,
      status = if (complete) "ok" else "selected",
      note = if (complete) NA_character_ else
        sprintf("%d of %d generated data sets contribute; the control mean is known only over all of them, so no adjustment",
                n, n_gen),
      stringsAsFactors = FALSE)
  }
  out <- bind_rows(out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}
