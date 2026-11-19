# Bias decomposition ------------------------------------------------------

#' Decompose bias against a secondary reference value
#'
#' The truth of an estimand is the causal (or target) quantity; a method
#' may instead converge to a different value, for example the identifying
#' expression it estimates (Morris, "causal bias vs statistical bias"). With
#' a second reference value per scenario, the total bias splits into
#' `estimate - reference` (statistical bias: does the method estimate what
#' it identifies?) and `reference - truth` (identification bias: does what
#' it identifies equal the estimand?). Both parts come with Monte Carlo
#' standard errors; the identification part is a property of the scenario
#' and has none when the reference is exact.
#'
#' @param run A run or results.
#' @param reference A function of the scenario returning the reference
#'   value, or a [sim_truth_exact()] or [sim_truth_numeric()] provider.
#' @param estimand Estimand; defaults to the first.
#' @param methods Methods; defaults to all targeting the estimand.
#' @param label Name of the reference (shown in the result).
#' @return A data frame with one row per scenario and method: `truth`,
#'   `reference`, `total_bias`, `statistical_bias`, `identification_bias`,
#'   the MCSE of the first two, and `n`.
#' @examples
#' # the sample mean of log-normal data identifies exp(mu + 1/2), not exp(mu)
#' study <- sim_study("lognormal", scenarios = sim_scenarios(n = 50L, mu = 0),
#'   generate = function(scenario, context) exp(rnorm(scenario$n, scenario$mu)),
#'   estimands = list(sim_estimand("median", sim_truth_exact(function(scenario) exp(scenario$mu)))),
#'   methods = list(mean = function(data) sim_estimate("median", mean(data))))
#' run <- sim_run(study, replications = 100, seed = 1, path = tempfile("adempR-run-"))
#' sim_bias_decomposition(run, function(scenario) exp(scenario$mu + 0.5), label = "population mean")
#' @export
sim_bias_decomposition <- function(run, reference, estimand = NULL, methods = NULL,
                                   label = "reference") {
  res <- sim_results(run)
  study <- res$manifest$study
  estimand <- estimand %||% names(study$estimands)[1]
  if (!estimand %in% names(study$estimands)) stopf("Unknown estimand `%s`.", estimand)
  fun <- if (inherits(reference, "adempr_truth")) {
    if (!reference$kind %in% c("exact", "numeric")) stopf("`reference` must be exact or numeric.")
    reference$fun
  } else reference
  check_function(fun, "reference", "scenario")
  methods <- methods %||% names(study$methods)[vapply(study$methods, function(m) {
    estimand %in% m$estimands
  }, logical(1))]
  est <- res$estimates
  est <- est[est$estimand == estimand & est$method %in% methods & est$valid &
               is.finite(est$estimate) & !is.na(est$truth), , drop = FALSE]
  out <- list()
  for (key in study$scenarios$keys) {
    idx <- scenario_index(study$scenarios, key)
    ref <- call_flexible(fun, list(scenario = study$scenarios$values[[idx]]))
    if (!is_scalar_number(ref)) stopf("`reference` must return one number per scenario.")
    for (mn in methods) {
      r <- est[est$scenario_key == key & est$method == mn, , drop = FALSE]
      n <- nrow(r)
      if (!n) next
      d_total <- r$estimate - r$truth
      d_stat <- r$estimate - ref
      mcse <- function(x) if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = study$scenarios$labels[idx], method = mn, estimand = estimand,
        reference_label = label, truth = mean(r$truth), reference = ref,
        total_bias = mean(d_total), total_mcse = mcse(d_total),
        statistical_bias = mean(d_stat), statistical_mcse = mcse(d_stat),
        identification_bias = ref - mean(r$truth), n = n, stringsAsFactors = FALSE)
    }
  }
  out <- bind_rows(out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

#' Check an algebraic identity between estimands
#'
#' Some estimands are linked by an identity that holds by construction, for
#' example `itt = pc * cace` under one-sided noncompliance, or
#' `log_or = log(or)`. A truth or an estimate that violates such an
#' identity reveals a bug in a truth function or a method. The identity is
#' a function of the estimand values (named as the estimands) and of the
#' scenario, returning a residual that should be zero.
#'
#' @param run A run or results.
#' @param identity Function of the estimand values (as named arguments, or a
#'   list `values`) and optionally `scenario`, returning one number.
#' @param target `"truth"` checks the truths per scenario; `"estimate"`
#'   checks every method that reports all the estimands involved, per
#'   replication.
#' @param tolerance Absolute residual above which a row is flagged.
#' @return A data frame with the residual and `flag` per scenario (and per
#'   method and replication for estimates).
#' @examples
#' study <- sim_study("ratio", scenarios = sim_scenarios(n = 40L, mu = 1),
#'   generate = function(scenario, context) rexp(scenario$n, 1 / scenario$mu),
#'   estimands = list(
#'     sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu)),
#'     sim_estimand("log_mean", sim_truth_exact(function(scenario) log(scenario$mu)))),
#'   methods = list(mle = function(data) list(sim_estimate("mean", mean(data)),
#'                                         sim_estimate("log_mean", log(mean(data))))))
#' run <- sim_run(study, replications = 20, seed = 1, path = tempfile("adempR-run-"))
#' sim_check_identity(run, function(mean, log_mean) log(mean) - log_mean, target = "estimate")
#' @export
sim_check_identity <- function(run, identity, target = c("truth", "estimate"), tolerance = 1e-8) {
  target <- match.arg(target)
  check_function(identity, "identity")
  res <- sim_results(run)
  study <- res$manifest$study
  fa <- names(formals(identity))
  needed <- intersect(fa, names(study$estimands))
  if (!length(needed) && !"values" %in% fa) {
    stopf("`identity` must take estimand names (or `values`) as arguments.")
  }
  if (!length(needed)) needed <- names(study$estimands)
  apply_id <- function(vals, scenario) {
    args <- c(as.list(vals), list(values = as.list(vals), scenario = scenario))
    v <- call_flexible(identity, args)
    if (!is_scalar_number(v) && !(length(v) == 1L && is.na(v))) stopf("`identity` must return one number.")
    as.numeric(v)
  }
  out <- list()
  if (target == "truth") {
    tr <- res$truth
    for (key in study$scenarios$keys) {
      idx <- scenario_index(study$scenarios, key)
      t <- tr[tr$scenario_key == key & tr$truth_status == "ok", , drop = FALSE]
      vals <- stats::setNames(t$truth[match(needed, t$estimand)], needed)
      if (anyNA(vals)) next
      out[[length(out) + 1L]] <- data.frame(scenario_id = study$scenarios$labels[idx],
        residual = apply_id(vals, study$scenarios$values[[idx]]), stringsAsFactors = FALSE)
    }
  } else {
    est <- res$estimates
    est <- est[est$estimand %in% needed & is.finite(est$estimate), , drop = FALSE]
    for (key in study$scenarios$keys) {
      idx <- scenario_index(study$scenarios, key)
      sc <- study$scenarios$values[[idx]]
      e <- est[est$scenario_key == key, , drop = FALSE]
      for (grp in split(e, list(e$method, e$replication), drop = TRUE)) {
        vals <- stats::setNames(grp$estimate[match(needed, grp$estimand)], needed)
        if (anyNA(vals)) next
        out[[length(out) + 1L]] <- data.frame(scenario_id = study$scenarios$labels[idx],
          method = grp$method[1], replication = grp$replication[1],
          residual = apply_id(vals, sc), stringsAsFactors = FALSE)
      }
    }
  }
  out <- bind_rows(out)
  if (is.null(out)) stopf("No scenario has values for all of: %s.", paste(needed, collapse = ", "))
  out$flag <- is.na(out$residual) | abs(out$residual) > tolerance
  rownames(out) <- NULL
  out
}
