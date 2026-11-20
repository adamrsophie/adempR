# Bias decomposition ------------------------------------------------------

#' Decompose bias against a secondary reference value
#'
#' The truth of an estimand is the causal (or target) quantity; a method
#' may instead converge to a different value, for example the identifying
#' expression it estimates (Morris, "causal bias vs statistical bias"). With
#' a second reference value per scenario, the total bias splits into
#' `estimate - reference` (statistical bias: does the method estimate what
#' it identifies?) and `reference - truth` (identification bias: does what
#' it identifies equal the estimand?).
#'
#' The reference is evaluated like a truth: a [sim_truth_numeric()]
#' provider may return `list(value, error)`, and the error bound is
#' reported in `reference_error`. Each part has a Monte Carlo standard
#' error. `total_mcse` includes the standard error of a shared Monte Carlo
#' truth. The identification part is a property of the scenario: with an
#' exact truth its MCSE is 0, with a Monte Carlo truth it is the truth's
#' standard error, and with a data set truth it is averaged over the
#' replications in which the method's fit is usable, so it can differ
#' between methods and has the MCSE of that average.
#'
#' @param run A run or results.
#' @param reference A function of the scenario returning the reference
#'   value, or a [sim_truth_exact()] or [sim_truth_numeric()] provider.
#' @param estimand Estimand; defaults to the first.
#' @param methods Methods; defaults to all targeting the estimand.
#' @param label Name of the reference (shown in the result).
#' @return A data frame with one row per scenario and method: `truth`,
#'   `reference`, `reference_error`, `total_bias`, `statistical_bias` and
#'   `identification_bias`, each with its MCSE, `n` and `n_units` (the
#'   independent units, antithetic pairs or replications, over which the
#'   MCSEs are taken).
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
  if (inherits(reference, "adempr_truth")) {
    if (!reference$kind %in% c("exact", "numeric")) stopf("`reference` must be exact or numeric.")
  } else {
    check_function(reference, "reference", "scenario")
    reference <- sim_truth_exact(reference)
  }
  methods <- methods %||% names(study$methods)[vapply(study$methods, function(m) {
    estimand %in% m$estimands
  }, logical(1))]
  est <- res$estimates
  est <- est[est$estimand == estimand & est$method %in% methods & est$valid &
               is.finite(est$estimate) & !is.na(est$truth), , drop = FALSE]
  out <- list()
  for (key in study$scenarios$keys) {
    idx <- scenario_index(study$scenarios, key)
    ref <- tryCatch(compute_truth_one(reference, study$scenarios$values[[idx]]),
      error = function(e) {
        stopf("`reference` must return one number per scenario (scenario %s): %s",
              study$scenarios$labels[idx], conditionMessage(e))
      })
    for (mn in methods) {
      r <- est[est$scenario_key == key & est$method == mn, , drop = FALSE]
      n <- nrow(r)
      if (!n) next
      # MCSEs over the independent units (antithetic pairs or replications)
      units <- rep_units(study, r$replication)
      mcse <- function(x) linear_unit_result(kernel_result(value = mean(x), n = n), x, units,
                                             study$level)$mcse
      # a shared truth estimated by simulation adds its own error to every
      # replication alike; a data set truth varies and is in the spread
      shared <- !all(r$truth_kind == "dataset")
      tse <- if (shared) max(r$truth_se[1], 0, na.rm = TRUE) else 0
      d_total <- r$estimate - r$truth
      d_stat <- r$estimate - ref$truth
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = study$scenarios$labels[idx], method = mn, estimand = estimand,
        reference_label = label, truth = mean(r$truth), reference = ref$truth,
        reference_error = ref$truth_error,
        total_bias = mean(d_total), total_mcse = sqrt(mcse(d_total)^2 + tse^2),
        statistical_bias = mean(d_stat), statistical_mcse = mcse(d_stat),
        identification_bias = ref$truth - mean(r$truth),
        identification_mcse = if (shared) tse else mcse(r$truth),
        n = n, n_units = length(unique(units)), stringsAsFactors = FALSE)
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
#' A residual fails when its absolute value exceeds `tolerance` times
#' `scale`, plus the residual's certified numerical error bound, plus `z`
#' standard errors of the residual. All three terms are in the units of the
#' residual: `scale` (1 by default) is a reference size for the residual,
#' given as a number or as a function of the estimand values (and the
#' scenario), for example the size of one side of the identity when its
#' values are large; it is never taken from the raw inputs, because a
#' residual such as `log(or) - log_or` does not share their scale. The
#' standard error comes from the standard errors of Monte Carlo truths,
#' and the error bound from the error bounds of numerical truths, each
#' propagated through the identity by numerical differentiation; exact
#' truths and estimates have neither. Every checked
#' unit gets a row: `status` is `pass`, `fail` or `unavailable` (a value
#' is missing, not usable or not computed), with the reason in `note`.
#'
#' @param run A run or results.
#' @param identity Function of the estimand values (as named arguments, or a
#'   list `values`) and optionally `scenario`, returning one number.
#' @param target `"truth"` checks the truths per scenario (per replication
#'   when one of them is a data set truth); `"estimate"` checks every
#'   method that targets all the estimands involved, per replication, using
#'   usable estimates only.
#' @param tolerance Residual allowed, as a multiple of `scale`.
#' @param estimands Estimands the identity uses. Defaults to the
#'   estimands named among its arguments, or all when it takes `values`.
#' @param z Number of residual standard errors added to the tolerance.
#' @param scale Reference size of the residual, in its units: a positive
#'   number, or a function of the estimand values (as named arguments or
#'   `values`) and optionally `scenario` returning one.
#' @return A data frame with one row per scenario (and per method or
#'   replication where applicable): `residual`, `residual_se`, `status`,
#'   `residual_error` (certified bound from numerical truths), `allowed`,
#'   `flag` (`TRUE` for a failure, `NA` when unavailable) and `note`.
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
sim_check_identity <- function(run, identity, target = c("truth", "estimate"), tolerance = 1e-8,
                               estimands = NULL, z = 3, scale = 1) {
  target <- match.arg(target)
  check_function(identity, "identity")
  if (!is_scalar_number(tolerance) || tolerance < 0) stopf("`tolerance` must be a nonnegative number.")
  if (!is_scalar_number(z) || z < 0) stopf("`z` must be a nonnegative number.")
  if (!is.function(scale) && (!is_scalar_number(scale) || scale <= 0)) {
    stopf("`scale` must be a positive number or a function of the estimand values.")
  }
  res <- sim_results(run)
  study <- res$manifest$study
  fa <- names(formals(identity))
  needed <- estimands %||% intersect(fa, names(study$estimands))
  unknown <- setdiff(needed, names(study$estimands))
  if (length(unknown)) stopf("Unknown estimand(s): %s.", paste(unknown, collapse = ", "))
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
  # gradient of the identity at the values, where a value has a standard
  # error or an error bound (numerical differentiation)
  grad <- function(vals, w, scenario, r0) {
    vapply(seq_along(vals), function(i) {
      if (w[i] == 0) return(0)
      h <- 1e-6 * max(1, abs(vals[i]))
      v <- vals
      v[i] <- v[i] + h
      (apply_id(v, scenario) - r0) / h
    }, numeric(1))
  }
  scale_of <- function(vals, scenario) {
    if (!is.function(scale)) return(scale)
    s <- call_flexible(scale, c(as.list(vals), list(values = as.list(vals), scenario = scenario)))
    if (!is_scalar_number(s) || !is.finite(s) || s <= 0) stopf("`scale` must return one positive number.")
    s
  }
  check_one <- function(vals, se, scenario, err = 0) {
    miss <- needed[is.na(vals)]
    if (length(miss)) {
      return(list(residual = NA_real_, residual_se = NA_real_, residual_error = NA_real_,
                  allowed = NA_real_, status = "unavailable",
                  note = paste("no usable value for", paste(miss, collapse = ", "))))
    }
    r0 <- apply_id(vals, scenario)
    if (!is.finite(r0)) {
      return(list(residual = r0, residual_se = NA_real_, residual_error = NA_real_,
                  allowed = NA_real_, status = "fail",
                  note = "the identity is not finite at these values"))
    }
    se[is.na(se)] <- 0
    err <- rep_len(err, length(vals))
    unstated <- needed[is.na(err)]
    err[is.na(err)] <- 0
    g <- grad(vals, as.numeric(se > 0 | err > 0), scenario, r0)
    rse <- sqrt(sum((g * se)^2))
    rerr <- sum(abs(g) * err)
    allowed <- tolerance * scale_of(vals, scenario) + rerr + z * rse
    list(residual = r0, residual_se = rse, residual_error = rerr, allowed = allowed,
         status = if (abs(r0) > allowed) "fail" else "pass",
         note = if (length(unstated)) paste("numerical truth without a stated error bound:",
                                            paste(unstated, collapse = ", ")) else NA_character_)
  }
  row <- function(lab, chk, ...) {
    data.frame(scenario_id = lab, ..., residual = chk$residual, residual_se = chk$residual_se,
               residual_error = chk$residual_error, allowed = chk$allowed,
               status = chk$status, note = chk$note, stringsAsFactors = FALSE)
  }
  out <- list()
  if (target == "truth") {
    tr <- res$truth
    ds <- needed[vapply(needed, function(id) study$estimands[[id]]$truth$kind == "dataset", logical(1))]
    dtr <- res$dataset_truth
    for (key in study$scenarios$keys) {
      idx <- scenario_index(study$scenarios, key)
      sc <- study$scenarios$values[[idx]]
      lab <- study$scenarios$labels[idx]
      t <- tr[tr$scenario_key == key & tr$truth_status == "ok", , drop = FALSE]
      vals <- stats::setNames(t$truth[match(needed, t$estimand)], needed)
      se <- stats::setNames(t$truth_se[match(needed, t$estimand)], needed)
      # error bounds of numerical truths (NA when not stated); other kinds
      # have none
      kinds <- t$truth_kind[match(needed, t$estimand)]
      err <- ifelse(kinds %in% "numeric", t$truth_error[match(needed, t$estimand)], 0)
      if (!length(ds)) {
        out[[length(out) + 1L]] <- row(lab, check_one(vals, se, sc, err))
        next
      }
      d <- if (!is.null(dtr)) dtr[dtr$scenario_key == key & dtr$estimand %in% ds, , drop = FALSE] else NULL
      reps <- sort(unique(res$replications$replication[res$replications$scenario_key == key &
                                                         res$replications$generation_status == "ok"]))
      for (rp in reps) {
        dr <- d[d$replication == rp, , drop = FALSE]
        v <- vals
        s <- se
        v[ds] <- dr$truth[match(ds, dr$estimand)]
        s[ds] <- 0
        e2 <- err
        e2[match(ds, needed)] <- 0
        out[[length(out) + 1L]] <- row(lab, check_one(v, s, sc, e2), replication = rp)
      }
    }
  } else {
    meths <- names(study$methods)[vapply(study$methods, function(m) all(needed %in% m$estimands),
                                         logical(1))]
    if (!length(meths)) stopf("No method targets all of: %s.", paste(needed, collapse = ", "))
    est <- res$estimates
    est <- est[est$estimand %in% needed & est$method %in% meths, , drop = FALSE]
    ok <- est$valid & is.finite(est$estimate)
    est$estimate[!ok] <- NA_real_
    att <- res$attempts
    for (key in study$scenarios$keys) {
      idx <- scenario_index(study$scenarios, key)
      sc <- study$scenarios$values[[idx]]
      lab <- study$scenarios$labels[idx]
      for (mn in meths) {
        reps <- sort(unique(att$replication[att$scenario_key == key & att$method == mn]))
        e <- est[est$scenario_key == key & est$method == mn, , drop = FALSE]
        for (rp in reps) {
          grp <- e[e$replication == rp, , drop = FALSE]
          vals <- stats::setNames(grp$estimate[match(needed, grp$estimand)], needed)
          out[[length(out) + 1L]] <- row(lab, check_one(vals, rep(0, length(needed)), sc),
                                         method = mn, replication = rp)
        }
      }
    }
  }
  out <- bind_rows(out)
  if (is.null(out)) stopf("Nothing to check for: %s.", paste(needed, collapse = ", "))
  out$flag <- ifelse(out$status == "unavailable", NA, out$status == "fail")
  out <- out[c(setdiff(names(out), c("flag", "note")), "flag", "note")]
  rownames(out) <- NULL
  out
}
