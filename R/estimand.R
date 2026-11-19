# Estimands and truth ------------------------------------------------------

EFFECT_TYPES <- c("unspecified", "marginal", "conditional", "subject_specific",
                  "population_average", "finite_population", "other")

#' Define an estimand
#'
#' An estimand records what is being estimated: the quantity, the target
#' population, the conditioning, the scale, the units and its truth. Methods
#' declare which estimands they target, and declared metadata is checked for
#' compatibility (for example, a conditional log odds ratio cannot be
#' evaluated against a marginal estimand). Validation can only check
#' declared metadata; it cannot prove that a user supplied definition is
#' scientifically correct.
#'
#' @param id Identifier used by methods in [sim_estimate()].
#' @param truth A truth provider: [sim_truth_exact()],
#'   [sim_truth_numeric()], [sim_truth_mc()], [sim_truth_dataset()] or
#'   [sim_truth_unknown()].
#' @param quantity Description of the quantity, for example `"mean"` or
#'   `"log odds ratio of treatment"`.
#' @param population Target population.
#' @param conditioning What the estimand conditions on (for example
#'   `"none"` or `"baseline covariates X"`).
#' @param effect One of `"unspecified"`, `"marginal"`, `"conditional"`,
#'   `"subject_specific"`, `"population_average"`, `"finite_population"`,
#'   `"other"`.
#' @param scale Scale of the estimand, for example `"identity"`,
#'   `"log_odds_ratio"`, `"risk_difference"`.
#' @param units Units of measurement.
#' @param horizon Time horizon for time to event estimands.
#' @param type Output type: `"scalar"`, or `"categorical"` for a class
#'   label coded as a number (a selected model, a number of components),
#'   evaluated with [sim_measures_categorical()]. Vector and function
#'   valued estimands are represented by one scalar component per point
#'   with [sim_estimand_grid()]; other types are rejected with an
#'   informative error.
#' @param null Value under the null hypothesis, used to label rejection
#'   rates as type I error or power.
#' @param description Free text.
#' @return An object of class `adempr_estimand`.
#' @examples
#' sim_estimand("mean", truth = sim_truth_exact(function(scenario) scenario$mu),
#'              population = "normal DGP", scale = "identity")
#' @export
sim_estimand <- function(id, truth = sim_truth_unknown(), quantity = NULL,
                         population = NULL, conditioning = NULL,
                         effect = "unspecified", scale = "identity", units = NULL,
                         horizon = NULL, type = "scalar", null = NULL,
                         description = NULL) {
  check_id(id, "id")
  if (!inherits(truth, "adempr_truth")) {
    stopf("`truth` for estimand `%s` must be a truth provider such as sim_truth_exact().", id)
  }
  effect <- match.arg(effect, EFFECT_TYPES)
  if (!is_scalar_character(scale)) stopf("`scale` must be a nonempty string.")
  if (!is_scalar_character(type) || !type %in% c("scalar", "categorical")) {
    stopf("Estimand type `%s` is not supported in adempR %s; only scalar and categorical estimands are implemented. Represent each component as its own scalar estimand (see sim_estimand_grid()).",
          paste(type, collapse = ""), utils::packageVersion("adempR"))
  }
  if (!is.null(null) && !is_scalar_number(null)) stopf("`null` must be a single number.")
  if (!is.null(horizon) && !is_scalar_number(horizon)) {
    stopf("`horizon` must be a single number.")
  }
  structure(list(
    id = id, truth = truth, quantity = quantity %||% id,
    population = population, conditioning = conditioning, effect = effect,
    scale = scale, units = units, horizon = horizon, type = type, null = null,
    description = description
  ), class = "adempr_estimand")
}

#' @export
print.adempr_estimand <- function(x, ...) {
  cat(sprintf("<adempR estimand `%s`: %s on the %s scale, %s effect; truth: %s>\n",
              x$id, x$quantity, x$scale, x$effect, x$truth$kind))
  invisible(x)
}

new_truth <- function(kind, ..., version = NULL, depends = list()) {
  structure(list(kind = kind, ..., version = version, depends = depends),
            class = c(paste0("adempr_truth_", kind), "adempr_truth"))
}

#' Truth providers
#'
#' Truth providers state how the true value of an estimand is obtained.
#'
#' * `sim_truth_exact()`: an analytic value computed from the scenario.
#' * `sim_truth_numeric()`: a deterministic numerical approximation (for
#'   example by [stats::integrate()]) with a reported error bound.
#' * `sim_truth_mc()`: a Monte Carlo approximation computed in its own
#'   stage, with its own random number stream, and stored with its
#'   standard error. Its uncertainty does not shrink when more outer
#'   replications are run.
#' * `sim_truth_dataset()`: a truth specific to each generated data set
#'   (for example a finite population or a sample average effect), returned
#'   by the generator through [sim_data()] and never shown to methods.
#' * `sim_truth_unknown()`: no truth. Metrics that need truth are reported
#'   as unavailable; precision, runtime and failure summaries remain.
#'
#' @param fun For `sim_truth_exact()` and `sim_truth_numeric()`, a function
#'   of `scenario` returning the value (`sim_truth_numeric()` may return a
#'   list with `value` and `error`). For `sim_truth_mc()`, a function of
#'   `scenario` and `n` that uses the random number generator and returns
#'   either a numeric vector of `n` independent contributions whose mean is
#'   the truth, or the result of [sim_truth_value()].
#' @param n Monte Carlo sample size passed to `fun`.
#' @param error For `sim_truth_numeric()`, the absolute error bound when
#'   `fun` returns a bare number.
#' @param version Optional version label; changing it invalidates cached
#'   truth.
#' @param depends Optional list of objects the truth depends on (data,
#'   constants), included in its fingerprint.
#' @param factors For `sim_truth_mc()`, optional names of the scenario
#'   fields the truth depends on. Scenarios with the same values of these
#'   fields share one Monte Carlo truth (computed once, from one stream),
#'   for example when the truth does not depend on the sample size. `fun`
#'   then must use only these fields.
#' @return An object of class `adempr_truth`.
#' @examples
#' sim_truth_exact(function(scenario) scenario$mu)
#' sim_truth_numeric(function(scenario) {
#'   r <- integrate(function(x) plogis(scenario$b0 + x) * dnorm(x), -Inf, Inf)
#'   list(value = r$value, error = r$abs.error)
#' })
#' sim_truth_mc(function(scenario, n) {
#'   x <- rnorm(n)
#'   plogis(scenario$b0 + x)
#' }, n = 1e5)
#' @export
sim_truth_exact <- function(fun, version = NULL, depends = list()) {
  check_function(fun, "fun", "scenario")
  new_truth("exact", fun = fun, version = version, depends = depends)
}

#' @rdname sim_truth_exact
#' @export
sim_truth_numeric <- function(fun, error = NA_real_, version = NULL,
                              depends = list()) {
  check_function(fun, "fun", "scenario")
  if (!identical(error, NA_real_) && (!is_scalar_number(error) || error < 0)) {
    stopf("`error` must be a single nonnegative number (or NA).")
  }
  new_truth("numeric", fun = fun, error = error, version = version,
            depends = depends)
}

#' @rdname sim_truth_exact
#' @export
sim_truth_mc <- function(fun, n = 1e5, version = NULL, depends = list(), factors = NULL) {
  check_function(fun, "fun", c("scenario", "n"))
  if (!is_count(n, 2)) stopf("`n` must be a whole number of at least 2.")
  if (!is.null(factors) && (!is.character(factors) || !length(factors) || anyNA(factors))) {
    stopf("`factors` must name the scenario fields the truth depends on.")
  }
  tr <- new_truth("mc", fun = fun, n = as.integer(n), version = version,
                  depends = depends)
  if (!is.null(factors)) tr$factors <- sort(unique(factors))
  tr
}

# Identity of the truth computation for one scenario: the scenario key, or,
# when the truth depends only on some factors, a key of their values, so
# scenarios that share them share one truth value and one stream.
truth_key <- function(truth, scenario, key) {
  if (is.null(truth$factors)) return(key)
  miss <- setdiff(truth$factors, names(scenario))
  if (length(miss)) stopf("The truth depends on unknown scenario field(s): %s.", paste(miss, collapse = ", "))
  paste0("factors-", stream_fingerprint(truth$factors, scenario[truth$factors]))
}

#' @rdname sim_truth_exact
#' @export
sim_truth_dataset <- function() new_truth("dataset")

#' @rdname sim_truth_exact
#' @export
sim_truth_unknown <- function() new_truth("unknown")

#' A truth value with its uncertainty
#'
#' Used by a [sim_truth_mc()] function that computes a nonlinear truth
#' itself (for example a marginal log odds ratio from simulated risks) and
#' reports its own standard error, for instance by the delta method.
#'
#' @param value The approximated truth.
#' @param se Its Monte Carlo standard error.
#' @param method Text describing how `se` was obtained.
#' @return An object of class `adempr_truth_value`.
#' @examples
#' sim_truth_value(0.5, se = 0.001, method = "delta method")
#' @export
sim_truth_value <- function(value, se, method = "user supplied") {
  if (!is_scalar_number(value)) stopf("`value` must be a single number.")
  if (!is_scalar_number(se) || se < 0) stopf("`se` must be a single nonnegative number.")
  structure(list(value = value, se = se, method = method),
            class = "adempr_truth_value")
}

truth_fingerprint <- function(truth) {
  if (is.null(truth$factors)) {
    return(fingerprint(truth$kind, truth$fun, truth$n, truth$error, truth$version,
                       truth$depends))
  }
  fingerprint(truth$kind, truth$fun, truth$n, truth$error, truth$version,
              truth$depends, truth$factors)
}

#' @export
print.adempr_truth <- function(x, ...) {
  cat(sprintf("<adempR truth provider: %s>\n", x$kind))
  invisible(x)
}

# Compute truth for one scenario. Returns a one-row data frame.
compute_truth_one <- function(truth, scenario, seed = NULL) {
  out <- data.frame(truth_kind = truth$kind, truth = NA_real_, truth_se = NA_real_,
                    truth_error = NA_real_, truth_n = NA_integer_,
                    truth_method = NA_character_, truth_status = "ok",
                    stringsAsFactors = FALSE)
  val <- switch(truth$kind,
    exact = {
      v <- deterministic_call("The exact truth function", truth$fun, scenario)
      if (!is_scalar_number(v)) stopf("An exact truth function must return a single number.")
      out$truth <- v
      out$truth_se <- 0
      out$truth_method <- "analytic"
    },
    numeric = {
      v <- deterministic_call("The numerical truth function", truth$fun, scenario)
      if (is.list(v)) {
        if (!is_scalar_number(v$value)) {
          stopf("A numeric truth function returning a list must contain a single number `value`.")
        }
        out$truth <- v$value
        out$truth_error <- v$error %||% truth$error
        if (!is.null(v$error) && (!is_scalar_number(v$error) || v$error < 0)) {
          stopf("The error bound returned by a numeric truth function must be a single nonnegative number.")
        }
      } else {
        if (!is_scalar_number(v)) stopf("A numeric truth function must return a single number.")
        out$truth <- v
        out$truth_error <- truth$error
      }
      out$truth_se <- 0
      out$truth_method <- "deterministic numerical approximation"
    },
    mc = {
      if (is.null(seed)) stopf("Monte Carlo truth needs an assigned random number stream.")
      v <- with_rng_state(seed, truth$fun(scenario, truth$n))
      if (inherits(v, "adempr_truth_value")) {
        out$truth <- v$value
        out$truth_se <- v$se
        out$truth_method <- paste("Monte Carlo;", v$method)
      } else {
        if (!is.numeric(v) || length(v) < 2L || anyNA(v) || any(!is.finite(v))) {
          stopf("A Monte Carlo truth function must return finite numeric contributions or sim_truth_value().")
        }
        out$truth <- mean(v)
        out$truth_se <- stats::sd(v) / sqrt(length(v))
        out$truth_method <- "Monte Carlo mean of independent contributions"
        if (length(v) != truth$n) {
          out$truth_method <- sprintf("%s (%d contributions returned; %d requested)",
                                      out$truth_method, length(v), truth$n)
        }
      }
      # the count actually used, not only the one requested
      out$truth_n <- if (inherits(v, "adempr_truth_value")) truth$n else length(v)
    },
    dataset = {
      out$truth_status <- "dataset"
      out$truth_method <- "specific to each generated data set"
    },
    unknown = {
      out$truth_status <- "unknown"
      out$truth_method <- "unknown"
    }
  )
  out
}

#' Vector and function valued estimands on a grid
#'
#' A function valued estimand (a survival curve, a dose response curve) or
#' a vector valued one is represented by one scalar component per grid
#' point, so every measure, comparison and diagnostic applies to each point.
#' `sim_estimand_grid()` creates the components, named `<id>_1`, `<id>_2`,
#' ...; `sim_estimate_grid()` returns a method's estimates for all of them;
#' `sim_summarise_grid()` collects a measure along the grid and adds its
#' average over the grid points (for example the integrated squared error
#' from `"mse"`), with an MCSE that accounts for the correlation between
#' points.
#'
#' @param id Identifier stem.
#' @param grid Numeric grid points.
#' @param truth Function of `scenario` and the grid point `at` returning the
#'   true value; wrapped in [sim_truth_exact()].
#' @param ... Passed to [sim_estimand()] for every component.
#' @return `sim_estimand_grid()`: a list of estimands with attributes
#'   recording the grid. `sim_estimate_grid()`: a list of estimates.
#'   `sim_summarise_grid()`: a data frame.
#' @examples
#' grid <- c(0.5, 1, 2)
#' est <- sim_estimand_grid("surv", grid,
#'   truth = function(scenario, at) exp(-scenario$rate * at))
#' study <- sim_study("curve", scenarios = sim_scenarios(n = 100L, rate = 0.5),
#'   generate = function(scenario, context) rexp(scenario$n, scenario$rate),
#'   estimands = est,
#'   methods = list(km = function(data) {
#'     sim_estimate_grid("surv", grid, vapply(grid, function(t) mean(data > t), numeric(1)))
#'   }))
#' run <- sim_run(study, replications = 30, seed = 1, path = tempfile("adempR-run-"))
#' sim_summarise_grid(run, "surv", "bias")
#' @export
sim_estimand_grid <- function(id, grid, truth, ...) {
  check_id(id, "id")
  if (!is.numeric(grid) || !length(grid) || anyNA(grid) || anyDuplicated(grid)) {
    stopf("`grid` must be distinct numbers.")
  }
  check_function(truth, "truth")
  lapply(seq_along(grid), function(i) {
    at <- grid[i]
    f <- local({
      at <- at
      truth <- truth
      function(scenario) truth(scenario = scenario, at = at)
    })
    # the grid point is part of the truth version, so each component's
    # truth has its own fingerprint even though the code is shared
    e <- sim_estimand(paste0(id, "_", i), truth = sim_truth_exact(f, version = paste0("at=", format(at, digits = 17))), ...)
    e$grid <- list(of = id, index = i, at = at)
    e
  })
}

#' @rdname sim_estimand_grid
#' @param values Estimates, one per grid point.
#' @param se,lower,upper Optional vectors of the same length.
#' @param level Interval level.
#' @export
sim_estimate_grid <- function(id, grid, values, se = NULL, lower = NULL, upper = NULL,
                              level = NA_real_) {
  k <- length(grid)
  if (length(values) != k) stopf("`values` must have one entry per grid point.")
  pick <- function(v, i) if (is.null(v)) NA_real_ else v[i]
  lapply(seq_len(k), function(i) {
    sim_estimate(paste0(id, "_", i), values[i], se = pick(se, i), lower = pick(lower, i),
                 upper = pick(upper, i), level = level)
  })
}

#' @rdname sim_estimand_grid
#' @param run A run or results.
#' @param measure One of `"bias"`, `"mse"` or `"coverage"` for the grid
#'   average; any measure for the pointwise table.
#' @export
sim_summarise_grid <- function(run, id, measure = "bias") {
  res <- sim_results(run)
  study <- res$manifest$study
  comp <- Filter(function(e) identical(e$grid$of, id), study$estimands)
  if (!length(comp)) stopf("No grid estimand `%s` in this study.", id)
  ids <- vapply(comp, `[[`, character(1), "id")
  at <- vapply(comp, function(e) e$grid$at, numeric(1))
  s <- as.data.frame(sim_summarise(res, measures = measure))
  s <- s[s$estimand %in% ids & s$measure == measure, , drop = FALSE]
  s$at <- at[match(s$estimand, ids)]
  point <- s[c("scenario_id", "method", "estimand", "at", "measure", "value", "mcse")]
  # grid average: mean over points of each replication's contribution, so
  # its MCSE uses the replication-level averages and is valid under any
  # correlation between points
  avg <- NULL
  if (measure %in% c("bias", "mse", "coverage")) {
    e <- res$estimates
    e <- e[e$estimand %in% ids & !is.na(e$truth), , drop = FALSE]
    e <- if (measure == "coverage") e[e$interval_status == "valid", , drop = FALSE] else
      e[e$valid & is.finite(e$estimate), , drop = FALSE]
    e$contrib <- switch(measure, bias = e$estimate - e$truth, mse = (e$estimate - e$truth)^2,
                        coverage = as.numeric(e$lower <= e$truth & e$truth <= e$upper))
    rows <- list()
    for (grp in split(e, list(e$scenario_key, e$method), drop = TRUE)) {
      per <- tapply(grp$contrib, grp$replication, function(v) if (length(v) == length(ids)) mean(v) else NA_real_)
      per <- per[!is.na(per)]
      n <- length(per)
      if (!n) next
      rows[[length(rows) + 1L]] <- data.frame(
        scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
        method = grp$method[1], estimand = id, at = NA_real_,
        measure = paste0("grid_mean_", measure), value = mean(per),
        mcse = if (n > 1) stats::sd(per) / sqrt(n) else NA_real_, stringsAsFactors = FALSE)
    }
    avg <- bind_rows(rows)
  }
  out <- bind_rows(list(point, avg))
  rownames(out) <- NULL
  out
}
