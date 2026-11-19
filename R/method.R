# Methods, generators and method outputs -----------------------------------

#' Define a data-generating mechanism
#'
#' Wraps a function `generate(scenario, context)` that simulates one data
#' set. The execution layer assigns the random number stream before calling
#' it; the function must not call [set.seed()]. It may return the analysis
#' data directly or [sim_data()] to attach protected metadata such as a
#' data set specific truth.
#'
#' Supported data are atomic vectors, matrices, data frames, data.table
#' objects and lists of these. Environments, external pointers and
#' reference class objects are rejected, because they cannot be isolated
#' between methods by copying.
#'
#' @param fun Function with arguments `scenario` and `context` (or
#'   `...`). `scenario` is the full scenario including privileged
#'   parameters; `context` holds `replication`, `scenario_key` and `stage`.
#' @param check Optional function of `scenario` that stops with an error
#'   when the scenario is an impossible configuration for this mechanism
#'   (for example a correlation matrix that is not positive semidefinite).
#'   Called by [sim_validate()] before anything runs.
#' @param version Optional version label; part of the fingerprint.
#' @param depends Optional list of objects the mechanism depends on (data
#'   sets, constants, helper functions), part of the fingerprint. Closure
#'   variables are fingerprinted automatically, but code in files, options
#'   and package implementations are not; declare what matters here.
#' @param packages Packages the mechanism uses. They are checked by
#'   [sim_validate()] and attached before study code runs, in the main
#'   process for the sequential backend and on every worker for the future
#'   backend, so unqualified calls resolve the same way in both.
#' @param description Free text for the protocol.
#' @return An object of class `adempr_generator`.
#' @examples
#' sim_generator(function(scenario, context) {
#'   data.frame(y = rnorm(scenario$n, scenario$mu))
#' }, description = "independent normal outcomes")
#' @export
sim_generator <- function(fun, check = NULL, version = NULL, depends = list(),
                          packages = character(), description = NULL) {
  check_function(fun, "fun", c("scenario", "context"))
  if (!is.null(check)) check_function(check, "check", "scenario")
  structure(list(fun = fun, check = check, version = version, depends = depends,
                 packages = packages, description = description),
            class = "adempr_generator")
}

generator_fingerprint <- function(g) {
  fingerprint(g$fun, g$version, g$depends)
}

#' Generated data with protected metadata
#'
#' A generator can return `sim_data()` to separate the analysis data, which
#' methods see, from metadata that methods must not see, such as a data set
#' specific truth, latent variables or future outcomes.
#'
#' @param data The analysis data passed to methods.
#' @param truth Named list of data set specific truth values, one per
#'   estimand whose truth provider is [sim_truth_dataset()].
#' @param meta Named list of other protected metadata, available to oracle
#'   methods and to diagnostics.
#' @return An object of class `adempr_data`.
#' @examples
#' sim_data(data.frame(y = rnorm(5)), truth = list(sate = 0.3))
#' @export
sim_data <- function(data, truth = list(), meta = list()) {
  if (!is.list(truth) || (length(truth) && is.null(names(truth)))) {
    stopf("`truth` must be a named list.")
  }
  for (n in names(truth)) {
    if (!is_scalar_number(truth[[n]])) {
      stopf("Data set specific truth `%s` must be a single number.", n)
    }
  }
  structure(list(data = data, truth = truth, meta = meta), class = "adempr_data")
}

# Detect objects that cannot be isolated by copying.
contains_reference <- function(x, depth = 0L) {
  if (depth > 50L) return(FALSE)
  if (is.environment(x) || typeof(x) %in% c("externalptr", "weakref", "bytecode")) {
    return(TRUE)
  }
  if (is.function(x)) return(FALSE)
  if (is.list(x)) {
    for (el in x) if (contains_reference(el, depth + 1L)) return(TRUE)
  }
  FALSE
}

check_generated <- function(x) {
  if (is.null(x)) stopf("The generator returned NULL.")
  if (isS4(x)) stopf("The generator returned an S4 object; supported outputs are vectors, matrices, data frames and lists of these.")
  if (contains_reference(x)) {
    stopf("The generator returned an object containing an environment or external pointer (class %s). Such objects cannot be isolated between methods; return plain data instead.",
          paste(class(x), collapse = "/"))
  }
  invisible(TRUE)
}

# A private copy of the data for one method, so that mutation by reference
# (for example data.table `:=` or setattr()) cannot leak into later methods.
isolate_data <- function(x) {
  if (inherits(x, "data.table") && requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::copy(x))
  }
  unserialize(serialize(x, NULL))
}

#' Define a method under evaluation
#'
#' Wraps an ordinary R function that analyzes one data set. The function is
#' called as `fun(data, scenario, context)`; it may omit `scenario` and
#' `context` from its arguments. `scenario` is a *public view* containing
#' only the fields named in `scenario_fields`, and `context` contains only
#' the replication number, the method identifier and an opaque scenario
#' key. Privileged data-generating parameters and truth are visible only to
#' methods declared as `oracle = TRUE`.
#'
#' Separating fitting from extraction: if `extract` is supplied, `fun`
#' returns a fitted object and `extract(fit, data, scenario, context)`
#' returns the estimates. Errors in the two steps are recorded as different
#' statuses.
#'
#' Each method has its own random number stream per replication, so adding
#' or reordering methods does not change any other method's randomness.
#' Methods listed in the same `coupling` group share one stream; this is an
#' explicit policy for sharing, for example, bootstrap samples or folds.
#'
#' @param fun Function analyzing one data set.
#' @param estimands Character vector of estimand identifiers the method
#'   targets. `NULL` means all estimands of the study.
#' @param scenario_fields Scenario fields the method may legitimately know
#'   (for example the sample size of a design). Defaults to none.
#' @param oracle If `TRUE`, the method receives the full scenario and the
#'   protected metadata of [sim_data()] in `context$meta` and
#'   `context$truth`. Use only for oracle benchmarks.
#' @param extract Optional extraction function.
#' @param declares Optional named list, by estimand, of declared metadata
#'   (`effect`, `scale`) that must match the estimand.
#' @param version Optional version label; part of the fingerprint.
#' @param depends Optional list of objects the method depends on; part of
#'   the fingerprint.
#' @param packages Packages the method uses; checked by [sim_validate()]
#'   and attached in both backends (see [sim_generator()]).
#' @param timeout Optional elapsed time limit in seconds. Enforced with
#'   [setTimeLimit()], which is cooperative: it interrupts R code at safe
#'   points but may not interrupt long running compiled code. This is a best
#'   effort limit, not a hard timeout.
#' @param coupling Optional name of a group of methods sharing one random
#'   number stream.
#' @param description Free text for the protocol.
#' @return An object of class `adempr_method`.
#' @examples
#' sim_method(function(data) {
#'   sim_estimate("mean", estimate = mean(data$y),
#'                se = sd(data$y) / sqrt(nrow(data)))
#' })
#' @export
sim_method <- function(fun, estimands = NULL, scenario_fields = character(),
                       oracle = FALSE, extract = NULL, declares = list(),
                       version = NULL, depends = list(), packages = character(),
                       timeout = NULL, coupling = NULL, description = NULL) {
  check_function(fun, "fun")
  fa <- names(formals(fun))
  if (!length(fa)) stopf("A method function must accept at least the `data` argument.")
  if (!is.null(extract)) check_function(extract, "extract")
  if (!is.null(estimands) && (!is.character(estimands) || anyNA(estimands))) {
    stopf("`estimands` must be a character vector of estimand identifiers.")
  }
  if (!is.character(scenario_fields)) stopf("`scenario_fields` must be a character vector.")
  if (!is.logical(oracle) || length(oracle) != 1L || is.na(oracle)) {
    stopf("`oracle` must be TRUE or FALSE.")
  }
  if (!is.null(timeout) && (!is_scalar_number(timeout) || timeout <= 0)) {
    stopf("`timeout` must be a positive number of seconds.")
  }
  if (!is.null(coupling)) check_id(coupling, "coupling")
  if (!is.list(declares)) stopf("`declares` must be a named list.")
  structure(list(fun = fun, estimands = estimands, scenario_fields = scenario_fields,
                 oracle = oracle, extract = extract, declares = declares,
                 version = version, depends = depends, packages = packages,
                 timeout = timeout, coupling = coupling, description = description),
            class = "adempr_method")
}

method_fingerprint <- function(m) {
  fingerprint(m$fun, m$extract, m$estimands, m$scenario_fields, m$oracle,
              m$declares, m$version, m$depends, m$coupling)
}

call_flexible <- function(f, args) {
  fa <- names(formals(f))
  if ("..." %in% fa) return(do.call(f, args))
  do.call(f, args[intersect(names(args), fa)])
}

#' Report the output of a method for one estimand
#'
#' The typed record a method returns for each estimand. Intervals and tests
#' are stored with their own metadata and are never derived from the
#' standard error behind the user's back: a confidence interval is used for
#' coverage only when `lower` and `upper` are supplied.
#'
#' @param estimand Estimand identifier.
#' @param estimate Point estimate. `Inf` and `-Inf` are legitimate extended
#'   real results; `NA` and `NaN` mean no estimate.
#' @param se Estimated standard error.
#' @param lower,upper Confidence interval limits. Infinite limits are
#'   allowed (one sided or unbounded intervals).
#' @param level Confidence level of the interval.
#' @param interval_type Text label such as `"wald"`, `"t"`, `"profile"`.
#' @param p_value P value of the test of `null`.
#' @param null Null value of the test.
#' @param alternative `"two.sided"`, `"less"` or `"greater"`.
#' @param alpha Significance level for rejection; if missing, rejection
#'   uses the study's `alpha`.
#' @param reject Optional explicit logical rejection decision, for tests
#'   without a p value.
#' @param df Degrees of freedom where relevant.
#' @param converged Optional logical convergence flag.
#' @param boundary Optional logical flag for estimates at a parameter space
#'   boundary.
#' @param effect,scale Optional declared metadata, checked against the
#'   estimand.
#' @param note Optional text.
#' @return An object of class `adempr_estimate`.
#' @examples
#' sim_estimate("mean", estimate = 0.1, se = 0.05, lower = 0.002,
#'              upper = 0.198, level = 0.95, interval_type = "wald")
#' @export
sim_estimate <- function(estimand, estimate, se = NA_real_, lower = NA_real_,
                         upper = NA_real_, level = NA_real_,
                         interval_type = NA_character_, p_value = NA_real_,
                         null = NA_real_, alternative = "two.sided",
                         alpha = NA_real_, reject = NA, df = NA_real_,
                         converged = NA, boundary = NA, effect = NULL,
                         scale = NULL, note = NA_character_) {
  if (!is_scalar_character(estimand)) stopf("`estimand` must be a single string.")
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  num1 <- function(x, nm) {
    if (is.null(x)) return(NA_real_)
    if (length(x) != 1L || !(is.numeric(x) || is.logical(x))) {
      stopf("`%s` must be a single number (estimand `%s`).", nm, estimand)
    }
    as.double(x)
  }
  lgl1 <- function(x, nm) {
    if (is.null(x)) return(NA)
    if (length(x) != 1L || !(is.logical(x) || is.numeric(x))) {
      stopf("`%s` must be a single logical value (estimand `%s`).", nm, estimand)
    }
    as.logical(x)
  }
  structure(list(
    estimand = estimand, estimate = num1(estimate, "estimate"),
    se = num1(se, "se"), lower = num1(lower, "lower"), upper = num1(upper, "upper"),
    level = num1(level, "level"),
    interval_type = as.character(interval_type %||% NA_character_)[1],
    p_value = num1(p_value, "p_value"), null = num1(null, "null"),
    alternative = alternative, alpha = num1(alpha, "alpha"),
    reject = lgl1(reject, "reject"), df = num1(df, "df"),
    converged = lgl1(converged, "converged"), boundary = lgl1(boundary, "boundary"),
    effect = effect, scale = scale, note = as.character(note %||% NA)[1]
  ), class = "adempr_estimate")
}

#' @export
print.adempr_estimate <- function(x, ...) {
  cat(sprintf("<adempR estimate `%s`: %s (SE %s)>\n", x$estimand,
              fmt_num(x$estimate), fmt_num(x$se)))
  invisible(x)
}

#' Failure and validity policy
#'
#' States how execution problems and questionable outputs are classified.
#' Every outcome is recorded either way; the policy decides which outputs
#' count as valid for performance measures.
#'
#' @param nonconverged `"invalid"` (default) treats estimates flagged
#'   `converged = FALSE` as unusable; `"valid"` keeps them.
#' @param warnings `"record"` (default) keeps outputs produced with warnings
#'   and records the warnings; `"invalid"` makes them unusable.
#' @param boundary `"record"` (default) or `"invalid"` for estimates flagged
#'   at a boundary.
#' @param infinite `"include"` (default) keeps infinite estimates as
#'   legitimate extended real results, so error based metrics can be
#'   infinite; `"invalid"` makes them unusable.
#' @param retries Number of infrastructure retries after an error. A retry
#'   reuses the same data set and random number state; it is not a new
#'   replication. Use it only for transient infrastructure faults: a
#'   statistical method's own restarts belong inside the method.
#' @return An object of class `adempr_failure_policy`.
#' @examples
#' sim_failure_policy(nonconverged = "invalid", warnings = "record")
#' @export
sim_failure_policy <- function(nonconverged = c("invalid", "valid"),
                               warnings = c("record", "invalid"),
                               boundary = c("record", "invalid"),
                               infinite = c("include", "invalid"),
                               retries = 0L) {
  if (!is_count(retries)) stopf("`retries` must be a nonnegative whole number.")
  structure(list(nonconverged = match.arg(nonconverged), warnings = match.arg(warnings),
                 boundary = match.arg(boundary), infinite = match.arg(infinite),
                 retries = as.integer(retries)),
            class = "adempr_failure_policy")
}

RAW_ESTIMATE_COLUMNS <- c("estimand", "estimate", "se", "lower", "upper", "level",
                          "interval_type", "p_value", "null", "alternative",
                          "alpha", "reject", "df", "converged", "boundary",
                          "declared_effect", "declared_scale", "note")

raw_row <- function(e) {
  data.frame(estimand = e$estimand, estimate = e$estimate, se = e$se,
             lower = e$lower, upper = e$upper, level = e$level,
             interval_type = e$interval_type, p_value = e$p_value, null = e$null,
             alternative = e$alternative, alpha = e$alpha, reject = e$reject,
             df = e$df, converged = e$converged, boundary = e$boundary,
             declared_effect = if (is.null(e$effect)) NA_character_ else as.character(e$effect)[1],
             declared_scale = if (is.null(e$scale)) NA_character_ else as.character(e$scale)[1],
             note = e$note, stringsAsFactors = FALSE)
}

missing_row <- function(estimand) {
  raw_row(sim_estimate(estimand, estimate = NA_real_, alternative = "two.sided",
                       note = "not returned by the method"))
}

# Convert a method's output to raw estimate rows. Errors are informative and
# are recorded by the caller as extraction failures. Validity is classified
# later, at evaluation, so that a change of failure policy does not require
# refitting.
outputs_to_rows <- function(out, targets) {
  if (inherits(out, "adempr_estimate")) out <- list(out)
  if (!is.list(out) || !length(out) ||
      !all(vapply(out, inherits, logical(1), "adempr_estimate"))) {
    stopf("The method returned an object of class %s; it must return sim_estimate() or a list of sim_estimate() results (or supply `extract`).",
          paste(class(out), collapse = "/"))
  }
  ids <- vapply(out, `[[`, character(1), "estimand")
  if (anyDuplicated(ids)) {
    stopf("The method returned several results for estimand(s) %s in one replication. Nested or repeated outputs (imputations, bootstrap draws, chains) must be combined within the method; the outer replication is the unit of performance evaluation.",
          paste(unique(ids[duplicated(ids)]), collapse = ", "))
  }
  unknown <- setdiff(ids, targets)
  if (length(unknown)) {
    stopf("The method returned undeclared estimand(s): %s.", paste(unknown, collapse = ", "))
  }
  rows <- lapply(targets, function(t) {
    j <- match(t, ids)
    if (is.na(j)) missing_row(t) else raw_row(out[[j]])
  })
  do.call(rbind, rows)
}

# Vectorized validity classification of estimate rows (joined with the
# attempt status and warning flag) under a failure policy and the
# estimand metadata.
classify_rows <- function(est, estimands, policy) {
  n <- nrow(est)
  if (!n) return(est)
  e <- est$estimate
  eff <- vapply(est$estimand, function(id) estimands[[id]]$effect, character(1))
  scl <- vapply(est$estimand, function(id) estimands[[id]]$scale, character(1))
  incompatible <- (!is.na(est$declared_effect) & est$declared_effect != eff) |
    (!is.na(est$declared_scale) & est$declared_scale != scl)
  # Checks of the fit and of the estimand apply to every output of a row
  # (point estimate, SE, interval, test), so they come first: an interval
  # only result from a nonconverged fit must not count toward coverage.
  status <- rep(NA_character_, n)
  set <- function(cond, value) status[is.na(status) & cond] <<- value
  set(est$attempt_status != "ok", est$attempt_status[is.na(status) & est$attempt_status != "ok"])
  set(incompatible, "incompatible")
  if (policy$nonconverged == "invalid") set(est$converged %in% FALSE, "nonconverged")
  if (policy$boundary == "invalid") set(est$boundary %in% TRUE, "boundary")
  if (policy$warnings == "invalid") set(est$had_warning %in% TRUE, "warning")
  set(is.na(e), "missing")
  set(is.infinite(e), if (policy$infinite == "include") "infinite" else "invalid_infinite")
  set(TRUE, "valid")
  est$estimate_status <- status
  est$valid <- status %in% c("valid", "infinite")
  est$se_status <- ifelse(is.na(est$se), "missing",
                   ifelse(!is.finite(est$se) | est$se < 0, "invalid", "valid"))
  lo <- est$lower
  up <- est$upper
  est$interval_status <- ifelse(is.na(lo) & is.na(up), "missing",
    ifelse(is.na(lo) | is.na(up) | lo > up | (is.infinite(lo) & lo > 0) |
             (is.infinite(up) & up < 0), "invalid", "valid"))
  est$test_status <- ifelse(!is.na(est$reject), "valid",
    ifelse(is.na(est$p_value), "missing",
    ifelse(est$p_value < 0 | est$p_value > 1, "invalid", "valid")))
  est
}
