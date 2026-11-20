# Performance measures -----------------------------------------------------
#
# Reference kernels. Each takes the eligible per replication values of one
# cell (scenario x method x estimand) and returns a list with value, mcse,
# lower, upper, n, status and note. Formulas and their assumptions are
# documented in vignette("performance-measures").

MEASURE_INFO <- data.frame(
  measure = c("bias", "relbias", "mae", "mse", "rmse", "empse", "mean_se",
              "rmvar", "se_ratio", "coverage", "noncoverage_lower",
              "noncoverage_upper", "be_coverage", "width", "rejection", "usable",
              "usable_coverage"),
  needs_truth = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                  TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
  eligibility = c(
    "usable estimate and truth",
    "usable estimate and nonzero truth",
    "usable estimate and truth",
    "usable estimate and truth",
    "usable estimate and truth",
    "usable estimate",
    "usable estimate with valid SE",
    "usable estimate with valid SE",
    "usable estimate with valid SE",
    "valid interval",
    "valid interval",
    "valid interval",
    "valid interval and usable estimate",
    "valid interval",
    "valid test result",
    "generated data sets",
    "generated data sets"
  ),
  direction = c("closer to 0", "closer to 0", "lower", "lower", "lower",
                "lower", "none", "none", "closer to 0", "closer to nominal",
                "closer to (1 - nominal) / 2", "closer to (1 - nominal) / 2", "closer to nominal", "lower",
                "depends on hypothesis", "higher", "higher"),
  mcse_method = c(
    "sd(e) / sqrt(R)",
    "MCSE(bias) / |truth|",
    "sd(|e|) / sqrt(R)",
    "sd(e^2) / sqrt(R)",
    "delta method: MCSE(MSE) / (2 RMSE)",
    "delta method with sample kurtosis: sqrt(Var(s^2)) / (2 s)",
    "sd(SE) / sqrt(R)",
    "delta method: sd(SE^2) / sqrt(R) / (2 sqrt(mean(SE^2)))",
    "Morris et al. (2019) Table 6 (normal theory for EmpSE)",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval",
    "binomial, conditional on the mean estimate (Morris et al. 2019)",
    "sd(width) / sqrt(R)",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval",
    "binomial: sqrt(p (1 - p) / R) with Clopper-Pearson interval"
  ),
  default = c(TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE,
              FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

#' Performance measures
#'
#' `sim_measures()` selects the performance measures to report.
#' `sim_measure_info()` lists the built in measures with their eligibility
#' rule, preferred direction and Monte Carlo standard error method.
#' Probability of a usable result and failure summaries are always reported,
#' whether or not they are selected.
#'
#' Built in measures: `bias`, `relbias` (relative bias, unavailable at zero
#' truth), `mae`, `mse`, `rmse`, `empse` (empirical standard error),
#' `mean_se` (mean reported standard error), `rmvar` (root mean reported
#' variance, the model based SE of Morris et al. 2019), `se_ratio` (relative
#' error in the model based SE, in percent), `coverage`,
#' `noncoverage_lower` (proportion of intervals entirely above the truth),
#' `noncoverage_upper` (entirely below), `be_coverage` (bias eliminated
#' coverage: intervals containing the replication's truth plus the mean
#' error, which for a constant truth is the mean estimate; it separates
#' miscalibrated width from bias; its MCSE comes from a jackknife that
#' re-estimates the mean error), `width`, `rejection` (labeled type
#' I error or power from the declared null and truth), `usable` and the
#' opt in composite `usable_coverage` (a usable interval that contains the
#' truth, over all generated data sets; this is not conditional coverage).
#'
#' @param ... Measure names or custom measures from [sim_measure()].
#' @return For `sim_measures()`, an object of class `adempr_measures`.
#' @examples
#' sim_measures("bias", "mse", "coverage")
#' sim_measure_info()
#' @export
sim_measures <- function(...) {
  ms <- list(...)
  # a plain list of custom measures (from the measure factories) is spliced
  ms <- do.call(c, lapply(ms, function(m) {
    if (is.list(m) && !inherits(m, "adempr_measure") && length(m) &&
        all(vapply(m, inherits, logical(1), "adempr_measure"))) unname(m) else list(m)
  }))
  if (!length(ms)) ms <- as.list(MEASURE_INFO$measure[MEASURE_INFO$default])
  out <- list()
  for (m in ms) {
    if (is.character(m)) {
      bad <- setdiff(m, MEASURE_INFO$measure)
      if (length(bad)) {
        stopf("Unknown performance measure(s): %s. See sim_measure_info().",
              paste(bad, collapse = ", "))
      }
      for (mm in m) out[[mm]] <- mm
    } else if (inherits(m, "adempr_measure")) {
      if (m$id %in% MEASURE_INFO$measure) {
        stopf("Custom measure id `%s` clashes with a built in measure.", m$id)
      }
      out[[m$id]] <- m
    } else {
      stopf("Measures must be names or sim_measure() objects.")
    }
  }
  structure(out, class = "adempr_measures")
}

#' @rdname sim_measures
#' @export
sim_measure_info <- function() MEASURE_INFO

#' Define a custom performance measure
#'
#' A custom measure computes one contribution per eligible replication and
#' aggregates them. With `aggregate = "mean"` the Monte Carlo standard
#' error is `sd / sqrt(R)`, valid for a mean of independent contributions.
#' For any other aggregation you must supply the uncertainty yourself:
#' `aggregate` must then return a list with `value` and `mcse` (which may be
#' `NA`); a generic `sd / sqrt(R)` is never attached to a nonlinear
#' statistic.
#'
#' @param id Identifier.
#' @param contribution Function of a data frame of eligible rows (columns
#'   `estimate`, `se`, `lower`, `upper`, `truth`, `p_value`, ...) returning
#'   one number per row. Return a logical vector for a proportion: it then
#'   gets an exact interval. Missing contributions give status
#'   `missing_contributions` and no value.
#' @param aggregate `"mean"` or a function of the contributions returning a
#'   list with `value` and `mcse`.
#' @param requires Character vector of required inputs among `"estimate"`,
#'   `"truth"`, `"se"`, `"interval"`, `"test"`.
#' @param direction Preferred direction, for display.
#' @param description Text.
#' @param estimand_type `NULL` (applies to every estimand), `"scalar"` or
#'   `"categorical"`: the type of estimand the measure applies to. Cells of
#'   other estimands get status `not_applicable`. Built in measures apply
#'   to scalar estimands; only `usable` also applies to categorical ones.
#' @return An object of class `adempr_measure`.
#' @examples
#' # probability that the absolute error exceeds 0.2
#' sim_measure("large_error",
#'   contribution = function(rows) abs(rows$estimate - rows$truth) > 0.2,
#'   requires = c("estimate", "truth"))
#' @export
sim_measure <- function(id, contribution, aggregate = "mean",
                        requires = c("estimate", "truth"), direction = "none",
                        description = NULL, estimand_type = NULL) {
  check_id(id, "id")
  if (!is.null(estimand_type) && !(is_scalar_character(estimand_type) &&
                                   estimand_type %in% c("scalar", "categorical"))) {
    stopf("`estimand_type` must be NULL, \"scalar\" or \"categorical\".")
  }
  check_function(contribution, "contribution")
  if (!identical(aggregate, "mean")) check_function(aggregate, "aggregate")
  requires <- match.arg(requires, c("estimate", "truth", "se", "interval", "test"),
                        several.ok = TRUE)
  structure(list(id = id, contribution = contribution, aggregate = aggregate,
                 requires = requires, direction = direction,
                 description = description, estimand_type = estimand_type),
            class = "adempr_measure")
}

# Whether a measure applies to an estimand of the given type. Built in
# measures treat the estimate as a number on a scale; for a class code only
# `usable` is meaningful.
measure_applies <- function(m, type) {
  type <- type %||% "scalar"
  if (inherits(m, "adempr_measure")) return(is.null(m$estimand_type) || identical(m$estimand_type, type))
  identical(type, "scalar") || identical(m, "usable")
}

kernel_result <- function(value = NA_real_, mcse = NA_real_, lower = NA_real_,
                          upper = NA_real_, n = 0L, status = "ok", note = NA_character_) {
  list(value = value, mcse = mcse, lower = lower, upper = upper, n = as.integer(n),
       status = status, note = note)
}

# Mean of independent contributions (bias, MAE, MSE, mean SE, width).
k_mean <- function(x, level = 0.95) {
  fin <- is.finite(x)
  k_mean_moments(n = length(x), mean = if (all(fin)) mean(x) else mean(x[fin]),
                 var = if (sum(fin) > 1) stats::var(x[fin]) else NA_real_,
                 pos_inf = sum(x == Inf, na.rm = TRUE), neg_inf = sum(x == -Inf, na.rm = TRUE),
                 level = level)
}

# The mean kernel from moments of the finite contributions and counts of
# infinite ones, so batch by batch accumulation gives the same result.
k_mean_moments <- function(n, mean, var, pos_inf = 0, neg_inf = 0, level = 0.95) {
  if (n == 0L) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  if (pos_inf + neg_inf > 0) {
    v <- if (pos_inf && neg_inf) NaN else if (pos_inf) Inf else -Inf
    return(kernel_result(value = v, n = n, status = "infinite",
                         note = sprintf("%d infinite contribution(s)", pos_inf + neg_inf)))
  }
  if (n == 1L) {
    return(kernel_result(value = mean, n = n, status = "insufficient",
                         note = "MCSE needs at least 2 replications"))
  }
  s <- sqrt(var / n)
  z <- stats::qnorm(1 - (1 - level) / 2)
  kernel_result(value = mean, mcse = s, lower = mean - z * s, upper = mean + z * s, n = n)
}

k_rmse <- function(e, level = 0.95) {
  k_root_mean_square(e, level, "RMSE")
}

# Root mean square of x, computed on a scaled copy so that squares cannot
# overflow or underflow when the answer itself is representable.
k_root_mean_square <- function(x, level, what) {
  n <- length(x)
  if (n == 0L) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  if (anyNA(x)) {
    return(kernel_result(n = n, status = "undefined", note = "undefined values (for example Inf - Inf)"))
  }
  if (any(is.infinite(x))) {
    return(kernel_result(value = Inf, n = n, status = "infinite", note = "infinite values"))
  }
  s <- max(abs(x))
  if (s == 0) {
    return(kernel_result(value = 0, n = n, status = "degenerate",
                         note = sprintf("%s is 0; the delta method MCSE is undefined", what)))
  }
  r <- k_mean((x / s)^2, level)
  out <- kernel_result(value = s * sqrt(r$value), n = n, status = r$status, note = r$note)
  if (!is.na(r$mcse)) {
    out$mcse <- s * r$mcse / (2 * sqrt(r$value))
    out$lower <- s * sqrt(max(r$lower, 0))
    out$upper <- s * sqrt(r$upper)
  }
  out
}

# Empirical SE. Var(s^2) = (mu4 - sigma^4 (n - 3) / (n - 1)) / n, then the
# delta method gives MCSE(s) = sqrt(Var(s^2)) / (2 s). Under normality this
# reduces approximately to s / sqrt(2 (n - 1)).
k_empse <- function(x, level = 0.95) {
  fin <- is.finite(x)
  xf <- x[fin]
  k_empse_moments(n = length(x), var = if (length(xf) > 1) stats::var(xf) else NA_real_,
                  m4 = if (length(xf)) mean((xf - mean(xf))^4) else NA_real_,
                  nonfinite = sum(!fin), level = level)
}

# The empirical SE kernel from the variance and the fourth central moment
# (mean of (x - mean)^4) of the contributions.
k_empse_moments <- function(n, var, m4, nonfinite = 0, level = 0.95) {
  if (n == 0L) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  if (nonfinite > 0) {
    return(kernel_result(value = NaN, n = n, status = "infinite",
                         note = "nonfinite estimates; the empirical SE is undefined"))
  }
  if (n < 2L) return(kernel_result(n = n, status = "insufficient", note = "needs at least 2 replications"))
  s <- sqrt(var)
  if (n < 4L) {
    return(kernel_result(value = s, n = n, status = "insufficient",
                         note = "MCSE needs at least 4 replications"))
  }
  if (s == 0) {
    return(kernel_result(value = 0, n = n, status = "degenerate",
                         note = "all estimates equal; MCSE undefined"))
  }
  v_s2 <- (m4 - s^4 * (n - 3) / (n - 1)) / n
  note <- NA_character_
  if (!is.finite(v_s2) || v_s2 <= 0) {
    v_s2 <- 2 * s^4 / (n - 1)
    note <- "kurtosis estimate unusable; normal theory MCSE used"
  }
  mc <- sqrt(v_s2) / (2 * s)
  z <- stats::qnorm(1 - (1 - level) / 2)
  kernel_result(value = s, mcse = mc, lower = max(s - z * mc, 0), upper = s + z * mc,
                n = n, note = note)
}

k_rmvar <- function(se, level = 0.95) {
  k_root_mean_square(se, level, "RMVAR")
}

# Relative error in model based SE, percent (Morris, White and Crowther
# 2019, Table 6). ModSE = sqrt(mean(se^2)).
k_se_ratio <- function(est, se, level = 0.95) {
  n <- length(est)
  if (n == 0L) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  if (any(!is.finite(est)) || any(!is.finite(se))) {
    return(kernel_result(value = NaN, n = n, status = "infinite",
                         note = "nonfinite estimates or SEs"))
  }
  if (n < 4L) return(kernel_result(n = n, status = "insufficient", note = "needs at least 4 replications"))
  # Work on a common scale so that squares cannot overflow or underflow.
  x <- est - mean(est)
  s <- max(abs(x), se)
  if (s == 0) return(kernel_result(n = n, status = "degenerate", note = "all estimates and SEs are equal to 0"))
  x <- x / s
  v <- (se / s)^2
  emp <- stats::sd(x)
  mod <- sqrt(mean(v))
  if (emp == 0) {
    return(kernel_result(n = n, status = "degenerate", note = "empirical SE is 0"))
  }
  if (mod == 0) {
    return(kernel_result(value = -100, n = n, status = "degenerate",
                         note = "every model SE is 0"))
  }
  val <- 100 * (mod / emp - 1)
  # Jackknife MCSE: leave-one-out values from running sums, so the
  # covariance between reported variances and the empirical variance is
  # included without assuming normal estimates.
  sx <- sum(x)
  sxx <- sum(x^2)
  sv <- sum(v)
  emp_i <- sqrt(((sxx - x^2) - (sx - x)^2 / (n - 1)) / (n - 2))
  mod_i <- sqrt((sv - v) / (n - 1))
  th <- 100 * (mod_i / emp_i - 1)
  mc <- if (all(is.finite(th))) sqrt((n - 1) / n * sum((th - mean(th))^2)) else NA_real_
  normal <- 100 * (mod / emp) * sqrt(stats::var(v) / (4 * n * mod^4) + 1 / (2 * (n - 1)))
  z <- stats::qnorm(1 - (1 - level) / 2)
  kernel_result(value = val, mcse = mc, lower = val - z * mc, upper = val + z * mc, n = n,
                note = sprintf("jackknife MCSE; the normal theory MCSE of Morris et al. (2019), which ignores the covariance of SEs and estimates, is %s",
                               format(signif(normal, 3))))
}

# Proportion with plug-in MCSE and an exact Clopper-Pearson interval, which
# stays informative when 0 or all events are observed.
k_prop <- function(x, level = 0.95) {
  x <- as.logical(x)
  k_prop_counts(sum(x), length(x), level)
}

# The proportion kernel from the number of events and of replications.
k_prop_counts <- function(k, n, level = 0.95) {
  if (n == 0L) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  p <- k / n
  a <- 1 - level
  lo <- if (k == 0) 0 else stats::qbeta(a / 2, k, n - k + 1)
  hi <- if (k == n) 1 else stats::qbeta(1 - a / 2, k + 1, n - k)
  status <- "ok"
  note <- NA_character_
  if (k == 0 || k == n) {
    status <- "boundary"
    note <- sprintf("%s events observed: the plug-in MCSE is 0 but the uncertainty is not; use the exact interval",
                    if (k == 0) "no" else "only")
  }
  r <- kernel_result(value = p, mcse = sqrt(p * (1 - p) / n), lower = lo, upper = hi,
                     n = n, status = status, note = note)
  r$proportion <- TRUE
  r
}

# Bias eliminated coverage: the proportion of intervals that contain the
# replication's truth shifted by the mean error (for a constant truth, the
# mean estimate). The same replications estimate that shift, so the
# indicators are dependent and the binomial interval does not apply. The
# MCSE and a normal interval come from a delete-a-group jackknife (at most
# 50 groups of consecutive rows) that re-estimates the shift each time.
k_be_coverage <- function(rows, sel, center_rows, truth, level = 0.95) {
  if (!any(sel)) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  err <- rows$estimate - truth
  if (!all(is.finite(err[center_rows]))) {
    return(kernel_result(n = sum(sel), status = "infinite",
                         note = "the mean error is not finite, so the bias cannot be removed"))
  }
  idx <- which(sel)
  cen <- which(center_rows)
  prop <- function(drop = integer()) {
    i <- setdiff(idx, drop)
    k <- setdiff(cen, drop)
    if (!length(i) || !length(k)) return(NA_real_)
    shift <- mean(err[k])
    mean(rows$lower[i] <= truth[i] + shift & truth[i] + shift <= rows$upper[i])
  }
  p <- prop()
  n <- length(idx)
  # The jackknife deletes whole replications: every row that enters the
  # coverage or the mean error, so both are recomputed without them.
  all <- sort(union(idx, cen))
  G <- min(length(all), 50L)
  mcse <- NA_real_
  if (G >= 2L) {
    grp <- split(all, cut(seq_along(all), G, labels = FALSE))
    th <- vapply(grp, prop, numeric(1))
    if (!anyNA(th)) mcse <- sqrt((G - 1) / G * sum((th - mean(th))^2))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  kernel_result(value = p, mcse = mcse, lower = max(0, p - z * mcse), upper = min(1, p + z * mcse),
                n = n, status = if (p %in% c(0, 1)) "boundary" else "ok",
                note = "intervals containing truth plus the mean error; jackknife MCSE re-estimates the mean error")
}

# Generic sd/sqrt(n) is attached only to the mean aggregate.
k_custom <- function(measure, rows, level = 0.95) {
  if (!nrow(rows)) return(kernel_result(status = "no_eligible", note = "no eligible replications"))
  x <- measure$contribution(rows)
  if (!(is.numeric(x) || is.logical(x)) || length(x) != nrow(rows)) {
    stopf("Custom measure `%s` must return one number per eligible row.", measure$id)
  }
  if (anyNA(x)) {
    return(kernel_result(n = length(x), status = "missing_contributions",
                         note = sprintf("%d of %d contributions are missing", sum(is.na(x)),
                                        length(x))))
  }
  # A proportion is declared by returning logical contributions; numbers
  # that happen to be 0 or 1 are averaged as numbers.
  binary <- is.logical(x)
  x <- as.double(x)
  if (identical(measure$aggregate, "mean")) {
    if (binary) return(k_prop(x == 1, level))
    return(k_mean(x, level))
  }
  r <- measure$aggregate(x)
  num_or_na <- function(v) is.null(v) || (length(v) == 1L && (is.numeric(v) || is.na(v)))
  if (!is.list(r) || !is_scalar_number(r$value) && !identical(r$value, NA_real_) ||
      !num_or_na(r$mcse) || !num_or_na(r$lower) || !num_or_na(r$upper) ||
      isTRUE(r$mcse < 0) || isTRUE(r$lower > r$upper)) {
    stopf("The aggregate of custom measure `%s` must return a list with a single numeric `value`, a nonnegative `mcse` (or NA) and optional `lower` <= `upper`.",
          measure$id)
  }
  mc <- r$mcse %||% NA_real_
  kernel_result(value = r$value, mcse = mc, lower = r$lower %||% NA_real_,
                upper = r$upper %||% NA_real_, n = length(x),
                status = if (is.na(mc)) "mcse_unavailable" else "ok",
                note = if (is.na(mc)) "no MCSE supplied for a nonlinear aggregate" else NA_character_)
}

# Does H0 hold at truth `theta` for null `h0` and alternative `alt`?
null_true <- function(theta, h0, alt) {
  ifelse(alt == "two.sided", abs(theta - h0) <= 1e-12 * pmax(1, abs(h0)),
  ifelse(alt == "less", theta >= h0, theta <= h0))
}

#' Measures for categorical estimands
#'
#' For an estimand of type `"categorical"` (a class coded as a number, such
#' as a selected model or an estimated number of components), returns
#' custom measures: `correct`, the proportion of replications in which the
#' estimate is the true class, `select_<class>`, the proportion in which
#' each declared class is chosen, and `select_other`, the proportion of
#' estimates outside the declared classes (an invalid output, not a
#' misclassification), so the selection proportions add up to one. Each
#' proportion has an exact interval. The measures apply to categorical
#' estimands only; for scalar estimands their cells are `not_applicable`.
#' A replication whose truth is not a declared class gives a missing
#' contribution to `correct`. Pass the measures to [sim_measures()] with
#' any other measures.
#'
#' @param classes Numeric codes of the classes.
#' @param labels Labels used in the measure ids; default the codes.
#' @return A list of [sim_measure()] objects.
#' @examples
#' sim_measures("usable", sim_measures_categorical(1:3))
#' @export
sim_measures_categorical <- function(classes, labels = classes) {
  if (!is.numeric(classes) || !length(classes) || anyNA(classes) || anyDuplicated(classes)) {
    stopf("`classes` must be distinct numbers.")
  }
  if (length(labels) != length(classes)) stopf("`labels` must have one entry per class.")
  if ("other" %in% as.character(labels)) stopf("`other` is reserved for estimates outside the classes.")
  out <- list(sim_measure("correct", local({
    classes <- classes
    function(rows) {
      ok <- rows$truth %in% classes
      ifelse(ok, rows$estimate %in% classes & rows$estimate == rows$truth, NA)
    }
  }), requires = c("estimate", "truth"), direction = "higher",
  description = "proportion of replications choosing the true class",
  estimand_type = "categorical"))
  for (i in seq_along(classes)) {
    cl <- classes[i]
    out[[i + 1L]] <- sim_measure(paste0("select_", labels[i]),
      local({
        cl <- cl
        function(rows) rows$estimate == cl
      }), requires = "estimate", description = sprintf("proportion choosing class %s", labels[i]),
      estimand_type = "categorical")
  }
  out[[length(out) + 1L]] <- sim_measure("select_other", local({
    classes <- classes
    function(rows) !(rows$estimate %in% classes)
  }), requires = "estimate", direction = "lower",
  description = "proportion of estimates outside the declared classes",
  estimand_type = "categorical")
  out
}

#' Interval score
#'
#' The interval score of Gneiting and Raftery (2007) for a central interval
#' at level `1 - alpha` rewards narrow intervals and penalizes misses in
#' proportion to their distance: `(u - l) + 2 / alpha * (l - y) 1(y < l) +
#' 2 / alpha * (y - u) 1(y > u)`, with `y` the truth and `alpha` taken from
#' each row's interval level. Lower is better. It is a mean of independent
#' contributions, so its MCSE is the usual one. Intervals are taken as
#' central intervals at their level; a one sided or unbounded interval has
#' an infinite width and therefore an infinite score, and only the side the
#' truth falls outside of is penalized. A row whose level is not strictly
#' between 0 and 1 has a missing contribution. With `scale`, each
#' contribution is divided by a per row scale (for example `"truth"` for a
#' relative score, or a number), so scenarios on different scales can be
#' compared; a row whose scale is zero, negative or not finite has a
#' missing contribution.
#'
#' @param id Measure id.
#' @param scale `NULL`, a column name of the estimates (for example
#'   `"truth"` or `"se"`), or a positive number.
#' @return A [sim_measure()].
#' @references Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
#'   prediction, and estimation. Journal of the American Statistical
#'   Association 102(477):359-378.
#' @examples
#' sim_measures("coverage", sim_measure_interval_score())
#' @export
sim_measure_interval_score <- function(id = "interval_score", scale = NULL) {
  if (!is.null(scale) && !(is_scalar_character(scale) ||
                           (is_scalar_number(scale) && is.finite(scale) && scale > 0))) {
    stopf("`scale` must be NULL, a column name or a positive finite number.")
  }
  sim_measure(id, function(rows) {
    a <- 1 - rows$level
    y <- rows$truth
    lo <- rows$lower
    up <- rows$upper
    # a penalty only on the side that misses, so an unbounded side adds
    # nothing and an unbounded interval scores Inf (its width)
    s <- (up - lo) + ifelse(y < lo, 2 / a * (lo - y), 0) + ifelse(y > up, 2 / a * (y - up), 0)
    s[!(is.finite(a) & a > 0 & a < 1)] <- NA_real_
    if (is.null(scale)) return(s)
    if (is.character(scale) && is.null(rows[[scale]])) {
      stopf("Interval score scale column `%s` is not a column of the estimates.", scale)
    }
    d <- rep_len(if (is.character(scale)) abs(rows[[scale]]) else scale, length(s))
    ifelse(is.finite(d) & d > 0, s / d, NA_real_)
  }, requires = c("interval", "truth"), direction = "lower",
  description = if (is.null(scale)) "interval score" else paste("interval score scaled by", scale))
}
