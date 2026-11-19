# Paired method comparisons ------------------------------------------------

COMPARE_LINEAR <- c("bias", "mae", "mse", "coverage", "width", "rejection", "usable")
COMPARE_NONLINEAR <- c("relative_precision", "empse_ratio", "rmse_ratio")

#' Compare methods using the pairing of data sets
#'
#' All methods are applied to the same data sets, so performance
#' differences are estimated from per replication differences, which keeps
#' their covariance. For a linear measure with per replication
#' contributions `a_r` and `b_r`, the difference is `mean(D)` with
#' `D_r = a_r - b_r` and Monte Carlo standard error `sd(D) / sqrt(R)`.
#' Nonlinear comparisons (`relative_precision`, `empse_ratio`,
#' `rmse_ratio`) use a paired bootstrap over replications that resamples
#' both methods together; `relative_precision` also reports the analytic
#' standard error of Morris et al. (2019).
#'
#' The comparison population is explicit: only replications where both
#' methods meet the measure's eligibility rule are used (the common valid
#' population), and the counts of excluded replications are reported. Conditioning
#' on joint success is not the original data-generating distribution when
#' failures are informative; read it together with the `usable`
#' comparison, which uses all generated data sets.
#'
#' A common truth cancels from differences in bias, so its uncertainty does
#' not enter; for other measures the truth does not cancel and a
#' sensitivity range is not computed here (see [sim_summarise()]).
#'
#' @param run A run or results.
#' @param methods Two or more methods. The first is the reference; every
#'   other method is compared with it.
#' @param estimand Estimand identifier (defaults to the first).
#' @param measures Measures among `"bias"`, `"mae"`, `"mse"`,
#'   `"coverage"`, `"width"`, `"rejection"`, `"usable"`,
#'   `"relative_precision"`, `"empse_ratio"`, `"rmse_ratio"`.
#' @param margin Optional scientifically meaningful tolerance; the output
#'   flags whether the confidence interval of the difference lies within
#'   `(-margin, margin)`.
#' @param bootstrap Number of paired bootstrap resamples for nonlinear
#'   comparisons.
#' @param seed Seed for the bootstrap (uses its own L'Ecuyer-CMRG stream and
#'   leaves the caller's random number state unchanged).
#' @return An object of class `adempr_comparison` (a data frame) with one
#'   row per scenario, comparison and measure.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rt(scenario$n, df = 3),
#'   estimands = list(sim_estimand("center", sim_truth_exact(function(scenario) 0))),
#'   methods = list(
#'     mean = function(data) sim_estimate("center", mean(data)),
#'     median = function(data) sim_estimate("center", median(data))))
#' run <- sim_run(study, replications = 100, seed = 1, path = tempfile("adempR-run-"))
#' sim_compare(run, c("mean", "median"), measures = c("mse", "relative_precision"),
#'             bootstrap = 200)
#' @export
sim_compare <- function(run, methods, estimand = NULL,
                        measures = c("bias", "mse", "coverage", "usable"),
                        margin = NULL, bootstrap = 1000L, seed = 1L) {
  res <- sim_results(run)
  study <- res$manifest$study
  if (length(methods) < 2L || !all(methods %in% names(study$methods))) {
    stopf("`methods` must name at least two methods of the study.")
  }
  estimand <- estimand %||% names(study$estimands)[1]
  if (!(estimand %in% names(study$estimands))) stopf("Unknown estimand %s.", estimand)
  bad <- setdiff(measures, c(COMPARE_LINEAR, COMPARE_NONLINEAR))
  if (length(bad)) stopf("Unsupported comparison measure(s): %s.", paste(bad, collapse = ", "))
  est <- res$estimates
  reps <- res$replications
  level <- study$level
  z <- stats::qnorm(1 - (1 - level) / 2)
  boot_state <- parallel::nextRNGStream(rng_base_seed(seed))
  out <- list()
  ref <- methods[1]
  for (key in study$scenarios$keys) {
    gen_ok <- reps$replication[reps$scenario_key == key & reps$generation_status == "ok"]
    get <- function(mn) {
      r <- est[est$scenario_key == key & est$method == mn & est$estimand == estimand, ,
               drop = FALSE]
      r <- r[r$replication %in% gen_ok, , drop = FALSE]
      r[match(gen_ok, r$replication), , drop = FALSE]
    }
    a <- get(ref)
    for (other in methods[-1]) {
      b <- get(other)
      for (ms in measures) {
        rr <- compare_one(ms, b, a, study, z, bootstrap, boot_state)
        within <- if (!is.null(margin) && !is.na(rr$lower)) {
          rr$lower > -margin && rr$upper < margin
        } else NA
        out[[length(out) + 1L]] <- data.frame(
          scenario_key = key,
          scenario_id = study$scenarios$labels[match(key, study$scenarios$keys)],
          estimand = estimand, method = other, reference = ref, measure = ms,
          contrast = rr$contrast, value = rr$value, mcse = rr$mcse, lower = rr$lower,
          upper = rr$upper, n_common = rr$n, n_valid_method = rr$n_b,
          n_valid_reference = rr$n_a, n_generated = length(gen_ok),
          correlation = rr$rho, within_margin = within, status = rr$status,
          note = rr$note, stringsAsFactors = FALSE)
      }
    }
  }
  structure(bind_rows(out), class = c("adempr_comparison", "data.frame"))
}

eligible_compare <- function(r, ms) {
  has_truth <- !is.na(r$truth)
  int_ok <- r$estimate_status %in% USABLE_FOR_INTERVAL & r$interval_status %in% "valid"
  test_ok <- r$estimate_status %in% USABLE_FOR_INTERVAL & r$test_status %in% "valid"
  usable <- r$valid %in% TRUE
  switch(ms,
    bias = , mae = , mse = , rmse_ratio = usable & has_truth & is.finite(r$estimate),
    relative_precision = , empse_ratio = usable & is.finite(r$estimate),
    coverage = int_ok & has_truth,
    width = int_ok & is.finite(r$upper - r$lower),
    rejection = test_ok,
    usable = rep(TRUE, nrow(r)))
}

contribution_compare <- function(r, ms, alpha) {
  e <- r$estimate - r$truth
  switch(ms,
    bias = r$estimate,
    mae = abs(e),
    mse = e^2,
    coverage = as.numeric(r$lower <= r$truth & r$truth <= r$upper),
    width = r$upper - r$lower,
    rejection = as.numeric(reject_indicator(r, alpha)),
    usable = as.numeric(r$valid %in% TRUE))
}

compare_one <- function(ms, b, a, study, z, B, boot_state) {
  ea <- eligible_compare(a, ms)
  eb <- eligible_compare(b, ms)
  ea[is.na(ea)] <- FALSE
  eb[is.na(eb)] <- FALSE
  common <- ea & eb
  n <- sum(common)
  base <- list(n = n, n_a = sum(ea), n_b = sum(eb), rho = NA_real_, mcse = NA_real_,
               lower = NA_real_, upper = NA_real_, value = NA_real_, status = "ok",
               note = NA_character_)
  if (n < 2L) {
    base$status <- "insufficient"
    base$contrast <- "difference"
    base$note <- "fewer than 2 common valid replications"
    return(base)
  }
  if (ms %in% COMPARE_LINEAR) {
    da <- contribution_compare(a[common, , drop = FALSE], ms, study$alpha)
    db <- contribution_compare(b[common, , drop = FALSE], ms, study$alpha)
    d <- db - da
    base$contrast <- "difference (method - reference)"
    base$value <- mean(d)
    base$mcse <- stats::sd(d) / sqrt(n)
    base$lower <- base$value - z * base$mcse
    base$upper <- base$value + z * base$mcse
    base$rho <- suppressWarnings(stats::cor(da, db))
    if (ms == "bias") base$note <- "common truth cancels: difference of mean estimates"
    if (ms == "usable") base$note <- "all generated data sets; not conditional on success"
    return(base)
  }
  xa <- a$estimate[common]
  xb <- b$estimate[common]
  ta <- a$truth[common]
  stat <- function(ia) {
    switch(ms,
      relative_precision = 100 * (stats::var(xa[ia]) / stats::var(xb[ia]) - 1),
      empse_ratio = stats::sd(xb[ia]) / stats::sd(xa[ia]),
      rmse_ratio = sqrt(mean((xb[ia] - ta[ia])^2) / mean((xa[ia] - ta[ia])^2)))
  }
  base$value <- stat(seq_len(n))
  base$rho <- suppressWarnings(stats::cor(xa, xb))
  base$contrast <- switch(ms,
    relative_precision = "percent increase in precision of method vs reference",
    empse_ratio = "ratio EmpSE(method) / EmpSE(reference)",
    rmse_ratio = "ratio RMSE(method) / RMSE(reference)")
  boots <- with_rng_state(boot_state, vapply(seq_len(B), function(i) {
    stat(sample.int(n, n, replace = TRUE))
  }, numeric(1)))
  boots <- boots[is.finite(boots)]
  if (length(boots) >= 10L) {
    q <- stats::quantile(boots, c((1 - study$level) / 2, 1 - (1 - study$level) / 2),
                         names = FALSE)
    base$lower <- q[1]
    base$upper <- q[2]
    base$mcse <- stats::sd(boots)
    base$note <- sprintf("paired bootstrap (%d resamples), percentile interval", length(boots))
  } else {
    base$status <- "bootstrap_failed"
  }
  if (ms == "relative_precision" && is.finite(base$rho) && n > 2L) {
    ratio <- stats::var(xa) / stats::var(xb)
    base$mcse <- 200 * ratio * sqrt((1 - base$rho^2) / (n - 1))
    base$note <- paste(base$note, "; MCSE from Morris et al. (2019), interval from the bootstrap")
  }
  if (!is.finite(base$value)) base$status <- "undefined"
  base
}

#' @export
print.adempr_comparison <- function(x, digits = 4, ...) {
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  need <- c("scenario_id", "method", "reference", "measure", "value", "mcse", "lower",
            "upper", "n_common", "status")
  if (!all(need %in% names(df))) return(print(df, ...))
  cat("<adempR paired comparison>\n")
  show <- data.frame(scenario = df$scenario_id, method = df$method,
                     reference = df$reference, measure = df$measure,
                     value = fmt_num(df$value, digits), mcse = fmt_num(df$mcse, 2),
                     interval = paste0("[", fmt_num(df$lower, 3), ", ", fmt_num(df$upper, 3), "]"),
                     n_common = df$n_common, status = df$status, stringsAsFactors = FALSE)
  print(show, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Pareto view across several performance criteria
#'
#' For each scenario and estimand, marks methods that are not dominated on
#' the chosen criteria: a method is dominated when another method is at
#' least as good on every criterion and strictly better on one. This shows
#' trade offs (for example accuracy against reliability and cost) without
#' declaring a universal winner or averaging ranks across scenarios.
#' Monte Carlo uncertainty is not taken into account; compare with the
#' standard errors in the summary before interpreting small differences.
#'
#' @param summary A summary from [sim_summarise()].
#' @param criteria Named character vector: names are measures (or
#'   `"runtime"`) and values are `"min"` or `"max"`.
#' @return A data frame with one row per scenario, estimand and method,
#'   the criteria values and a logical `pareto` column.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rt(scenario$n, df = 3),
#'   estimands = list(sim_estimand("center", sim_truth_exact(function(scenario) 0))),
#'   methods = list(
#'     mean = function(data) sim_estimate("center", mean(data)),
#'     median = function(data) sim_estimate("center", median(data))))
#' run <- sim_run(study, replications = 50, seed = 1, path = tempfile("adempR-run-"))
#' sim_pareto(sim_summarise(run, c("rmse", "usable")),
#'            c(rmse = "min", usable = "max", runtime = "min"))
#' @export
sim_pareto <- function(summary, criteria = c(rmse = "min", usable = "max", runtime = "min")) {
  stopifnot(inherits(summary, "adempr_summary"))
  if (is.null(names(criteria)) || !all(criteria %in% c("min", "max"))) {
    stopf("`criteria` must be a named vector of \"min\" or \"max\".")
  }
  s <- as.data.frame(summary)
  rt <- attr(summary, "runtime")
  cells <- unique(s[s$scenario_key != "weighted", c("scenario_key", "scenario_id", "estimand", "method")])
  vals <- sapply(names(criteria), function(cr) {
    vapply(seq_len(nrow(cells)), function(i) {
      if (cr == "runtime") {
        if (is.null(rt)) return(NA_real_)
        v <- rt$mean_seconds[rt$scenario_key == cells$scenario_key[i] & rt$method == cells$method[i]]
        return(if (length(v)) v[1] else NA_real_)
      }
      v <- s$value[s$scenario_key == cells$scenario_key[i] & s$method == cells$method[i] &
                     s$estimand == cells$estimand[i] & s$measure == cr]
      if (length(v)) v[1] else NA_real_
    }, numeric(1))
  })
  if (is.null(dim(vals))) vals <- matrix(vals, nrow = nrow(cells), dimnames = list(NULL, names(criteria)))
  out <- cbind(cells, as.data.frame(vals))
  sign <- ifelse(criteria == "min", 1, -1)
  out$pareto <- NA
  grp <- paste(cells$scenario_key, cells$estimand)
  for (g in unique(grp)) {
    idx <- which(grp == g)
    m <- sweep(as.matrix(vals[idx, , drop = FALSE]), 2, sign, `*`)
    for (j in seq_along(idx)) {
      if (anyNA(m[j, ])) next
      dominated <- FALSE
      for (k in seq_along(idx)) {
        if (k == j || anyNA(m[k, ])) next
        if (all(m[k, ] <= m[j, ]) && any(m[k, ] < m[j, ])) {
          dominated <- TRUE
          break
        }
      }
      out$pareto[idx[j]] <- !dominated
    }
  }
  rownames(out) <- NULL
  out
}
