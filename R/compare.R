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
#' Only methods that target the estimand can be compared. `usable` uses the
#' replications attempted by both methods, so work that has not run yet is
#' not counted as failure. Infinite estimates excluded from a comparison
#' are counted in the note.
#'
#' A common truth cancels from differences in bias, so its uncertainty does
#' not enter. For MAE, MSE and coverage differences against a Monte Carlo
#' truth it does not cancel: the interval is conditional on the truth, and
#' `truth_sens_lower` and `truth_sens_upper` give the range of the
#' difference (for `rmse_ratio`, of the ratio) while the shared truth moves
#' by two of its standard errors for both methods together. Replications
#' are paired by number, and every scenario, comparison and measure has
#' its own bootstrap stream derived from `seed`, so the results do not
#' depend on the order of batches, scenarios or measures. Infinite
#' interval widths excluded from `width` are counted in the note. Bootstrap comparisons report how many
#' resamples were undefined and withhold the interval (status `unstable`)
#' when more than 1 percent were.
#'
#' @param run A run or results.
#' @param methods Two or more methods. The first is the reference; every
#'   other method is compared with it.
#' @param estimand Estimand identifier (defaults to the first).
#' @param measures Measures among `"bias"`, `"mae"`, `"mse"`,
#'   `"coverage"`, `"width"`, `"rejection"`, `"usable"`,
#'   `"relative_precision"`, `"empse_ratio"`, `"rmse_ratio"`.
#' @param margin Optional scientifically meaningful tolerance; the output
#'   flags whether the confidence interval lies within `(-margin, margin)`
#'   for differences and the percent change in precision, and within
#'   `(1 / (1 + margin), 1 + margin)` for the ratios `empse_ratio` and
#'   `rmse_ratio`. No equivalence is claimed (`within_margin = NA`) when
#'   the status is not `ok`, for example when every paired difference is
#'   equal (status `boundary`).
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
  untargeted <- methods[!vapply(methods, function(mn) estimand %in% study$methods[[mn]]$estimands,
                                logical(1))]
  if (length(untargeted)) {
    stopf("Method(s) %s do not target estimand %s and cannot be compared on it.",
          paste(untargeted, collapse = ", "), estimand)
  }
  if (!is.null(margin) && (!is_scalar_number(margin) || margin <= 0)) {
    stopf("`margin` must be a positive number.")
  }
  bad <- setdiff(measures, c(COMPARE_LINEAR, COMPARE_NONLINEAR))
  if (length(bad)) stopf("Unsupported comparison measure(s): %s.", paste(bad, collapse = ", "))
  est <- res$estimates
  reps <- res$replications
  level <- study$level
  z <- stats::qnorm(1 - (1 - level) / 2)
  if (!is_count(seed)) stopf("`seed` must be a single nonnegative whole number.")
  out <- list()
  ref <- methods[1]
  for (key in study$scenarios$keys) {
    # sorted, so the pairing and the bootstrap do not depend on the order in
    # which batches were read
    gen_ok <- sort(unique(reps$replication[reps$scenario_key == key & reps$generation_status == "ok"]))
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
        # each cell has its own bootstrap stream, keyed by what it compares
        boot_state <- rng_base_seed(seed, paste0("compare:", key, ":", estimand, ":", ref, ":",
                                                 other, ":", ms))
        rr <- compare_one(ms, b, a, study, z, bootstrap, boot_state)
        # Ratios are equal at 1: a margin m means the interval lies within
        # (1 / (1 + m), 1 + m). Differences (and the percent change in
        # precision) are equal at 0. No equivalence claim is made when the
        # interval is degenerate.
        within <- if (!is.null(margin) && !is.na(rr$lower) && rr$status == "ok") {
          if (ms %in% c("empse_ratio", "rmse_ratio")) {
            rr$lower > 1 / (1 + margin) && rr$upper < 1 + margin
          } else rr$lower > -margin && rr$upper < margin
        } else NA
        out[[length(out) + 1L]] <- data.frame(
          scenario_key = key,
          scenario_id = study$scenarios$labels[match(key, study$scenarios$keys)],
          estimand = estimand, method = other, reference = ref, measure = ms,
          contrast = rr$contrast, value = rr$value, mcse = rr$mcse, lower = rr$lower,
          upper = rr$upper, n_common = rr$n, n_valid_method = rr$n_b,
          n_valid_reference = rr$n_a, n_generated = length(gen_ok),
          correlation = rr$rho, within_margin = within,
          truth_sens_lower = rr$sens[1], truth_sens_upper = rr$sens[2], status = rr$status,
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
    # a replication the method has not attempted (no row yet) is not a failure
    usable = !is.na(r$method))
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
  notes <- character()
  if (ms %in% c("bias", "mae", "mse", "rmse_ratio", "relative_precision", "empse_ratio")) {
    nonfinite <- sum((a$valid %in% TRUE & !is.finite(a$estimate)) |
                       (b$valid %in% TRUE & !is.finite(b$estimate)))
    if (nonfinite) notes <- sprintf("%d replication(s) with an infinite estimate excluded", nonfinite)
  }
  if (ms == "width") {
    int_ok <- function(r) r$estimate_status %in% USABLE_FOR_INTERVAL & r$interval_status %in% "valid"
    infinite <- sum((int_ok(a) & !is.finite(a$upper - a$lower)) | (int_ok(b) & !is.finite(b$upper - b$lower)))
    if (infinite) notes <- sprintf("%d replication(s) with an infinite interval width excluded", infinite)
  }
  tse <- a$truth_se[common][1]
  mc_truth <- ms %in% c("mae", "mse", "coverage", "rmse_ratio") && any(common) &&
    identical(a$truth_kind[common][1], "mc") && is.finite(tse) && tse > 0
  if (mc_truth) {
    notes <- c(notes, sprintf("conditional on the Monte Carlo truth (SE %s), which does not cancel; see truth_sens_lower/upper",
                              format(signif(tse, 3))))
  }
  base <- list(n = n, n_a = sum(ea), n_b = sum(eb), rho = NA_real_, mcse = NA_real_,
               lower = NA_real_, upper = NA_real_, value = NA_real_, status = "ok",
               note = NA_character_, sens = c(NA_real_, NA_real_))
  finish <- function(base) {
    base$note <- paste(c(if (!is.na(base$note)) base$note, notes), collapse = "; ")
    if (!nzchar(base$note)) base$note <- NA_character_
    base
  }
  if (n < 2L) {
    base$status <- "insufficient"
    base$contrast <- "difference"
    base$note <- "fewer than 2 common valid replications"
    return(finish(base))
  }
  if (ms %in% COMPARE_LINEAR) {
    ac <- a[common, , drop = FALSE]
    bc <- b[common, , drop = FALSE]
    diff_at <- function(shift) {
      ac$truth <- ac$truth + shift
      bc$truth <- bc$truth + shift
      contribution_compare(bc, ms, study$alpha) - contribution_compare(ac, ms, study$alpha)
    }
    da <- contribution_compare(ac, ms, study$alpha)
    db <- contribution_compare(bc, ms, study$alpha)
    d <- db - da
    base$contrast <- "difference (method - reference)"
    base$value <- mean(d)
    units <- rep_units(study, a$replication[common])
    base$mcse <- if (anyDuplicated(units)) cluster_mcse(d, units) else stats::sd(d) / sqrt(n)
    if (anyDuplicated(units)) notes <- c(notes, "MCSE from antithetic pairs")
    base$lower <- base$value - z * base$mcse
    base$upper <- base$value + z * base$mcse
    base$rho <- suppressWarnings(stats::cor(da, db))
    if (all(d == d[1])) {
      # Every paired difference is equal (for proportions: no discordant
      # pair). The plug-in MCSE of 0 does not show that the methods agree.
      base$status <- "boundary"
      notes <- c(notes, "all paired differences are equal; the MCSE of 0 is not evidence of equivalence")
    }
    if (mc_truth) {
      # the shared truth moves both methods together
      v <- vapply(seq(-2, 2, length.out = 41L) * tse, function(s) mean(diff_at(s)), numeric(1))
      base$sens <- range(v)
    }
    if (ms == "bias") base$note <- "common truth cancels: difference of mean estimates"
    if (ms == "usable") base$note <- "replications attempted by both methods; not conditional on success"
    return(finish(base))
  }
  xa <- a$estimate[common]
  xb <- b$estimate[common]
  ta <- a$truth[common]
  # with a data set truth the precision comparisons use the errors, as the
  # empirical SE of sim_summarise() does
  ds <- identical(study$estimands[[a$estimand[1]]]$truth$kind, "dataset")
  pa <- if (ds) xa - ta else xa
  pb <- if (ds) xb - b$truth[common] else xb
  stat <- function(ia) {
    switch(ms,
      relative_precision = 100 * (stats::var(pa[ia]) / stats::var(pb[ia]) - 1),
      empse_ratio = stats::sd(pb[ia]) / stats::sd(pa[ia]),
      rmse_ratio = sqrt(mean((xb[ia] - ta[ia])^2) / mean((xa[ia] - ta[ia])^2)))
  }
  base$value <- stat(seq_len(n))
  base$rho <- suppressWarnings(stats::cor(pa, pb))
  base$contrast <- switch(ms,
    relative_precision = "percent increase in precision of method vs reference",
    empse_ratio = "ratio EmpSE(method) / EmpSE(reference)",
    rmse_ratio = "ratio RMSE(method) / RMSE(reference)")
  # antithetic pairs are resampled whole
  units <- rep_units(study, a$replication[common])
  boots <- with_rng_state(boot_state, vapply(seq_len(B), function(i) {
    stat(boot_index(units))
  }, numeric(1)))
  failed <- sum(!is.finite(boots))
  boots <- boots[is.finite(boots)]
  # Dropping undefined resamples conditions the bootstrap distribution; with
  # more than 1 percent of them the interval is withheld.
  if (length(boots) >= 10L && failed <= 0.01 * B) {
    q <- stats::quantile(boots, c((1 - study$level) / 2, 1 - (1 - study$level) / 2),
                         names = FALSE)
    base$lower <- q[1]
    base$upper <- q[2]
    base$mcse <- stats::sd(boots)
    base$note <- sprintf("paired bootstrap (%d resamples, %d undefined), percentile interval; MCSE is the bootstrap SD",
                         B, failed)
  } else {
    base$status <- if (failed > 0.01 * B) "unstable" else "bootstrap_failed"
    base$note <- sprintf("%d of %d bootstrap resamples undefined (zero denominator); no interval",
                         failed, B)
  }
  if (ms == "relative_precision" && is.finite(base$rho) && n > 2L) {
    ratio <- stats::var(xa) / stats::var(xb)
    analytic <- 200 * ratio * sqrt((1 - base$rho^2) / (n - 1))
    notes <- c(notes, sprintf("analytic MCSE %s (Morris et al. 2019; assumes bivariate normal estimates)",
                              format(signif(analytic, 3))))
  }
  if (mc_truth && ms == "rmse_ratio") {
    # the shared truth moves both methods together
    ta0 <- ta
    v <- vapply(seq(-2, 2, length.out = 41L) * tse, function(s) {
      ta <<- ta0 + s
      stat(seq_len(n))
    }, numeric(1))
    ta <- ta0
    base$sens <- range(v[is.finite(v)])
  }
  if (!is.finite(base$value)) base$status <- "undefined"
  finish(base)
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
