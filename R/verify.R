# Reproducibility and seed sensitivity ------------------------------------

#' Compare two runs of the same study
#'
#' With two runs that used the same seed (1 and 1L are the same seed),
#' `sim_verify()` checks that every stored outcome is reproduced: estimate
#' rows are matched on scenario, estimand, replication and method, and
#' every estimate column (point, SE, interval, test, flags), the attempt
#' status, the warning flag, the truth and the validity status must agree;
#' generated data sets must have the same outcome and checksum. Missing
#' values equal missing values and equal infinities are equal. Rows present
#' in one run only make the check fail. This is the check for a rerun on
#' another machine or with another backend.
#'
#' With two runs that used different seeds, it is a seed sensitivity check:
#' each performance measure is compared between the runs with
#' `z = (a - b) / sqrt(mcse_a^2 + mcse_b^2)`, using `mcse_total` when the
#' truth is uncertain and, for proportions, the pooled binomial variance
#' with one pseudo event each way, so a proportion of 0 or 1 is not treated
#' as certain. The runs are independent, so about 5% of cells exceed `|z| >
#' 1.96` by chance (the result gives the expected number of flags); cells
#' beyond `threshold` suggest that the conclusions depend on the seed (too
#' few replications, a heavy tailed estimator, or a generator that does
#' not use the stream it was given). Cells that cannot be compared make
#' `ok` `NA` (indeterminate) when no other cell is flagged.
#'
#' @param a,b Runs, results or run directories.
#' @param measures Measures compared in the seed sensitivity mode; defaults
#'   to all measures of the study.
#' @param threshold Absolute z above which a cell is flagged.
#' @param tolerance Absolute difference under which estimates count as
#'   identical in the reproduction mode.
#' @return A list of class `adempr_verify` with `mode` (`"reproduce"` or
#'   `"seed_sensitivity"`), `table` and `ok` (`TRUE`, `FALSE` or, for an
#'   indeterminate seed comparison, `NA`).
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data))),
#'   measures = c("bias", "empse"))
#' a <- sim_run(study, replications = 30, seed = 1, path = tempfile("adempR-a-"))
#' b <- sim_run(study, replications = 30, seed = 2, path = tempfile("adempR-b-"))
#' sim_verify(a, b)
#' @export
sim_verify <- function(a, b, measures = NULL, threshold = 3, tolerance = 0) {
  if (!is_scalar_number(threshold) || threshold <= 0) stopf("`threshold` must be a positive number.")
  if (!is_scalar_number(tolerance) || tolerance < 0) stopf("`tolerance` must be a nonnegative number.")
  if (is.character(a)) a <- sim_open(a)
  if (is.character(b)) b <- sim_open(b)
  ra <- sim_results(a)
  rb <- sim_results(b)
  if (!identical(ra$manifest$study_fingerprint, rb$manifest$study_fingerprint)) {
    stopf("The runs belong to different studies (study fingerprints differ).")
  }
  # 1 and 1L are the same seed
  same_seed <- isTRUE(all.equal(as.numeric(ra$manifest$seed), as.numeric(rb$manifest$seed),
                                tolerance = 0)) &&
    identical(ra$manifest$rng_domain %||% "", rb$manifest$rng_domain %||% "")
  out <- if (same_seed) verify_reproduce(ra, rb, tolerance) else
    verify_seeds(ra, rb, measures, threshold)
  out$seeds <- c(a = ra$manifest$seed, b = rb$manifest$seed)
  class(out) <- "adempr_verify"
  out
}

# Equality of two columns: missing equals missing, equal infinities are
# equal, and finite numbers agree within the tolerance.
same_values <- function(x, y, tolerance) {
  if (is.numeric(x) && is.numeric(y)) {
    both_na <- is.na(x) & is.na(y)
    eq <- !is.na(x) & !is.na(y) & (x == y | (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance))
    return(both_na | eq)
  }
  x <- as.character(x)
  y <- as.character(y)
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

# Every stored outcome is compared on the complete set of keys: estimate
# columns, statuses and the truth of each estimate row, and the outcome
# and data checksum of each generated data set. Rows present in one run
# only are reported and make the check fail.
verify_reproduce <- function(ra, rb, tolerance) {
  keys <- c("scenario_key", "estimand", "replication", "method")
  cols <- c(setdiff(RAW_ESTIMATE_COLUMNS, "estimand"), "fallback", "attempt_status",
            "had_warning", "truth", "estimate_status")
  cmp <- function(xa, xb, keys, cols) {
    cols <- intersect(cols, intersect(names(xa), names(xb)))
    xa$.in_a <- TRUE
    xb$.in_b <- TRUE
    m <- merge(xa[c(keys, cols, ".in_a")], xb[c(keys, cols, ".in_b")], by = keys,
               suffixes = c("_a", "_b"), all = TRUE)
    both <- m$.in_a %in% TRUE & m$.in_b %in% TRUE
    same <- rep(TRUE, nrow(m))
    diff <- 0
    for (cn in cols) {
      x <- m[[paste0(cn, "_a")]]
      y <- m[[paste0(cn, "_b")]]
      same <- same & same_values(x, y, tolerance)
      if (is.numeric(x) && is.numeric(y)) {
        d <- abs(x - y)
        d[same_values(x, y, 0)] <- 0
        d[xor(is.na(x), is.na(y))] <- Inf
        diff <- max(c(diff, d[both]), na.rm = TRUE)
      }
    }
    list(m = m, both = both, same = same & both, diff = diff)
  }
  e <- cmp(ra$estimates, rb$estimates, keys, cols)
  g <- cmp(ra$replications, rb$replications, c("scenario_key", "replication"),
           c("generation_status", "data_md5"))
  tab <- data.frame(n_a = nrow(ra$estimates), n_b = nrow(rb$estimates),
                    n_compared = sum(e$both), n_identical = sum(e$same),
                    only_a = sum(!e$m$.in_b %in% TRUE), only_b = sum(!e$m$.in_a %in% TRUE),
                    datasets_differing = sum(!g$same), max_abs_diff = e$diff)
  mism <- e$m[!e$same, keys, drop = FALSE]
  list(mode = "reproduce", table = tab, mismatches = mism,
       ok = nrow(e$m) > 0 && all(e$same) && all(g$same))
}

# Seed sensitivity: each cell's difference against its Monte Carlo error.
# The truth error enters through mcse_total; a proportion uses the pooled
# binomial variance with one pseudo event in each direction, so a cell at 0
# or 1 is never treated as certain. Cells that cannot be compared (no
# MCSE, or present in one run only) make the verdict indeterminate (NA)
# unless another cell is flagged.
verify_seeds <- function(ra, rb, measures, threshold) {
  sa <- as.data.frame(sim_summarise(ra, measures = measures))
  sb <- as.data.frame(sim_summarise(rb, measures = measures))
  keys <- c("scenario_key", "scenario_id", "method", "estimand", "measure")
  take <- function(x) {
    x$se <- if (is.null(x$mcse_total)) x$mcse else ifelse(is.na(x$mcse_total), x$mcse, x$mcse_total)
    # with antithetic pairs the independent units are the pairs; a mean of
    # G independent values in [0, 1] with mean p has variance at most
    # p (1 - p) / G, so the pooled binomial variance below uses G
    if (!is.null(x$n_units)) x$n_eligible <- ifelse(is.na(x$n_units), x$n_eligible, x$n_units)
    x[c(keys, "value", "se", "n_eligible")]
  }
  m <- merge(take(sa), take(sb), by = keys, suffixes = c("_a", "_b"), all = TRUE)
  dv <- m$value_a - m$value_b
  se <- sqrt(m$se_a^2 + m$se_b^2)
  prop <- m$measure %in% PROPORTION_MEASURES & !is.na(dv) & m$n_eligible_a > 0 & m$n_eligible_b > 0
  if (any(prop)) {
    na <- m$n_eligible_a[prop]
    nb <- m$n_eligible_b[prop]
    pbar <- (m$value_a[prop] * na + m$value_b[prop] * nb + 1) / (na + nb + 2)
    se[prop] <- sqrt(pbar * (1 - pbar) * (1 / na + 1 / nb))
  }
  m$z <- ifelse(se > 0, dv / se, ifelse(dv == 0, 0, sign(dv) * Inf))
  m$flag <- !is.na(m$z) & abs(m$z) > threshold
  m$flag[is.na(m$z)] <- NA
  n_ok <- sum(!is.na(m$z))
  ok <- if (any(m$flag %in% TRUE)) FALSE else if (anyNA(m$flag)) NA else TRUE
  list(mode = "seed_sensitivity", table = m, ok = ok, n_indeterminate = sum(is.na(m$z)),
       expected_flags = n_ok * 2 * stats::pnorm(-threshold))
}

#' @export
print.adempr_verify <- function(x, ...) {
  if (x$mode == "reproduce") {
    t <- x$table
    cat(sprintf("Reproduction check (seed %s): %d of %d matched estimates identical, max difference %g.\n",
                x$seeds[[1]], t$n_identical, t$n_compared, t$max_abs_diff))
    if (t$only_a || t$only_b) cat(sprintf("Rows in one run only: %d in a, %d in b.\n", t$only_a, t$only_b))
    if (t$datasets_differing) cat(sprintf("%d generated data sets differ.\n", t$datasets_differing))
  } else {
    t <- x$table
    cat(sprintf("Seed sensitivity (seeds %s and %s): %d cells, %d beyond the threshold (%.2g expected by chance), %d not comparable, %.0f%% with |z| > 1.96.\n",
                x$seeds[[1]], x$seeds[[2]], nrow(t), sum(t$flag %in% TRUE), x$expected_flags,
                x$n_indeterminate, 100 * mean(abs(t$z) > 1.96, na.rm = TRUE)))
    if (any(t$flag %in% TRUE)) print(t[t$flag %in% TRUE, c("scenario_id", "method", "measure", "value_a", "value_b", "z")],
                                     row.names = FALSE)
  }
  invisible(x)
}
