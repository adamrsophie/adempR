# Reproducibility and seed sensitivity ------------------------------------

#' Compare two runs of the same study
#'
#' With two runs that used the same seed, `sim_verify()` checks that every
#' stored estimate is reproduced: rows are matched on scenario, estimand,
#' replication and method, and the result reports how many are identical
#' and the largest absolute difference. This is the check for a rerun on
#' another machine or with another backend.
#'
#' With two runs that used different seeds, it is a seed sensitivity check:
#' each performance measure is compared between the runs with
#' `z = (a - b) / sqrt(mcse_a^2 + mcse_b^2)`. The runs are independent, so
#' about 5% of cells exceed `|z| > 1.96` by chance; cells beyond
#' `threshold` suggest that the conclusions depend on the seed (too few
#' replications, a heavy tailed estimator, or a generator that does not use
#' the stream it was given).
#'
#' @param a,b Runs, results or run directories.
#' @param measures Measures compared in the seed sensitivity mode; defaults
#'   to all measures of the study.
#' @param threshold Absolute z above which a cell is flagged.
#' @param tolerance Absolute difference under which estimates count as
#'   identical in the reproduction mode.
#' @return A list of class `adempr_verify` with `mode` (`"reproduce"` or
#'   `"seed_sensitivity"`), `table` and `ok`.
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
  if (is.character(a)) a <- sim_open(a)
  if (is.character(b)) b <- sim_open(b)
  ra <- sim_results(a)
  rb <- sim_results(b)
  if (!identical(ra$manifest$study_fingerprint, rb$manifest$study_fingerprint)) {
    stopf("The runs belong to different studies (study fingerprints differ).")
  }
  if (identical(ra$manifest$seed, rb$manifest$seed)) {
    keys <- c("scenario_key", "estimand", "replication", "method")
    ea <- ra$estimates
    eb <- rb$estimates
    m <- merge(ea[, c(keys, "estimate", "se", "lower", "upper")],
               eb[, c(keys, "estimate", "se", "lower", "upper")],
               by = keys, suffixes = c("_a", "_b"))
    d <- function(col) {
      x <- m[[paste0(col, "_a")]]
      y <- m[[paste0(col, "_b")]]
      same <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) <= tolerance)
      diff <- ifelse(is.na(x) | is.na(y), NA_real_, abs(x - y))
      list(same = same, diff = diff)
    }
    parts <- lapply(c("estimate", "se", "lower", "upper"), d)
    same <- Reduce(`&`, lapply(parts, `[[`, "same"))
    maxdiff <- suppressWarnings(max(unlist(lapply(parts, `[[`, "diff")), na.rm = TRUE))
    tab <- data.frame(n_a = nrow(ea), n_b = nrow(eb), n_compared = nrow(m),
                      n_identical = sum(same),
                      max_abs_diff = if (is.finite(maxdiff)) maxdiff else 0)
    out <- list(mode = "reproduce", table = tab, mismatches = m[!same, keys, drop = FALSE],
                ok = nrow(m) > 0 && all(same))
  } else {
    sa <- as.data.frame(sim_summarise(a, measures = measures))
    sb <- as.data.frame(sim_summarise(b, measures = measures))
    keys <- c("scenario_key", "scenario_id", "method", "estimand", "measure")
    m <- merge(sa[, c(keys, "value", "mcse")], sb[, c(keys, "value", "mcse")],
               by = keys, suffixes = c("_a", "_b"))
    dv <- m$value_a - m$value_b
    m$z <- dv / sqrt(m$mcse_a^2 + m$mcse_b^2)
    # Both MCSEs zero (e.g. every replication usable in both runs): equal
    # values agree exactly, different values disagree without doubt.
    zero <- !is.na(dv) & m$mcse_a %in% 0 & m$mcse_b %in% 0
    m$z[zero] <- ifelse(dv[zero] == 0, 0, sign(dv[zero]) * Inf)
    m$flag <- !is.na(m$z) & abs(m$z) > threshold
    out <- list(mode = "seed_sensitivity", table = m, ok = !any(m$flag))
  }
  out$seeds <- c(a = ra$manifest$seed, b = rb$manifest$seed)
  class(out) <- "adempr_verify"
  out
}

#' @export
print.adempr_verify <- function(x, ...) {
  if (x$mode == "reproduce") {
    t <- x$table
    cat(sprintf("Reproduction check (seed %s): %d of %d matched estimates identical, max difference %g.\n",
                x$seeds[[1]], t$n_identical, t$n_compared, t$max_abs_diff))
    if (t$n_a != t$n_b) cat(sprintf("Row counts differ: %d vs %d.\n", t$n_a, t$n_b))
  } else {
    t <- x$table
    cat(sprintf("Seed sensitivity (seeds %s and %s): %d cells, %d beyond the threshold, %.0f%% with |z| > 1.96.\n",
                x$seeds[[1]], x$seeds[[2]], nrow(t), sum(t$flag),
                100 * mean(abs(t$z) > 1.96, na.rm = TRUE)))
    if (any(t$flag)) print(t[t$flag, c("scenario_id", "method", "measure", "value_a", "value_b", "z")],
                           row.names = FALSE)
  }
  invisible(x)
}
