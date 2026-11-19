# Precision-driven stopping with predeclared looks -------------------------
#
# Construction: K predeclared looks at committed replication prefixes
# n_1 < ... < n_K and C monitored Bernoulli cells. At every look each cell
# gets an exact Clopper-Pearson interval at level 1 - alpha / (K C).
# By the union bound, with probability at least 1 - alpha all intervals at
# all looks contain their true proportions simultaneously, whatever the
# stopping rule. Stopping when all half widths are below the target
# therefore keeps the reported intervals valid. The guarantee is
# conservative and covers Bernoulli measures of independent replications
# only.

SEQUENTIAL_MEASURES <- c("coverage", "rejection", "usable", "noncoverage_lower",
                         "noncoverage_upper", "usable_coverage")

cp_interval <- function(k, n, conf) {
  a <- 1 - conf
  lo <- ifelse(k == 0, 0, stats::qbeta(a / 2, k, n - k + 1))
  hi <- ifelse(k == n, 1, stats::qbeta(1 - a / 2, k + 1, n - k))
  cbind(lower = lo, upper = hi)
}

#' Run a study with precision-driven stopping
#'
#' Runs the study in stages up to predeclared numbers of replications
#' (`looks`), and stops at the first look where every monitored Bernoulli
#' measure (coverage, rejection, usable result, noncoverage) has an
#' interval half width below `half_width` in every scenario and method.
#' Each look uses only the committed prefix of replications `1..n_k`, never
#' whichever jobs finished first.
#'
#' Intervals are exact Clopper-Pearson intervals at level
#' `1 - alpha / (K C)` for `K` looks and `C` monitored cells (an alpha
#' spending split by the union bound). With probability at least
#' `1 - alpha`, all intervals at all looks contain the true proportions
#' simultaneously, so the reported intervals remain valid after optional
#' stopping. Ordinary fixed sample intervals do not have this property. The
#' guarantee assumes independent replications and Bernoulli measures; it
#' does not cover continuous measures such as bias or MSE, Monte Carlo
#' truth, or adaptive changes to the scenarios, which are refused.
#' Stopping the Monte Carlo study is unrelated to stopping rules inside a
#' simulated clinical trial.
#'
#' @param study A [sim_study()].
#' @param looks Increasing replication counts at which to look.
#' @param half_width Target half width of every monitored interval.
#' @param seed,path Passed to [sim_run()].
#' @param measures Monitored measures (Bernoulli only).
#' @param alpha Familywise error rate across looks and cells.
#' @param ... Passed to [sim_run()].
#' @return An object of class `adempr_sequential` with the final `run`,
#'   the `looks` table (one row per look and cell) and `stopped_at`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(t = function(data) {
#'     tt <- t.test(data)
#'     sim_estimate("mean", mean(data), lower = tt$conf.int[1], upper = tt$conf.int[2],
#'                  level = 0.95)
#'   }))
#' sq <- sim_run_sequential(study, looks = c(100, 200, 400), half_width = 0.05,
#'                          seed = 1, path = tempfile("adempR-run-"))
#' sq
#' @export
sim_run_sequential <- function(study, looks, half_width, seed, path,
                               measures = "coverage", alpha = 0.05, ...) {
  stopifnot(inherits(study, "adempr_study"))
  if (!is.numeric(looks) || length(looks) < 1L || any(diff(looks) <= 0) ||
      any(looks < 2) || any(looks != round(looks))) {
    stopf("`looks` must be increasing whole numbers of at least 2.")
  }
  bad <- setdiff(measures, SEQUENTIAL_MEASURES)
  if (length(bad)) {
    stopf("Sequential stopping is implemented only for Bernoulli measures (%s); not for %s.",
          paste(SEQUENTIAL_MEASURES, collapse = ", "), paste(bad, collapse = ", "))
  }
  if (any(vapply(study$estimands, function(e) e$truth$kind == "mc", logical(1))) &&
      any(measures %in% c("coverage", "noncoverage_lower", "noncoverage_upper", "usable_coverage"))) {
    stopf("Sequential stopping with Monte Carlo truth is not supported: the truth's uncertainty is not part of the guarantee.")
  }
  if (!is_scalar_number(half_width) || half_width <= 0 || half_width >= 0.5) {
    stopf("`half_width` must be in (0, 0.5).")
  }
  K <- length(looks)
  cells <- sum(vapply(study$methods, function(m) length(m$estimands), integer(1))) *
    length(study$scenarios) * length(measures)
  conf <- 1 - alpha / (K * cells)
  tab <- list()
  stopped <- NA_integer_
  run <- NULL
  for (k in seq_len(K)) {
    run <- sim_run(study, replications = looks[k], seed = seed, path = path, ...)
    s <- as.data.frame(sim_summarise(run, measures = do.call(sim_measures, as.list(measures))))
    s <- s[s$measure %in% measures, , drop = FALSE]
    kk <- round(s$value * s$n_eligible)
    ci <- cp_interval(kk, s$n_eligible, conf)
    hw <- pmax(s$value - ci[, "lower"], ci[, "upper"] - s$value)
    hw[s$n_eligible == 0] <- Inf
    tab[[k]] <- data.frame(look = k, replications = looks[k], scenario_id = s$scenario_id,
                           method = s$method, estimand = s$estimand, measure = s$measure,
                           value = s$value, n_eligible = s$n_eligible,
                           lower = ci[, "lower"], upper = ci[, "upper"], half_width = hw,
                           stringsAsFactors = FALSE)
    if (all(hw <= half_width)) {
      stopped <- k
      break
    }
  }
  structure(list(run = run, looks = bind_rows(tab), stopped_at = stopped,
                 conf_level_per_interval = conf, alpha = alpha, cells = cells,
                 half_width = half_width),
            class = "adempr_sequential")
}

#' @export
print.adempr_sequential <- function(x, ...) {
  cat("<adempR sequential run>\n")
  if (is.na(x$stopped_at)) {
    cat(sprintf("Target half width %g not reached by the last look (%d replications).\n",
                x$half_width, max(x$looks$replications)))
  } else {
    cat(sprintf("Stopped at look %d (%d replications).\n", x$stopped_at,
                x$looks$replications[x$looks$look == x$stopped_at][1]))
  }
  cat(sprintf("Each interval at level %.6f; all %d cells at all looks hold jointly with probability >= %.3f.\n",
              x$conf_level_per_interval, x$cells, 1 - x$alpha))
  last <- x$looks[x$looks$look == max(x$looks$look), , drop = FALSE]
  print(last[c("scenario_id", "method", "measure", "value", "lower", "upper", "half_width")],
        row.names = FALSE)
  invisible(x)
}
