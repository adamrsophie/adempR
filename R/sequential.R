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

# Estimands whose truth is only known approximately: a Monte Carlo truth,
# or a numerical truth whose error bound is not zero (or not stated). A
# guarantee for a truth dependent measure computed against them holds for
# the approximate reference value, not for the parameter.
uncertain_truths <- function(study, estimands = names(study$estimands)) {
  bad <- vapply(study$estimands[estimands], function(e) {
    k <- e$truth$kind
    identical(k, "mc") || (identical(k, "numeric") && !isTRUE(e$truth$error == 0))
  }, logical(1))
  names(bad)[bad]
}

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
#' The design (looks, target, measures, alpha and number of cells) and the
#' looks table are stored in the manifest after every look. A run directory
#' not started by `sim_run_sequential()`, or started with a different
#' design, is refused, and a target that cannot be reached even at the last
#' look is refused before any work. Results, summaries and reports of the
#' directory use the replications up to the analyzed look only, and
#' summaries and reports state that the run was sequential and where the
#' valid simultaneous intervals are.
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
  bad <- setdiff(measures, PROPORTION_MEASURES)
  if (length(bad)) {
    stopf("Sequential stopping is implemented only for Bernoulli measures (%s); not for %s.",
          paste(PROPORTION_MEASURES, collapse = ", "), paste(bad, collapse = ", "))
  }
  unc <- uncertain_truths(study)
  if (length(unc) && any(measures %in% c("coverage", "noncoverage_lower", "noncoverage_upper",
                                         "usable_coverage"))) {
    stopf("Sequential stopping of coverage needs an exact truth: estimand(s) %s have a Monte Carlo truth or a numerical truth whose error bound is not zero. The intervals at the looks would hold for coverage of the approximate reference value, not of the parameter.",
          paste(unc, collapse = ", "))
  }
  if (isTRUE(study$generator$antithetic)) {
    stopf("Sequential stopping needs independent replications; antithetic pairs are dependent, so the binomial intervals at the looks would not hold.")
  }
  if (!is_scalar_number(half_width) || half_width <= 0 || half_width >= 0.5) {
    stopf("`half_width` must be in (0, 0.5).")
  }
  if (!is_scalar_number(alpha) || alpha <= 0 || alpha >= 1) {
    stopf("`alpha` must be in (0, 1).")
  }
  looks <- as.integer(looks)
  K <- length(looks)
  cells <- sum(vapply(study$methods, function(m) length(m$estimands), integer(1))) *
    length(study$scenarios$keys) * length(measures)
  conf <- 1 - alpha / (K * cells)
  # The narrowest possible interval at the last look (all or no successes):
  # if even that is wider than the target, the design cannot stop early.
  best <- min(cp_half_widths(looks[K], conf))
  if (best > half_width) {
    stopf("The target half width %g cannot be reached even at the last look (%d replications): the narrowest interval there has half width %.4f. Add later looks or relax the target.",
          half_width, looks[K], best)
  }
  # the monitored cells are part of the design: the same count of other
  # cells would be a different design
  cell_ids <- sort(unlist(lapply(names(study$methods), function(mn) {
    outer(study$scenarios$keys, study$methods[[mn]]$estimands, paste, mn, sep = "|")
  })))
  design <- list(looks = looks, half_width = half_width, measures = measures,
                 alpha = alpha, cells = cells, conf_level_per_interval = conf,
                 cell_ids = cell_ids)
  check_sequential_directory(path, design)
  # the design is stored before the first replication runs, so a run
  # interrupted during its first look resumes as itself
  attr(study, "adempr_sequential") <- design
  tab <- list()
  stopped <- NA_integer_
  run <- NULL
  for (k in seq_len(K)) {
    run <- sim_run(study, replications = looks[k], seed = seed, path = path, ...)
    # Each look analyzes the committed prefix 1..looks[k] only.
    run$manifest$sequential <- list(design = design, analyzed = looks[k])
    s <- as.data.frame(sim_summarise(run, measures = do.call(sim_measures, as.list(measures))))
    s <- s[s$measure %in% measures, , drop = FALSE]
    kk <- round(s$value * s$n_eligible)
    ci <- cp_interval(kk, s$n_eligible, conf)
    hw <- pmax(s$value - ci[, "lower"], ci[, "upper"] - s$value)
    # an empty or indeterminate cell (for example usable coverage with
    # missing data set truths) is not yet precise
    hw[s$n_eligible == 0 | is.na(hw)] <- Inf
    tab[[k]] <- data.frame(look = k, replications = looks[k], scenario_id = s$scenario_id,
                           method = s$method, estimand = s$estimand, measure = s$measure,
                           value = s$value, n_eligible = s$n_eligible,
                           lower = ci[, "lower"], upper = ci[, "upper"], half_width = hw,
                           stringsAsFactors = FALSE)
    if (all(hw <= half_width)) stopped <- k
    # Recorded after every look, so an interrupted design resumes as itself.
    run$manifest$sequential <- list(design = design, analyzed = looks[k], stopped_at = stopped,
                                    looks = bind_rows(tab))
    update_manifest(path, function(man) {
      man$sequential <- run$manifest$sequential
      man
    })
    if (!is.na(stopped)) break
  }
  structure(list(run = run, looks = bind_rows(tab), stopped_at = stopped,
                 conf_level_per_interval = conf, alpha = alpha, cells = cells,
                 half_width = half_width),
            class = "adempr_sequential")
}

cp_half_widths <- function(n, conf) {
  k <- 0:n
  ci <- cp_interval(k, n, conf)
  p <- k / n
  pmax(p - ci[, "lower"], ci[, "upper"] - p)
}

# A sequential design is fixed when the directory is first used. A directory
# started by plain sim_run(), or one with a different design, is refused:
# analyzing it at the looks would use replications beyond the look, and a
# changed design would spend alpha again.
check_sequential_directory <- function(path, design) {
  if (!file.exists(run_paths(path)$manifest)) return(invisible())
  man <- read_manifest(path)
  old <- man$sequential$design
  if (is.null(old)) {
    stopf("%s was not started by sim_run_sequential(); a sequential design needs a new run directory.", path)
  }
  if (!identical(old, design)) {
    stopf("%s was started with a different sequential design (looks, target, measures, alpha or cells). Changing the design after looking spends alpha again; use a new run directory.", path)
  }
  invisible()
}

sequential_note <- function(sq) {
  d <- sq$design
  paste0("Sequential run: analyzed at ", sq$analyzed, " replications (",
         if (is.na(sq$stopped_at %||% NA)) "target not reached by the last look" else
           paste0("stopped at look ", sq$stopped_at, " of ", length(d$looks)),
         "). The intervals above are fixed sample intervals and do not keep their coverage after optional stopping; the valid simultaneous Clopper-Pearson intervals at level ",
         format(d$conf_level_per_interval, digits = 6), " per cell are in the sequential looks table.")
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
