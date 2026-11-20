# Calibration of method settings -------------------------------------------

#' Tune a method setting on calibration streams
#'
#' Choosing a tuning parameter (a penalty, a bandwidth, a number of
#' iterations) on the same data sets that are later used to evaluate the
#' method favors it. `sim_calibrate_method()` evaluates a grid of settings
#' on data sets drawn from a separate stream family (the `"calibration"`
#' domain), so the production run never sees the data used for tuning. The
#' chosen setting and the grid are returned for the protocol (`settings_source`
#' of [sim_method()]).
#'
#' @param study A study; its methods are replaced by one method per setting.
#' @param method A function of the settings (named arguments) returning a
#'   method function or a [sim_method()].
#' @param grid A data frame (or a named list expanded with
#'   `expand.grid()`) of settings, one row per candidate.
#' @param measure Measure that decides: `"bias"` (absolute), `"mae"`,
#'   `"mse"`, `"rmse"`, `"empse"`, `"coverage"` (distance from the nominal
#'   level) or `"width"`; smaller loss is better.
#' @param estimand Estimand; defaults to the first.
#' @param replications Calibration replications per scenario.
#' @param seed Calibration seed.
#' @param path Directory for the calibration run.
#' @param ... Passed to [sim_run()].
#' Each setting is evaluated on its own usable replications in each
#' scenario (as in [sim_summarise()]); its usable proportion is reported
#' beside the loss, so a setting that often fails is visible. The loss is
#' averaged over scenarios with the design weights when the scenarios have
#' them, else equally. A setting without a loss in a scenario (it never
#' produced a usable result there) cannot be chosen; a scenario where no
#' setting has a loss (for example an unknown truth) is left out and
#' named in `excluded`. The MCSE of the mean loss treats scenarios as
#' independent; differences between settings evaluated on the same data
#' sets are usually more precise than the two MCSEs suggest (see
#' [sim_compare()] on the calibration run).
#'
#' @return An object of class `adempr_calibration` with `table` (one row
#'   per setting and scenario: value, MCSE, loss, usable proportion, number
#'   of contributions, status), `overall` (per setting: mean loss, its MCSE,
#'   mean usable proportion and the number of scenarios without a loss,
#'   ordered by loss), `best` (the chosen settings, a one row data frame),
#'   `excluded` (scenarios not used) and `run`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 30L),
#'   generate = function(scenario, context) rt(scenario$n, df = 3),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data))))
#' cal <- sim_calibrate_method(study,
#'   method = function(trim) function(data) sim_estimate("mean", mean(data, trim = trim)),
#'   grid = list(trim = c(0, 0.1, 0.2)), measure = "mse", replications = 50, seed = 7)
#' cal$best
#' @export
sim_calibrate_method <- function(study, method, grid, measure = "mse", estimand = NULL,
                                 replications = 100L, seed,
                                 path = tempfile("adempR-calibration-"), ...) {
  stopifnot(inherits(study, "adempr_study"))
  if (!is.function(method)) stopf("`method` must be a function of the settings.")
  measure <- match.arg(measure, RANK_MEASURES)
  if (is.list(grid) && !is.data.frame(grid)) grid <- expand.grid(grid, stringsAsFactors = FALSE)
  if (!is.data.frame(grid) || nrow(grid) < 2L || is.null(names(grid))) {
    stopf("`grid` must be a data frame or a named list giving at least two settings.")
  }
  estimand <- estimand %||% names(study$estimands)[1]
  if (!estimand %in% names(study$estimands)) stopf("Unknown estimand `%s`.", estimand)
  ids <- sprintf("setting_%d", seq_len(nrow(grid)))
  methods <- lapply(seq_len(nrow(grid)), function(i) {
    m <- do.call(method, as.list(grid[i, , drop = FALSE]))
    if (is.function(m)) m <- sim_method(m, estimands = estimand)
    if (!inherits(m, "adempr_method")) stopf("`method` must return a function or a sim_method().")
    m
  })
  names(methods) <- ids
  cal <- study
  cal$methods <- lapply(methods, function(m) {
    if (is.null(m$estimands)) m$estimands <- names(study$estimands)
    m
  })
  attr(cal, "adempr_derived") <- list(kind = "calibration", parent = derived_parent(study))
  run <- sim_run(cal, replications = replications, seed = seed, path = path,
                 label = "calibration", ...)
  res <- sim_results(run)
  sm <- as.data.frame(sim_summarise(res, measures = unique(c(measure, "usable"))))
  sm <- sm[sm$estimand == estimand & sm$scenario_key %in% cal$scenarios$keys, , drop = FALSE]
  est <- res$estimates
  keys <- cal$scenarios$keys
  rows <- list()
  for (k in keys) for (id in ids) {
    cell <- sm[sm$scenario_key == k & sm$method == id, , drop = FALSE]
    v <- cell[cell$measure == measure, , drop = FALSE]
    u <- cell$value[cell$measure == "usable"]
    nominal <- NA_real_
    if (measure == "coverage") {
      lv <- unique(est$level[est$scenario_key == k & est$method == id & est$estimand == estimand &
                               !is.na(est$level)])
      nominal <- if (length(lv) == 1L) lv else NA_real_
    }
    val <- if (nrow(v)) v$value[1] else NA_real_
    loss <- switch(measure, bias = abs(val), coverage = abs(val - nominal), val)
    rows[[length(rows) + 1L]] <- data.frame(
      setting = id, scenario_key = k,
      scenario_id = cal$scenarios$labels[match(k, keys)],
      value = val, mcse = if (nrow(v)) v$mcse[1] else NA_real_, loss = loss,
      usable = if (length(u)) u[1] else NA_real_, n = if (nrow(v)) v$n_eligible[1] else 0L,
      status = if (nrow(v)) v$status[1] else "no_results",
      stringsAsFactors = FALSE)
  }
  tab <- bind_rows(rows)
  # A scenario where no setting has a loss (for example an unknown truth)
  # cannot inform the choice and is reported. Elsewhere a setting without a
  # loss failed there and cannot be chosen, so a setting that fails in a
  # hard scenario is not rescued by dropping that scenario.
  informative <- vapply(keys, function(k) any(is.finite(tab$loss[tab$scenario_key == k])), logical(1))
  excluded <- cal$scenarios$labels[!informative]
  w <- cal$scenarios$weights %||% rep(1, length(keys))
  w <- stats::setNames(w / sum(w[informative]), keys)
  use <- tab$scenario_key %in% keys[informative]
  overall <- do.call(rbind, lapply(ids, function(id) {
    t <- tab[use & tab$setting == id, , drop = FALSE]
    wk <- w[t$scenario_key]
    data.frame(setting = id,
               loss = if (all(is.finite(t$loss))) sum(wk * t$loss) else Inf,
               mcse = sqrt(sum(wk^2 * t$mcse^2)),
               usable = sum(wk * t$usable),
               scenarios_without_loss = sum(!is.finite(t$loss)),
               stringsAsFactors = FALSE)
  }))
  if (!any(is.finite(overall$loss))) {
    stopf("No setting has a %s loss in every informative scenario; nothing can be chosen.", measure)
  }
  overall <- overall[order(overall$loss), , drop = FALSE]
  rownames(overall) <- NULL
  best_id <- overall$setting[1]
  best <- grid[match(best_id, ids), , drop = FALSE]
  rownames(best) <- NULL
  tab <- cbind(grid[match(tab$setting, ids), , drop = FALSE], tab)
  rownames(tab) <- NULL
  structure(list(table = tab, overall = overall, best = best, measure = measure,
                 estimand = estimand, excluded = excluded,
                 weighted = !is.null(cal$scenarios$weights), run = run),
            class = "adempr_calibration")
}

#' @export
print.adempr_calibration <- function(x, ...) {
  cat(sprintf("Calibration on %s (%s) over %d settings, on calibration streams.\n",
              x$measure, x$estimand, nrow(x$overall)))
  cat("Chosen: ", paste(names(x$best), vapply(x$best, format, character(1)), sep = " = ",
                        collapse = ", "), "\n", sep = "")
  o <- x$overall
  cat(sprintf("Mean loss (%s): %s.\n", if (isTRUE(x$weighted)) "design weights" else "equal weights",
              paste(sprintf("%s %s (MCSE %s, usable %s)", o$setting, format(signif(o$loss, 4)),
                            format(signif(o$mcse, 2)), format(round(o$usable, 3))), collapse = "; ")))
  if (length(x$excluded)) {
    cat(sprintf("Scenarios without a loss for any setting, not used: %s.\n",
                paste(x$excluded, collapse = ", ")))
  }
  invisible(x)
}
