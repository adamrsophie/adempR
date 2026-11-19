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
#' @return An object of class `adempr_calibration` with `table` (one row
#'   per setting and scenario: value, MCSE, loss), `overall` (mean loss per
#'   setting across scenarios), `best` (the chosen settings, a one row data
#'   frame), and `run`.
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
  r <- sim_rank(run, measure = measure, estimand = estimand, B = 10L, seed = 1L)
  tab <- data.frame(setting = r$method, r[c("scenario_id", "value", "loss", "rank", "n_common")],
                    stringsAsFactors = FALSE)
  tab <- cbind(grid[match(tab$setting, ids), , drop = FALSE], tab)
  rownames(tab) <- NULL
  overall <- stats::aggregate(loss ~ setting, data = tab, FUN = mean)
  overall <- overall[order(overall$loss), , drop = FALSE]
  best_id <- overall$setting[1]
  best <- grid[match(best_id, ids), , drop = FALSE]
  rownames(best) <- NULL
  structure(list(table = tab, overall = overall, best = best, measure = measure,
                 estimand = estimand, run = run), class = "adempr_calibration")
}

#' @export
print.adempr_calibration <- function(x, ...) {
  cat(sprintf("Calibration on %s (%s) over %d settings, on calibration streams.\n",
              x$measure, x$estimand, nrow(x$overall)))
  cat("Chosen: ", paste(names(x$best), vapply(x$best, format, character(1)), sep = " = ",
                        collapse = ", "), "\n", sep = "")
  invisible(x)
}
