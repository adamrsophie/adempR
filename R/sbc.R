# Simulation-based calibration ---------------------------------------------

#' Simulation-based calibration of a Bayesian method
#'
#' Simulation-based calibration (Talts et al. 2018) checks a Bayesian
#' implementation: draw the parameter from its prior, simulate data, fit
#' the posterior; if the implementation is correct, the rank of the drawn
#' parameter among the posterior draws is uniform. In adempR the generator
#' draws the parameter and returns it as a data set specific truth
#' (`sim_data(truth = )` with [sim_truth_dataset()]); the method is wrapped
#' with `sim_method_sbc()`, which turns a function returning posterior draws
#' into estimates (posterior mean, SD and central interval) and stores
#' `n_quantiles` posterior quantiles per fit as diagnostics. `sim_sbc()`
#' computes the normalized rank of each truth among those quantiles, bins
#' them, and tests uniformity (chi-squared over bins) per scenario.
#'
#' @param fun Function of the data (and optionally `scenario`, `context`)
#'   returning a numeric vector of posterior draws.
#' @param estimand Estimand id.
#' @param level Level of the reported central interval.
#' @param n_quantiles Number of stored posterior quantiles; the resolution
#'   of the rank.
#' @param ... Passed to [sim_method()].
#' @return `sim_method_sbc()`: a method. `sim_sbc()`: an object of class
#'   `adempr_sbc` with `ranks` (one per fit), `bins` (counts per scenario)
#'   and `tests` (chi-squared statistic and p value per scenario).
#' @references Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A
#'   (2018). Validating Bayesian inference algorithms with simulation-based
#'   calibration. arXiv:1804.06788.
#' @examples
#' # conjugate normal model: exact posterior draws are calibrated
#' study <- sim_study("sbc", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) {
#'     mu <- rnorm(1)
#'     sim_data(rnorm(scenario$n, mu), truth = list(mu = mu))
#'   },
#'   estimands = list(sim_estimand("mu", sim_truth_dataset())),
#'   methods = list(exact = sim_method_sbc(function(data) {
#'     n <- length(data)
#'     rnorm(500, sum(data) / (n + 1), sqrt(1 / (n + 1)))
#'   }, "mu")))
#' run <- sim_run(study, replications = 100, seed = 1, path = tempfile("adempR-run-"))
#' sim_sbc(run, "exact")
#' @export
sim_method_sbc <- function(fun, estimand, level = 0.95, n_quantiles = 99L, ...) {
  check_function(fun, "fun")
  check_id(estimand, "estimand")
  if (!is_count(n_quantiles, 9)) stopf("`n_quantiles` must be a whole number of at least 9.")
  probs <- seq_len(n_quantiles) / (n_quantiles + 1)
  a <- (1 - level) / 2
  sim_method(fun,
    extract = function(fit) {
      if (!is.numeric(fit) || length(fit) < 2L || anyNA(fit)) stop("The method must return posterior draws.")
      q <- stats::quantile(fit, c(a, 1 - a), names = FALSE)
      sim_estimate(estimand, mean(fit), se = stats::sd(fit), lower = q[1], upper = q[2],
                   level = level, interval_type = "posterior")
    },
    diagnostics = function(fit) {
      q <- stats::quantile(fit, probs, names = FALSE, type = 7)
      stats::setNames(q, sprintf("sbc_q%04d", seq_along(q)))
    },
    estimands = estimand, ...)
}

#' @rdname sim_method_sbc
#' @param run A run or results.
#' @param method Method wrapped with `sim_method_sbc()`.
#' @param bins Number of rank bins.
#' @export
sim_sbc <- function(run, method, estimand = NULL, bins = 20L) {
  res <- sim_results(run)
  study <- res$manifest$study
  if (!method %in% names(study$methods)) stopf("Unknown method `%s`.", method)
  estimand <- estimand %||% study$methods[[method]]$estimands[1]
  dg <- res$diagnostics
  if (!is.null(dg)) dg <- dg[dg$method == method & grepl("^sbc_q", dg$name), , drop = FALSE]
  if (is.null(dg) || !nrow(dg)) stopf("Method `%s` stored no posterior quantiles; wrap it with sim_method_sbc().", method)
  est <- res$estimates
  est <- est[est$method == method & est$estimand == estimand & !is.na(est$truth), , drop = FALSE]
  if (!is_count(bins, 2)) stopf("`bins` must be a whole number of at least 2.")
  rk <- lapply(split(dg, paste(dg$scenario_key, dg$replication, sep = "\r")), function(d) {
    tr <- est$truth[est$scenario_key == d$scenario_key[1] & est$replication == d$replication[1]]
    if (length(tr) != 1L) return(NULL)
    q <- d$value[order(d$name)]
    # normalized rank of the truth among L quantiles, in (0, 1)
    data.frame(scenario_key = d$scenario_key[1], replication = d$replication[1],
               rank = (sum(q < tr) + 0.5) / (length(q) + 1), stringsAsFactors = FALSE)
  })
  ranks <- bind_rows(rk)
  if (is.null(ranks)) stopf("No fit has both posterior quantiles and a truth.")
  ranks$scenario_id <- study$scenarios$labels[match(ranks$scenario_key, study$scenarios$keys)]
  ranks$bin <- pmin(bins, floor(ranks$rank * bins) + 1L)
  bt <- list()
  tests <- list()
  for (k in unique(ranks$scenario_key)) {
    r <- ranks[ranks$scenario_key == k, , drop = FALSE]
    counts <- tabulate(r$bin, bins)
    exp_n <- nrow(r) / bins
    chi <- sum((counts - exp_n)^2 / exp_n)
    lab <- r$scenario_id[1]
    bt[[length(bt) + 1L]] <- data.frame(scenario_id = lab, bin = seq_len(bins), count = counts,
                                        expected = exp_n, stringsAsFactors = FALSE)
    tests[[length(tests) + 1L]] <- data.frame(scenario_id = lab, n = nrow(r), chisq = chi,
      df = bins - 1L, p_value = stats::pchisq(chi, bins - 1L, lower.tail = FALSE),
      note = if (exp_n < 5) "fewer than 5 expected per bin; use fewer bins" else NA_character_,
      stringsAsFactors = FALSE)
  }
  structure(list(ranks = ranks, bins = bind_rows(bt), tests = bind_rows(tests), method = method,
                 estimand = estimand), class = "adempr_sbc")
}

#' @export
print.adempr_sbc <- function(x, ...) {
  cat(sprintf("Simulation-based calibration of %s (%s): %d fits.\n", x$method, x$estimand,
              nrow(x$ranks)))
  print(x$tests, row.names = FALSE, digits = 3)
  invisible(x)
}
