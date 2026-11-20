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
#' `n_draws` equally spaced posterior draws per fit as diagnostics.
#' `sim_sbc()` ranks each truth among those draws (the number of draws
#' below it, 0 to `n_draws`), breaking ties at random with a draw fixed by
#' the scenario, replication and method, so a point mass posterior still
#' gives uniform ranks. Under a correct posterior and a generator that
#' draws the parameter from the prior, the rank is uniform on its
#' `n_draws + 1` values when the stored draws are independent; thin
#' correlated (MCMC) draws before returning them. The ranks are binned and
#' uniformity is tested per scenario with a chi-squared test whose expected
#' counts are the exact probabilities of the discrete rank in each bin.
#' Every scenario of the study is reported, starting from the
#' replications planned for it: those without a rank are counted as
#' excluded and split into failed generations, replications not generated,
#' failed fits and fits without a truth or a full set of draws. A scenario
#' without any rank keeps its row with no test (`p_value` is `NA`).
#'
#' @param fun Function of the data (and optionally `scenario`, `context`)
#'   returning a numeric vector of posterior draws.
#' @param estimand Estimand id.
#' @param level Level of the reported central interval.
#' @param n_draws Number of stored posterior draws (equally spaced over the
#'   returned draws); the rank takes `n_draws + 1` values. The method must
#'   return at least this many draws.
#' @param ... Passed to [sim_method()].
#' @return `sim_method_sbc()`: a method. `sim_sbc()`: an object of class
#'   `adempr_sbc` with `ranks` (one per fit), `bins` (counts and exact
#'   expected counts per scenario) and `tests` (per scenario: ranks `n`,
#'   `planned`, `excluded` and its parts `generation_failed`,
#'   `not_generated`, `fit_failed` and `no_truth_or_draws`, the chi-squared
#'   statistic, its p value and a note).
#' @references Talts S, Betancourt M, Simpson D, Vehtari A, Gelman A
#'   (2018). Validating Bayesian inference algorithms with simulation-based
#'   calibration. arXiv:1804.06788. Modrak M, Moon AH, Kim S, Burkner P,
#'   Huurre N, Faltejskova K, Gelman A, Vehtari A (2023). Simulation-based
#'   calibration checking for Bayesian computation: the choice of test
#'   quantities shapes sensitivity. Bayesian Analysis.
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
sim_method_sbc <- function(fun, estimand, level = 0.95, n_draws = 100L, ...) {
  check_function(fun, "fun")
  check_id(estimand, "estimand")
  if (!is_count(n_draws, 9)) stopf("`n_draws` must be a whole number of at least 9.")
  n_draws <- as.integer(n_draws)
  a <- (1 - level) / 2
  m <- sim_method(fun,
    extract = function(fit) {
      if (!is.numeric(fit) || length(fit) < 2L || anyNA(fit)) stop("The method must return posterior draws.")
      q <- stats::quantile(fit, c(a, 1 - a), names = FALSE)
      sim_estimate(estimand, mean(fit), se = stats::sd(fit), lower = q[1], upper = q[2],
                   level = level, interval_type = "posterior")
    },
    diagnostics = function(fit) {
      if (length(fit) < n_draws) {
        stop(sprintf("The method returned %d posterior draws; sim_method_sbc() stores %d.", length(fit), n_draws))
      }
      d <- fit[round(seq(1, length(fit), length.out = n_draws))]
      stats::setNames(d, sprintf("sbc_d%05d", seq_along(d)))
    },
    estimands = estimand, ...)
  # the number of stored draws, so that a run whose fits all failed still
  # knows the rank scale
  m$sbc_draws <- n_draws
  m
}

#' @rdname sim_method_sbc
#' @param run A run or results.
#' @param method Method wrapped with `sim_method_sbc()`.
#' @param bins Number of rank bins; at most `n_draws + 1`.
#' @export
sim_sbc <- function(run, method, estimand = NULL, bins = 20L) {
  res <- sim_results(run)
  study <- res$manifest$study
  if (!method %in% names(study$methods)) stopf("Unknown method `%s`.", method)
  estimand <- estimand %||% study$methods[[method]]$estimands[1]
  if (!is_count(bins, 2)) stopf("`bins` must be a whole number of at least 2.")
  dg <- res$diagnostics
  if (!is.null(dg)) dg <- dg[dg$method == method & grepl("^sbc_d", dg$name), , drop = FALSE]
  est <- res$estimates
  est <- est[est$method == method & est$estimand == estimand, , drop = FALSE]
  if ((is.null(dg) || !nrow(dg)) && any(est$attempt_status == "ok")) {
    stopf("Method `%s` stored no posterior draws; wrap it with sim_method_sbc().", method)
  }
  if (is.null(dg)) dg <- data.frame(scenario_key = character(), replication = integer(),
                                    name = character(), value = numeric())
  fits <- split(dg, paste(dg$scenario_key, dg$replication, sep = "\r"))
  L <- if (length(fits)) max(vapply(fits, nrow, integer(1))) else
    as.integer(study$methods[[method]]$sbc_draws %||% NA_integer_)
  if (!is.na(L) && bins > L + 1L) stopf("`bins` (%d) exceeds the %d possible ranks.", bins, L + 1L)
  no_truth <- character()
  rk <- lapply(fits, function(d) {
    tr <- est$truth[est$scenario_key == d$scenario_key[1] & est$replication == d$replication[1]]
    if (length(tr) != 1L || is.na(tr) || nrow(d) != L) {
      no_truth <<- c(no_truth, d$scenario_key[1])
      return(NULL)
    }
    x <- d$value
    below <- sum(x < tr)
    ties <- sum(x == tr)
    # ties are broken uniformly at random, with a draw fixed by the fit
    u <- strtoi(substr(md5_string(c("sbc", d$scenario_key[1], d$replication[1], method)), 1, 7), 16L) / 16^7
    data.frame(scenario_key = d$scenario_key[1], replication = d$replication[1],
               rank = below + floor(u * (ties + 1)), stringsAsFactors = FALSE)
  })
  ranks <- bind_rows(rk)
  if (is.null(ranks)) ranks <- data.frame(scenario_key = character(), replication = integer(),
                                          rank = integer(), stringsAsFactors = FALSE)
  ranks$scenario_id <- study$scenarios$labels[match(ranks$scenario_key, study$scenarios$keys)]
  # bin of each of the L + 1 rank values and the exact bin probabilities
  if (!is.na(L)) {
    bin_of <- floor((0:L) * bins / (L + 1)) + 1L
    prob <- tabulate(bin_of, bins) / (L + 1)
    ranks$bin <- bin_of[ranks$rank + 1L]
  } else {
    prob <- rep(1 / bins, bins)
    ranks$bin <- integer()
  }
  reps <- res$replications
  bt <- list()
  tests <- list()
  # every scenario of the study is reported, with the replications planned
  # for it and why each one gave no rank, so that failures stay visible
  for (k in study$scenarios$keys) {
    r <- ranks[ranks$scenario_key == k, , drop = FALSE]
    lab <- study$scenarios$labels[match(k, study$scenarios$keys)]
    planned <- unname(planned_replications(res$manifest, k))
    gk <- if (is.null(reps)) reps else reps[reps$scenario_key == k, , drop = FALSE]
    n_gen_failed <- if (is.null(gk)) 0L else sum(gk$generation_status != "ok")
    n_generated <- if (is.null(gk)) 0L else sum(gk$generation_status == "ok")
    with_draws <- length(unique(dg$replication[dg$scenario_key == k]))
    n_fit_failed <- max(length(unique(est$replication[est$scenario_key == k])) - with_draws, 0L)
    n_no_truth <- sum(no_truth == k)
    counts <- tabulate(r$bin, bins)
    exp_n <- nrow(r) * prob
    test_ok <- nrow(r) > 0L
    chi <- if (test_ok) sum((counts - exp_n)^2 / exp_n) else NA_real_
    bt[[length(bt) + 1L]] <- data.frame(scenario_id = lab, bin = seq_len(bins), count = counts,
                                        expected = exp_n, stringsAsFactors = FALSE)
    tests[[length(tests) + 1L]] <- data.frame(scenario_id = lab, n = nrow(r),
      planned = planned, excluded = max(planned - nrow(r), 0L),
      generation_failed = n_gen_failed, not_generated = max(planned - n_gen_failed - n_generated, 0L),
      fit_failed = n_fit_failed, no_truth_or_draws = n_no_truth,
      chisq = chi, df = bins - 1L,
      p_value = if (test_ok) stats::pchisq(chi, bins - 1L, lower.tail = FALSE) else NA_real_,
      note = if (!test_ok) "no usable ranks; no test" else
        if (min(exp_n) < 5) "fewer than 5 expected in some bin; use fewer bins" else NA_character_,
      stringsAsFactors = FALSE)
  }
  structure(list(ranks = ranks, bins = bind_rows(bt), tests = bind_rows(tests), method = method,
                 estimand = estimand, n_draws = L), class = "adempr_sbc")
}

#' @export
print.adempr_sbc <- function(x, ...) {
  cat(sprintf("Simulation-based calibration of %s (%s): %d fits.\n", x$method, x$estimand,
              nrow(x$ranks)))
  print(x$tests, row.names = FALSE, digits = 3)
  invisible(x)
}
