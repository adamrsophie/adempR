# Regressions for statistical invariants identified in the read-only contract.
# These tests intentionally fail adempR commit
# 814be78113bf0c52dfd7217550d75c2b881779d2. No source patch is applied.
# Run with testthat::test_file(...), after loading the installed adempR.
# M06 in findings.md is a scaling defect measured by a separate benchmark;
# this file deliberately contains no machine-dependent timing assertion.

library(testthat)
library(adempR)

.contract_stats_with_git <- function(fun) {
  old <- options(adempR.git_dir = tempdir())
  on.exit(options(old), add = TRUE)
  fun()
}

.contract_stats_study <- function(methods,
                              truth = sim_truth_exact(function(scenario) 0),
                              generate = function(scenario, context) rnorm(1),
                              measures = c("bias", "coverage")) {
  sim_study("contract_statistics", scenarios = sim_scenarios(n = 1L),
    generate = generate, estimands = list(sim_estimand("mu", truth)),
    methods = methods, measures = measures)
}

.contract_stats_run <- function(study, replications = 20, seed = 21, ...) {
  .contract_stats_with_git(function() {
    sim_run(study, replications = replications, seed = seed,
            path = tempfile("contract-stats-run-"), ...)
  })
}

.contract_stats_anti_generator <- function() {
  sim_generator(function(scenario, context) {
    z <- rnorm(2)
    if (context$antithetic) z <- -z
    data.frame(y = z[1], x = z[2])
  }, antithetic = TRUE, record = function(data) list(control = data$x))
}

.contract_stats_anti_study <- function() {
  .contract_stats_study(
    methods = list(m = function(data) {
      sim_estimate("mu", data$y, lower = data$y - 1, upper = data$y + 1,
                   level = .95)
    }), generate = .contract_stats_anti_generator(),
    measures = c("bias", "mse", "coverage"))
}

.contract_stats_pair_mcse <- function(x, replication) {
  pair_means <- tapply(x, (replication + 1L) %/% 2L, mean)
  stats::sd(pair_means) / sqrt(length(pair_means))
}

test_that("M01 streaming antithetic summaries use independent pairs across batches", {
  run <- .contract_stats_run(.contract_stats_anti_study(), chunk_size = 7)
  e <- sim_results(run)$estimates
  err <- e$estimate - e$truth
  expected <- c(
    bias = .contract_stats_pair_mcse(err, e$replication),
    mse = .contract_stats_pair_mcse(err^2, e$replication),
    coverage = .contract_stats_pair_mcse(
      as.numeric(e$lower <= e$truth & e$truth <= e$upper), e$replication))
  got <- sim_summarise_stream(run, c("bias", "mse", "coverage"))
  expect_equal(got$mcse, unname(expected[got$measure]), tolerance = 1e-12)
  expect_equal(expected[["bias"]], 0, tolerance = 1e-12)
  expect_gt(expected[["mse"]], 0)
})

test_that("M01 control variates do not invent variance in exactly mirrored pair means", {
  run <- .contract_stats_run(.contract_stats_anti_study(), chunk_size = 7)
  got <- sim_control_variate(run, "control", control_mean = 0, measure = "bias")
  expect_equal(got$value, 0, tolerance = 1e-12)
  expect_equal(got$value_adjusted, 0, tolerance = 1e-12)
  # Every pair mean of both outcome and control is identically zero,
  # whatever finite control coefficient is selected.
  expect_equal(got$mcse, 0, tolerance = 1e-12)
  expect_equal(got$mcse_adjusted, 0, tolerance = 1e-12)
})

test_that("M01 a grid of identical errors retains the pair-based MCSE", {
  grid <- c(0, 1)
  study <- sim_study("contract_grid_pairs", scenarios = sim_scenarios(n = 1L),
    generate = .contract_stats_anti_generator(),
    estimands = sim_estimand_grid("mu", grid,
      truth = function(scenario, at) 0),
    methods = list(m = function(data) {
      sim_estimate_grid("mu", grid, rep(data$y, length(grid)))
    }), measures = c("bias", "mse"))
  run <- .contract_stats_run(study, seed = 27, chunk_size = 7)
  e <- sim_results(run)$estimates
  e <- e[e$estimand == "mu_1", ]
  for (measure in c("bias", "mse")) {
    x <- e$estimate - e$truth
    if (measure == "mse") x <- x^2
    expected <- .contract_stats_pair_mcse(x, e$replication)
    got <- sim_summarise_grid(run, "mu", measure)
    avg <- got[got$measure == paste0("grid_mean_", measure), ]
    expect_equal(avg$mcse, expected, tolerance = 1e-12)
  }
})

test_that("M01 EXT03 bias decomposition respects antithetic independent units", {
  run <- .contract_stats_run(.contract_stats_anti_study())
  got <- sim_bias_decomposition(run, function(scenario) 0)
  e <- sim_results(run)$estimates
  expected <- .contract_stats_pair_mcse(e$estimate - e$truth, e$replication)
  expect_equal(got$total_mcse, expected, tolerance = 1e-12)
})

test_that("M02 all-events and no-events antithetic uncertainty is nondegenerate or unavailable", {
  # Fixed observed sample, not an assumption that its generating law is
  # degenerate. Either a valid bounded-pair interval or an explicit refusal
  # is acceptable; a raw binomial interval for four dependent rows is not
  # established by this test.
  for (covered in c(TRUE, FALSE)) {
    study <- .contract_stats_study(
      generate = .contract_stats_anti_generator(),
      methods = list(m = function(data, context) {
        x <- c(-.5, .5, -.25, .25)[context$replication] + if (covered) 0 else 10
        sim_estimate("mu", x, lower = x - 1, upper = x + 1, level = .95)
      }), measures = "coverage")
    got <- tryCatch(as.data.frame(sim_summarise(.contract_stats_run(study, 4))),
                    error = identity)
    if (inherits(got, "error")) {
      expect_match(conditionMessage(got), "antithetic|pair|boundary|interval|independent",
                   ignore.case = TRUE)
    } else {
      cov <- got[got$measure == "coverage", ]
      expect_equal(cov$value, as.numeric(covered))
      finite_interval <- is.finite(cov$lower) & is.finite(cov$upper)
      expect_true(all(!finite_interval | cov$upper > cov$lower),
        info = "A finite boundary sample must not be reported as a certainty interval [0,0] or [1,1].")
    }
  }
})

test_that("M03 equivalent coverage and usable coverage have equal positive pilot requirements", {
  study <- .contract_stats_study(
    methods = list(m = function(data) sim_estimate("mu", 10, lower = 9, upper = 11)),
    measures = c("coverage", "usable_coverage"))
  pilot <- .contract_stats_with_git(function() {
    sim_pilot(study, replications = 20, seed = 22,
      targets = c(coverage = .01, usable_coverage = .01),
      path = tempfile("contract-stats-pilot-"))
  })
  p <- pilot$precision
  cov <- p[p$measure == "coverage", ]
  usable_cov <- p[p$measure == "usable_coverage", ]
  expect_gt(usable_cov$required_eligible, 0)
  expect_equal(usable_cov$required_eligible, cov$required_eligible)
  expect_equal(usable_cov$required_replications, cov$required_replications)
})

test_that("M03 an all-zero usable coverage sample does not satisfy an arbitrarily small MCSE target", {
  study <- .contract_stats_study(
    methods = list(m = function(data) sim_estimate("mu", 10, lower = 9, upper = 11)),
    measures = "usable_coverage")
  run <- .contract_stats_with_git(function() {
    sim_run_adaptive(study, "usable_coverage", target = .001, seed = 23,
      path = tempfile("contract-stats-adaptive-"), initial = 20,
      max_replications = 100, max_rounds = 1)
  })
  first <- run$allocation[1, ]
  expect_gt(first$precision, 0)
  expect_false(first$target_met)
})

test_that("M04 unchanged run and results handles cannot bypass frozen CS tuning", {
  study <- .contract_stats_study(methods = list(m = function(data) {
    sim_estimate("mu", data, lower = data - 1, upper = data + 1, level = .95)
  }))
  run <- .contract_stats_run(study, 30, seed = 24)
  old_results <- sim_results(run)
  first <- sim_confidence_sequence(run, "coverage", t_opt = 20)
  expect_equal(first$t_opt, 20)
  expect_error(sim_confidence_sequence(run, "coverage", t_opt = 80),
               "tuned|tuning|fixed|before|same", ignore.case = TRUE)
  expect_error(sim_confidence_sequence(run, "coverage", t_opt = 20, alpha = .1),
               "alpha|tuning|fixed|same", ignore.case = TRUE)
  expect_error(sim_confidence_sequence(old_results, "coverage", t_opt = 80),
               "tuned|tuning|fixed|before|same", ignore.case = TRUE)
})

.contract_stats_approximate_truth_study <- function() {
  .contract_stats_study(
    truth = sim_truth_numeric(function(scenario) .5, error = .5),
    generate = function(scenario, context) 0,
    methods = list(m = function(data) {
      sim_estimate("mu", .5, lower = .4, upper = .6, level = .95)
    }), measures = "coverage")
}

.contract_stats_expect_truth_contract <- function(result) {
  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), "truth|reference|uncertain|numeric|approxim",
                 ignore.case = TRUE)
    return(invisible(NULL))
  }
  tab <- if (inherits(result, "adempr_sequential")) result$looks else as.data.frame(result)
  # Accept a safe refusal, a robust uncertainty envelope, or explicit
  # conditional-reference labeling. This deliberately does not prescribe
  # a new column name or require one particular statistical repair.
  robust_interval <- all(tab$lower <= 0 & tab$upper >= 1, na.rm = FALSE)
  public_text <- paste(c(capture.output(print(result)), capture.output(print(tab))), collapse = " ")
  explicit_target <- grepl("conditional|reference value|approximat|truth sensitiv|truth uncertaint",
                            public_text, ignore.case = TRUE)
  sens_lo <- grep("truth.*sens.*lower|sens.*truth.*lower", names(tab), value = TRUE)
  sens_hi <- grep("truth.*sens.*upper|sens.*truth.*upper", names(tab), value = TRUE)
  sensitivity_envelope <- length(sens_lo) > 0 && length(sens_hi) > 0 &&
    all(tab[[sens_lo[1]]] <= 0 & tab[[sens_hi[1]]] >= 1)
  expect_true(isTRUE(robust_interval) || explicit_target || isTRUE(sensitivity_envelope),
    info = paste("Coverage of numerical reference .5 is 1, but allowed true parameter 0 has coverage 0.",
                 "A guarantee must be refused, truth-sensitive, or explicitly conditional on the reference."))
}

test_that("M05 exact sequential coverage handles nonzero numerical truth uncertainty explicitly", {
  study <- .contract_stats_approximate_truth_study()
  result <- tryCatch(.contract_stats_with_git(function() {
    sim_run_sequential(study, looks = 100, half_width = .05, seed = 25,
      path = tempfile("contract-stats-sequential-"), measures = "coverage")
  }), error = identity)
  .contract_stats_expect_truth_contract(result)
})

test_that("M05 confidence sequence coverage handles nonzero numerical truth uncertainty explicitly", {
  run <- .contract_stats_run(.contract_stats_approximate_truth_study(), 100, seed = 25)
  result <- tryCatch(sim_confidence_sequence(run, "coverage"), error = identity)
  .contract_stats_expect_truth_contract(result)
})

test_that("S01 interval-only and test-only results appear in their own validity counts", {
  study <- .contract_stats_study(methods = list(m = function(data) {
    sim_estimate("mu", NA_real_, lower = -1, upper = 1, p_value = .2, level = .95)
  }), measures = c("coverage", "rejection"))
  summary <- sim_summarise(.contract_stats_run(study, 6, seed = 26))
  failures <- attr(summary, "failures")
  s <- as.data.frame(summary)
  expect_equal(s$n_eligible[s$measure == "coverage"], 6L)
  expect_equal(s$n_eligible[s$measure == "rejection"], 6L)
  expect_equal(failures$interval_valid, 6)
  expect_equal(failures$test_valid, 6)
  expect_equal(failures$usable, 0)
})

test_that("S02 finite representable empirical SD does not silently overflow as status ok", {
  study <- .contract_stats_study(methods = list(m = function(data, context) {
    sim_estimate("mu", c(-1, 1, -2, 2)[context$replication] * 1e200)
  }), measures = "empse")
  s <- as.data.frame(sim_summarise(.contract_stats_run(study, 4, seed = 28)))
  e <- s[s$measure == "empse", ]
  # Scaling gives a representable exact numerical reference. A clearly
  # signaled numerical failure is safer than reporting Inf/NaN as ok, but
  # full support should recover this finite result.
  finite_correct <- is.finite(e$value) &&
    abs(e$value / 1e200 - stats::sd(c(-1, 1, -2, 2))) < 1e-12
  flagged <- !e$status %in% c("ok", "boundary", "unstable") &&
    grepl("numer|overflow|finite", paste(e$status, e$note), ignore.case = TRUE)
  expect_true(finite_correct || flagged)
})
