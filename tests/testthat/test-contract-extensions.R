# Focused regression tests for extension-contract findings at adempR commit
# 814be78113bf0c52dfd7217550d75c2b881779d2.
# These assert the intended behavior and intentionally fail on that commit.
library(adempR)

contract_ext_run_path <- function() tempfile("adempr_contract_extension_run_")

contract_ext_make_mean_study <- function(truth = sim_truth_exact(function(scenario) 0)) {
  sim_study("contract_extension_mean", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("mean", truth)),
    methods = list(avg = function(data) sim_estimate("mean", mean(data$y))),
    measures = c("bias", "mse"))
}

testthat::test_that("EXT01 calibration cannot select a setting without informative losses", {
  study <- contract_ext_make_mean_study(sim_truth_unknown())
  answer <- tryCatch(
    sim_calibrate_method(study,
      method = function(offset) {
        function(data) sim_estimate("mean", mean(data$y) + offset)
      },
      grid = list(offset = c(0, 1)), measure = "mse", replications = 3,
      seed = 12, path = contract_ext_run_path()),
    error = function(e) e)

  # Either an explicit unavailable/no-selection result or a relevant error
  # is acceptable. A miscellaneous setup error is not counted as success.
  if (inherits(answer, "error")) {
    testthat::expect_match(tolower(conditionMessage(answer)),
      "informative|loss|truth|select|choose")
  } else {
    no_choice <- is.null(answer$best) ||
      (is.data.frame(answer$best) && nrow(answer$best) == 0L)
    testthat::expect_true(no_choice,
      info = "All losses are unavailable: the first grid row is not a valid best setting.")
  }
})

testthat::test_that("EXT02 Weibull PH inversion avoids intermediate overflow and underflow", {
  for (linear_predictor in c(1000, -1000)) {
    set.seed(7)
    actual <- sim_rweibull_ph(rep(linear_predictor, 4), lambda = 1,
                              shape = 1000)$event_time
    set.seed(7)
    uniform <- runif(4)
    # Independent algebraic oracle: work on the log-time scale and
    # exponentiate only after the final representable value is obtained.
    expected <- exp((log(-log(uniform)) - linear_predictor) / 1000)
    testthat::expect_equal(actual, expected, tolerance = 1e-12,
      info = paste("Finite lp", linear_predictor, "has finite nonzero event times."))
  }
})

testthat::test_that("EXT03 bias decomposition uses independent antithetic pairs", {
  generator <- sim_generator(function(scenario, context) {
    z <- rnorm(scenario$n)
    if (context$antithetic) z <- -z
    data.frame(y = z)
  }, antithetic = TRUE)
  study <- sim_study("contract_extension_antithetic", scenarios = sim_scenarios(n = 5L),
    generate = generator,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(avg = function(data) sim_estimate("mean", mean(data$y))),
    measures = "bias")
  run <- sim_run(study, replications = 10, seed = 6, path = contract_ext_run_path())
  result <- sim_bias_decomposition(run, function(scenario) 0)
  estimates <- sim_results(run)$estimates
  pair <- (estimates$replication + 1L) %/% 2L
  pair_errors <- tapply(estimates$estimate - estimates$truth, pair, mean)
  expected <- sd(pair_errors) / sqrt(length(pair_errors))
  testthat::expect_equal(expected, 0, tolerance = 1e-12)
  testthat::expect_equal(result$total_mcse, expected, tolerance = 1e-12)
  testthat::expect_equal(result$statistical_mcse, expected, tolerance = 1e-12)
})

testthat::test_that("EXT04 SBC retains a wholly failed scenario and its denominator", {
  study <- sim_study("contract_extension_sbc", scenarios = sim_scenarios(fail = c(FALSE, TRUE)),
    generate = function(scenario, context) {
      mu <- rnorm(1)
      sim_data(data.frame(y = rnorm(3, mu), fail = scenario$fail), truth = list(mu = mu))
    },
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = list(m = sim_method_sbc(function(data) {
      if (data$fail[1]) stop("deliberate full-scenario failure")
      rnorm(20, sum(data$y) / 4, 0.5)
    }, "mu", n_draws = 9)), measures = "bias")
  run <- sim_run(study, replications = 6, seed = 3, path = contract_ext_run_path())
  result <- sim_sbc(run, "m", bins = 2)
  testthat::expect_setequal(result$tests$scenario_id, study$scenarios$labels)
  failed <- result$tests[result$tests$scenario_id == "fail=TRUE", , drop = FALSE]
  testthat::expect_equal(nrow(failed), 1L)
  if (nrow(failed) == 1L) {
    testthat::expect_equal(failed$n, 0)
    testthat::expect_equal(failed$excluded, 6)
    testthat::expect_true(is.na(failed$p_value))
  }
})

testthat::test_that("EXT05 a false nonlinear identity fails despite large raw inputs", {
  study <- sim_study("contract_extension_identity", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) data.frame(y = 0),
    estimands = list(
      sim_estimand("or", sim_truth_exact(function(scenario) 1e10)),
      sim_estimand("log_or", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) {
      list(sim_estimate("or", 1e10), sim_estimate("log_or", 0))
    }))
  run <- sim_run(study, replications = 2, seed = 1, path = contract_ext_run_path())
  result <- sim_check_identity(run, function(or, log_or) log(or) - log_or)
  testthat::expect_equal(result$residual, log(1e10))
  testthat::expect_identical(result$status, "fail")
  testthat::expect_true(result$flag)
})

testthat::test_that("EXT06 an aliased path cannot put the blinding key inside its export", {
  run <- sim_run(contract_ext_make_mean_study(), replications = 2, seed = 5,
                 path = contract_ext_run_path())
  export_dir <- tempfile("adempr_contract_blinded_")
  inside_alias <- file.path(dirname(export_dir), ".", basename(export_dir), "key.rds")
  testthat::expect_error(
    sim_export_blinded(run, dir = export_dir, key = inside_alias, seed = 1),
    regexp = "outside|inside|contain|apart|separate")
  testthat::expect_false(file.exists(file.path(export_dir, "key.rds")))
})
