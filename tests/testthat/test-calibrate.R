test_that("method settings are tuned on calibration streams, not production streams", {
  st <- sim_study("cal", scenarios = sim_scenarios(n = 30L),
    generate = function(scenario, context) data.frame(y = rt(scenario$n, df = 2)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))))
  cal <- sim_calibrate_method(st,
    method = function(trim) function(data) sim_estimate("mean", mean(data$y, trim = trim)),
    grid = list(trim = c(0, 0.25)), measure = "mse", replications = 60, seed = 3,
    path = run_path())
  expect_equal(cal$best$trim, 0.25)
  expect_equal(nrow(cal$table), 2L)
  expect_equal(cal$run$manifest$rng_domain, "calibration")
  # calibration data differ from production data with the same seed
  prod <- sim_run(st, replications = 5, seed = 3, path = run_path())
  pe <- sim_results(prod)$estimates
  ce <- sim_results(cal$run)$estimates
  ce <- ce[ce$method == "setting_1" & ce$replication <= 5, ]
  expect_false(isTRUE(all.equal(sort(pe$estimate), sort(ce$estimate))))
  expect_output(print(cal), "Chosen: trim = 0.25")
  expect_error(sim_calibrate_method(st, method = 1, grid = list(a = 1), seed = 1), "function")
})

test_that("calibration penalizes settings that fail, uses weights and reports uncertainty", {
  st <- sim_study("cal2", scenarios = sim_scenarios(n = c(8L, 40L), .weights = c(1, 3)),
    generate = function(scenario, context) data.frame(y = rt(scenario$n, df = 2)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))))
  # the heavy trim is best where it runs but fails in the small scenario
  cal <- sim_calibrate_method(st,
    method = function(trim) function(data) {
      if (trim > 0.3 && nrow(data) < 10) stop("too few observations")
      sim_estimate("mean", mean(data$y, trim = trim))
    },
    grid = list(trim = c(0, 0.2, 0.4)), measure = "mse", replications = 40, seed = 4,
    path = run_path())
  o <- cal$overall
  expect_true(all(c("loss", "mcse", "usable", "scenarios_without_loss") %in% names(o)))
  expect_identical(o$loss[o$setting == "setting_3"], Inf)
  expect_equal(o$scenarios_without_loss[o$setting == "setting_3"], 1)
  expect_lt(o$usable[o$setting == "setting_3"], 1)
  expect_false(identical(cal$best$trim, 0.4))
  # weighted mean of the per scenario losses
  t2 <- cal$table[cal$table$setting == "setting_2", ]
  expect_equal(o$loss[o$setting == "setting_2"], sum(c(0.25, 0.75) * t2$loss))
  expect_equal(o$mcse[o$setting == "setting_2"], sqrt(sum(c(0.25, 0.75)^2 * t2$mcse^2)))
  expect_true(all(is.finite(cal$table$mcse[cal$table$setting != "setting_3"])))
  expect_length(cal$excluded, 0L)
  expect_output(print(cal), "design weights")
})

test_that("calibration chooses nothing without an informative loss", {
  mk <- function(truth, method = function(data) sim_estimate("mean", mean(data$y)), w = NULL) {
    sc <- sim_scenarios(n = c(5L, 6L), .weights = w)
    sim_study("cal_none", scenarios = sc,
      generate = function(scenario, context) data.frame(y = stats::rnorm(scenario$n)),
      estimands = list(sim_estimand("mean", truth)),
      methods = list(avg = method), measures = "mse")
  }
  cal <- function(st, f = function(offset) function(data) sim_estimate("mean", mean(data$y) + offset))
    sim_calibrate_method(st, method = f, grid = list(offset = c(0, 1)), measure = "mse",
                         replications = 3, seed = 1, path = run_path())
  # every truth unknown
  expect_error(cal(mk(sim_truth_unknown())), "no setting can be chosen")
  # every setting fails
  expect_error(cal(mk(sim_truth_exact(function(scenario) 0)),
                   function(offset) function(data) stop("fails")), "can be chosen")
  # the only scenario with a loss has weight 0
  expect_error(cal(mk(sim_truth_exact(function(scenario) 0), w = c(0, 1)),
                   function(offset) function(data) {
                     if (nrow(data) == 6L) stop("fails")
                     sim_estimate("mean", mean(data$y) + offset)
                   }), "no setting can be chosen")
})
