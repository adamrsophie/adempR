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
