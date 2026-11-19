test_that("factor influence attributes variation to the factor that drives it", {
  st <- sim_study("fi", scenarios = sim_scenarios(n = c(10L, 40L), sd = c(1, 4), shift = c(0, 1)),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, scenario$shift, scenario$sd)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$shift))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))),
    measures = "empse")
  run <- sim_run(st, replications = 60, seed = 1, path = run_path())
  fi <- sim_factor_influence(run, "empse")
  sh <- stats::setNames(fi$anova$share, fi$anova$term)
  expect_gt(sh[["sd"]], sh[["shift"]])
  expect_equal(sum(fi$anova$share, na.rm = TRUE), 1)
  expect_equal(nrow(fi$cells), 8L)
  expect_output(print(fi), "does not predict")
  fi2 <- sim_factor_influence(run, "empse", factors = c("n", "sd"), interactions = TRUE)
  expect_true("n:sd" %in% fi2$anova$term)
  expect_error(sim_factor_influence(run, "empse", factors = "nope"), "factors")
})
