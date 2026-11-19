sbc_study <- function(methods) {
  sim_study("sbc", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) {
      mu <- rnorm(1)
      sim_data(data.frame(y = rnorm(scenario$n, mu)), truth = list(mu = mu))
    },
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = methods, measures = c("bias", "coverage"))
}

test_that("SBC accepts an exact posterior and rejects an overconfident one", {
  exact <- sim_method_sbc(function(data) {
    n <- nrow(data)
    rnorm(400, sum(data$y) / (n + 1), sqrt(1 / (n + 1)))
  }, "mu")
  narrow <- sim_method_sbc(function(data) {
    n <- nrow(data)
    rnorm(400, sum(data$y) / (n + 1), 0.3 * sqrt(1 / (n + 1)))
  }, "mu")
  run <- sim_run(sbc_study(list(exact = exact, narrow = narrow)), replications = 200, seed = 4,
                 path = run_path())
  good <- sim_sbc(run, "exact", bins = 10)
  bad <- sim_sbc(run, "narrow", bins = 10)
  expect_gt(good$tests$p_value, 0.01)
  expect_lt(bad$tests$p_value, 1e-6)
  expect_equal(sum(good$bins$count), 200L)
  expect_true(all(good$ranks$rank > 0 & good$ranks$rank < 1))
  # the estimates are posterior summaries
  s <- as.data.frame(sim_summarise(run))
  expect_lt(s$value[s$method == "narrow" & s$measure == "coverage"], 0.6)
  expect_output(print(good), "Simulation-based calibration")
  expect_error(sim_sbc(run, "nope"), "Unknown method")
})
