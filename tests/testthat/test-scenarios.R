
test_that("every scenario constructor rejects reserved and duplicated factor names", {
  for (nm in c("scenario_key", "scenario_id", "replication")) {
    l <- stats::setNames(list(1L), nm)
    expect_error(do.call(sim_scenarios, l), "reserved")
    expect_error(sim_scenarios_table(stats::setNames(data.frame(1L), nm)), "reserved")
    expect_error(sim_scenarios_sample(function() c(l, list(n = 1L)), 1, seed = 1,
                                      distribution = "point"), "reserved")
  }
  expect_error(sim_scenarios_sample(function() list(n = 1L, n = 2L), 1, seed = 1,
                                    distribution = "point"), "Duplicated")
  expect_error(sim_scenarios(weight = 1:2, .weights = c(1, 1)), "weight")
  # sampled scenarios keep the content keys in every public table
  sc <- sim_scenarios_sample(function() list(n = sample(5:9, 1)), 3, seed = 2, distribution = "uniform")
  st <- sim_study("keys", scenarios = sc, generate = function(scenario, context) stats::rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("mean", mean(data))))
  res <- sim_results(sim_run(st, 2, 1, run_path(), progress = FALSE))
  expect_setequal(unique(res$estimates$scenario_key), res$scenarios$scenario_key)
  expect_equal(as.data.frame(sc)$scenario_key, sc$keys)
})
