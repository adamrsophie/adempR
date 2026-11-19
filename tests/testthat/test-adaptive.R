test_that("adaptive allocation gives replications where the Monte Carlo error is", {
  st <- sim_study("ad", scenarios = sim_scenarios(sd = c(1, 4)),
    generate = function(scenario, context) data.frame(y = rnorm(20, 0, scenario$sd)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))),
    measures = "bias")
  p <- run_path()
  run <- sim_run_adaptive(st, "bias", target = 0.06, seed = 1, path = p, initial = 20)
  s <- as.data.frame(sim_summarise(run))
  n <- stats::setNames(s$n_planned, s$scenario_id)
  expect_equal(n[["sd=1"]], 20L)
  expect_gt(n[["sd=4"]], 100L)
  expect_true(all(s$mcse <= 0.06 + 1e-12))
  expect_true(all(run$allocation$target_met[run$allocation$round == max(run$allocation$round)]))
  # per scenario counts survive reopening, and the extra data are what a
  # plain run of that many replications would give for that scenario
  re <- sim_results(sim_open(p))$estimates
  big <- sim_run(st, replications = n[["sd=4"]], seed = 1, path = run_path())
  be <- sim_results(big)$estimates
  k <- st$scenarios$keys[2]
  expect_equal(re$estimate[re$scenario_key == k], be$estimate[be$scenario_key == k])
  expect_equal(nrow(sim_flow(run)), 2L)
  expect_equal(sim_flow(run)$planned, unname(n[sim_flow(run)$scenario_id]))
})
