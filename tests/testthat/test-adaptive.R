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

test_that("adaptive allocation never accepts a zero plug-in MCSE and records every call", {
  st <- sim_study("ad0", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(wide = function(data) sim_estimate("mean", mean(data$y), lower = -100, upper = 100)),
    measures = "coverage")
  p <- run_path()
  # coverage 1 has MCSE 0; the exact interval [0.83, 1] for 20 hits is far
  # wider than the target, so the scenario is extended to the cap
  run <- sim_run_adaptive(st, "coverage", target = 0.001, seed = 1, path = p, initial = 20,
                          max_replications = 60, max_rounds = 2)
  a <- run$allocation
  expect_false(a$target_met[1])
  hw <- (1 - stats::qbeta(0.025, 20, 1)) / (2 * stats::qnorm(0.975))
  expect_equal(a$precision[1], hw, tolerance = 1e-6)
  expect_equal(max(a$replications), 60L)
  expect_match(sim_open(p)$manifest$adaptive$inference, "Heuristic")
  # a second call appends to the history
  run2 <- sim_run_adaptive(st, "coverage", target = 0.5, seed = 1, path = p, initial = 20)
  h <- sim_open(p)$manifest$adaptive$allocation
  expect_equal(sort(unique(h$call)), 1:2)
  expect_equal(nrow(h), nrow(a) + 1L)
})

test_that("a scenario with an undefined MCSE is unresolved and is reported", {
  st <- sim_study("adna", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(none = function(data) sim_estimate("mean", NA_real_)),
    measures = "bias")
  expect_warning(run <- sim_run_adaptive(st, "bias", target = 0.1, seed = 1, path = run_path(),
                                         initial = 5), "undefined")
  expect_true(is.na(run$allocation$target_met[1]))
  expect_equal(run$allocation$replications, 5L)
})
