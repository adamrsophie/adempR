test_that("bias splits into statistical and identification parts", {
  st <- sim_study("lognormal", scenarios = sim_scenarios(n = 50L, mu = c(0, 1)),
    generate = function(scenario, context) data.frame(y = exp(rnorm(scenario$n, scenario$mu))),
    estimands = list(sim_estimand("median", sim_truth_exact(function(scenario) exp(scenario$mu)))),
    methods = list(mean = function(data) sim_estimate("median", mean(data$y))))
  run <- sim_run(st, replications = 200, seed = 1, path = run_path())
  d <- sim_bias_decomposition(run, sim_truth_exact(function(scenario) exp(scenario$mu + 0.5)))
  expect_equal(nrow(d), 2L)
  expect_equal(d$total_bias, d$statistical_bias + d$identification_bias)
  expect_equal(d$identification_bias, exp(c(0, 1) + 0.5) - exp(c(0, 1)))
  # the sample mean is unbiased for what it identifies
  expect_true(all(abs(d$statistical_bias) < 3 * d$statistical_mcse))
  expect_true(all(abs(d$total_bias) > 3 * d$total_mcse))
  expect_error(sim_bias_decomposition(run, function(scenario) c(1, 2)), "one number")
})

test_that("identities between estimands are checked on truths and estimates", {
  st <- sim_study("ratio", scenarios = sim_scenarios(n = 40L, mu = c(1, 2)),
    generate = function(scenario, context) data.frame(y = rexp(scenario$n, 1 / scenario$mu)),
    estimands = list(
      sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu)),
      sim_estimand("log_mean", sim_truth_exact(function(scenario) log(scenario$mu) + (scenario$mu == 2)))),
    methods = list(
      good = function(data) list(sim_estimate("mean", mean(data$y)), sim_estimate("log_mean", log(mean(data$y)))),
      bad = function(data) list(sim_estimate("mean", mean(data$y)), sim_estimate("log_mean", mean(log(data$y))))))
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  tr <- sim_check_identity(run, function(mean, log_mean) log(mean) - log_mean)
  expect_equal(tr$flag, c(FALSE, TRUE))
  es <- sim_check_identity(run, function(mean, log_mean) log(mean) - log_mean, target = "estimate")
  expect_false(any(es$flag[es$method == "good"]))
  expect_true(all(es$flag[es$method == "bad"]))
  expect_error(sim_check_identity(run, function(a) a), "estimand names")
})
