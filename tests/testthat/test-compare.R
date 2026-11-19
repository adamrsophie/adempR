# Scientific test 2: paired comparisons retain the covariance.

test_that("paired differences agree with direct calculations", {
  st <- sim_study("pair",
    scenarios = sim_scenarios(n = 15L),
    generate = function(scenario, context) rt(scenario$n, df = 3),
    estimands = list(sim_estimand("center", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) {
      tt <- t.test(data)
      sim_estimate("center", mean(data), lower = tt$conf.int[1], upper = tt$conf.int[2],
                   level = 0.95)
    }, median = function(data) {
      w <- suppressWarnings(wilcox.test(data, conf.int = TRUE))
      sim_estimate("center", median(data), lower = w$conf.int[1], upper = w$conf.int[2],
                   level = 0.95)
    }))
  run <- sim_run(st, replications = 300, seed = 4, path = run_path())
  e <- sim_results(run)$estimates
  a <- e[e$method == "mean", ]
  b <- e[e$method == "median", ]
  stopifnot(identical(a$replication, b$replication))
  cmp <- sim_compare(run, c("mean", "median"), measures = c("mse", "coverage", "bias"))
  d <- b$estimate^2 - a$estimate^2
  r <- cmp[cmp$measure == "mse", ]
  expect_equal(r$value, mean(d))
  expect_equal(r$mcse, sd(d) / sqrt(300))
  # covariance retained: much smaller than treating the methods as independent
  expect_lt(r$mcse, sqrt(var(a$estimate^2) / 300 + var(b$estimate^2) / 300))
  cv <- cmp[cmp$measure == "coverage", ]
  ca <- as.numeric(a$lower <= 0 & 0 <= a$upper)
  cb <- as.numeric(b$lower <= 0 & 0 <= b$upper)
  expect_equal(cv$value, mean(cb - ca))
  expect_equal(cv$mcse, sd(cb - ca) / sqrt(300))
  bi <- cmp[cmp$measure == "bias", ]
  expect_equal(bi$value, mean(b$estimate - a$estimate))
})

test_that("relative precision matches Morris et al. and the bootstrap is reproducible", {
  st <- sim_study("rp",
    scenarios = sim_scenarios(n = 20L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("center", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) sim_estimate("center", mean(data)),
                   median = function(data) sim_estimate("center", median(data))))
  run <- sim_run(st, replications = 200, seed = 5, path = run_path())
  e <- sim_results(run)$estimates
  xa <- e$estimate[e$method == "mean"]
  xb <- e$estimate[e$method == "median"]
  c1 <- sim_compare(run, c("mean", "median"), measures = "relative_precision", bootstrap = 300)
  expect_equal(c1$value, 100 * (var(xa) / var(xb) - 1))
  rho <- cor(xa, xb)
  expect_equal(c1$mcse, 200 * (var(xa) / var(xb)) * sqrt((1 - rho^2) / 199))
  # the median is less precise than the mean for normal data
  expect_lt(c1$upper, 0)
  c2 <- sim_compare(run, c("mean", "median"), measures = "relative_precision", bootstrap = 300)
  expect_identical(c1$lower, c2$lower)
})

test_that("Pareto view marks dominated methods without ranking scenarios", {
  st <- sim_study("pareto",
    scenarios = sim_scenarios(n = 20L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("center", sim_truth_exact(function(scenario) 0))),
    methods = list(mean = function(data) sim_estimate("center", mean(data)),
                   worse = function(data) sim_estimate("center", mean(data) + 1)))
  run <- sim_run(st, replications = 20, seed = 5, path = run_path())
  pv <- sim_pareto(sim_summarise(run, c("rmse", "usable")), c(rmse = "min", usable = "max"))
  expect_true(pv$pareto[pv$method == "mean"])
  expect_false(pv$pareto[pv$method == "worse"])
})
