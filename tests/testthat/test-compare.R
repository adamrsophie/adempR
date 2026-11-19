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
  analytic <- 200 * (var(xa) / var(xb)) * sqrt((1 - rho^2) / 199)
  expect_match(c1$note, paste("analytic MCSE", format(signif(analytic, 3))), fixed = TRUE)
  expect_match(c1$note, "bootstrap SD")
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

test_that("comparison edge cases are reported, not hidden", {
  m <- list(
    a = function(data) sim_estimate("mean", mean(data$y), lower = -10, upper = 10, level = 0.95),
    b = function(data) sim_estimate("mean", median(data$y), lower = -10, upper = 10, level = 0.95))
  st <- normal_study(n = 10L, mu = 0, methods = m)
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  cp <- as.data.frame(sim_compare(run, c("a", "b"), measures = c("coverage", "mse"),
                                  margin = 0.001))
  cv <- cp[cp$measure == "coverage", ]
  expect_identical(cv$status, "boundary")
  expect_true(is.na(cv$within_margin))
  expect_error(sim_compare(run, c("a", "b"), margin = -1), "margin")
  # a method that does not target the estimand is refused
  st2 <- sim_study("two", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("m", sim_truth_exact(function(scenario) 0)),
                     sim_estimand("v", sim_truth_exact(function(scenario) 1))),
    methods = list(a = sim_method(function(data) sim_estimate("m", mean(data)), estimands = "m"),
                   b = sim_method(function(data) sim_estimate("v", var(data)), estimands = "v")))
  r2 <- sim_run(st2, replications = 3, seed = 1, path = run_path())
  expect_error(sim_compare(r2, c("a", "b"), estimand = "m"), "do not target")
})

test_that("usable comparisons ignore replications a method has not attempted", {
  st <- normal_study(n = 10L, mu = 0)
  p <- run_path()
  sim_run(st, replications = 4, seed = 1, path = p)
  m2 <- c(st$methods, list(other = sim_method(function(data) sim_estimate("mean", mean(data$y)),
                                              estimands = "mean")))
  st2 <- normal_study(n = 10L, mu = 0, methods = m2)
  run <- sim_run(st2, replications = 4, seed = 1, path = p)
  # drop the new method's results for two replications
  res <- sim_results(run)
  keep <- !(res$estimates$method == "other" & res$estimates$replication > 2)
  res$estimates <- res$estimates[keep, ]
  cp <- as.data.frame(sim_compare(res, c("t_interval", "other"), measures = "usable"))
  expect_equal(cp$value, 0)
  expect_equal(cp$n_common, 2)
})
