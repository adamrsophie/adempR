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

test_that("comparisons do not depend on row order and report excluded widths", {
  st <- sim_study("order",
    scenarios = sim_scenarios(n = c(10L, 20L)),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mu", sim_truth_exact(function(scenario) 0))),
    methods = list(
      a = function(data) sim_estimate("mu", mean(data), lower = mean(data) - 1, upper = mean(data) + 1,
                                      level = 0.95),
      b = function(data) {
        m <- median(data)
        # an unbounded interval now and then
        sim_estimate("mu", m, lower = m - 1, upper = if (data[1] > 1) Inf else m + 1, level = 0.95)
      }))
  run <- sim_run(st, replications = 40, seed = 3, path = run_path())
  res <- sim_results(run)
  cmp <- sim_compare(res, c("a", "b"), measures = c("empse_ratio", "width"), bootstrap = 200)
  rev_res <- res
  rev_res$replications <- res$replications[rev(seq_len(nrow(res$replications))), ]
  rev_res$estimates <- res$estimates[rev(seq_len(nrow(res$estimates))), ]
  expect_identical(sim_compare(rev_res, c("a", "b"), measures = c("empse_ratio", "width"),
                               bootstrap = 200), cmp)
  # each cell has its own stream: measures in another order give the same rows
  cmp2 <- sim_compare(res, c("a", "b"), measures = c("width", "empse_ratio"), bootstrap = 200)
  key <- function(x) x[order(x$scenario_key, x$measure), c("value", "lower", "upper")]
  expect_equal(key(cmp2), key(cmp), ignore_attr = TRUE)
  w <- cmp[cmp$measure == "width", ]
  n_inf <- vapply(st$scenarios$keys, function(k) {
    e <- res$estimates[res$estimates$scenario_key == k & res$estimates$method == "b", ]
    sum(!is.finite(e$upper))
  }, numeric(1))
  expect_true(any(n_inf > 0))
  for (i in seq_len(nrow(w))) {
    k <- n_inf[[w$scenario_key[i]]]
    if (k) expect_match(w$note[i], sprintf("%d replication(s) with an infinite interval width excluded", k),
                        fixed = TRUE)
  }
})

test_that("the RMSE ratio reports its sensitivity to a Monte Carlo truth", {
  mk <- function(x) data.frame(replication = seq_along(x), method = "m", estimate = x, truth = 0,
                               truth_se = 0.1, truth_kind = "mc", valid = TRUE,
                               estimate_status = "ok", interval_status = NA, test_status = NA,
                               lower = NA, upper = NA)
  set.seed(1)
  a <- mk(rnorm(50, 0.05))
  b <- mk(rnorm(50, 0.3, 2))
  study <- list(generator = list(), alpha = 0.05, level = 0.95)
  rr <- adempR:::compare_one("rmse_ratio", b, a, study, 1.96, 100, adempR:::rng_base_seed(1, "t"))
  grid <- vapply(seq(-0.2, 0.2, length.out = 41L), function(s) {
    sqrt(mean((b$estimate - s)^2) / mean((a$estimate - s)^2))
  }, numeric(1))
  expect_equal(rr$sens, range(grid))
  expect_true(rr$sens[1] < rr$value && rr$value < rr$sens[2])
})

test_that("the Pareto view leaves out weighted rows formed within factor levels", {
  st <- normal_study(n = c(10L, 20L), mu = c(0, 1), measures = c("bias", "mse"))
  st$scenarios$weights <- rep(0.25, 4)
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  pv <- sim_pareto(sim_summarise(run, c("mse", "usable"), weighted = TRUE, by = "mu"),
                   c(mse = "min", usable = "max"))
  expect_false(any(startsWith(pv$scenario_key, "weighted")))
  expect_setequal(unique(pv$scenario_key), st$scenarios$keys)
})
