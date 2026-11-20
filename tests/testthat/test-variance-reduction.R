vr_study <- function(gen, scen = sim_scenarios(n = 20L, mu = c(0, 1))) {
  sim_study("vr", scenarios = scen, generate = gen,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))),
    measures = c("bias", "empse"))
}

test_that("common random numbers give scenarios that differ only in a shift the same draws", {
  st <- vr_study(sim_generator(function(scenario, context) data.frame(y = scenario$mu + rnorm(scenario$n)),
                               common = "n"))
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  e <- sim_results(run)$estimates
  k <- st$scenarios$keys
  expect_equal(e$estimate[e$scenario_key == k[2]] - e$estimate[e$scenario_key == k[1]], rep(1, 5))
  # without common random numbers the scenarios are independent
  st0 <- vr_study(function(scenario, context) data.frame(y = scenario$mu + rnorm(scenario$n)))
  e0 <- sim_results(sim_run(st0, replications = 5, seed = 1, path = run_path()))$estimates
  expect_false(isTRUE(all.equal(e0$estimate[e0$scenario_key == k[2]] - e0$estimate[e0$scenario_key == k[1]], rep(1, 5))))
  # the generator fingerprint changes, other generators keep theirs
  expect_false(identical(adempR:::generator_fingerprint(st$generator), adempR:::generator_fingerprint(st0$generator)))
  bad <- vr_study(sim_generator(function(scenario, context) data.frame(y = rnorm(scenario$n)), common = "nope"))
  expect_error(sim_run(bad, replications = 1, seed = 1, path = run_path()), "unknown scenario factor")
})

test_that("antithetic pairs share a stream and their summary uses pair means", {
  st <- vr_study(sim_generator(function(scenario, context) {
    z <- rnorm(scenario$n)
    if (context$antithetic) z <- -z
    data.frame(y = scenario$mu + z)
  }, antithetic = TRUE), scen = sim_scenarios(n = 20L, mu = 0))
  run <- sim_run(st, replications = 10, seed = 2, path = run_path())
  e <- sim_results(run)$estimates
  e <- e[order(e$replication), ]
  expect_equal(e$estimate[c(1, 3, 5)], -e$estimate[c(2, 4, 6)])
  a <- sim_summarise_antithetic(run, "bias")
  expect_equal(a$value, 0)
  expect_equal(a$mcse, 0)
  expect_equal(a$pairs, 5L)
  expect_error(sim_summarise_antithetic(sim_run(vr_study(function(scenario, context) data.frame(y = rnorm(5))),
                                                replications = 2, seed = 1, path = run_path())), "antithetic")
})

test_that("a control variate with known mean reduces the Monte Carlo error", {
  st <- vr_study(sim_generator(function(scenario, context) {
    x <- rnorm(scenario$n)
    data.frame(y = scenario$mu + x + rnorm(scenario$n, 0, 0.2), x = x)
  }, record = function(data) c(xbar = mean(data$x))), scen = sim_scenarios(n = 20L, mu = 0))
  run <- sim_run(st, replications = 80, seed = 3, path = run_path())
  cv <- sim_control_variate(run, "xbar", 0, "bias")
  expect_lt(cv$mcse_adjusted, cv$mcse / 3)
  expect_gt(cv$variance_ratio, 9)
  expect_equal(cv$value, as.data.frame(sim_summarise(run))$value[1])
  expect_error(sim_control_variate(run, "nope", 0), "No recorded quantity")
})

test_that("selected contributions block the control variate adjustment and are counted", {
  # the method fails when the first observation is large, so the fits that
  # contribute are a selection whose control mean is not the known one
  st <- sim_study("vr", scenarios = sim_scenarios(n = 20L, mu = 0),
    generate = sim_generator(function(scenario, context) {
      x <- rnorm(scenario$n)
      data.frame(y = scenario$mu + x + rnorm(scenario$n, 0, 0.2), x = x)
    }, record = function(data) c(xbar = mean(data$x))),
    estimands = list(sim_estimand("mean", sim_truth_unknown())),
    methods = list(mean = function(data) {
      if (data$x[1] > 1) stop("large first value")
      sim_estimate("mean", mean(data$y))
    }), measures = c("bias"))
  run <- sim_run(st, replications = 40, seed = 4, path = run_path())
  e <- sim_results(run)$estimates
  n_ok <- sum(e$valid & is.finite(e$estimate))
  expect_lt(n_ok, 40)
  cv <- sim_control_variate(run, "xbar", 0, "mean")
  # the mean needs no truth: every successful fit contributes
  expect_equal(cv$n, n_ok)
  expect_equal(cv$n_generated, 40L)
  expect_identical(cv$status, "selected")
  expect_true(is.na(cv$value_adjusted) && is.na(cv$mcse_adjusted) && is.na(cv$coefficient))
  expect_match(cv$note, sprintf("%d of 40", n_ok))
  expect_equal(cv$value, mean(e$estimate[e$valid & is.finite(e$estimate)]))
})

test_that("antithetic summaries count the pairs left incomplete by a failed member", {
  st <- sim_study("vr", scenarios = sim_scenarios(n = 20L, mu = 0),
    generate = sim_generator(function(scenario, context) {
      z <- rnorm(scenario$n)
      if (context$antithetic) z <- -z
      data.frame(y = scenario$mu + z)
    }, antithetic = TRUE),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu))),
    methods = list(mean = function(data) {
      if (data$y[1] > 1) stop("large first value")
      sim_estimate("mean", mean(data$y))
    }), measures = c("bias"))
  run <- sim_run(st, replications = 20, seed = 5, path = run_path())
  e <- sim_results(run)$estimates
  ok <- e$replication[e$valid & is.finite(e$estimate)]
  pair <- (ok + 1L) %/% 2L
  a <- sim_summarise_antithetic(run, "bias")
  expect_equal(a$pairs, sum(table(pair) == 2L))
  expect_equal(a$incomplete_pairs, length(unique(pair)) - a$pairs)
  expect_gt(a$incomplete_pairs, 0L)
})

test_that("every summary, comparison and ranking of antithetic runs uses the pairs", {
  gen <- sim_generator(function(scenario, context) {
    z <- rnorm(scenario$n)
    if (context$antithetic) z <- -z
    data.frame(y = scenario$mu + z)
  }, antithetic = TRUE)
  st <- sim_study("anti", scenarios = sim_scenarios(n = 20L, mu = 0), generate = gen,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu))),
    methods = list(
      mean = function(data) {
        m <- mean(data$y); s <- sd(data$y) / sqrt(nrow(data))
        sim_estimate("mean", m, se = s, lower = m - 1.96 * s, upper = m + 1.96 * s, level = 0.95)
      },
      trim = function(data) sim_estimate("mean", mean(data$y, trim = 0.2))),
    measures = c("bias", "mse", "coverage", "empse"))
  run <- sim_run(st, replications = 40, seed = 3, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  e <- sim_results(run)$estimates
  em <- e[e$method == "mean", ]
  em <- em[order(em$replication), ]
  pair <- (em$replication + 1) %/% 2
  # linear measures: the delete-one-pair jackknife equals sd(pair means) / sqrt(pairs)
  pm_mse <- tapply((em$estimate - em$truth)^2, pair, mean)
  get <- function(ms, mt = "mean") s[s$method == mt & s$measure == ms, ]
  expect_equal(get("mse")$mcse, sd(pm_mse) / sqrt(20), tolerance = 1e-10)
  # the estimates of a pair are mirror images: every pair mean of the error is 0
  expect_equal(get("bias")$mcse, 0, tolerance = 1e-12)
  cov_x <- as.numeric(em$lower <= em$truth & em$truth <= em$upper)
  expect_equal(get("coverage")$mcse, sd(tapply(cov_x, pair, mean)) / sqrt(20), tolerance = 1e-10)
  expect_true(is.finite(get("empse")$mcse))
  expect_match(get("empse")$note, "antithetic pairs")
  ref <- sim_summarise_antithetic(run, "mse")
  expect_equal(get("mse")$mcse, ref$mcse[ref$method == "mean"], tolerance = 1e-10)
  # paired comparison: MCSE from pair means of the differences
  cmp <- as.data.frame(sim_compare(run, c("mean", "trim"), measures = "mse", bootstrap = 50))
  et <- e[e$method == "trim", ]
  et <- et[order(et$replication), ]
  d <- (et$estimate - et$truth)^2 - (em$estimate - em$truth)^2
  expect_equal(cmp$mcse[cmp$measure == "mse"], sd(tapply(d, pair, mean)) / sqrt(20), tolerance = 1e-10)
  # ranking resamples pairs, confidence sequences refuse
  expect_s3_class(sim_rank(run, "mse", B = 20), "data.frame")
  expect_error(sim_confidence_sequence(run, "coverage"), "independent replications")
})

test_that("antithetic summaries use the spread of nonzero pair means", {
  # the pair members are not exact mirrors, so the pair means vary
  st <- vr_study(sim_generator(function(scenario, context) {
    z <- rnorm(scenario$n)
    if (context$antithetic) z <- -z
    data.frame(y = scenario$mu + z + rnorm(scenario$n, 0, 0.5))
  }, antithetic = TRUE), scen = sim_scenarios(n = 20L, mu = 0))
  run <- sim_run(st, replications = 20, seed = 6, path = run_path())
  e <- sim_results(run)$estimates
  e <- e[order(e$replication), ]
  pm <- tapply(e$estimate - e$truth, (e$replication + 1L) %/% 2L, mean)
  a <- sim_summarise_antithetic(run, "bias")
  expect_gt(sd(pm), 0)
  expect_equal(a$value, mean(pm))
  expect_equal(a$mcse, sd(pm) / sqrt(10))
})
