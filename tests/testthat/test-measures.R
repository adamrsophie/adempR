# Reference kernels against hand calculations (expected values are written
# out from the formulas, never produced by the functions under test).

test_that("bias and MSE kernels match hand calculations", {
  e <- c(-1, 0, 1, 2)
  b <- adempR:::k_mean(e)
  expect_equal(b$value, 0.5)
  expect_equal(b$mcse, sqrt(5 / 3) / 2, tolerance = 1e-12)
  r <- adempR:::k_rmse(e)
  expect_equal(r$value, sqrt(1.5))
  # MCSE(MSE) = sd(e^2) / sqrt(R) = sqrt(3) / 2; delta method divides by 2 RMSE
  expect_equal(r$mcse, (sqrt(3) / 2) / (2 * sqrt(1.5)), tolerance = 1e-12)
})

test_that("empirical SE MCSE uses the kurtosis based variance of s^2", {
  x <- c(1, 2, 3, 4, 10)
  s2 <- 50 / 4
  m4 <- (81 + 16 + 1 + 0 + 1296) / 5
  v_s2 <- (m4 - s2^2 * (5 - 3) / (5 - 1)) / 5
  r <- adempR:::k_empse(x)
  expect_equal(r$value, sqrt(s2))
  expect_equal(r$mcse, sqrt(v_s2) / (2 * sqrt(s2)), tolerance = 1e-12)
  # for normal data it approaches the normal theory s / sqrt(2 (n - 1))
  set.seed(11)
  z <- rnorm(5000)
  rz <- adempR:::k_empse(z)
  expect_equal(rz$mcse, sd(z) / sqrt(2 * 4999), tolerance = 0.1)
})

test_that("relative error in model SE follows Morris et al. Table 6", {
  est <- c(1, 2, 3, 4, 10)
  se <- c(1, 1, 2, 2, 3)
  mod <- sqrt(mean(c(1, 1, 4, 4, 9)))
  emp <- sqrt(12.5)
  v_se2 <- var(c(1, 1, 4, 4, 9))
  r <- adempR:::k_se_ratio(est, se)
  expect_equal(r$value, 100 * (mod / emp - 1))
  normal <- 100 * (mod / emp) * sqrt(v_se2 / (4 * 5 * mod^4) + 1 / (2 * 4))
  expect_match(r$note, format(signif(normal, 3)), fixed = TRUE)
  # the MCSE is the jackknife over replications
  th <- vapply(1:5, function(i) 100 * (sqrt(mean(se[-i]^2)) / sd(est[-i]) - 1), numeric(1))
  expect_equal(r$mcse, sqrt(4 / 5 * sum((th - mean(th))^2)))
  # zero model SEs and huge values do not give NaN with status ok
  expect_identical(adempR:::k_se_ratio(1:4, rep(0, 4))$status, "degenerate")
  big <- adempR:::k_se_ratio(c(1, 2, 3, 4, 10) * 1e200, se * 1e200)
  expect_equal(big$value, r$value)
})

test_that("RMVAR differs from mean SE and uses the delta method", {
  se <- c(1, 2, 3, 4)
  expect_equal(adempR:::k_mean(se)$value, 2.5)
  r <- adempR:::k_rmvar(se)
  expect_equal(r$value, sqrt(30 / 4))
  expect_equal(r$mcse, (sd(se^2) / 2) / (2 * sqrt(7.5)), tolerance = 1e-12)
})

test_that("proportions keep nonzero uncertainty at the boundaries", {
  p <- adempR:::k_prop(c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(p$value, 0.75)
  expect_equal(p$mcse, sqrt(0.75 * 0.25 / 4))
  z <- adempR:::k_prop(rep(FALSE, 10))
  expect_equal(z$value, 0)
  expect_equal(z$mcse, 0)
  expect_equal(z$upper, 1 - 0.025^(1 / 10), tolerance = 1e-10)
  expect_identical(z$status, "boundary")
  a <- adempR:::k_prop(rep(TRUE, 10))
  expect_equal(a$lower, 0.025^(1 / 10), tolerance = 1e-10)
})

test_that("empty, single and degenerate inputs return honest statuses", {
  expect_identical(adempR:::k_mean(numeric())$status, "no_eligible")
  one <- adempR:::k_mean(3)
  expect_identical(one$status, "insufficient")
  expect_equal(one$value, 3)
  expect_true(is.na(one$mcse))
  expect_identical(adempR:::k_prop(logical())$status, "no_eligible")
  expect_identical(adempR:::k_rmse(c(0, 0, 0))$status, "degenerate")
  expect_identical(adempR:::k_empse(c(2, 2, 2, 2, 2))$status, "degenerate")
  inf <- adempR:::k_mean(c(1, Inf, 2))
  expect_identical(inf$status, "infinite")
  expect_equal(inf$value, Inf)
})

test_that("custom measures never get a generic MCSE for nonlinear aggregates", {
  m <- sim_measure("med_abs", contribution = function(rows) abs(rows$estimate - rows$truth),
                   aggregate = function(x) list(value = median(x)))
  r <- adempR:::k_custom(m, data.frame(estimate = c(1, 2, 3), truth = 0))
  expect_equal(r$value, 2)
  expect_true(is.na(r$mcse))
  expect_identical(r$status, "mcse_unavailable")
  lin <- sim_measure("abs_err", contribution = function(rows) abs(rows$estimate - rows$truth))
  r2 <- adempR:::k_custom(lin, data.frame(estimate = c(1, 2, 3), truth = 0))
  expect_equal(r2$mcse, sd(c(1, 2, 3)) / sqrt(3))
  expect_error(sim_measures(sim_measure("bias", function(rows) rows$estimate)), "clashes")
})

test_that("mergeable accumulators match batch calculations", {
  set.seed(5)
  x <- rnorm(3001, 1e6, 3)
  y <- x + rnorm(3001)
  parts <- split(seq_along(x), rep(1:7, length.out = length(x)))
  acc <- Reduce(sim_merge, lapply(parts, function(i) {
    sim_accumulate(sim_accumulator(), x[i], y[i])
  }))
  m <- sim_accumulated(acc)
  expect_equal(m$n, 3001)
  expect_equal(m$mean_x, mean(x), tolerance = 1e-14)
  expect_equal(m$var_x, var(x), tolerance = 1e-10)
  expect_equal(m$cov_xy, cov(x, y), tolerance = 1e-10)
  expect_equal(m$mcse_mean_diff, sd(x - y) / sqrt(3001), tolerance = 1e-8)
})

test_that("replication planning reproduces the textbook numbers", {
  expect_equal(as.integer(sim_plan_replications("proportion", p = 0.95, mcse = 0.005)), 1900L)
  expect_equal(as.integer(sim_plan_replications("proportion", mcse = 0.005)), 10000L)
  expect_equal(as.integer(sim_plan_replications("mean", sd = 0.2, mcse = 0.01)), 400L)
  hw <- sim_plan_replications("proportion", p = 0.95, half_width = 0.01)
  expect_equal(as.integer(hw), ceiling(qnorm(0.975)^2 * 0.95 * 0.05 / 0.01^2))
  # a pilot that saw no failures is not taken at face value
  expect_gt(adempR:::conservative_pq(0, 20), 0)
})

test_that("rejection uses p <= alpha, validated metadata and the eligible rows", {
  expect_error(sim_estimate("mean", 1, alpha = -1), "alpha")
  expect_error(sim_estimate("mean", 1, level = 2), "level")
  m <- list(
    edge = sim_method(function(data) sim_estimate("mean", 0, p_value = 0.05, null = 0,
                                                  lower = -1, upper = 1, level = 0.9),
                      estimands = "mean"),
    mixed = sim_method(function(data, context) {
      sim_estimate("mean", 0, lower = -1, upper = 1,
                   level = if (context$replication == 1) 0.5 else 0.95)
    }, estimands = "mean"))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = c("coverage", "rejection"))
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 3, seed = 1, path = run_path())))
  rj <- s[s$method == "edge" & s$measure == "rejection", ]
  expect_equal(rj$value, 1)
  expect_identical(rj$label, "type_i_error")
  expect_match(s$note[s$method == "edge" & s$measure == "coverage"], "nominal level 0.9")
  expect_match(s$note[s$method == "mixed" & s$measure == "coverage"], "several nominal levels")
  # an approximate truth close to the null gets the neutral label
  rows <- data.frame(truth = 0.01, truth_se = 0.2, null = 0, alternative = "two.sided")
  expect_identical(adempR:::rejection_label(rows, list(null = 0)), "rejection")
})

test_that("truth sensitivity ranges include interior extrema", {
  rows <- data.frame(estimate = c(0, 0), truth = c(0, 0), lower = c(-0.5, -0.5),
                     upper = c(0.5, 0.5), valid = TRUE, estimate_status = "valid",
                     interval_status = "valid")
  expect_equal(adempR:::truth_sensitivity("coverage", rows, list(level = 0.95), NULL, 2), c(0, 1))
  st <- list(level = 0.95)
  expect_equal(adempR:::truth_sensitivity("mse", rows, st, NULL, 2), c(0, 4))
})

test_that("bias eliminated coverage handles varying truth, infinity and dependence", {
  rows <- data.frame(estimate = c(0, 10, 0, 10), lower = c(-0.1, 9.9, -0.1, 9.9),
                     upper = c(0.1, 10.1, 0.1, 10.1))
  r <- adempR:::k_be_coverage(rows, rep(TRUE, 4), rep(TRUE, 4), c(0, 10, 0, 10))
  expect_equal(r$value, 1)
  rows$estimate <- c(Inf, -Inf, 0, 0)
  r <- adempR:::k_be_coverage(rows, rep(TRUE, 4), rep(TRUE, 4), c(0, 10, 0, 10))
  expect_identical(r$status, "infinite")
  set.seed(2)
  x <- rnorm(200)
  rows <- data.frame(estimate = x, lower = x - 1.5, upper = x + 1.5)
  r <- adempR:::k_be_coverage(rows, rep(TRUE, 200), rep(TRUE, 200), rep(0, 200))
  expect_true(r$mcse > 0 && r$lower < r$value && r$value < r$upper)
})

test_that("weighted summaries keep boundary uncertainty and ignore zero weights", {
  m <- list(wide = function(data) sim_estimate("mean", 0, lower = -100, upper = 100, level = 0.95))
  st <- normal_study(n = c(10L, 20L, 30L), mu = 0, measures = c("coverage", "rmse"), methods = m)
  st$scenarios$weights <- c(0.5, 0.5, 0)
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 2, seed = 1, path = run_path()),
                                   weighted = TRUE))
  wc <- s[s$scenario_key == "weighted" & s$measure == "coverage", ]
  expect_equal(wc$value, 1)
  expect_identical(wc$status, "boundary")
  expect_lt(wc$lower, 1)
  wr <- s[s$scenario_key == "weighted" & s$measure == "rmse", ]
  expect_equal(wr$value, 0)
})

test_that("custom measures refuse missing contributions and invalid aggregates", {
  rows <- data.frame(estimate = c(1, 2, NA))
  m <- sim_measure("m", contribution = function(rows) rows$estimate)
  expect_identical(adempR:::k_custom(m, rows)$status, "missing_contributions")
  num01 <- sim_measure("n", contribution = function(rows) c(0, 1, 1))
  r <- adempR:::k_custom(num01, rows)
  expect_true(is.na(r$lower) || r$status != "boundary")
  expect_equal(r$mcse, sd(c(0, 1, 1)) / sqrt(3))
  bad <- sim_measure("b", contribution = function(rows) c(1, 2, 3),
                     aggregate = function(x) list(value = "x", mcse = 1))
  expect_error(adempR:::k_custom(bad, rows), "single numeric")
})

test_that("paired accumulators keep small differences and large counts", {
  x <- 1e8 * c(1, -1, 1, -1)
  y <- x + c(0, 1, 0, -1)
  a <- sim_accumulate(sim_accumulator(), x[1:2], y[1:2])
  b <- sim_accumulate(sim_accumulator(), x[3:4], y[3:4])
  m <- sim_accumulated(sim_merge(a, b))
  expect_equal(m$mcse_mean_diff, sd(x - y) / 2)
  acc <- sim_accumulate(sim_accumulator(), 1)
  for (i in 1:40) acc <- sim_merge(acc, acc)
  expect_equal(sim_accumulated(acc)$n, 2^40)
  expect_error(sim_accumulate(sim_accumulator(), c(1, Inf)), "finite")
})

test_that("root mean squares, relative bias and planning handle extreme values", {
  expect_equal(adempR:::k_rmse(c(1e200, 1e200))$value, 1e200)
  expect_equal(adempR:::k_rmse(c(1e-200, 2e-200))$value, sqrt(mean(c(1, 4))) * 1e-200)
  expect_equal(adempR:::k_rmvar(c(1e155, 2e155))$value, sqrt(2.5) * 1e155)
  expect_error(sim_plan_replications("proportion", mcse = 1e-6), "relax")
  # truth near zero relative to the spread is flagged
  st <- normal_study(n = 10L, mu = 1e-6, measures = c("relbias"))
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 5, seed = 1, path = run_path())))
  expect_identical(s$status[s$measure == "relbias"], "unstable")
})

test_that("the robust z score uses mad() once", {
  # replication 6 has robust z of about 6.1 with mad(), 4.1 with a double scaling
  vals <- c(-1, -0.5, 0, 0.5, 1, 7)
  m <- list(a = function(data, context) sim_estimate("mean", vals[context$replication]))
  st <- normal_study(n = 5L, mu = 0, methods = m)
  run <- sim_run(st, replications = 6, seed = 1, path = run_path())
  d <- sim_diagnose(run, outlier_z = 5)
  expect_true(any(d$outliers$quantity == "estimate"))
})

test_that("diagnostics flag fixed data, test and interval disagreement, degenerate estimates", {
  fixed <- sim_study("fixed", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) data.frame(y = seq_len(scenario$n)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 5.5)),
                     sim_estimand("median", sim_truth_exact(function(scenario) 5.5))),
    methods = list(m = function(data) {
      list(sim_estimate("mean", mean(data$y)), sim_estimate("median", mean(data$y)))
    }))
  d <- sim_diagnose(sim_run(fixed, replications = 3, seed = 1, path = run_path()))
  expect_equal(d$repeated_data$distinct, 1L)
  expect_true(any(grepl("all 3 replications", d$degenerate$problem)))
  expect_true(any(d$degenerate$estimand == "mean = median"))
  # a method whose p value comes from a different statistic than its interval
  odd <- normal_study(n = 10L, mu = 0, methods = list(odd = function(data) {
    m <- mean(data$y); se <- sd(data$y) / sqrt(10)
    sim_estimate("mean", m, se = se, lower = m - 1.96 * se, upper = m + 1.96 * se,
                 level = 0.95, p_value = 0.01, null = 0)
  }))
  d2 <- sim_diagnose(sim_run(odd, replications = 6, seed = 1, path = run_path()))
  expect_gt(d2$test_interval_disagreement$n_disagree, 0)
  expect_null(sim_diagnose(sim_run(normal_study(n = 10L, mu = 0), replications = 6, seed = 1,
                                   path = run_path()))$repeated_data)
})

test_that("weighted averages can be formed within levels of a factor and show their range", {
  st <- normal_study(n = c(10L, 20L), mu = c(0, 1), measures = "bias")
  st$scenarios$weights <- rep(0.25, 4)
  run <- sim_run(st, replications = 20, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run, weighted = TRUE, by = "mu"))
  w <- s[grepl("^weighted", s$scenario_key) & s$measure == "bias", ]
  expect_equal(nrow(w), 2L)
  expect_true(all(grepl("2 scenarios, values from", w$note)))
  expect_error(sim_summarise(run, weighted = TRUE, by = "nope"), "Unknown")
})

test_that("weighted groups are typed tuples, dependent scenarios get a valid bound, counts are totals", {
  # factor values whose concatenations collide stay separate groups
  sc <- sim_scenarios(a = c("x, y", "x"), b = c("z", "y, z"), .weights = rep(1, 4))
  st <- sim_study("grp", scenarios = sc,
    generate = function(scenario, context) data.frame(y = rnorm(5) + nchar(scenario$a)),
    estimands = list(sim_estimand("m", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("m", mean(data$y))), measures = "bias")
  run <- sim_run(st, replications = 4, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run, weighted = TRUE, by = c("a", "b")))
  # four distinct (a, b) tuples; ("x, y", "z") and ("x", "y, z") read the
  # same when pasted together but are different groups
  expect_equal(length(unique(s$scenario_key[startsWith(s$scenario_key, "weighted: ")])), 4L)
  # common random numbers: two scenarios with identical errors
  gen <- sim_generator(function(scenario, context) data.frame(y = scenario$mu + rnorm(20)),
                       common = "n")
  sc2 <- sim_scenarios(n = 20L, mu = c(0, 1), .weights = c(1, 1))
  st2 <- sim_study("crn", scenarios = sc2,
    generate = gen,
    estimands = list(sim_estimand("m", sim_truth_exact(function(scenario) scenario$mu))),
    methods = list(m = function(data) sim_estimate("m", mean(data$y))), measures = "bias")
  s2 <- as.data.frame(sim_summarise(sim_run(st2, replications = 30, seed = 2, path = run_path()),
                                    weighted = TRUE))
  per <- s2$mcse[s2$scenario_key != "weighted" & s2$measure == "bias"]
  expect_equal(per[1], per[2], tolerance = 1e-12)
  wrow <- s2[s2$scenario_key == "weighted" & s2$measure == "bias", ]
  # the weighted mean of two identical error vectors has the per scenario MCSE
  expect_equal(wrow$mcse, per[1], tolerance = 1e-12)
  expect_match(wrow$note, "share random numbers")
  expect_equal(wrow$n_planned, 60L)
  expect_equal(wrow$n_eligible, 60L)
})

test_that("weighted boundary intervals cover usable_coverage and custom logical proportions", {
  hit <- sim_measure("hit", function(rows) rep(TRUE, nrow(rows)), requires = "estimate")
  st <- normal_study(n = c(10L, 20L), mu = 0,
                     measures = sim_measures("usable_coverage", hit))
  st$scenarios$weights <- c(0.5, 0.5)
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 5, seed = 1, path = run_path()),
                                   weighted = TRUE))
  w <- s[s$scenario_key == "weighted" & s$measure == "hit", ]
  expect_equal(w$status, "boundary")
  expect_true(is.finite(w$lower) && w$lower < 1)
})

test_that("truth sensitivity finds the exact extrema of squared error and coverage", {
  # Monte Carlo truth 0 with standard error 0.577; a constant estimate of 1/3
  mc <- sim_truth_mc(function(scenario, n) c(-1, 1, -1, 1), n = 4)
  st <- sim_study("sens", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(odd = context$replication %% 2 == 1),
    estimands = list(sim_estimand("m", mc)),
    methods = list(c = function(data) {
      # intervals [-1, 0] and [1e-10, 1]: no interval covers a truth inside the gap
      if (data$odd[1]) sim_estimate("m", 1 / 3, lower = -1, upper = 0, level = 0.95)
      else sim_estimate("m", 1 / 3, lower = 1e-10, upper = 1, level = 0.95)
    }), measures = c("mse", "coverage"))
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 4, seed = 1, path = run_path())))
  mse <- s[s$measure == "mse", ]
  expect_equal(mse$truth_sens_lower, 0)
  cov <- s[s$measure == "coverage", ]
  expect_equal(cov$value, 0.5)
  expect_equal(cov$truth_sens_lower, 0)
  expect_equal(cov$truth_sens_upper, 0.5)
})

test_that("a numerical truth error bound gives a sensitivity range", {
  st <- sim_study("numsens", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("m", sim_truth_numeric(function(scenario) list(value = 0, error = 0.25)))),
    methods = list(c = function(data) sim_estimate("m", mean(data$y))), measures = "bias")
  run <- sim_run(st, replications = 10, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  b <- s[s$measure == "bias", ]
  expect_equal(c(b$truth_sens_lower, b$truth_sens_upper), b$value + c(-0.25, 0.25))
  expect_true(is.na(b$mcse_total))
})

test_that("bias eliminated coverage deletes whole replications and handles a fixed unknown truth", {
  # rows 1 to 4 have intervals, rows 5 and 6 only estimates that move the center
  rows <- data.frame(estimate = c(0, 0.2, -0.2, 0.1, 6, -6), lower = c(-1, -1, -1, -1, NA, NA),
                     upper = c(1, 1, 1, 1, NA, NA))
  sel <- c(rep(TRUE, 4), FALSE, FALSE)
  r <- adempR:::k_be_coverage(rows, sel, rep(TRUE, 6), rep(0, 6))
  # deleting replication 5 or 6 moves the center by 1, out of every interval
  th <- c(rep(1, 4), 0, 0)
  expect_equal(r$mcse, sqrt(5 / 6 * sum((th - mean(th))^2)))
  # with a fixed unknown truth the target is the mean estimate
  st <- sim_study("unk", scenarios = sim_scenarios(n = 20L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, 3)),
    estimands = list(sim_estimand("m", sim_truth_unknown())),
    methods = list(t = function(data) {
      ci <- t.test(data$y)$conf.int
      sim_estimate("m", mean(data$y), lower = ci[1], upper = ci[2], level = 0.95)
    }), measures = "be_coverage")
  run <- sim_run(st, replications = 30, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  b <- s[s$measure == "be_coverage", ]
  e <- sim_results(run)$estimates
  expect_equal(b$value, mean(e$lower <= mean(e$estimate) & mean(e$estimate) <= e$upper))
  expect_match(b$note, "mean estimate")
})

test_that("an undefined leave one out SE ratio is not reported as ok", {
  r <- adempR:::k_se_ratio(c(0, 0, 0, 1), rep(1, 4))
  expect_identical(r$status, "mcse_unavailable")
  expect_true(is.na(r$mcse))
  expect_equal(r$value, 100 * (1 / sd(c(0, 0, 0, 1)) - 1))
  expect_match(r$note, "jackknife MCSE is undefined")
})
