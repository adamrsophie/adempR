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
  expect_equal(adempR:::truth_sensitivity("coverage", rows, list(level = 0.95), NULL, 1), c(0, 1))
  st <- list(level = 0.95)
  expect_equal(adempR:::truth_sensitivity("mse", rows, st, NULL, 1), c(0, 4))
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
