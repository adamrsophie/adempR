# Scientific test 1: the t interval for a normal mean has known performance.

test_that("normal mean: bias, empirical SE, coverage and rejection agree with theory", {
  study <- normal_study(n = 20L, mu = c(0, 1))
  run <- sim_run(study, replications = 1000, seed = 2026, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  get <- function(mu, m) s[s$measure == m & grepl(paste0("mu=", mu, ","), s$scenario_id), ]
  for (mu in c(0, 1)) {
    b <- get(mu, "bias")
    expect_lt(abs(b$value), 4 * b$mcse)
    es <- get(mu, "empse")
    expect_lt(abs(es$value - 1 / sqrt(20)), 4 * es$mcse)
    # E[s] = c4 sigma, with c4 = sqrt(2 / (n - 1)) Gamma(n / 2) / Gamma((n - 1) / 2)
    c4 <- sqrt(2 / 19) * exp(lgamma(10) - lgamma(9.5))
    ms <- get(mu, "mean_se")
    expect_lt(abs(ms$value - c4 / sqrt(20)), 4 * ms$mcse)
    cv <- get(mu, "coverage")
    expect_lt(abs(cv$value - 0.95), 4 * sqrt(0.95 * 0.05 / 1000))
    expect_equal(cv$n_eligible, 1000L)
  }
  t1 <- get(0, "rejection")
  expect_identical(t1$label, "type_i_error")
  expect_lt(abs(t1$value - 0.05), 4 * sqrt(0.05 * 0.95 / 1000))
  pw <- get(1, "rejection")
  expect_identical(pw$label, "power")
  true_power <- power.t.test(n = 20, delta = 1, sd = 1, type = "one.sample")$power
  expect_lt(abs(pw$value - true_power), 4 * sqrt(true_power * (1 - true_power) / 1000) + 0.002)
  expect_true(all(s$value[s$measure == "usable"] == 1))
})

test_that("results agree with rsimsum where definitions align", {
  skip_if_not_installed("rsimsum")
  study <- normal_study(n = 15L, mu = 0.3)
  run <- sim_run(study, replications = 300, seed = 7, path = run_path())
  est <- sim_results(run)$estimates
  s <- as.data.frame(sim_summarise(run))
  ss <- rsimsum::simsum(data = est, estvarname = "estimate", se = "se", true = 0.3)
  rs <- rsimsum::tidy(summary(ss))
  pick <- function(stat) rs[rs$stat == stat, ]
  expect_equal(s$value[s$measure == "bias"], pick("bias")$est, tolerance = 1e-10)
  expect_equal(s$mcse[s$measure == "bias"], pick("bias")$mcse, tolerance = 1e-10)
  expect_equal(s$value[s$measure == "empse"], pick("empse")$est, tolerance = 1e-10)
  expect_equal(s$value[s$measure == "mse"], pick("mse")$est, tolerance = 1e-10)
  expect_equal(s$value[s$measure == "se_ratio"], pick("relerror")$est, tolerance = 1e-8)
  expect_equal(s$mcse[s$measure == "se_ratio"], pick("relerror")$mcse, tolerance = 1e-8)
})
