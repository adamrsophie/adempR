test_that("censoring calibration reaches the target proportion", {
  sampler <- function(n) sim_rweibull_ph(rnorm(n), lambda = 0.1, shape = 1.5)$event_time
  cal <- sim_calibrate_censoring(0.4, sampler, horizon = 8, n = 2e4, seed = 1)
  expect_equal(cal$achieved, 0.4, tolerance = 1e-6)
  set.seed(5)
  d <- sim_rweibull_ph(rnorm(4e4), lambda = 0.1, shape = 1.5, censoring_rate = cal$rate, horizon = 8)
  expect_equal(1 - mean(d$status), 0.4, tolerance = 0.015)
  expect_error(sim_calibrate_censoring(0.01, sampler, horizon = 0.1, n = 1e3), "above the target")
})

test_that("Gaussian copula draws have the requested margins and rank correlation", {
  set.seed(2)
  d <- sim_rcopula(2e4, list(x = function(p) qexp(p, 2), g = function(p) qgamma(p, 3)),
                   matrix(c(1, 0.6, 0.6, 1), 2))
  expect_named(d, c("x", "g"))
  expect_equal(mean(d$x), 0.5, tolerance = 0.03)
  expect_equal(cor(d$x, d$g, method = "spearman"), 6 / pi * asin(0.3), tolerance = 0.02)
  expect_error(sim_rcopula(10, list(x = qnorm), diag(2)), "1 by 1")
})

test_that("plasmode draws resample real rows and add simulated columns", {
  set.seed(3)
  pm <- sim_plasmode(mtcars, n = 40, outcome = function(d) d$hp * 0 + 1)
  expect_equal(nrow(pm), 40L)
  expect_true(all(pm$y == 1))
  expect_equal(pm$mpg, mtcars$mpg[attr(pm, "rows")])
  sub <- sim_plasmode(mtcars, n = 10, replace = FALSE, outcome = function(d) list(a = d$am, b = 2))
  expect_false(anyDuplicated(attr(sub, "rows")) > 0)
  expect_true(all(c("a", "b") %in% names(sub)))
  expect_error(sim_plasmode(mtcars, n = 40, replace = FALSE), "cannot exceed")
})
