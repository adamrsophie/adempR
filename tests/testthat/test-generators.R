test_that("matrix checks never repair", {
  expect_silent(sim_check_matrix(diag(3)))
  bad <- matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3)
  expect_error(sim_check_matrix(bad), "not positive semidefinite")
  expect_error(sim_check_matrix(matrix(c(1, 0.5, 0.4, 1), 2)), "symmetric")
  expect_error(sim_check_matrix(matrix(c(2, 0, 0, 1), 2)), "unit diagonal")
  expect_error(sim_rmvnorm(10, c(0, 0), matrix(c(1, 2, 2, 1), 2)), "positive semidefinite")
})

test_that("multivariate normal draws reach their target correlation", {
  set.seed(9)
  x <- sim_rmvnorm(20000, c(1, -1), matrix(c(1, 0.6, 0.6, 2), 2))
  expect_equal(cor(x)[1, 2], 0.6 / sqrt(2), tolerance = 0.03)
  expect_equal(colMeans(x), c(1, -1), tolerance = 0.05)
})

test_that("intercept calibration achieves the marginal prevalence", {
  cal <- sim_calibrate_intercept(0.2, function(n) 1.5 * rnorm(n), n = 2e5, seed = 4)
  expect_equal(cal$achieved, 0.2, tolerance = 1e-8)
  # independent check on fresh draws
  set.seed(10)
  fresh <- mean(plogis(cal$intercept + 1.5 * rnorm(2e5)))
  expect_lt(abs(fresh - 0.2), 4 * sqrt(0.2 * 0.8 / 2e5) + 0.002)
  expect_error(sim_calibrate_intercept(1.2, function(n) rnorm(n)), "probabilities")
})

test_that("survival and missingness generators hit their targets", {
  set.seed(11)
  d <- sim_rweibull_ph(rep(0, 20000), lambda = 0.2, shape = 1, horizon = 3)
  # exponential with rate 0.2 censored at 3: P(event) = 1 - exp(-0.6)
  expect_equal(mean(d$status), 1 - exp(-0.6), tolerance = 0.02)
  expect_true(all(d$time <= 3))
  dd <- data.frame(x = rnorm(5000), y = rnorm(5000))
  m <- sim_add_missing(dd, "y", 0.3, "MAR", predictors = "x", weights = 1.5)
  expect_equal(mean(is.na(m$y)), 0.3, tolerance = 0.03)
  expect_gt(mean(m$x[is.na(m$y)]), mean(m$x[!is.na(m$y)]))
  expect_error(sim_add_missing(dd, "y", 0.3, "MAR", predictors = "y", weights = 1), "MNAR")
  cl <- sim_rclustered(400, 10, icc = 0.3)
  expect_equal(var(tapply(cl$u, cl$cluster, mean)), 0.3, tolerance = 0.06)
  expect_error(sim_rclustered(10, 5, icc = 1), "icc")
})

test_that("DGP checks distinguish fluctuation from a wrong mechanism", {
  st <- sim_study("chk", scenarios = sim_scenarios(n = 200L, p = c(0.1, 0.5)),
    generate = function(scenario, context) rbinom(scenario$n, 1, scenario$p),
    estimands = list(sim_estimand("p")),
    methods = list(a = function(data) sim_estimate("p", mean(data))))
  ok <- sim_check_dgp(st, function(d) c(prev = mean(d)), function(scenario) c(prev = scenario$p),
                      n = 200, seed = 1)
  expect_false(any(ok$flag))
  wrong <- sim_check_dgp(st, function(d) c(prev = mean(d)),
                         function(scenario) c(prev = scenario$p + 0.05), n = 200, seed = 1)
  expect_true(all(wrong$flag))
})

test_that("an indefinite block is found next to a large variance", {
  m <- diag(c(1e12, 1, 1))
  m[2, 3] <- m[3, 2] <- 10
  expect_error(sim_check_matrix(m, "covariance"), "not positive semidefinite")
  expect_error(sim_rmvnorm(5, c(0, 0, 0), m), "not positive semidefinite")
})

test_that("symmetry is checked on the correlation scale", {
  m <- diag(c(1e14, 1, 1))
  m[2, 3] <- 2
  m[3, 2] <- 0
  expect_error(sim_check_matrix(m, "covariance"), "not symmetric")
  # a large but symmetric covariance passes
  ok <- diag(c(1e14, 1, 1))
  ok[2, 3] <- ok[3, 2] <- 0.5
  expect_silent(sim_check_matrix(ok, "covariance"))
})
