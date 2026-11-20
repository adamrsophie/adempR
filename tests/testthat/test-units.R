anti_study <- function(method, record = function(data) list(control = data$x)) {
  sim_study("units", scenarios = sim_scenarios(n = 1L),
    generate = sim_generator(function(scenario, context) {
      z <- stats::rnorm(2)
      if (context$antithetic) z <- -z
      data.frame(y = z[1], x = z[2])
    }, antithetic = TRUE, record = record),
    estimands = list(sim_estimand("mu", sim_truth_exact(function(scenario) 0))),
    methods = list(m = method), measures = c("bias", "mse", "coverage"))
}

# hand computed delete-one-pair jackknife of a mean over rows
jk_mean <- function(x, rep) {
  u <- (rep + 1L) %/% 2L
  th <- vapply(unique(u), function(k) mean(x[u != k]), numeric(1))
  G <- length(th)
  sqrt((G - 1) / G * sum((th - mean(th))^2))
}

test_that("every summary path takes the MCSE from antithetic pairs, also across chunks and with gaps", {
  # nonzero varying pair means: the estimate is y plus a replication shift
  # that differs within pairs; replication 5 fails, so pair 3 has one
  # member, and the run ends on an odd replication
  method <- function(data, context) {
    if (context$replication == 5L) stop("fit failed")
    x <- data$y + context$replication / 10
    sim_estimate("mu", x, lower = x - 1, upper = x + 1, level = 0.95)
  }
  st <- anti_study(method)
  run <- sim_run(st, 9, 3, run_path(), chunk_size = 3, progress = FALSE)
  e <- sim_results(run)$estimates
  e <- e[e$valid, ]
  err <- e$estimate - e$truth
  cov <- as.numeric(e$lower <= 0 & 0 <= e$upper)
  s <- as.data.frame(sim_summarise(run))
  get <- function(d, m, col = "mcse") d[d$measure == m, col]
  expect_equal(get(s, "bias"), jk_mean(err, e$replication), tolerance = 1e-12)
  expect_equal(get(s, "mse"), jk_mean(err^2, e$replication), tolerance = 1e-12)
  expect_equal(get(s, "coverage"), jk_mean(cov, e$replication), tolerance = 1e-12)
  expect_equal(get(s, "bias", "n_units"), 5L)
  expect_equal(get(s, "bias", "n_eligible"), 8L)
  ss <- sim_summarise_stream(run, c("bias", "mse", "coverage"))
  for (m in c("bias", "mse", "coverage")) {
    expect_equal(get(ss, m), get(s, m), tolerance = 1e-10)
    expect_equal(get(ss, m, "value"), get(s, m, "value"), tolerance = 1e-12)
  }
  d <- sim_bias_decomposition(run, function(scenario) 0)
  expect_equal(d$total_mcse, get(s, "bias"), tolerance = 1e-12)
  expect_equal(d$n_units, 5L)
})

test_that("the closed form agrees with the generic pair jackknife", {
  method <- function(data, context) sim_estimate("mu", data$y + context$replication / 7)
  st <- anti_study(method)
  run <- sim_run(st, 7, 4, run_path(), progress = FALSE)
  res <- sim_results(run)
  rows <- res$estimates
  r <- adempR:::measure_cell("mse", rows, st, st$estimands$mu)
  fast <- adempR:::pair_jackknife("mse", rows, st, st$estimands$mu, r)
  # force the generic path through a custom measure with the same value
  m <- sim_measure("mse2", function(rows) (rows$estimate - rows$truth)^2)
  r2 <- adempR:::measure_cell(m, rows, st, st$estimands$mu)
  slow <- adempR:::pair_jackknife(m, rows, st, st$estimands$mu, r2)
  expect_equal(fast$mcse, slow$mcse, tolerance = 1e-12)
})

test_that("antithetic proportions with equal pair means are not reported as certain", {
  for (shift in c(0, 10)) {
    method <- function(data, context) {
      x <- c(-0.5, 0.5, -0.25, 0.25, -0.1, 0.1)[context$replication] + shift
      sim_estimate("mu", x, lower = x - 1, upper = x + 1, level = 0.95)
    }
    run <- sim_run(anti_study(method), 6, 1, run_path(), progress = FALSE)
    cov <- as.data.frame(sim_summarise(run))
    cov <- cov[cov$measure == "coverage", ]
    expect_equal(cov$status, "boundary")
    # exact limits over 3 independent pairs, not over 6 replications
    if (shift == 0) {
      expect_equal(cov$lower, unname(adempR:::cp_interval(3, 3, 0.95)[1, "lower"]))
      expect_equal(cov$upper, 1)
    } else {
      expect_equal(cov$lower, 0)
      expect_equal(cov$upper, unname(adempR:::cp_interval(0, 3, 0.95)[1, "upper"]))
    }
  }
  # pair means of 0.5 (one member covers): the bound is the whole range
  method <- function(data, context) {
    x <- if (context$replication %% 2L) 0 else 10
    sim_estimate("mu", x, lower = x - 1, upper = x + 1, level = 0.95)
  }
  cov <- as.data.frame(sim_summarise(sim_run(anti_study(method), 6, 1, run_path(), progress = FALSE)))
  cov <- cov[cov$measure == "coverage", ]
  expect_equal(cov$value, 0.5)
  expect_equal(c(cov$lower, cov$upper), c(0, 1))
})

test_that("control variates use pair means for the coefficient and the MCSEs", {
  # a control whose pair means vary (x^2 is the same in both members)
  method <- function(data, context) sim_estimate("mu", data$y + data$x^2 / 2 + context$replication / 10)
  st <- anti_study(method, record = function(data) list(control = data$x^2))
  run <- sim_run(st, 10, 2, run_path(), progress = FALSE)
  res <- sim_results(run)
  cv <- sim_control_variate(run, "control", control_mean = 1, measure = "bias")
  e <- res$estimates
  rec <- res$records
  cc <- rec$value[match(e$replication, rec$replication)]
  x <- e$estimate - e$truth
  u <- (e$replication + 1L) %/% 2L
  xb <- tapply(x, u, mean)
  cb <- tapply(cc, u, mean)
  b <- stats::cov(xb, cb) / stats::var(cb)
  expect_equal(cv$coefficient, b, tolerance = 1e-12)
  expect_equal(cv$mcse, stats::sd(xb) / sqrt(5), tolerance = 1e-12)
  expect_equal(cv$mcse_adjusted, stats::sd(xb - b * cb) / sqrt(5), tolerance = 1e-12)
  expect_equal(cv$value_adjusted, mean(x - b * (cc - 1)), tolerance = 1e-12)
  expect_equal(cv$n_units, 5L)
})

test_that("grid averages keep eligible infinite estimates", {
  mk <- function(vals) sim_study("grid_inf", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) 0,
    estimands = sim_estimand_grid("g", c(1, 2), function(scenario, at) at),
    methods = list(m = function(data) sim_estimate_grid("g", c(1, 2), vals)),
    measures = "bias")
  avg <- function(vals) {
    out <- sim_summarise_grid(sim_run(mk(vals), 3, 1, run_path(), progress = FALSE), "g", "bias")
    out[out$measure == "grid_mean_bias", ]
  }
  a <- avg(c(-Inf, 2))
  expect_identical(a$value, -Inf)
  expect_equal(a$status, "infinite")
  expect_equal(a$n_eligible, 3L)
  b <- avg(c(Inf, -Inf))
  expect_true(is.na(b$value))
  expect_equal(b$status, "undefined")
  expect_equal(b$n_eligible, 3L)
  f <- avg(c(1, 2))
  expect_equal(f$value, 0)
  expect_equal(f$status, "ok")
})

test_that("reference checks of antithetic proportions use the pairs", {
  f <- function(data, context) {
    x <- c(-0.5, 0.5, -0.25, 0.25, -0.1, 0.1)[context$replication]
    sim_estimate("mu", x, lower = x - 1, upper = x + 1, level = 0.95)
  }
  st <- anti_study(sim_method(f, reference = c(coverage = 0.4)))
  chk <- sim_check_reference(sim_run(st, 6, 1, run_path(), progress = FALSE))
  chk <- chk[chk$measure == "coverage", ]
  # three independent pairs, all covering: 0.4 is inside the bounded pair
  # interval, although a binomial interval over six replications excludes it
  expect_equal(chk$check, "bounded pair interval")
  expect_false(chk$flag)
  expect_gt(unname(adempR:::cp_interval(6, 6, 0.95)[1, "lower"]), 0.4)
})

test_that("pair jackknife matches an explicit delete-one-pair oracle", {
  set.seed(11)
  st <- list(level = 0.95, generator = list(antithetic = TRUE))
  est <- list(truth = list(kind = "exact"))
  oracle <- function(f, x, rep) {
    u <- (rep + 1L) %/% 2L
    th <- vapply(unique(u), function(k) f(x[u != k]), numeric(1))
    G <- length(th)
    sqrt((G - 1) / G * sum((th - mean(th))^2))
  }
  for (n in c(10L, 11L)) {
    z <- rnorm(n)
    rows <- data.frame(estimate = z, truth = 0, valid = TRUE, estimate_status = "valid",
                       interval_status = "missing", se_status = "missing",
                       test_status = "missing", se = NA_real_, lower = NA_real_,
                       upper = NA_real_, level = 0.95, replication = seq_len(n))
    # linear (closed form), complete and incomplete pairs
    r <- adempR:::measure_cell("mse", rows, st, est)
    expect_equal(adempR:::pair_jackknife("mse", rows, st, est, r)$mcse,
                 oracle(function(x) mean(x^2), z, rows$replication), tolerance = 1e-10)
    # nonlinear (generic path)
    r <- adempR:::measure_cell("rmse", rows, st, est)
    expect_equal(adempR:::pair_jackknife("rmse", rows, st, est, r)$mcse,
                 oracle(function(x) sqrt(mean(x^2)), z, rows$replication), tolerance = 1e-10)
  }
})
