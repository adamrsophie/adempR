test_that("confidence sequences are wider than fixed intervals and cover the truth uniformly", {
  st <- normal_study(n = 15L, mu = 0, measures = c("coverage", "bias"))
  run <- sim_run(st, replications = 60, seed = 2, path = run_path())
  cs <- sim_confidence_sequence(run, "coverage")
  s <- as.data.frame(sim_summarise(run))
  ci <- s[s$measure == "coverage", ]
  expect_equal(cs$estimate, ci$value)
  expect_true(cs$lower <= ci$lower && cs$upper >= ci$upper)
  expect_equal(cs$type, "exact")
  p <- sim_confidence_sequence(run, "bias", path = TRUE)
  expect_equal(nrow(p), 60L)
  expect_equal(p$type[1], "asymptotic")
  # simulation check of time-uniform validity for Bernoulli(0.95) streams
  set.seed(1)
  bern_cs <- function(x) {
    t <- seq_along(x); m <- cumsum(x) / t
    h <- adempR:::robbins_halfwidth(t, 0.05, 200)
    any(0.95 < m - h | 0.95 > m + h)
  }
  miss <- mean(replicate(400, bern_cs(rbinom(500, 1, 0.95))))
  expect_lt(miss, 0.05)
  expect_error(sim_confidence_sequence(run, alpha = 2), "alpha")
})

test_that("confidence sequence tuning is fixed at the first look and survives extensions", {
  st <- normal_study(n = 15L, mu = 0, measures = c("coverage", "bias"))
  p <- run_path()
  run <- sim_run(st, replications = 20, seed = 2, path = p)
  first <- sim_confidence_sequence(run, "bias")
  expect_equal(first$t_opt, 20)
  run <- sim_run(st, replications = 80, seed = 2, path = p)
  later <- sim_confidence_sequence(run, "bias", path = TRUE)
  expect_equal(unique(later$t_opt), 20)
  # the first 20 looks are unchanged by the extension
  expect_equal(later$upper[20], first$upper, tolerance = 1e-12)
  expect_error(sim_confidence_sequence(run, "bias", t_opt = 80), "tuned with t_opt = 20")
  expect_error(sim_confidence_sequence(run, "bias", alpha = 0.1), "same alpha")
  # the running SD matches the direct computation
  x <- c(3, 1, 4, 1, 5, 9, 2, 6)
  expect_equal(adempR:::running_sd(x)[-1], vapply(2:8, function(k) sd(x[1:k]), numeric(1)))
})

test_that("rejection sequences use the shared decision and eligibility; zero SD gives no interval", {
  flag_m <- function(data, context) {
    r <- context$replication
    if (r == 1) sim_estimate("mean", 0, reject = TRUE)                     # decision only
    else if (r == 2) sim_estimate("mean", 0, reject = TRUE, p_value = 0.9) # explicit wins
    else sim_estimate("mean", 0, p_value = 0.01, converged = FALSE)       # excluded fit
  }
  st <- normal_study(n = 5L, mu = 0, methods = list(m = flag_m), measures = "rejection")
  run <- sim_run(st, replications = 3, seed = 1, path = run_path())
  cs <- sim_confidence_sequence(run, "rejection", path = TRUE)
  s <- as.data.frame(sim_summarise(run))
  expect_equal(nrow(cs), 2L)
  expect_equal(cs$estimate[2], s$value[s$measure %in% c("rejection", "type_i_error", "power")])
  # a constant error: running SD 0, so no asymptotic interval
  const <- normal_study(n = 5L, mu = 0, methods = list(c = function(data) sim_estimate("mean", 1)),
                        measures = "bias")
  b <- sim_confidence_sequence(sim_run(const, replications = 4, seed = 1, path = run_path()), "bias")
  expect_true(is.na(b$lower) && is.na(b$upper))
})

test_that("the tuning on disk decides, whatever handle is reused", {
  st <- normal_study(n = 5L, mu = 0, measures = "coverage")
  p <- run_path()
  run <- sim_run(st, 20, 1, p, progress = FALSE)
  res <- sim_results(run)
  a <- sim_confidence_sequence(run, "coverage", t_opt = 20)
  expect_true(all(a$tuning_recorded))
  # same tuning is idempotent, from the run and from older results
  expect_equal(sim_confidence_sequence(res, "coverage", t_opt = 20)$upper, a$upper)
  expect_error(sim_confidence_sequence(res, "coverage", t_opt = 80), "tuned")
  expect_error(sim_confidence_sequence(run, "coverage", alpha = 0.1), "alpha")
  # extension and reopening keep the sequence
  run2 <- sim_run(st, 30, 1, p, progress = FALSE)
  expect_equal(unique(sim_confidence_sequence(run2, "coverage")$t_opt), 20)
  expect_error(sim_confidence_sequence(sim_open(p), "coverage", t_opt = 30), "tuned")
})

test_that("monitoring against an approximate truth is refused or widened", {
  mk <- function(truth) sim_study("cs_truth", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) 0,
    estimands = list(sim_estimand("mu", truth)),
    methods = list(m = function(data, context) {
      x <- 0.5 + (context$replication %% 2) / 100
      sim_estimate("mu", x, lower = x - 0.1, upper = x + 0.1)
    }), measures = c("bias", "coverage"))
  approx <- mk(sim_truth_numeric(function(scenario) 0.5, error = 0.5))
  unstated <- mk(sim_truth_numeric(function(scenario) 0.5))
  exact0 <- mk(sim_truth_numeric(function(scenario) 0.5, error = 0))
  for (st in list(approx, unstated)) {
    expect_error(sim_confidence_sequence(sim_run(st, 10, 1, run_path(), progress = FALSE), "coverage"),
                 "exact truth")
    expect_error(sim_run_sequential(st, looks = 10, half_width = 0.05, seed = 1, path = run_path(),
                                    measures = "coverage"), "exact truth")
  }
  # bias against a stated error bound: widened by the bound
  b <- sim_confidence_sequence(sim_run(approx, 10, 1, run_path(), progress = FALSE), "bias")
  expect_equal(b$truth_error, 0.5)
  expect_lte(b$lower, b$estimate - 0.5)
  expect_gte(b$upper, b$estimate + 0.5)
  expect_error(sim_confidence_sequence(sim_run(unstated, 10, 1, run_path(), progress = FALSE), "bias"),
               "exact truth")
  # an exact or zero error truth is monitored as before
  cs <- sim_confidence_sequence(sim_run(exact0, 10, 1, run_path(), progress = FALSE), "coverage")
  expect_equal(cs$estimate, 1)
  expect_equal(cs$truth_error, 0)
})
