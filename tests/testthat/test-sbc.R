sbc_study <- function(methods) {
  sim_study("sbc", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) {
      mu <- rnorm(1)
      sim_data(data.frame(y = rnorm(scenario$n, mu)), truth = list(mu = mu))
    },
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = methods, measures = c("bias", "coverage"))
}

test_that("SBC accepts an exact posterior and rejects an overconfident one", {
  exact <- sim_method_sbc(function(data) {
    n <- nrow(data)
    rnorm(400, sum(data$y) / (n + 1), sqrt(1 / (n + 1)))
  }, "mu")
  narrow <- sim_method_sbc(function(data) {
    n <- nrow(data)
    rnorm(400, sum(data$y) / (n + 1), 0.3 * sqrt(1 / (n + 1)))
  }, "mu")
  run <- sim_run(sbc_study(list(exact = exact, narrow = narrow)), replications = 200, seed = 4,
                 path = run_path())
  good <- sim_sbc(run, "exact", bins = 10)
  bad <- sim_sbc(run, "narrow", bins = 10)
  expect_gt(good$tests$p_value, 0.01)
  expect_lt(bad$tests$p_value, 1e-6)
  expect_equal(sum(good$bins$count), 200L)
  expect_true(all(good$ranks$rank >= 0 & good$ranks$rank <= 100))
  # the estimates are posterior summaries
  s <- as.data.frame(sim_summarise(run))
  expect_lt(s$value[s$method == "narrow" & s$measure == "coverage"], 0.6)
  expect_output(print(good), "Simulation-based calibration")
  expect_error(sim_sbc(run, "nope"), "Unknown method")
})

test_that("SBC ranks are uniform for point masses and small draw counts, and exclusions are counted", {
  # point mass prior and exact point mass posterior: every draw ties
  pm <- sim_study("pm", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) sim_data(data.frame(y = 0), truth = list(mu = 0)),
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = list(point = sim_method_sbc(function(data) rep(0, 9), "mu", n_draws = 9)),
    measures = "bias")
  r <- sim_sbc(sim_run(pm, replications = 200, seed = 1, path = run_path()), "point", bins = 2)
  expect_gt(r$tests$p_value, 0.001)
  expect_equal(r$bins$expected, c(100, 100))
  # nine draws give ten ranks; 20 bins are refused
  expect_error(sim_sbc(sim_run(pm, replications = 2, seed = 1, path = run_path()), "point", bins = 20),
               "exceeds the 10 possible ranks")
  # an exact posterior with few draws: the rank is uniform on 0..9, and
  # 3 bins have exact probabilities 4/10, 3/10, 3/10
  few <- sim_method_sbc(function(data) {
    n <- nrow(data)
    rnorm(9, sum(data$y) / (n + 1), sqrt(1 / (n + 1)))
  }, "mu", n_draws = 9)
  fails <- sim_method_sbc(function(data, context) {
    if (context$replication %% 10 == 0) stop("divergent")
    rnorm(20)
  }, "mu", n_draws = 9)
  run <- sim_run(sbc_study(list(few = few, fails = fails)), replications = 600, seed = 2,
                 path = run_path())
  s <- sim_sbc(run, "few", bins = 3)
  expect_equal(s$bins$expected, 600 * c(4, 3, 3) / 10)
  expect_gt(s$tests$p_value, 0.001)
  expect_equal(sim_sbc(run, "fails", bins = 3)$tests$excluded, 60L)
  expect_error(sim_method_sbc(function(data) 1, "mu", n_draws = 2), "n_draws")
})

test_that("SBC reports every planned scenario and why fits gave no rank", {
  mk <- function(fail_gen = FALSE) sim_study("sbc_fail", scenarios = sim_scenarios(fail = c(FALSE, TRUE)),
    generate = function(scenario, context) {
      if (fail_gen && context$replication == 1L) stop("no data")
      mu <- stats::rnorm(1)
      sim_data(data.frame(y = stats::rnorm(3, mu), fail = scenario$fail), truth = list(mu = mu))
    },
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = list(m = sim_method_sbc(function(data) {
      if (data$fail[1]) stop("fails")
      stats::rnorm(20, sum(data$y) / 4, 0.5)
    }, "mu", n_draws = 9)), measures = "bias")
  out <- sim_sbc(sim_run(mk(TRUE), 6, 3, run_path(), progress = FALSE), "m", bins = 2)
  t <- out$tests
  expect_equal(nrow(t), 2L)
  ok <- t[t$scenario_id == "fail=FALSE", ]
  bad <- t[t$scenario_id == "fail=TRUE", ]
  expect_equal(c(ok$n, ok$generation_failed, ok$excluded), c(5L, 1L, 1L))
  expect_equal(c(bad$n, bad$generation_failed, bad$fit_failed, bad$excluded), c(0L, 1L, 5L, 6L))
  expect_true(is.na(bad$p_value))
  expect_false(is.na(ok$p_value))
  expect_equal(t$excluded, t$generation_failed + t$not_generated + t$fit_failed + t$no_truth_or_draws)
})

test_that("SBC with every fit failing reports counts and no tests", {
  st <- sim_study("sbc_allfail", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) sim_data(data.frame(y = 1), truth = list(mu = 0)),
    estimands = list(sim_estimand("mu", sim_truth_dataset())),
    methods = list(m = sim_method_sbc(function(data) stop("fails"), "mu", n_draws = 9)),
    measures = "bias")
  out <- sim_sbc(sim_run(st, 4, 1, run_path(), progress = FALSE), "m", bins = 2)
  expect_equal(out$tests$fit_failed, 4L)
  expect_equal(out$tests$n, 0L)
  expect_true(is.na(out$tests$p_value))
})
