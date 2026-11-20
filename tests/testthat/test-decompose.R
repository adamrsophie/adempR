test_that("bias splits into statistical and identification parts", {
  st <- sim_study("lognormal", scenarios = sim_scenarios(n = 50L, mu = c(0, 1)),
    generate = function(scenario, context) data.frame(y = exp(rnorm(scenario$n, scenario$mu))),
    estimands = list(sim_estimand("median", sim_truth_exact(function(scenario) exp(scenario$mu)))),
    methods = list(mean = function(data) sim_estimate("median", mean(data$y))))
  run <- sim_run(st, replications = 200, seed = 1, path = run_path())
  d <- sim_bias_decomposition(run, sim_truth_exact(function(scenario) exp(scenario$mu + 0.5)))
  expect_equal(nrow(d), 2L)
  expect_equal(d$total_bias, d$statistical_bias + d$identification_bias)
  expect_equal(d$identification_bias, exp(c(0, 1) + 0.5) - exp(c(0, 1)))
  # the sample mean is unbiased for what it identifies
  expect_true(all(abs(d$statistical_bias) < 3 * d$statistical_mcse))
  expect_true(all(abs(d$total_bias) > 3 * d$total_mcse))
  expect_error(sim_bias_decomposition(run, function(scenario) c(1, 2)), "one number")
})

test_that("identities between estimands are checked on truths and estimates", {
  st <- sim_study("ratio", scenarios = sim_scenarios(n = 40L, mu = c(1, 2)),
    generate = function(scenario, context) data.frame(y = rexp(scenario$n, 1 / scenario$mu)),
    estimands = list(
      sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu)),
      sim_estimand("log_mean", sim_truth_exact(function(scenario) log(scenario$mu) + (scenario$mu == 2)))),
    methods = list(
      good = function(data) list(sim_estimate("mean", mean(data$y)), sim_estimate("log_mean", log(mean(data$y)))),
      bad = function(data) list(sim_estimate("mean", mean(data$y)), sim_estimate("log_mean", mean(log(data$y))))))
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  tr <- sim_check_identity(run, function(mean, log_mean) log(mean) - log_mean)
  expect_equal(tr$flag, c(FALSE, TRUE))
  es <- sim_check_identity(run, function(mean, log_mean) log(mean) - log_mean, target = "estimate")
  expect_false(any(es$flag[es$method == "good"]))
  expect_true(all(es$flag[es$method == "bad"]))
  expect_error(sim_check_identity(run, function(a) a), "estimand names")
})

test_that("the reference follows the truth contract and shared truth error enters the MCSE", {
  st <- sim_study("mcmean", scenarios = sim_scenarios(n = 30L, mu = 1),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, scenario$mu)),
    estimands = list(sim_estimand("mean", sim_truth_mc(function(scenario, n) rnorm(n, scenario$mu), n = 400))),
    methods = list(m = function(data) sim_estimate("mean", mean(data$y))))
  run <- sim_run(st, replications = 30, seed = 2, path = run_path())
  d <- sim_bias_decomposition(run, sim_truth_numeric(function(scenario) list(value = scenario$mu, error = 0.1)))
  expect_equal(d$reference, 1)
  expect_equal(d$reference_error, 0.1)
  tse <- sim_truth(run)$truth_se[1]
  e <- sim_results(run)$estimates
  expect_equal(d$total_mcse, sqrt(var(e$estimate - e$truth) / 30 + tse^2))
  expect_equal(d$identification_mcse, tse)
})

test_that("a data set truth gives the identification part an MCSE", {
  st <- sim_study("sate", scenarios = sim_scenarios(n = 30L),
    generate = function(scenario, context) {
      y <- rnorm(scenario$n)
      sim_data(data.frame(y = y), truth = list(sate = mean(y) / 2))
    },
    estimands = list(sim_estimand("sate", sim_truth_dataset())),
    methods = list(m = function(data) sim_estimate("sate", mean(data$y))))
  run <- sim_run(st, replications = 20, seed = 3, path = run_path())
  d <- sim_bias_decomposition(run, function(scenario) 0)
  e <- sim_results(run)$estimates
  expect_equal(d$identification_bias, -mean(e$truth))
  expect_equal(d$identification_mcse, sd(e$truth) / sqrt(20))
  # identity on data set truths is checked per replication
  chk <- sim_check_identity(run, function(sate) sate - sate)
  expect_equal(nrow(chk), 20L)
  expect_true(all(chk$status == "pass"))
})

test_that("identity checks use truth errors, usable estimates and report unavailable units", {
  st <- sim_study("mcid", scenarios = sim_scenarios(n = 30L, mu = c(1, 2)),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, scenario$mu)),
    estimands = list(
      sim_estimand("a", sim_truth_exact(function(scenario) scenario$mu)),
      sim_estimand("b", sim_truth_mc(function(scenario, n) rnorm(n, scenario$mu), n = 400)),
      sim_estimand("c", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) list(
      sim_estimate("a", mean(data$y)),
      sim_estimate("b", if (mean(data$y) > 1.9) NA_real_ else mean(data$y)),
      sim_estimate("c", 0))))
  run <- sim_run(st, replications = 10, seed = 4, path = run_path())
  tr <- sim_check_identity(run, function(a, b) a - b)
  # the truths differ by Monte Carlo error only: an absolute 1e-8 would fail
  expect_true(all(abs(tr$residual) > 1e-8))
  expect_true(all(tr$residual_se > 0))
  expect_identical(tr$status, c("pass", "pass"))
  es <- sim_check_identity(run, function(a, b) a - b, target = "estimate")
  e <- sim_results(run)$estimates
  nb <- sum(is.na(e$estimate[e$estimand == "b"]))
  expect_gt(nb, 0L)
  expect_equal(sum(es$status == "unavailable"), nb)
  expect_true(all(is.na(es$flag[es$status == "unavailable"])))
  expect_match(es$note[es$status == "unavailable"][1], "no usable value for b")
  expect_equal(nrow(es), 20L)
  # `values` with only the estimands the identity uses
  v <- sim_check_identity(run, function(values) values$a - values$a, estimands = "a")
  expect_identical(v$status, c("pass", "pass"))
  expect_error(sim_check_identity(run, function(values) 0, estimands = "zz"), "Unknown estimand")
  # the tolerance is in residual units: 1 by default, or a declared scale
  expect_identical(sim_check_identity(run, function(a) 0.015, tolerance = 0.01)$status,
                   c("fail", "fail"))
  rel <- sim_check_identity(run, function(a) 0.015, tolerance = 0.01,
                            scale = function(a) max(1, abs(a)))
  expect_identical(rel$status, c("fail", "pass"))
})

test_that("identity tolerances are in the units of the residual", {
  mk <- function(or, log_or, est_or = or, est_log = log_or, numeric_error = NULL) {
    tr_log <- if (is.null(numeric_error)) sim_truth_exact(function(scenario) log_or) else
      sim_truth_numeric(function(scenario) log_or, error = numeric_error)
    sim_study("ident_units", scenarios = sim_scenarios(n = 1L),
      generate = function(scenario, context) 0,
      estimands = list(sim_estimand("or", sim_truth_exact(function(scenario) or)),
                       sim_estimand("log_or", tr_log)),
      methods = list(m = function(data) list(sim_estimate("or", est_or), sim_estimate("log_or", est_log))))
  }
  chk <- function(st, ...) sim_check_identity(sim_run(st, 2, 1, run_path(), progress = FALSE),
                                              function(or, log_or) log(or) - log_or, ...)
  # an exactly false log identity with a huge raw input fails
  r <- chk(mk(1e10, 0))
  expect_identical(r$status, "fail")
  expect_equal(r$residual, log(1e10))
  # a true one with huge inputs passes
  expect_identical(chk(mk(1e10, log(1e10)))$status, "pass")
  # a numerical truth's error bound widens the allowance by its effect on
  # the residual, not by the raw value
  r <- chk(mk(exp(2), 2.001, numeric_error = 0.01))
  expect_identical(r$status, "pass")
  expect_equal(r$residual_error, 0.01, tolerance = 1e-6)
  expect_identical(chk(mk(exp(2), 2.1, numeric_error = 0.01))$status, "fail")
  # a product identity in large units passes with a declared scale
  st <- sim_study("ident_prod", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) 0,
    estimands = list(sim_estimand("itt", sim_truth_exact(function(scenario) 2e6 * 0.5 + 1e-3)),
                     sim_estimand("pc", sim_truth_exact(function(scenario) 0.5)),
                     sim_estimand("cace", sim_truth_exact(function(scenario) 2e6))),
    methods = list(m = function(data) list(sim_estimate("itt", 1), sim_estimate("pc", 1),
                                           sim_estimate("cace", 1))))
  run <- sim_run(st, 2, 1, run_path(), progress = FALSE)
  id <- function(itt, pc, cace) itt - pc * cace
  expect_identical(sim_check_identity(run, id, tolerance = 1e-12)$status, "fail")
  expect_identical(sim_check_identity(run, id, tolerance = 1e-8, scale = function(itt) abs(itt))$status,
                   "pass")
  expect_error(sim_check_identity(run, id, scale = -1), "positive")
})
