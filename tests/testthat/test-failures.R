# Scientific tests 3, 4 and 6: failures, partial outputs and edge cases.

selective_study <- function() {
  sim_study("selective",
    scenarios = sim_scenarios(n = 25L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      always = function(data) {
        tt <- t.test(data)
        sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
                     lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95)
      },
      # fails exactly when its estimate would be large: informative failure
      fragile = function(data) {
        if (mean(data) > 0.1) stop("did not converge")
        tt <- t.test(data)
        sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
                     lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95)
      }),
    measures = c("bias", "coverage", "usable_coverage"))
}

test_that("conditional, common-valid and composite performance are different quantities", {
  run <- sim_run(selective_study(), replications = 400, seed = 3, path = run_path())
  res <- sim_results(run)
  e <- res$estimates
  fr <- e[e$method == "fragile", ]
  ok <- fr$valid
  expect_true(any(!ok))
  expect_true(all(fr$estimate_status[!ok] == "error"))
  s <- as.data.frame(sim_summarise(run))
  fs <- s[s$method == "fragile", ]
  # conditional bias uses only the usable estimates and is selected downwards
  expect_equal(fs$value[fs$measure == "bias"], mean(fr$estimate[ok]))
  expect_equal(fs$n_eligible[fs$measure == "bias"], sum(ok))
  expect_lt(fs$value[fs$measure == "bias"], 0)
  # usable probability over all generated data sets
  expect_equal(fs$value[fs$measure == "usable"], mean(ok))
  expect_equal(fs$n_generated[fs$measure == "usable"], 400L)
  # composite: usable interval containing the truth, over all generated data sets
  covered <- fr$lower <= 0 & 0 <= fr$upper
  expect_equal(fs$value[fs$measure == "usable_coverage"], sum(covered[ok]) / 400)
  expect_equal(fs$value[fs$measure == "coverage"], mean(covered[ok]))
  expect_false(isTRUE(all.equal(fs$value[fs$measure == "usable_coverage"],
                                fs$value[fs$measure == "coverage"])))
  # paired comparison on the common valid population
  cmp <- sim_compare(run, c("always", "fragile"), measures = c("bias", "usable"))
  expect_equal(cmp$n_common[cmp$measure == "bias"], sum(ok))
  expect_equal(cmp$n_common[cmp$measure == "usable"], 400L)
  expect_equal(cmp$value[cmp$measure == "usable"], mean(ok) - 1)
  # failures are reported by default
  f <- attr(sim_summarise(run), "failures")
  expect_true(any(grepl("error:", f$statuses[f$method == "fragile"])))
})

test_that("estimates without SE count for bias but not for SE calibration; intervals stand alone", {
  study <- sim_study("partial",
    scenarios = sim_scenarios(n = 30L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      no_se = function(data, context) {
        sim_estimate("mean", mean(data),
                     se = if (context$replication %% 2 == 0) NA else sd(data) / sqrt(30))
      },
      ci_only = function(data) {
        tt <- t.test(data)
        sim_estimate("mean", NA, lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95)
      }),
    measures = c("bias", "mean_se", "coverage"))
  run <- sim_run(study, replications = 40, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  ns <- s[s$method == "no_se", ]
  expect_equal(ns$n_eligible[ns$measure == "bias"], 40L)
  expect_equal(ns$n_eligible[ns$measure == "mean_se"], 20L)
  ci <- s[s$method == "ci_only", ]
  expect_equal(ci$n_eligible[ci$measure == "coverage"], 40L)
  expect_identical(ci$status[ci$measure == "bias"], "no_eligible")
  expect_equal(ci$value[ci$measure == "usable"], 0)
})

test_that("zero truth, infinite estimates, boundaries and invalid intervals are handled deliberately", {
  study <- sim_study("edges",
    scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      plain = function(data) sim_estimate("mean", mean(data)),
      blowup = function(data, context) {
        sim_estimate("mean", if (context$replication == 2) Inf else mean(data),
                     lower = 1, upper = -1, boundary = context$replication == 3)
      }),
    measures = c("bias", "relbias", "mse", "coverage"))
  run <- sim_run(study, replications = 10, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  expect_identical(s$status[s$method == "plain" & s$measure == "relbias"], "undefined")
  bl <- s[s$method == "blowup", ]
  expect_identical(bl$status[bl$measure == "bias"], "infinite")
  expect_equal(bl$value[bl$measure == "mse"], Inf)
  expect_identical(bl$status[bl$measure == "coverage"], "no_eligible")
  est <- sim_results(run)$estimates
  expect_true(all(est$interval_status[est$method == "blowup"] == "invalid"))
  expect_identical(est$estimate_status[est$method == "blowup" & est$replication == 2], "infinite")
  # a stricter policy makes infinite and boundary estimates unusable
  strict <- sim_amend(study, failure = sim_failure_policy(infinite = "invalid", boundary = "invalid"),
                      reason = "stricter policy")
  run2 <- sim_run(strict, replications = 10, seed = 1, path = run$path)
  e2 <- sim_results(run2)$estimates
  b2 <- e2[e2$method == "blowup", ]
  expect_identical(b2$estimate_status[b2$replication == 2], "invalid_infinite")
  expect_identical(b2$estimate_status[b2$replication == 3], "boundary")
  # one-sided intervals with infinite limits are valid
  expect_identical(adempR:::classify_rows(data.frame(
    estimand = "mean", estimate = 1, se = NA, lower = -Inf, upper = 2, p_value = NA, reject = NA,
    converged = NA, boundary = NA, declared_effect = NA, declared_scale = NA,
    attempt_status = "ok", had_warning = FALSE),
    list(mean = sim_estimand("mean")), sim_failure_policy())$interval_status, "valid")
})

test_that("warnings, nonconvergence, timeouts and extraction errors get separate statuses", {
  study <- sim_study("statuses",
    scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      warns = function(data) {
        warning("numerical issue")
        sim_estimate("mean", mean(data))
      },
      nonconv = function(data) sim_estimate("mean", mean(data), converged = FALSE),
      badtype = function(data) list(a = 1),
      slow = sim_method(function(data) {
        repeat {
          Sys.sleep(0.05)
        }
      }, timeout = 0.3),
      split = sim_method(function(data) lm(data ~ 1),
                         extract = function(fit) sim_estimate("mean", unname(coef(fit)))))
  )
  run <- sim_run(study, replications = 2, seed = 1, path = run_path())
  res <- sim_results(run)
  st <- function(m) unique(res$estimates$estimate_status[res$estimates$method == m])
  expect_identical(st("warns"), "valid")
  expect_true(all(res$estimates$had_warning[res$estimates$method == "warns"]))
  expect_identical(st("nonconv"), "nonconverged")
  expect_identical(st("badtype"), "extraction_error")
  expect_identical(st("slow"), "timeout")
  expect_identical(st("split"), "valid")
  att <- res$attempts
  expect_true(all(att$elapsed[att$method == "slow"] >= 0.25))
  expect_identical(unique(att$failed_stage[att$method == "badtype"]), "extract")
})

test_that("generation failures are recorded and never redrawn", {
  study <- sim_study("genfail",
    scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) {
      if (context$replication == 4) stop("impossible draw")
      rnorm(scenario$n)
    },
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(avg = function(data) sim_estimate("mean", mean(data))))
  run <- sim_run(study, replications = 6, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  expect_equal(s$n_generated[1], 5L)
  expect_equal(s$n_generation_failed[1], 1L)
  expect_equal(s$n_planned[1], 6L)
  reps <- sim_results(run)$replications
  expect_identical(reps$generation_status[reps$replication == 4], "error")
})
