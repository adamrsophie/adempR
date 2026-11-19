# Operational test 15 and information boundaries.

test_that("a method mutating its data cannot affect later methods", {
  gen <- function(scenario, context) data.frame(y = rnorm(scenario$n))
  mk <- function(methods) sim_study("iso", scenarios = sim_scenarios(n = 20L), generate = gen,
    estimands = list(sim_estimand("total", sim_truth_exact(function(scenario) 0))),
    methods = methods)
  reader <- function(data) sim_estimate("total", sum(data$y))
  mutator <- function(data) {
    data$y[] <- 1e6
    attr(data, "touched") <- TRUE
    sim_estimate("total", 0)
  }
  r1 <- sim_run(mk(list(a_mutator = mutator, b_reader = reader)), 5, seed = 1, path = run_path())
  r2 <- sim_run(mk(list(b_reader = reader)), 5, seed = 1, path = run_path())
  e1 <- sim_results(r1)$estimates
  e2 <- sim_results(r2)$estimates
  expect_identical(e1$estimate[e1$method == "b_reader"], e2$estimate)
})

test_that("data.table mutation by reference is isolated", {
  skip_if_not_installed("data.table")
  st <- sim_study("dt", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) data.table::data.table(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("total", sim_truth_exact(function(scenario) 0))),
    methods = list(
      a = function(data) {
        data.table::set(data, j = "y", value = 0)
        sim_estimate("total", 0)
      },
      b = function(data) sim_estimate("total", sum(data$y))))
  run <- sim_run(st, 3, seed = 1, path = run_path())
  e <- sim_results(run)$estimates
  expect_true(all(e$estimate[e$method == "b"] != 0))
})

test_that("environments and external pointers are rejected as generated data", {
  st <- sim_study("env", scenarios = sim_scenarios(n = 3L),
    generate = function(scenario, context) list(x = rnorm(3), e = new.env()),
    estimands = list(sim_estimand("m")),
    methods = list(a = function(data) sim_estimate("m", 1)))
  run <- sim_run(st, 2, seed = 1, path = run_path())
  r <- sim_results(run)$replications
  expect_true(all(r$generation_status == "error"))
  expect_match(r$generation_message[1], "environment")
})

test_that("methods see only declared scenario fields; oracles see everything", {
  st <- sim_study("views", scenarios = sim_scenarios(n = 5L, secret = 3.14),
    generate = function(scenario, context) sim_data(rnorm(scenario$n), meta = list(latent = 1)),
    estimands = list(sim_estimand("m", sim_truth_exact(function(scenario) scenario$secret))),
    methods = list(
      public = sim_method(function(data, scenario, context) {
        sim_estimate("m", length(scenario) + 10 * is.null(context$meta) +
                       100 * ("n" %in% names(scenario)))
      }, scenario_fields = "n"),
      oracle = sim_method(function(data, scenario, context) {
        sim_estimate("m", scenario$secret + context$meta$latent)
      }, oracle = TRUE)))
  run <- sim_run(st, 1, seed = 1, path = run_path())
  e <- sim_results(run)$estimates
  expect_equal(e$estimate[e$method == "public"], 1 + 10 + 100)
  expect_equal(e$estimate[e$method == "oracle"], 4.14)
  v <- sim_validate(st, error = FALSE)
  expect_true(any(grepl("oracle", v$message)))
})

test_that("reference objects hidden in attributes are refused, data.table internals allowed", {
  skip_if_not_installed("data.table")
  shared <- new.env()
  assign("v", 1, envir = shared)
  mk <- function(gen) sim_study("attr", scenarios = sim_scenarios(n = 5L), generate = gen,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(a = function(data) {
      e <- attr(data, "box"); if (!is.null(e)) assign("v", 99, envir = e)
      sim_estimate("mean", mean(data$y))
    }, b = function(data) {
      e <- attr(data, "box")
      sim_estimate("mean", if (is.null(e)) mean(data$y) else get("v", envir = e))
    }))
  dt_env <- mk(function(scenario, context) {
    d <- data.table::data.table(y = rnorm(scenario$n))
    data.table::setattr(d, "box", shared)
    d
  })
  r <- sim_run(dt_env, replications = 1, seed = 1, path = run_path())
  expect_equal(sim_results(r)$replications$generation_status, "error")
  expect_match(sim_results(r)$replications$generation_message, "attribute")
  expect_equal(get("v", envir = shared), 1)
  # nested: an environment in an attribute of a list column
  df_env <- mk(function(scenario, context) {
    d <- data.frame(y = rnorm(scenario$n))
    attr(d$y, "box") <- shared
    d
  })
  r2 <- sim_run(df_env, replications = 1, seed = 1, path = run_path())
  expect_equal(sim_results(r2)$replications$generation_status, "error")
  # a plain data.table keeps working
  plain <- mk(function(scenario, context) data.table::data.table(y = c(1, 2, 3)))
  r3 <- sim_run(plain, replications = 1, seed = 1, path = run_path())
  expect_equal(sim_results(r3)$estimates$estimate, c(2, 2))
})
