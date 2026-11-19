test_that("fallback values are stored and used only under the policy", {
  m <- list(fragile = sim_method(function(data, context) {
    if (context$replication %% 2 == 0) stop("no convergence")
    sim_estimate("mean", mean(data$y))
  }, estimands = "mean", fallback = function(data) sim_estimate("mean", median(data$y))))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = "bias")
  path <- run_path()
  run <- sim_run(st, replications = 4, seed = 1, path = path)
  e <- sim_results(run)$estimates
  expect_equal(sum(e$fallback), 2L)
  expect_equal(sum(e$valid), 2L)
  expect_true(all(e$attempt_status[e$fallback] == "error"))
  st2 <- st
  st2$failure <- sim_failure_policy(fallback = "use")
  run2 <- sim_run(st2, replications = 4, seed = 1, path = path)
  e2 <- sim_results(run2)$estimates
  expect_equal(sum(e2$valid), 4L)
  expect_equal(sum(e2$used_fallback), 2L)
})

test_that("memory failures, fit diagnostics and reference checks are recorded", {
  m <- list(
    oom = sim_method(function(data) stop("cannot allocate vector of size 9 Gb"), estimands = "mean"),
    diag = sim_method(function(data) list(m = mean(data$y), it = 7L),
                      extract = function(fit) sim_estimate("mean", fit$m),
                      diagnostics = function(fit) c(iterations = fit$it),
                      estimands = "mean"),
    known = sim_method(function(data) sim_estimate("mean", mean(data$y)), estimands = "mean",
                       role = "benchmark", reference = c(bias = 0.5),
                       settings_source = "package defaults"))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = "bias")
  run <- sim_run(st, replications = 30, seed = 1, path = run_path())
  res <- sim_results(run)
  expect_true(all(res$attempts$status[res$attempts$method == "oom"] == "memory"))
  expect_equal(unique(res$diagnostics$value), 7)
  expect_equal(nrow(res$diagnostics), 30L)
  d <- sim_diagnose(run)
  expect_true(d$reference_check$flag[d$reference_check$method == "known"])
  rc <- sim_check_reference(run, threshold = 1e6)
  expect_false(any(rc$flag))
  expect_equal(rc$reference, 0.5)
  expect_null(sim_check_reference(sim_run(normal_study(n = 10L, mu = 0), replications = 2,
                                          seed = 1, path = run_path())))
  expect_error(sim_method(function(data) 1, reference = c(nope = 1)), "reference")
})

test_that("generators record realized properties of each data set", {
  g <- sim_generator(function(scenario, context) {
    y <- rnorm(scenario$n)
    data.frame(y = y, obs = rbinom(scenario$n, 1, 0.7))
  }, record = function(data) c(n = nrow(data), observed = mean(data$obs)))
  st <- sim_study("rec", scenarios = sim_scenarios(n = c(20L, 40L)), generate = g,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("mean", mean(data$y))))
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  r <- sim_results(run)$records
  expect_equal(nrow(r), 2 * 5 * 2)
  d <- sim_diagnose(run)
  expect_equal(sort(d$records$mean[d$records$name == "n"]), c(20, 40))
  expect_true(all(d$records$min[d$records$name == "observed"] >= 0))
})

test_that("validation warns about code that resets the seed", {
  helper <- function(x) {
    set.seed(42)
    mean(x) + rnorm(1, 0, 0.01)
  }
  m <- list(reseeds = sim_method(function(data) sim_estimate("mean", helper(data$y)),
                                 estimands = "mean"),
            fine = sim_method(function(data) sim_estimate("mean", mean(data$y) + rnorm(1, 0, 0.01)),
                              estimands = "mean"))
  st <- normal_study(n = 10L, mu = 0, methods = m)
  v <- sim_validate(st, smoke = TRUE, error = FALSE)
  expect_true(any(v$check == "rng:reseeds" & v$status == "warning"))
  expect_false(any(v$check == "rng:fine"))
  expect_false(any(v$check == "rng:generator"))
})
