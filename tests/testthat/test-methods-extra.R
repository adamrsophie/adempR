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

test_that("a changed record function or package list extends a run", {
  gen <- function(scenario, context) data.frame(y = rnorm(scenario$n))
  mk <- function(g) sim_study("rec2", scenarios = sim_scenarios(n = 10L), generate = g,
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("mean", mean(data$y))))
  path <- run_path()
  run <- sim_run(mk(sim_generator(gen, record = function(data) c(n = nrow(data)))),
                 replications = 3, seed = 1, path = path)
  est3 <- sim_results(run)$estimates$estimate
  # a new record function: the run extends, earlier estimates stay, and
  # only records of the new function are summarized
  run2 <- sim_run(mk(sim_generator(gen, record = function(data) c(sd = sd(data$y)))),
                  replications = 5, seed = 1, path = path)
  res <- sim_results(run2)
  expect_equal(res$estimates$estimate[1:3], est3)
  expect_equal(nrow(res$estimates), 5)
  expect_setequal(res$records$name, "sd")
  expect_setequal(res$records$replication, 4:5)
  expect_match(unlist(lapply(res$manifest$revisions, `[[`, "changes")),
               "record function of the generator changed", all = FALSE)
  expect_false("record_fp" %in% names(res$records))
  # declared packages: a revision, not a stop
  run3 <- sim_run(mk(sim_generator(gen, record = function(data) c(sd = sd(data$y)), packages = "stats")),
                  replications = 6, seed = 1, path = path)
  res3 <- sim_results(run3)
  expect_equal(res3$estimates$estimate[1:5], res$estimates$estimate)
  expect_match(unlist(lapply(res3$manifest$revisions, `[[`, "changes")),
               "packages declared by the generator changed (none to stats)", fixed = TRUE, all = FALSE)
  expect_setequal(res3$records$replication, 4:6)
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

test_that("reference checks flag zero MCSE mismatches, use exact intervals and respect scope", {
  m <- list(
    # a constant estimate of 1 against a known bias of 0: MCSE 0, a
    # decisive mismatch
    const = sim_method(function(data) sim_estimate("mean", 1), estimands = "mean",
                       reference = c(bias = 0)),
    # intervals that always cover against a reference coverage of 0.95
    wide = sim_method(function(data) sim_estimate("mean", mean(data$y), lower = -100, upper = 100),
                      estimands = "mean", reference = c(coverage = 0.95)),
    # an exact zero bias with MCSE 0 passes
    exact = sim_method(function(data) sim_estimate("mean", 0), estimands = "mean",
                       reference = c(bias = 0)),
    # a reference that applies to one scenario only
    scoped = sim_method(function(data) sim_estimate("mean", mean(data$y)), estimands = "mean",
                        reference = data.frame(measure = "bias", value = 0, scenario_id = "n=10, mu=0, sigma=1")))
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = m, measures = c("bias", "coverage"))
  run <- sim_run(st, replications = 200, seed = 1, path = run_path())
  rc <- sim_check_reference(run)
  expect_true(all(rc$flag[rc$method == "const"]))
  expect_true(all(is.infinite(rc$z[rc$method == "const"])))
  w <- rc[rc$method == "wide", ]
  expect_true(all(w$flag))
  expect_equal(w$check, rep("exact interval", 2))
  # with 200 hits the exact interval at the 3 sigma level excludes 0.95
  expect_equal(w$lower[1], stats::qbeta(2 * stats::pnorm(-3) / 2, 200, 1), tolerance = 1e-10)
  expect_false(any(rc$flag[rc$method == "exact"]))
  expect_equal(rc$scenario_id[rc$method == "scoped"], "n=10, mu=0, sigma=1")
  expect_error(sim_method(function(data) 1, reference = c(bias = Inf)), "reference")
  bad <- m["scoped"]
  bad$scoped$reference$scenario_id <- "nope"
  st2 <- normal_study(n = 10L, mu = 0, methods = bad, measures = "bias")
  expect_error(sim_check_reference(sim_run(st2, replications = 2, seed = 1, path = run_path())),
               "unknown estimand")
})
