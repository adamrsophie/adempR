# Operational tests 13 and 14: random number identity.

test_that("the caller's RNG kind and state are restored, including a missing seed", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  set.seed(42)
  before <- .Random.seed
  run <- sim_run(normal_study(n = 10L, mu = 0), replications = 3, seed = 1, path = run_path())
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1:3], c("Mersenne-Twister", "Inversion", "Rounding"))
  rm(".Random.seed", envir = globalenv())
  sim_run(normal_study(n = 10L, mu = 0), replications = 3, seed = 1, path = run_path())
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the registry is append only and ignores display order", {
  k1 <- data.frame(stage = "dgp", scenario_key = c("b", "a"), component = "data")
  r1 <- rng_registry(k1)
  r1b <- rng_registry(k1[2:1, ])
  expect_identical(r1, r1b)
  r2 <- rng_registry(rbind(k1, data.frame(stage = "dgp", scenario_key = "0", component = "data")), r1)
  expect_identical(r2[1:2, ], r1)
  expect_identical(r2$index[r2$scenario_key == "0"], 3L)
})

test_that("results do not depend on scenario display order, method order or chunk size", {
  m <- stochastic_methods()
  s1 <- normal_study(n = c(10L, 20L), mu = c(0, 1), methods = m)
  s2 <- normal_study(n = c(20L, 10L), mu = c(1, 0), methods = rev(m))
  r1 <- sim_run(s1, replications = 12, seed = 77, path = run_path(), chunk_size = 1)
  r2 <- sim_run(s2, replications = 12, seed = 77, path = run_path(), chunk_size = 7)
  expect_identical(est_signature(r1), est_signature(r2))
})

test_that("extending replications and adding a method preserve existing streams", {
  m <- stochastic_methods()
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = m)
  p <- run_path()
  r100 <- sim_run(st, replications = 10, seed = 5, path = p)
  first <- est_signature(r100)
  reg1 <- r100$manifest$rng_registry
  m3 <- c(m, list(third = function(data) sim_estimate("mean", mean(data$y) + runif(1, -0.01, 0.01))))
  st3 <- normal_study(n = c(10L, 20L), mu = 0, methods = m3)
  r200 <- sim_run(st3, replications = 20, seed = 5, path = p)
  reg2 <- r200$manifest$rng_registry
  expect_identical(reg2[seq_len(nrow(reg1)), ], reg1)
  ext <- est_signature(r200)
  old_rows <- ext[ext$method %in% names(m) & ext$replication <= 10, ]
  rownames(old_rows) <- rownames(first) <- NULL
  expect_identical(old_rows, first)
  # the same extension in a single run of 20 replications gives the same values
  # for the original methods
  r20 <- sim_run(st, replications = 20, seed = 5, path = run_path())
  a <- est_signature(r20)
  b <- ext[ext$method %in% names(m), ]
  rownames(a) <- rownames(b) <- NULL
  expect_identical(a, b)
  # the added method is reproducible by replay
  rp <- sim_replay(r200, scenario = 1, replication = 15, methods = "third")
  expect_true(rp$matches)
})

test_that("the future backend reproduces the sequential reference", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_dev_namespace()
  st <- normal_study(n = c(10L, 20L), mu = c(0, 1), methods = stochastic_methods())
  seq_run <- sim_run(st, replications = 15, seed = 31, path = run_path())
  par_run <- sim_run(st, replications = 15, seed = 31, path = run_path(), backend = "future",
                     workers = 2, chunk_size = 4)
  expect_identical(est_signature(seq_run), est_signature(par_run))
  op <- future::plan(future::sequential)
  on.exit(future::plan(op), add = TRUE)
  par1 <- sim_run(st, replications = 15, seed = 31, path = run_path(), backend = "future",
                  chunk_size = 9)
  expect_identical(est_signature(seq_run), est_signature(par1))
})

test_that("a stochastic extract uses the method stream, not the caller's", {
  st <- sim_study("x", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = sim_method(function(data) mean(data),
      extract = function(fit) sim_estimate("mean", fit + rnorm(1)))))
  set.seed(1)
  r1 <- sim_run(st, 3, seed = 1, path = tempfile("adempR-x-"), progress = FALSE)
  set.seed(2)
  before <- .Random.seed
  r2 <- sim_run(st, 3, seed = 1, path = tempfile("adempR-x-"), progress = FALSE)
  expect_identical(.Random.seed, before)
  expect_identical(sim_results(r1)$estimates$estimate, sim_results(r2)$estimates$estimate)
  expect_true(sim_replay(r1, scenario = 1, replication = 2)$matches)
})

test_that("exact truth functions and generator checks must be deterministic", {
  st <- sim_study("d", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) rnorm(1)))),
    methods = list(m = function(data) sim_estimate("mean", mean(data))))
  before <- .Random.seed
  expect_error(sim_run(st, 2, seed = 1, path = tempfile("adempR-d-"), progress = FALSE),
               "validation failed")
  expect_match(as.data.frame(sim_validate(st, error = FALSE))$message, "must be deterministic", all = FALSE)
  expect_identical(.Random.seed, before)
})

test_that("sim_data(NULL) is recorded as a generation failure", {
  st <- sim_study("n", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) sim_data(NULL),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("mean", mean(data))))
  run <- sim_run(st, 2, seed = 1, path = tempfile("adempR-n-"), progress = FALSE)
  expect_true(all(sim_results(run)$replications$generation_status == "error"))
})

test_that("both backends resolve global helpers and declared packages the same way", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_dev_namespace()
  # top-level helpers, as in a user script
  assign("adempr_test_scale", 2, envir = globalenv())
  assign("adempr_test_helper", eval(quote(function(y) mean(y) * adempr_test_scale), globalenv()),
         envir = globalenv())
  on.exit(rm(list = c("adempr_test_scale", "adempr_test_helper"), envir = globalenv()), add = TRUE)
  m <- list(
    helper = eval(quote(function(data) sim_estimate("mean", adempr_test_helper(data$y))),
                  globalenv()),
    pkg = sim_method(eval(quote(function(data) {
      sim_estimate("mean", mean(data$y) + nchar(file_path_sans_ext("a.txt")))
    }), globalenv()), estimands = "mean", packages = "tools"))
  st <- normal_study(n = 10L, mu = 0, methods = m)
  s <- sim_run(st, replications = 4, seed = 3, path = run_path())
  f <- sim_run(st, replications = 4, seed = 3, path = run_path(), backend = "future", workers = 2)
  expect_true(all(sim_results(s)$estimates$estimate_status == "valid"))
  expect_identical(est_signature(s), est_signature(f))
})

test_that("a missing declared package fails validation", {
  st <- normal_study(n = 10L, mu = 0, methods = list(m = sim_method(
    function(data) sim_estimate("mean", 1), estimands = "mean", packages = "notAPackageXYZ")))
  v <- sim_validate(st, error = FALSE)
  expect_true(any(v$check == "package:notAPackageXYZ" & v$status == "error"))
})

test_that("the future backend refuses closures that reach global objects", {
  assign("adempr_test_helper2", eval(quote(function(y) mean(y)), globalenv()), envir = globalenv())
  on.exit(rm("adempr_test_helper2", envir = globalenv()), add = TRUE)
  make <- eval(quote(function() function(data) sim_estimate("mean", adempr_test_helper2(data$y))),
               globalenv())
  st <- normal_study(n = 10L, mu = 0, methods = list(m = make()))
  expect_identical(adempR:::unshipped_globals(st), "adempr_test_helper2")
  expect_length(adempR:::unshipped_globals(normal_study(n = 10L, mu = 0,
                                                        methods = stochastic_methods())), 0L)
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  expect_error(sim_run(st, replications = 2, seed = 1, path = run_path(), backend = "future"),
               "adempr_test_helper2")
})

test_that("registry lookups refuse aliased indices and fractional replications", {
  reg <- rng_registry(data.frame(stage = "dgp", scenario_key = c("a", "b"), component = "data"))
  bad <- reg
  bad$index <- c(1L, 1L)
  expect_error(rng_task_seed(1, bad, "dgp", "a", "data", 1), "one index")
  expect_error(rng_registry(data.frame(stage = "dgp", scenario_key = "c", component = "data"), bad),
               "one index")
  expect_error(rng_task_seed(1, reg, "dgp", "a", "data", 1.5), "whole number")
  expect_error(rng_task_seed(1, reg, "dgp", "a", "data", 0), "positive")
  expect_error(rng_registry(data.frame(stage = "dgp", scenario_key = NA, component = "data")),
               "complete")
  expect_false(identical(rng_task_seed(1, reg, "dgp", "a", "data", 1),
                         rng_task_seed(1, reg, "dgp", "a", "data", 1, domain = "pilot")))
})

test_that("declared packages are attached only while study code runs", {
  skip_if_not_installed("tools")
  skip_if("package:tools" %in% search(), "tools already attached")
  m <- list(uses = sim_method(function(data) {
    sim_estimate("mean", nchar(file_ext("a.txt")))
  }, estimands = "mean", packages = "tools"))
  st <- normal_study(n = 10L, mu = 0, methods = m)
  before <- search()
  run <- sim_run(st, replications = 2, seed = 1, path = run_path())
  expect_equal(sim_results(run)$estimates$estimate, c(3, 3))
  expect_identical(search(), before)
  expect_true(sim_replay(run, scenario = 1, replication = 1)$matches)
  expect_identical(search(), before)
  sim_validate(st, smoke = TRUE)
  expect_identical(search(), before)
})
