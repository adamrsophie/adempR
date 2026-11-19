test_that("resume refuses new work when a stored method's globals changed", {
  assign("adempr_test_shift", 1, envir = globalenv())
  on.exit(rm("adempr_test_shift", envir = globalenv()), add = TRUE)
  m <- list(shifted = function(data) sim_estimate("mean", mean(data$y) + adempr_test_shift))
  st <- normal_study(n = 10L, mu = 0, methods = m)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  assign("adempr_test_shift", 7, envir = globalenv())
  expect_error(sim_resume(path, replications = 4), "differ from those of the committed results")
  # reading is still allowed and unchanged
  expect_equal(nrow(sim_results(sim_open(path))$estimates), 2L)
  # with the original global restored, the extension runs
  assign("adempr_test_shift", 1, envir = globalenv())
  expect_no_error(sim_resume(path, replications = 4))
})

test_that("rereading a completed run reports the truth evaluated with it", {
  assign("adempr_test_truth", 1, envir = globalenv())
  on.exit(rm("adempr_test_truth", envir = globalenv()), add = TRUE)
  st <- normal_study(n = 10L, mu = 0)
  st$estimands$mean$truth <- sim_truth_exact(function(scenario) adempr_test_truth)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  b1 <- sim_truth(sim_open(path))$truth
  assign("adempr_test_truth", 3, envir = globalenv())
  expect_equal(sim_truth(sim_open(path))$truth, b1)
  expect_equal(b1, 1)
})

test_that("formulas are identified by the values their environment supplies", {
  mk <- function(k) { force(k); y ~ I(x * k) }
  expect_false(identical(fingerprint(mk(1)), fingerprint(mk(2))))
  expect_identical(fingerprint(mk(1)), fingerprint(mk(1)))
  # data columns that the environment does not supply do not matter
  mk2 <- function() y ~ x
  expect_identical(fingerprint(mk2()), fingerprint(local(y ~ x, envir = new.env())))
})

test_that("environments are encoded with their parents and functions without a depth cutoff", {
  outer <- new.env(parent = globalenv())
  assign("a", 1, envir = outer)
  e1 <- new.env(parent = outer)
  outer2 <- new.env(parent = globalenv())
  assign("a", 2, envir = outer2)
  e2 <- new.env(parent = outer2)
  expect_false(identical(fingerprint(e1), fingerprint(e2)))
  # a chain of eight nested helpers: the innermost value still matters
  chain <- function(v) {
    f <- function() v
    for (i in 1:8) { g <- f; f <- local(function() g(), list(g = g)) }
    f
  }
  expect_false(identical(fingerprint(chain(1)), fingerprint(chain(2))))
  # mutually recursive functions terminate
  env <- new.env(parent = globalenv())
  local({ f <- function(x) if (x > 0) g(x - 1) else 0; g <- function(x) f(x) }, envir = env)
  expect_type(fingerprint(env$f), "character")
})

test_that("scenario values edited after sim_scenarios() are refused", {
  st <- normal_study(n = 1L, mu = 0)
  st$scenarios$values[[1]]$n <- 100L
  expect_error(expect_output(sim_validate(st)), "validation failed")
  expect_match(sim_validate(st, error = FALSE)$message, "edited after sim_scenarios", all = FALSE)
})

test_that("records, packages and the GPU flag enter fingerprints only when used", {
  g0 <- sim_generator(function(scenario, context) 1)
  g1 <- sim_generator(function(scenario, context) 1, record = function(data) c(a = 1))
  g2 <- sim_generator(function(scenario, context) 1, record = function(data) c(a = 2))
  g3 <- sim_generator(function(scenario, context) 1, packages = "stats")
  expect_identical(generator_fingerprint(g0),
                   fingerprint(g0$fun, g0$version, g0$depends))
  fps <- vapply(list(g0, g1, g2, g3), generator_fingerprint, character(1))
  expect_equal(length(unique(fps)), 4L)
  f <- function(data) sim_estimate("mean", 1)
  m0 <- sim_method(f)
  expect_false(identical(method_fingerprint(m0), method_fingerprint(sim_method(f, packages = "stats"))))
  expect_false(identical(method_fingerprint(m0), method_fingerprint(sim_method(f, gpu = TRUE))))
})

test_that("a retry never pairs fits from different data sets", {
  # external inputs (files) are invisible to fingerprints
  input <- tempfile(); flag <- tempfile()
  writeLines("1", input); writeLines("fail", flag)
  st <- sim_study("pairing", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) data.frame(y = as.numeric(readLines(input))),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      a = function(data) sim_estimate("mean", data$y),
      b = function(data) {
        if (file.exists(flag)) stop("license server unavailable")
        sim_estimate("mean", data$y)
      }))
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  unlink(flag)
  writeLines("2", input)
  expect_error(sim_run(st, replications = 2, seed = 1, path = path, retry_failed = TRUE),
               "generated different data than the committed data set")
  # with the input restored, the retry pairs b with the same data as a
  writeLines("1", input)
  r <- sim_run(st, replications = 2, seed = 1, path = path, retry_failed = TRUE)
  e <- sim_results(r)$estimates
  expect_equal(e$estimate[e$method == "a"], e$estimate[e$method == "b"])
})

test_that("diagnostics of an earlier round do not survive a retry that produced none", {
  flag <- tempfile(); writeLines("x", flag)
  m <- sim_method(function(data) {
    if (file.exists(flag)) stop("cannot allocate vector of size 9 Gb")
    sim_estimate("mean", mean(data$y))
  }, diagnostics = function(fit) NULL)
  st <- normal_study(n = 5L, mu = 0, methods = list(m = m))
  path <- run_path()
  sim_run(st, replications = 1, seed = 1, path = path)
  unlink(flag)
  r <- sim_run(st, replications = 1, seed = 1, path = path, retry_failed = TRUE)
  res <- sim_results(r)
  expect_equal(res$estimates$attempt_status, "ok")
  expect_equal(NROW(res$diagnostics), 0L)
})

test_that("a method's declared packages do not mask functions for other methods", {
  skip_if_not_installed("dplyr")
  ms <- list(
    uses_stats = function(data) sim_estimate("mean", sum(filter(1:3, 1))),
    declares_dplyr = sim_method(function(data) sim_estimate("mean", nrow(dplyr::tibble(a = 1))),
                                packages = "dplyr"))
  st <- normal_study(n = 5L, mu = 0, methods = ms)
  r <- sim_run(st, replications = 2, seed = 1, path = run_path())
  e <- sim_results(r)$estimates
  expect_equal(e$estimate[e$method == "uses_stats"], c(6, 6))
  expect_false("package:dplyr" %in% search())
})
