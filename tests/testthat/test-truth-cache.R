test_that("a lost or damaged truth cache is recomputed only from the original definition", {
  env <- new.env()
  env$offset <- 0
  truth <- local(function(scenario) offset, envir = env)
  st <- sim_study("truth_cache", scenarios = sim_scenarios(n = 3L),
                  generate = function(scenario, context) rep(0, scenario$n),
                  estimands = list(sim_estimand("mean", sim_truth_exact(truth))),
                  methods = list(a = function(data) sim_estimate("mean", 0)), measures = "bias")
  p <- run_path()
  run <- sim_run(st, 2, 1, p, progress = FALSE)
  t0 <- sim_truth(run)
  files <- function() list.files(file.path(p, "truth"), pattern = "\\.rds$", full.names = TRUE)
  # unchanged definition: a missing entry is recomputed to the same value
  unlink(files())
  expect_equal(sim_truth(run)$truth, t0$truth)
  expect_length(files(), 1L)
  # changed definition: an intact entry keeps the original value
  env$offset <- 100
  expect_equal(sim_truth(run)$truth, 0)
  # changed definition and a damaged entry: refused
  writeLines("damaged", files()[1])
  expect_error(sim_truth(run), "definition changed")
  # changed definition and a missing entry: refused, nothing published
  unlink(list.files(file.path(p, "truth"), full.names = TRUE))
  expect_error(sim_truth(run), "definition changed")
  expect_length(files(), 0L)
  # the original definition restored: recomputed under the stored identity
  env$offset <- 0
  expect_equal(sim_truth(run)$truth_fp, t0$truth_fp)
  expect_equal(sim_truth(run)$truth, 0)
})
