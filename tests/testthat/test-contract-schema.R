library(adempR)

testthat::test_that("A16 grid bias preserves the default extended-real eligibility policy", {
  old <- options(adempR.git_dir = tempdir())
  on.exit(options(old), add = TRUE)
  st <- sim_study("contract_grid_inf", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) 0,
    estimands = sim_estimand_grid("g", c(1, 2), function(scenario, at) at),
    methods = list(m = function(data) sim_estimate_grid("g", c(1, 2), c(Inf, Inf))),
    measures = "bias")
  path <- tempfile("contract-grid-inf-")
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  run <- sim_run(st, 4, 2, path, progress = FALSE)
  out <- sim_summarise_grid(run, "g", "bias")
  avg <- out[out$measure == "grid_mean_bias", , drop = FALSE]
  testthat::expect_equal(avg$n_eligible, 4)
  testthat::expect_identical(avg$value, Inf)
})

testthat::test_that("A17 sampled scenarios reject reserved scientific identity names", {
  testthat::expect_error(
    sim_scenarios_sample(function() list(scenario_key = "user_key", n = 1L),
                         1, seed = 3, distribution = "point mass"),
    "reserved|scenario_key|factor", ignore.case = TRUE)
  testthat::expect_error(
    sim_scenarios_sample(function() list(scenario_id = "user_label", n = 1L),
                         1, seed = 3, distribution = "point mass"),
    "reserved|scenario_id|factor", ignore.case = TRUE)
})
