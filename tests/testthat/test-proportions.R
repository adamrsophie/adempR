test_that("coverage and usable coverage share one event registry in planning", {
  expect_true(all(c("coverage", "usable_coverage", "rejection", "usable") %in%
                    adempR:::PROPORTION_MEASURES))
  mk <- function(method) sim_study("props", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) stats::rnorm(1),
    estimands = list(sim_estimand("mu", sim_truth_exact(function(scenario) 0))),
    methods = list(m = method), measures = c("coverage", "usable_coverage"))
  req <- function(method) {
    p <- sim_pilot(mk(method), replications = 20, seed = 2, path = run_path(),
                   targets = c(coverage = 0.02, usable_coverage = 0.02))$precision
    stats::setNames(p$required_eligible, p$measure)
  }
  # all results usable: equal requirements at 0, 1 and inside
  for (centre in list(function(d) 10, function(d) 0, function(d) d)) {
    r <- req(function(data) {
      x <- centre(data)
      sim_estimate("mu", x, lower = x - 1, upper = x + 1)
    })
    expect_equal(r[["coverage"]], r[["usable_coverage"]])
    expect_gt(r[["coverage"]], 0)
  }
  # failed fits are misses for usable coverage, but not eligible for coverage
  method <- function(data, context) {
    if (context$replication %% 4L == 0L) stop("no fit")
    sim_estimate("mu", 0, lower = -1, upper = 1)
  }
  s <- as.data.frame(sim_summarise(sim_run(mk(method), 20, 3, run_path(), progress = FALSE)))
  expect_equal(s$value[s$measure == "coverage"], 1)
  expect_equal(s$value[s$measure == "usable_coverage"], 0.75)
  expect_equal(s$n_eligible[s$measure == "usable_coverage"], 20L)
  expect_equal(s$n_eligible[s$measure == "coverage"], 15L)
})
