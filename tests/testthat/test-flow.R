test_that("the flow of replications adds up", {
  st <- sim_study("flow", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) {
      if (context$replication == 2) stop("no data")
      rnorm(scenario$n)
    },
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      ok = function(data) sim_estimate("mean", mean(data)),
      bad = function(data, context) {
        if (context$replication == 3) stop("boom")
        if (context$replication == 4) return(sim_estimate("mean", NaN))
        sim_estimate("mean", mean(data))
      }))
  run <- suppressWarnings(sim_run(st, replications = 6, seed = 1, path = run_path()))
  f <- sim_flow(run)
  parts <- c("generation_failed", "not_generated", "not_run", "fit_failed", "excluded", "usable")
  expect_equal(rowSums(f[parts]), f$planned)
  bad <- f[f$method == "bad", ]
  expect_equal(bad$generation_failed, 1L)
  expect_equal(bad$fit_failed, 1L)
  expect_equal(bad$excluded, 1L)
  expect_equal(bad$usable, 3L)
  expect_equal(f$usable[f$method == "ok"], 5L)
})

test_that("acceptability uses the Monte Carlo interval", {
  st <- normal_study(n = 10L, mu = c(0, 2), measures = "rejection")
  run <- sim_run(st, replications = 60, seed = 1, path = run_path())
  a <- sim_acceptable(run)
  expect_equal(a$range_lower[1], 0.025)
  expect_true(all(a$verdict %in% c("inside", "outside", "uncertain")))
  # power near one is far outside the type I error range
  expect_equal(a$verdict[a$scenario_id == st$scenarios$labels[2]], "outside")
  expect_error(sim_acceptable(run, "bias"), "range")
})

test_that("a large sample run flags a method aimed at a different quantity", {
  st <- normal_study(n = 10L, mu = 1, methods = list(
    good = function(data) sim_estimate("mean", mean(data$y), se = sd(data$y) / sqrt(nrow(data))),
    off = function(data) sim_estimate("mean", median(data$y) + 0.5,
                                      se = sd(data$y) / sqrt(nrow(data)))))
  v <- sim_validate(st, large_n = list(n = 20000L), error = FALSE)
  expect_equal(v$status[grepl("large_n:.*:good", v$check)], "ok")
  expect_equal(v$status[grepl("large_n:.*:off", v$check)], "warning")
  expect_error(sim_validate(st, large_n = 1e5), "named list")
})
