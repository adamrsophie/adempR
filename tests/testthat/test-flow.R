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
  expect_true(a$verdict[a$label == "type_i_error"] %in% c("inside", "outside", "uncertain"))
  # the type I error range is not applied to a power cell
  pw <- a$scenario_id == st$scenarios$labels[2]
  expect_identical(a$label[pw], "power")
  expect_true(is.na(a$verdict[pw]))
  expect_match(a$verdict_note[pw], "null hypothesis holds")
  # an explicit range applies everywhere
  b <- sim_acceptable(run, range = c(0.8, 1))
  expect_equal(b$verdict[pw], "inside")
  expect_error(sim_acceptable(run, "bias"), "range")
})

test_that("rejection decisions use the study alpha and numerical truth bounds", {
  st <- normal_study(n = 10L, mu = 0, measures = "rejection", alpha = 0.1)
  run <- sim_run(st, replications = 20, seed = 1, path = run_path())
  # the default range uses the study's alpha
  expect_equal(sim_acceptable(run)$range_lower, 0.05)
  # diagnostics compare a test at the study's alpha with an interval at
  # 1 - alpha; a 95% interval is not compared with a test at 0.1
  expect_null(sim_diagnose(run)$test_interval_disagreement)
  # a numerical truth whose error bound reaches the null is not labeled power
  st2 <- sim_study("numnull", scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("mean", sim_truth_numeric(function(scenario) list(value = 0.01, error = 0.1)))),
    methods = list(t = function(data) sim_estimate("mean", mean(data$y), p_value = t.test(data$y)$p.value, null = 0)),
    measures = "rejection")
  s <- as.data.frame(sim_summarise(sim_run(st2, replications = 5, seed = 1, path = run_path())))
  expect_identical(s$label[s$measure == "rejection"], "rejection")
})

test_that("a sequential run is judged with its monitoring intervals", {
  st <- normal_study(n = 10L, mu = 0, measures = "rejection")
  sq <- sim_run_sequential(st, looks = c(100, 200), half_width = 0.2, seed = 1,
                           path = run_path(), measures = "rejection")
  a <- sim_acceptable(sq$run)
  s <- as.data.frame(sim_summarise(sq$run, measures = "rejection"))
  s <- s[s$measure == "rejection", ]
  ci <- adempR:::cp_interval(round(s$value * s$n_eligible), s$n_eligible,
                             sq$conf_level_per_interval)
  expect_equal(a$lower, unname(ci[, "lower"]))
  expect_equal(a$upper, unname(ci[, "upper"]))
  expect_identical(a$interval, "sequential monitoring")
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
