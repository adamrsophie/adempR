test_that("categorical estimands get selection proportions and accuracy", {
  st <- sim_study("select", scenarios = sim_scenarios(n = c(20L, 200L)),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, 0.3)),
    estimands = list(sim_estimand("model", type = "categorical",
                                  truth = sim_truth_exact(function(scenario) 2))),
    methods = list(test = function(data) sim_estimate("model", if (t.test(data$y)$p.value < 0.05) 2 else 1)),
    measures = sim_measures("usable", sim_measures_categorical(1:2, c("null", "alt"))))
  run <- sim_run(st, replications = 40, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  e <- sim_results(run)$estimates
  for (k in unique(s$scenario_key)) {
    sk <- s[s$scenario_key == k, ]
    ek <- e[e$scenario_key == k, ]
    expect_equal(sk$value[sk$measure == "correct"], mean(ek$estimate == 2))
    expect_equal(sk$value[sk$measure == "select_null"] + sk$value[sk$measure == "select_alt"], 1)
  }
  expect_error(sim_estimand("x", type = "matrix"), "not supported")
})

test_that("categorical measures apply only to categorical estimands and count invalid classes", {
  st <- sim_study("mixed", scenarios = sim_scenarios(n = 30L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, 0.3)),
    estimands = list(
      sim_estimand("model", type = "categorical", truth = sim_truth_exact(function(scenario) 2)),
      sim_estimand("mean", truth = sim_truth_exact(function(scenario) 0.3))),
    methods = list(m = function(data) {
      # class 3 is not a declared class: an invalid output
      cl <- if (data$y[1] > 1.2) 3 else if (t.test(data$y)$p.value < 0.05) 2 else 1
      list(sim_estimate("model", cl), sim_estimate("mean", mean(data$y)))
    }),
    measures = sim_measures("bias", "usable", sim_measures_categorical(1:2)))
  run <- sim_run(st, replications = 60, seed = 5, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  e <- sim_results(run)$estimates
  em <- e[e$estimand == "model", ]
  expect_true(any(em$estimate == 3))
  sm <- s[s$estimand == "model", ]
  v <- stats::setNames(sm$value, sm$measure)
  expect_equal(v[["select_1"]] + v[["select_2"]] + v[["select_other"]], 1)
  expect_equal(v[["select_other"]], mean(em$estimate == 3))
  expect_equal(v[["correct"]], mean(em$estimate == 2))
  expect_identical(sm$status[sm$measure == "bias"], "not_applicable")
  expect_equal(sm$value[sm$measure == "usable"], 1)
  ss <- s[s$estimand == "mean", ]
  expect_true(all(ss$status[ss$measure %in% c("correct", "select_1", "select_2", "select_other")] ==
                    "not_applicable"))
  expect_equal(ss$value[ss$measure == "bias"], mean(e$estimate[e$estimand == "mean"]) - 0.3)
  expect_error(sim_measure("x", function(rows) 1, estimand_type = "vector"), "estimand_type")
  # the streaming summary scopes built in measures the same way
  st_s <- sim_summarise_stream(run, "bias")
  expect_identical(st_s$status[st_s$estimand == "model"], "not_applicable")
  expect_equal(st_s$value[st_s$estimand == "mean"], ss$value[ss$measure == "bias"])
})

test_that("the interval score matches a hand computation", {
  st <- normal_study(n = 20L, mu = 1, measures = sim_measures("coverage", sim_measure_interval_score(),
                                                              sim_measure_interval_score("is_rel", scale = "truth")))
  run <- sim_run(st, replications = 30, seed = 2, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  e <- sim_results(run)$estimates
  is <- (e$upper - e$lower) + 40 * (e$lower - 1) * (1 < e$lower) + 40 * (1 - e$upper) * (1 > e$upper)
  expect_equal(s$value[s$measure == "interval_score"], mean(is))
  expect_equal(s$mcse[s$measure == "interval_score"], sd(is) / sqrt(30))
  expect_equal(s$value[s$measure == "is_rel"], mean(is))
  expect_error(sim_measure_interval_score(scale = -1), "scale")
  expect_error(sim_measure_interval_score(scale = Inf), "scale")
  # contributions: unbounded sides, vector scale, invalid level and scale
  f <- sim_measure_interval_score()$contribution
  rows <- data.frame(level = c(0.95, 0.95, 0.95, 1), truth = c(0, 3, 0, 0),
                     lower = c(-Inf, -Inf, -1, -1), upper = c(Inf, 1, 1, 1), se = c(1, 2, 0, 1))
  v <- f(rows)
  expect_identical(v[1:2], c(Inf, Inf))
  expect_equal(v[3], 2)
  expect_true(is.na(v[4]))
  g <- sim_measure_interval_score("s", scale = 2)$contribution
  expect_equal(g(rows[3, ][c(1, 1), ]), c(1, 1))
  h <- sim_measure_interval_score("s", scale = "se")$contribution
  expect_true(is.na(h(rows)[3]))
  expect_error(sim_measure_interval_score("s", scale = "nope")$contribution(rows), "not a column")
  # a one sided interval that misses on its bounded side
  rows2 <- data.frame(level = 0.9, truth = 2, lower = -Inf, upper = 1)
  expect_identical(f(rows2), Inf)
})

test_that("grid estimands give pointwise and grid averaged performance", {
  grid <- c(0.5, 1, 2)
  st <- sim_study("curve", scenarios = sim_scenarios(n = 100L, rate = 0.5),
    generate = function(scenario, context) data.frame(t = rexp(scenario$n, scenario$rate)),
    estimands = sim_estimand_grid("surv", grid, truth = function(scenario, at) exp(-scenario$rate * at)),
    methods = list(km = function(data) {
      p <- vapply(grid, function(x) mean(data$t > x), numeric(1))
      se <- sqrt(p * (1 - p) / nrow(data))
      sim_estimate_grid("surv", grid, p, se = se, lower = p - 1.96 * se, upper = p + 1.96 * se,
                        level = 0.95)
    }), measures = c("bias", "mse", "coverage"))
  run <- sim_run(st, replications = 40, seed = 3, path = run_path())
  tr <- sim_results(run)$truth
  expect_equal(sort(tr$truth), sort(exp(-0.5 * grid)))
  g <- sim_summarise_grid(run, "surv", "mse")
  expect_equal(sum(!is.na(g$at)), 3L)
  avg <- g[g$measure == "grid_mean_mse", ]
  expect_equal(avg$value, mean(g$value[!is.na(g$at)]))
  # the MCSE from replication-level averages, computed by hand
  e <- sim_results(run)$estimates
  per <- tapply((e$estimate - e$truth)^2, e$replication, mean)
  expect_equal(avg$mcse, sd(per) / sqrt(length(per)), tolerance = 1e-12)
  expect_equal(avg$n_eligible, 40L)
  expect_equal(avg$n_incomplete, 0L)
  expect_error(sim_summarise_grid(run, "nope"), "No grid estimand")
})

test_that("grid estimates match by coordinate, check lengths and follow eligibility", {
  grid <- c(1, 2)
  expect_error(sim_estimate_grid("g", c(1, 2, 3), c(1, 2, 3), se = c(0.1, 0.2)),
               "one entry per grid point")
  expect_error(sim_estimate_grid("g", grid, 1), "one entry per grid point")
  mk <- function(method) sim_study("grid2", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) data.frame(y = 0),
    estimands = sim_estimand_grid("g", grid, truth = function(scenario, at) at),
    methods = list(m = method), measures = c("bias", "coverage"))
  # the method lists the grid in reverse order with matching values
  rev_m <- function(data) sim_estimate_grid("g", rev(grid), rev(grid))
  b <- sim_summarise_grid(sim_run(mk(rev_m), replications = 2, seed = 1, path = run_path()), "g")
  expect_equal(b$value[!is.na(b$at)], c(0, 0))
  # an estimate at a point the estimand does not have is refused
  bad <- function(data) sim_estimate_grid("g", c(1, 3), c(1, 3))
  r <- sim_run(mk(bad), replications = 1, seed = 1, path = run_path())
  expect_match(sim_results(r)$attempts$message[1], "not a point of that grid")
  # nonconverged fits on even replications have intervals that miss; the
  # grid coverage excludes them exactly like the pointwise coverage
  nc <- function(data, context) {
    odd <- context$replication %% 2 == 1
    lapply(sim_estimate_grid("g", grid, grid, lower = grid + if (odd) -1 else 5,
                             upper = grid + if (odd) 1 else 6), function(e) {
      e$converged <- odd
      e
    })
  }
  r <- sim_run(mk(nc), replications = 10, seed = 1, path = run_path())
  g <- sim_summarise_grid(r, "g", "coverage")
  expect_equal(g$value[!is.na(g$at)], c(1, 1))
  avg <- g[g$measure == "grid_mean_coverage", ]
  expect_equal(avg$value, 1)
  expect_equal(avg$n_eligible, 5L)
  expect_equal(avg$n_incomplete, 5L)
})
