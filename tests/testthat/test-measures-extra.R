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
  expect_true(is.finite(avg$mcse))
  expect_error(sim_summarise_grid(run, "nope"), "No grid estimand")
})
