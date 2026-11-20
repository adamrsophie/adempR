test_that("rankings come with bootstrap uncertainty and a Kendall distance", {
  m <- list(mean = function(data) sim_estimate("mean", mean(data$y)),
            noisy = function(data) sim_estimate("mean", mean(data$y) + rnorm(1, 0, 1)))
  st <- normal_study(n = 20L, mu = 0, methods = m, measures = "mse")
  run <- sim_run(st, replications = 40, seed = 1, path = run_path())
  set.seed(9); before <- .Random.seed
  r <- sim_rank(run, "mse", B = 100, reference = c("mean", "noisy"))
  expect_identical(.Random.seed, before)
  expect_equal(r$rank[r$method == "mean"], 1)
  expect_gt(r$p_first[r$method == "mean"], 0.9)
  expect_equal(sum(r$p_first), 1)
  expect_equal(unique(r$kendall), 0)
  # hand check of the value
  e <- sim_results(run)$estimates
  expect_equal(r$value[r$method == "noisy"], mean(e$estimate[e$method == "noisy"]^2))
  expect_equal(adempR:::kendall_distance(c(1, 2, 3), c(3, 2, 1)), 1)
  expect_identical(sim_rank(run, "mse", B = 50, seed = 2), sim_rank(run, "mse", B = 50, seed = 2))
})

test_that("tied methods share first place, scenarios have their own streams, skips are reported", {
  m <- list(a = function(data) sim_estimate("mean", mean(data$y)),
            b = function(data) sim_estimate("mean", mean(data$y)),
            c = function(data) sim_estimate("mean", mean(data$y) + 1))
  st <- normal_study(n = c(10L, 20L, 30L), mu = 0, methods = m, measures = "mse")
  run <- sim_run(st, replications = 20, seed = 1, path = run_path())
  r <- sim_rank(run, "mse", B = 50)
  s1 <- r[r$scenario_key == st$scenarios$keys[1], ]
  expect_equal(s1$first, c(0.5, 0.5, 0))
  expect_equal(s1$p_first, c(0.5, 0.5, 0))
  lb <- sim_leaderboard(run, "mse", B = 50)
  expect_equal(lb$first[lb$method %in% c("a", "b")], c(1.5, 1.5))
  # a scenario's result does not depend on the others
  noisy <- list(x = function(data) sim_estimate("mean", mean(data$y) + rnorm(1, 0, 0.3)),
                y = function(data) sim_estimate("mean", median(data$y)))
  st2 <- normal_study(n = c(10L, 20L), mu = 0, methods = noisy, measures = "mse")
  run2 <- sim_run(st2, replications = 20, seed = 1, path = run_path())
  full <- sim_rank(run2, "mse", B = 50)
  k <- st2$scenarios$keys[2]
  res <- sim_results(run2)
  res$estimates <- res$estimates[res$estimates$scenario_key == k, ]
  alone <- sim_rank(res, "mse", B = 50)
  expect_equal(alone$p_first[alone$scenario_key == k], full$p_first[full$scenario_key == k])
  # the first scenario now has no estimates: reported, not dropped
  expect_equal(unique(alone$status[alone$scenario_key != k]),
               "fewer than two replications usable by every method")
  expect_error(sim_rank(run2, methods = c("x", "nope")), "do not exist")
})

test_that("EmpSE ranking uses the errors when the truth varies", {
  st <- sim_study("dsr", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) sim_data(data.frame(y = context$replication),
                                                    truth = list(mean = context$replication)),
    estimands = list(sim_estimand("mean", sim_truth_dataset())),
    methods = list(exact = function(data) sim_estimate("mean", data$y),
                   off = function(data) sim_estimate("mean", data$y + rnorm(1)) ),
    measures = "empse")
  r <- sim_rank(sim_run(st, replications = 10, seed = 1, path = run_path()), "empse", B = 20)
  expect_equal(r$value[r$method == "exact"], 0)
  expect_equal(r$rank[r$method == "exact"], 1)
})

test_that("a ranking that an uncertain truth can reverse is labeled", {
  mk <- function(se) sim_study("rank_truth", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_mc(function(scenario, n) {
      sim_truth_value(0.6, se = se, n = n)
    }, n = 10))),
    methods = list(a = function(data) sim_estimate("mean", 0),
                   b = function(data) sim_estimate("mean", 1)))
  run <- sim_run(mk(1), replications = 10, seed = 1, path = run_path())
  r <- sim_rank(run, "bias", B = 20)
  # b wins for the reported truth 0.6, a wins for truths below 0.5
  expect_equal(r$p_first[r$method == "b"], 1)
  expect_true(all(r$truth_conditional))
  expect_equal(unique(r$truth_halfwidth), stats::qnorm(0.975) * sim_truth(run)$truth_se, tolerance = 1e-8)
  expect_true(all(r$first_possible))
  expect_true(all(r$status == "truth_sensitive"))
  # a precise truth cannot reverse it
  r2 <- sim_rank(sim_run(mk(0.01), replications = 10, seed = 1, path = run_path()), "bias", B = 20)
  expect_true(all(r2$status == "ok"))
  expect_equal(r2$first_possible, c(FALSE, TRUE))
  # the empirical SE does not use the truth
  r3 <- sim_rank(sim_run(mk(1), replications = 10, seed = 1, path = run_path()), "empse", B = 20)
  expect_false(any(r3$truth_conditional))
})
