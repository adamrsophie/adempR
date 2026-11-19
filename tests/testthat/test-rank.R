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
