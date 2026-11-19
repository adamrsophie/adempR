test_that("methods added to a living benchmark see the same data and do not change others", {
  st <- normal_study(n = c(10L, 30L), mu = 0, methods = stochastic_methods(), measures = c("bias", "mse"))
  p <- run_path()
  run <- sim_run(st, replications = 8, seed = 2, path = p)
  before <- est_signature(run)
  run2 <- sim_add_methods(p, list(median = function(data) sim_estimate("mean", median(data$y))))
  after <- est_signature(run2)
  old <- after[after$method %in% names(stochastic_methods()), ]
  rownames(old) <- rownames(before) <- NULL
  expect_identical(old, before)
  # the new method saw the stored data sets: replay agrees
  expect_true(sim_replay(run2, scenario = 1, replication = 3, methods = "median")$matches)
  expect_equal(length(run2$manifest$study$amendments), 1L)
  lb <- sim_leaderboard(run2, "mse", B = 20)
  expect_setequal(lb$method, c("jitter_a", "boot_b", "median"))
  expect_true(all(diff(lb$mean_rank) >= 0))
  rv <- Filter(function(r) "added method median" %in% r$changes, run2$manifest$revisions)
  expect_identical(lb$added[lb$method == "median"], rv[[1]]$time)
  expect_error(sim_add_methods(p, list(median = function(data) 1)), "already in the benchmark")
})
