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

test_that("additions from a stale handle keep earlier additions, and dates match names exactly", {
  st <- normal_study(n = 10L, mu = 0, measures = "mse")
  p <- run_path()
  run <- sim_run(st, replications = 3, seed = 2, path = p)
  sim_add_methods(run, list(a.b = function(data) sim_estimate("mean", median(data$y))))
  run3 <- sim_add_methods(run, list(axb = function(data) sim_estimate("mean", mean(data$y) + 1)))
  expect_setequal(names(run3$manifest$study$methods), c("t_interval", "a.b", "axb"))
  rv <- run3$manifest$revisions
  t_ab <- Filter(function(r) "a.b" %in% r$added_methods, rv)[[1]]$time
  t_axb <- Filter(function(r) "axb" %in% r$added_methods, rv)[[1]]$time
  lb <- sim_leaderboard(run3, "mse", B = 10)
  expect_identical(lb$added[lb$method == "a.b"], t_ab)
  expect_identical(lb$added[lb$method == "axb"], t_axb)
  # a concurrent addition committed between reading and locking is refused
  stale <- sim_open(p)$manifest
  new <- sim_amend(stale$study, methods = c(stale$study$methods,
                                            list(z = function(data) sim_estimate("mean", 0))),
                   reason = "z")
  attr(new, "adempr_base") <- "not-the-current-fingerprint"
  expect_error(sim_run(new, replications = 3, seed = 2, path = p), "changed after it was read")
})
