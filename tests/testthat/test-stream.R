test_that("batch by batch summaries agree with the in-memory summaries", {
  st <- normal_study(n = c(10L, 30L), mu = c(0, 1), methods = stochastic_methods(),
                     measures = c("bias", "mse", "empse", "mean_se"))
  p <- run_path()
  run <- sim_run(st, replications = 20, seed = 7, path = p, chunk_size = 7)
  run <- sim_run(st, replications = 30, seed = 7, path = p, chunk_size = 7)
  a <- as.data.frame(sim_summarise(run))
  b <- sim_summarise_stream(run, c("bias", "mse", "empse", "mean_se"))
  m <- merge(a[c("scenario_key", "method", "estimand", "measure", "value", "mcse", "n_eligible")],
             b, by = c("scenario_key", "method", "estimand", "measure"))
  expect_equal(nrow(m), nrow(b))
  expect_equal(m$value.x, m$value.y, tolerance = 1e-10)
  expect_equal(m$n_eligible, m$n)
  bm <- m[m$measure %in% c("bias", "mse"), ]
  expect_equal(bm$mcse.x, bm$mcse.y, tolerance = 1e-10)
  # coverage on a study with intervals
  st2 <- normal_study(n = 15L, mu = 0, measures = "coverage")
  r2 <- sim_run(st2, replications = 25, seed = 1, path = run_path(), chunk_size = 4)
  c1 <- as.data.frame(sim_summarise(r2))
  c2 <- sim_summarise_stream(r2, "coverage")
  expect_equal(c2$value, c1$value[c1$measure == "coverage"])
})
