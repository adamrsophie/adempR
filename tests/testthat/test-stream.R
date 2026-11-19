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
  expect_equal(nrow(b), nrow(a[a$measure %in% c("bias", "mse", "empse", "mean_se"), ]))
  expect_equal(m$value.x, m$value.y, tolerance = 1e-10)
  expect_equal(m$n_eligible, m$n)
  # MCSEs agree for every measure, including the kurtosis based EmpSE MCSE
  expect_equal(m$mcse.x, m$mcse.y, tolerance = 1e-10)
  # coverage on a study with intervals
  st2 <- normal_study(n = 15L, mu = 0, measures = "coverage")
  r2 <- sim_run(st2, replications = 25, seed = 1, path = run_path(), chunk_size = 4)
  c1 <- as.data.frame(sim_summarise(r2))
  c2 <- sim_summarise_stream(r2, "coverage")
  expect_equal(c2$value, c1$value[c1$measure == "coverage"])
})

test_that("streaming summaries agree on data set truths, boundaries, empty cells and conflicts", {
  st <- sim_study("dst", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) sim_data(data.frame(y = context$replication),
                                                    truth = list(mean = context$replication)),
    estimands = list(sim_estimand("mean", sim_truth_dataset())),
    methods = list(exact = function(data) sim_estimate("mean", data$y, lower = data$y - 1, upper = data$y + 1),
                   none = function(data) sim_estimate("mean", NA_real_)),
    measures = c("empse", "coverage", "bias"))
  run <- sim_run(st, replications = 6, seed = 1, path = run_path(), chunk_size = 2)
  a <- as.data.frame(sim_summarise(run))
  b <- sim_summarise_stream(run, c("bias", "empse", "coverage"))
  m <- merge(a[c("method", "measure", "value", "mcse", "status")], b,
             by = c("method", "measure"))
  expect_equal(nrow(m), 6L)
  # the estimator equals the data set truths: EmpSE of the errors is 0
  expect_equal(b$value[b$method == "exact" & b$measure == "empse"], 0)
  expect_equal(m$status.x, m$status.y)
  expect_equal(m$value.x, m$value.y)
  # coverage 1 is a boundary with an exact interval in both
  cv <- b[b$method == "exact" & b$measure == "coverage", ]
  expect_equal(cv$status, "boundary")
  expect_equal(cv$lower, stats::qbeta(0.025, 6, 1))
  # a method with no usable estimate still has its cells
  expect_true(all(b$status[b$method == "none"] == "no_eligible"))
  # a conflicting duplicate stops, as in sim_results()
  f <- list.files(file.path(run$path, "batches"), pattern = "rds$", full.names = TRUE)[1]
  o <- readRDS(f)
  o$estimates$estimate <- o$estimates$estimate + 1
  g <- file.path(dirname(f), "zz-conflict.rds")
  saveRDS(o, g)
  writeLines(adempR:::file_md5(g), sub("rds$", "md5", g))
  expect_error(sim_summarise_stream(run, "bias"), "Conflicting duplicate")
  expect_error(sim_results(run), "Conflicting duplicate")
})

test_that("accumulated fourth moments merge exactly", {
  x <- rexp(500)
  a <- sim_accumulate(sim_accumulator(), x[1:123])
  b <- sim_accumulate(sim_accumulator(), x[124:500])
  m <- sim_accumulated(sim_merge(a, b))
  expect_equal(m$m4_x, mean((x - mean(x))^4), tolerance = 1e-10)
})
