# Scientific tests 9 and 11.

test_that("inner resampling does not inflate the number of independent replications", {
  st <- sim_study("nested",
    scenarios = sim_scenarios(n = 30L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      boot = function(data) {
        b <- vapply(1:50, function(i) mean(sample(data, replace = TRUE)), numeric(1))
        sim_estimate("mean", mean(data), se = sd(b))
      },
      draws = function(data) {
        # repeated outputs for one estimand must be combined by the method
        lapply(1:3, function(i) sim_estimate("mean", mean(data)))
      }))
  run <- sim_run(st, replications = 40, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  b <- s[s$method == "boot" & s$measure == "bias", ]
  expect_equal(b$n_eligible, 40L)
  e <- sim_results(run)$estimates
  expect_equal(nrow(e[e$method == "boot", ]), 40L)
  expect_true(all(e$estimate_status[e$method == "draws"] == "extraction_error"))
  att <- sim_results(run)$attempts
  expect_match(att$message[att$method == "draws"][1], "outer replication is the unit")
})

test_that("alpha spending keeps simultaneous coverage across repeated looks", {
  looks <- c(50, 100, 200, 400)
  p <- 0.9
  alpha <- 0.05
  conf <- 1 - alpha / length(looks)
  set.seed(123)
  nsim <- 4000
  miss_adj <- miss_naive <- logical(nsim)
  for (i in seq_len(nsim)) {
    x <- rbinom(max(looks), 1, p)
    k <- cumsum(x)[looks]
    adj <- adempR:::cp_interval(k, looks, conf)
    naive <- adempR:::cp_interval(k, looks, 0.95)
    miss_adj[i] <- any(adj[, "lower"] > p | adj[, "upper"] < p)
    miss_naive[i] <- any(naive[, "lower"] > p | naive[, "upper"] < p)
  }
  tol <- 3 * sqrt(alpha * (1 - alpha) / nsim)
  expect_lte(mean(miss_adj), alpha + tol)
  # repeated looks at fixed sample intervals inflate the error
  expect_gt(mean(miss_naive), alpha)
})

test_that("sequential runs stop on committed prefixes and refuse unsupported measures", {
  st <- normal_study(n = 20L, mu = 0, measures = c("coverage"))
  p <- run_path()
  sq <- sim_run_sequential(st, looks = c(40, 80, 400), half_width = 0.08, seed = 1, path = p)
  expect_true(sq$stopped_at %in% 1:3)
  last <- sq$looks[sq$looks$look == sq$stopped_at, ]
  expect_true(all(last$half_width <= 0.08))
  expect_equal(sq$run$manifest$replications, unique(last$replications))
  expect_error(sim_run_sequential(st, looks = c(40, 80), half_width = 0.05, seed = 1,
                                  path = run_path(), measures = "bias"), "Bernoulli")
})
