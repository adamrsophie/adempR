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

test_that("sequential designs are stored, analyze only the prefix and cannot be changed", {
  st <- normal_study(n = c(10L, 20L), mu = 0, measures = c("coverage"))
  # a directory started by sim_run() holds replications beyond the first look
  p0 <- run_path()
  sim_run(st, replications = 12, seed = 1, path = p0)
  expect_error(sim_run_sequential(st, looks = c(4, 8), half_width = 0.45, seed = 1, path = p0),
               "not started by sim_run_sequential")
  # unattainable target and invalid alpha are refused before any work
  expect_error(sim_run_sequential(st, looks = c(4, 8), half_width = 0.05, seed = 1,
                                  path = run_path()), "cannot be reached")
  expect_error(sim_run_sequential(st, looks = c(4, 8), half_width = 0.45, seed = 1,
                                  path = run_path(), alpha = 1.5), "alpha")
  p <- run_path()
  sq <- sim_run_sequential(st, looks = c(30, 60, 400), half_width = 0.2, seed = 1, path = p)
  expect_false(is.na(sq$stopped_at))
  n_stop <- sq$looks$replications[sq$looks$look == sq$stopped_at][1]
  run <- sim_open(p)
  expect_identical(run$manifest$sequential$analyzed, n_stop)
  expect_identical(run$manifest$sequential$design$looks, c(30L, 60L, 400L))
  # extending the directory with sim_run() does not change what is analyzed
  sim_run(st, replications = n_stop + 10L, seed = 1, path = p)
  s <- as.data.frame(sim_summarise(sim_open(p)))
  expect_true(all(s$n_generated <= n_stop))
  expect_match(paste(capture.output(print(sim_summarise(sim_open(p)))), collapse = " "),
               "Sequential run")
  f <- tempfile(fileext = ".html")
  sim_report(sim_open(p), f, plots = FALSE)
  expect_match(paste(readLines(f), collapse = " "), "Sequential stopping")
  # rerunning with different looks is refused; the same design resumes
  expect_error(sim_run_sequential(st, looks = c(30, 90, 400), half_width = 0.2, seed = 1, path = p),
               "different sequential design")
  sq2 <- sim_run_sequential(st, looks = c(30, 60, 400), half_width = 0.2, seed = 1, path = p)
  expect_identical(sq2$stopped_at, sq$stopped_at)
})

test_that("the sequential runner splits alpha over looks and cells on prefix counts", {
  # two scenarios, two methods: four monitored cells; intervals checked
  # against binom.test() on coverage counts recomputed from the prefix
  m <- normal_study()$methods
  m2 <- list(t_a = m$t_interval, t_b = m$t_interval)
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = m2, measures = "coverage")
  looks <- c(20, 40)
  sq <- sim_run_sequential(st, looks = looks, half_width = 0.15, seed = 5, path = run_path(),
                           alpha = 0.1)
  expect_equal(sq$cells, 4L)
  conf <- 1 - 0.1 / (2 * 4)
  expect_equal(sq$conf_level_per_interval, conf)
  e <- sim_results(sq$run)$estimates
  lk <- sq$looks
  for (i in seq_len(nrow(lk))) {
    key <- st$scenarios$keys[match(lk$scenario_id[i], st$scenarios$labels)]
    r <- e[e$scenario_key == key & e$method == lk$method[i] &
             e$replication <= lk$replications[i], ]
    k <- sum(r$lower <= r$truth & r$truth <= r$upper)
    bt <- stats::binom.test(k, nrow(r), conf.level = conf)$conf.int
    expect_equal(lk$n_eligible[i], nrow(r))
    expect_equal(c(lk$lower[i], lk$upper[i]), as.numeric(bt), tolerance = 1e-10)
  }
})
