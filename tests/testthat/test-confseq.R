test_that("confidence sequences are wider than fixed intervals and cover the truth uniformly", {
  st <- normal_study(n = 15L, mu = 0, measures = c("coverage", "bias"))
  run <- sim_run(st, replications = 60, seed = 2, path = run_path())
  cs <- sim_confidence_sequence(run, "coverage")
  s <- as.data.frame(sim_summarise(run))
  ci <- s[s$measure == "coverage", ]
  expect_equal(cs$estimate, ci$value)
  expect_true(cs$lower <= ci$lower && cs$upper >= ci$upper)
  expect_equal(cs$type, "exact")
  p <- sim_confidence_sequence(run, "bias", path = TRUE)
  expect_equal(nrow(p), 60L)
  expect_equal(p$type[1], "asymptotic")
  # simulation check of time-uniform validity for Bernoulli(0.95) streams
  set.seed(1)
  bern_cs <- function(x) {
    t <- seq_along(x); m <- cumsum(x) / t
    h <- adempR:::robbins_halfwidth(t, 0.05, 200)
    any(0.95 < m - h | 0.95 > m + h)
  }
  miss <- mean(replicate(400, bern_cs(rbinom(500, 1, 0.95))))
  expect_lt(miss, 0.05)
  expect_error(sim_confidence_sequence(run, alpha = 2), "alpha")
})
