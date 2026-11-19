# Scientific tests 7 and 8: estimand compatibility and numerical truth.

logistic_study <- function(truth_n = 2e5) {
  b0 <- -1; bA <- 1; bX <- 2
  sim_study("marginal_vs_conditional",
    scenarios = sim_scenarios(n = 400L, b0 = b0, bA = bA, bX = bX),
    generate = function(scenario, context) {
      x <- rnorm(scenario$n)
      a <- rbinom(scenario$n, 1, 0.5)
      y <- rbinom(scenario$n, 1, plogis(scenario$b0 + scenario$bA * a + scenario$bX * x))
      data.frame(y = y, a = a, x = x)
    },
    estimands = list(
      sim_estimand("marginal_lor", effect = "marginal", scale = "log_odds_ratio",
        population = "all randomized", conditioning = "none",
        truth = sim_truth_mc(function(scenario, n) {
          x <- rnorm(n)
          p1 <- plogis(scenario$b0 + scenario$bA + scenario$bX * x)
          p0 <- plogis(scenario$b0 + scenario$bX * x)
          m1 <- mean(p1); m0 <- mean(p0)
          val <- qlogis(m1) - qlogis(m0)
          g <- cbind(p1 / (m1 * (1 - m1)), -p0 / (m0 * (1 - m0)))
          se <- sd(drop(g %*% c(1, 1))) / sqrt(n)
          sim_truth_value(val, se, "delta method")
        }, n = truth_n)),
      sim_estimand("conditional_lor", effect = "conditional", scale = "log_odds_ratio",
        population = "all randomized", conditioning = "x",
        truth = sim_truth_exact(function(scenario) scenario$bA))),
    methods = list(
      unadjusted = sim_method(function(data) {
        f <- glm(y ~ a, binomial, data = data)
        sim_estimate("marginal_lor", unname(coef(f)["a"]),
                     se = sqrt(vcov(f)["a", "a"]), effect = "marginal")
      }, estimands = "marginal_lor"),
      adjusted = sim_method(function(data) {
        f <- glm(y ~ a + x, binomial, data = data)
        sim_estimate("conditional_lor", unname(coef(f)["a"]), se = sqrt(vcov(f)["a", "a"]))
      }, estimands = "conditional_lor", declares = list(conditional_lor = list(effect = "conditional")))),
    measures = c("bias", "empse"))
}

test_that("a conditional estimate cannot be assessed against a marginal estimand", {
  st <- logistic_study()
  bad <- st$methods
  bad$adjusted$estimands <- "marginal_lor"
  bad$adjusted$declares <- list(marginal_lor = list(effect = "conditional"))
  expect_error(sim_amend(st, methods = bad, reason = "wrong pairing") |> sim_validate(),
               "validation failed")
  v <- sim_validate(sim_amend(st, methods = bad, reason = "wrong pairing"), error = FALSE)
  expect_true(any(grepl("does not estimate this estimand", v$message)))
  # runtime declaration is also checked
  rt <- st$methods
  rt$unadjusted <- sim_method(function(data) {
    sim_estimate("marginal_lor", 0.5, effect = "conditional")
  }, estimands = "marginal_lor")
  run <- sim_run(sim_amend(st, methods = rt, reason = "runtime mismatch"), replications = 2,
                 seed = 1, path = run_path())
  e <- sim_results(run)$estimates
  expect_true(all(e$estimate_status[e$method == "unadjusted"] == "incompatible"))
})

test_that("the generating coefficient is the wrong truth for the marginal log odds ratio", {
  st <- logistic_study()
  run <- sim_run(st, replications = 300, seed = 2026, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  tr <- sim_truth(run)
  marg <- tr$truth[tr$estimand == "marginal_lor"]
  # noncollapsibility: the marginal log odds ratio is attenuated
  expect_lt(marg, 1 - 0.2)
  expect_gt(tr$truth_se[tr$estimand == "marginal_lor"], 0)
  b_marg <- s[s$method == "unadjusted" & s$measure == "bias", ]
  expect_lt(abs(b_marg$value), 4 * b_marg$mcse_total)
  # the same estimates judged against the conditional coefficient look badly biased
  est <- sim_results(run)$estimates
  un <- est$estimate[est$method == "unadjusted"]
  wrong_bias <- mean(un) - 1
  expect_gt(abs(wrong_bias), 10 * b_marg$mcse)
  b_cond <- s[s$method == "adjusted" & s$measure == "bias", ]
  expect_lt(abs(b_cond$value), 4 * b_cond$mcse + 0.05)
})

test_that("Monte Carlo truth uncertainty does not shrink with more outer replications", {
  st <- logistic_study(truth_n = 5e4)
  p <- run_path()
  r1 <- sim_run(st, replications = 50, seed = 9, path = p)
  s1 <- as.data.frame(sim_summarise(r1))
  r2 <- sim_run(st, replications = 150, seed = 9, path = p)
  s2 <- as.data.frame(sim_summarise(r2))
  pick <- function(s) s[s$method == "unadjusted" & s$measure == "bias", ]
  expect_equal(pick(s1)$truth_se, pick(s2)$truth_se)
  expect_lt(pick(s2)$mcse, pick(s1)$mcse)
  expect_gt(pick(s2)$mcse_total, pick(s2)$truth_se)
  expect_equal(pick(s2)$mcse_total, sqrt(pick(s2)$mcse^2 + pick(s2)$truth_se^2))
  expect_true(pick(s2)$truth_sens_lower < pick(s2)$value && pick(s2)$value < pick(s2)$truth_sens_upper)
  # cached truth is reused, not recomputed
  expect_length(list.files(file.path(p, "truth"), pattern = "\\.rds$"), 1L)
})

test_that("unknown truth disables truth based measures only", {
  st <- sim_study("unknown",
    scenarios = sim_scenarios(n = 10L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean")),
    methods = list(avg = function(data) sim_estimate("mean", mean(data), se = 1)),
    measures = c("bias", "empse", "mean_se"))
  expect_warning(v <- sim_validate(st, error = FALSE), NA)
  run <- sim_run(st, replications = 10, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  expect_identical(s$status[s$measure == "bias"], "unavailable")
  expect_identical(s$status[s$measure == "empse"], "ok")
})

test_that("data set specific truth comes from the generator and stays hidden from methods", {
  st <- sim_study("dataset_truth",
    scenarios = sim_scenarios(n = 50L),
    generate = function(scenario, context) {
      y0 <- rnorm(scenario$n); y1 <- y0 + rnorm(scenario$n, 1)
      sim_data(data.frame(y = y1 - y0), truth = list(sate = mean(y1 - y0)))
    },
    estimands = list(sim_estimand("sate", truth = sim_truth_dataset(), effect = "finite_population")),
    methods = list(
      avg = function(data, scenario, context) {
        stopifnot(is.null(context$truth), length(scenario) == 0)
        sim_estimate("sate", mean(data$y))
      },
      oracle = sim_method(function(data, context) sim_estimate("sate", context$truth$sate),
                          oracle = TRUE)),
    measures = c("bias", "mse"))
  run <- sim_run(st, replications = 5, seed = 1, path = run_path())
  s <- as.data.frame(sim_summarise(run))
  expect_equal(s$value[s$method == "oracle" & s$measure == "mse"], 0)
  expect_equal(s$value[s$method == "avg" & s$measure == "mse"], 0)
  expect_identical(s$truth_kind[1], "dataset")
})

test_that("a usable interval with missing data set truth is indeterminate, not a miss", {
  st <- sim_study("partial_truth",
    scenarios = sim_scenarios(n = 20L),
    generate = function(scenario, context) {
      y <- rnorm(scenario$n)
      tr <- if (context$replication %% 2 == 0) list(m = 0) else list()
      sim_data(data.frame(y = y), truth = tr)
    },
    estimands = list(sim_estimand("m", truth = sim_truth_dataset())),
    methods = list(wide = function(data) {
      sim_estimate("m", mean(data$y), lower = -10, upper = 10, level = 0.95)
    }),
    measures = c("coverage", "usable_coverage"))
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 4, seed = 1, path = run_path())))
  cov <- s[s$measure == "coverage", ]
  expect_equal(cov$value, 1)
  expect_equal(cov$n_eligible, 2)
  uc <- s[s$measure == "usable_coverage", ]
  expect_identical(uc$status, "indeterminate")
  expect_true(is.na(uc$value))
  expect_lt(uc$lower, 0.5)
  expect_equal(uc$upper, 1)
})

test_that("empirical SE uses the errors when the truth varies between data sets", {
  st <- sim_study("vt", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) {
      t <- context$replication * 100
      sim_data(data.frame(y = t + c(-1, 1)[1 + context$replication %% 2]), truth = list(m = t))
    },
    estimands = list(sim_estimand("m", truth = sim_truth_dataset())),
    methods = list(a = function(data) sim_estimate("m", data$y[1], se = 1)),
    measures = c("empse", "se_ratio"))
  s <- as.data.frame(sim_summarise(sim_run(st, replications = 6, seed = 1, path = run_path())))
  expect_equal(s$value[s$measure == "empse"], sd(rep(c(1, -1), 3)))
  expect_match(s$note[s$measure == "empse"], "errors")
})

test_that("truth providers report their real precision", {
  expect_error(sim_truth_numeric(function(scenario) 1, error = -1), "nonnegative")
  tr <- sim_truth_mc(function(scenario, n) c(1, 3), n = 100)
  row <- adempR:::compute_truth_one(tr, list(), seed = adempR:::rng_base_seed(1))
  expect_equal(row$truth_n, 2)
  expect_match(row$truth_method, "2 contributions returned; 100 requested")
  bad <- sim_truth_numeric(function(scenario) list(value = 1, error = -1))
  expect_error(adempR:::compute_truth_one(bad, list()), "nonnegative")
})
