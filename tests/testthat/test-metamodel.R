test_that("factor influence attributes variation to the factor that drives it", {
  st <- sim_study("fi", scenarios = sim_scenarios(n = c(10L, 40L), sd = c(1, 4), shift = c(0, 1)),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, scenario$shift, scenario$sd)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$shift))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))),
    measures = "empse")
  run <- sim_run(st, replications = 60, seed = 1, path = run_path())
  fi <- sim_factor_influence(run, "empse")
  sh <- stats::setNames(fi$anova$share, fi$anova$term)
  expect_gt(sh[["sd"]], sh[["shift"]])
  expect_equal(sum(fi$anova$share, na.rm = TRUE), 1)
  expect_equal(nrow(fi$cells), 8L)
  expect_output(print(fi), "does not predict")
  fi2 <- sim_factor_influence(run, "empse", factors = c("n", "sd"), interactions = TRUE)
  expect_true("n:sd" %in% fi2$anova$term)
  expect_error(sim_factor_influence(run, "empse", factors = "nope"), "factors")
})

test_that("factor influence does not depend on factor order, names or the weight column", {
  sc <- sim_scenarios(n = c(10L, 40L), `dose level` = c(1, 4), shift = c(0, 1),
                      .constraints = list(sim_constraint(function(s) s$n == 40L && s$shift == 1, "not run")),
                      .weights = c(1, 2, 1, 2, 1, 2))
  st <- sim_study("fi2", scenarios = sc,
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n, scenario$shift,
                                                                 scenario[["dose level"]])),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$shift))),
    methods = list(mean = function(data) sim_estimate("mean", mean(data$y))),
    measures = "empse")
  run <- sim_run(st, replications = 40, seed = 2, path = run_path())
  fi <- sim_factor_influence(run, "empse")
  expect_setequal(fi$anova$term, c("n", "dose level", "shift", "Residuals"))
  rev <- sim_factor_influence(run, "empse", factors = c("shift", "dose level", "n"))
  share <- function(x) stats::setNames(x$anova$share, x$anova$term)[c("n", "dose level", "shift")]
  expect_equal(share(rev), share(fi))
  # each term last: the drop in weighted RSS when the term leaves the model
  d <- fi$cells
  w <- 1 / d$mcse^2
  full <- lm(value ~ factor(n) + factor(`dose level`) + factor(shift), d, weights = w)
  no_n <- lm(value ~ factor(`dose level`) + factor(shift), d, weights = w)
  expect_equal(fi$anova$sum_sq[fi$anova$term == "n"], deviance(no_n) - deviance(full))
  expect_identical(fi$weighting, "inverse squared MCSE")
  expect_output(print(fi), "type II")
  inter <- sim_factor_influence(run, "empse", factors = c("n", "shift"), interactions = TRUE)
  # the constrained design has no cell with n = 40 and shift = 1
  expect_equal(inter$anova$df[inter$anova$term == "n:shift"], 0)
})
