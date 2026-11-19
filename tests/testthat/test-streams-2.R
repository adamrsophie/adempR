test_that("smoke validation never draws the production data", {
  seen <- new.env()
  st <- sim_study("streams", scenarios = sim_scenarios(n = 3L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(m = function(data) sim_estimate("mean", mean(data$y))))
  path <- run_path()
  r <- sim_run(st, replications = 1, seed = 5, path = path)
  prod <- sim_results(r)$estimates$estimate
  v <- sim_validate(st, smoke = TRUE, seed = 5, error = FALSE)
  # replay the validation task to obtain its estimate
  ctx <- list(seed = 5, domain = "validation",
              registry = rng_registry(rng_keys(st, "validation", "validation_method")),
              cache = new.env(), dgp_stage = "validation", method_stage = "validation_method")
  val <- execute_task(st, ctx, st$scenarios$keys, 1L, "m")$estimates$estimate
  expect_false(isTRUE(all.equal(prod, val)))
})

test_that("smoke validation and DGP checks work with common random numbers and antithetic pairs", {
  st <- sim_study("crn", scenarios = sim_scenarios(n = 20L, shift = c(0, 1)),
    generate = sim_generator(function(scenario, context) {
      z <- rnorm(scenario$n)
      if (isTRUE(context$antithetic)) z <- -z
      data.frame(y = z + scenario$shift)
    }, common = "n", antithetic = TRUE),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$shift))),
    methods = list(m = function(data) sim_estimate("mean", mean(data$y))))
  v <- sim_validate(st, smoke = TRUE, error = FALSE)
  expect_false(any(v$status == "error"))
  set.seed(42); before <- runif(1); set.seed(42)
  d <- sim_check_dgp(st, summaries = function(data) c(m = mean(data$y), r = runif(1)),
                     targets = function(scenario) c(m = scenario$shift), n = 20, seed = 1)
  expect_equal(runif(1), before)
  # antithetic pairs cancel exactly: the pair means equal the shift
  expect_equal(d$mean, c(0, 1))
  expect_equal(d$mcse, c(0, 0))
})
