test_that("scenario keys depend on content, not order", {
  a <- sim_scenarios(n = c(10L, 20L), mu = c(0, 1))
  b <- sim_scenarios(mu = c(1, 0), n = c(20L, 10L))
  expect_setequal(a$keys, b$keys)
  expect_length(a, 4L)
  expect_error(sim_scenarios_table(data.frame(n = c(1, 1))), "Duplicate")
  expect_error(sim_scenarios(n = c(1, NA)), "missing")
  expect_error(sim_scenarios(replication = 1), "reserved")
  expect_error(sim_scenarios(n = 1:1000, m = 1:1000), "max_scenarios")
})

test_that("constraints separate deliberate exclusions from invalid configurations", {
  sc <- sim_scenarios(n = c(10L, 100L), p = c(2L, 20L),
    .constraints = list(sim_constraint(function(s) s$n < 5 * s$p, "too few observations per predictor")))
  expect_length(sc, 3L)
  expect_equal(nrow(sc$excluded), 1L)
  expect_match(sc$excluded$reason, "too few")
  expect_error(sim_scenarios(rho = c(0.5, 1.5),
    .constraints = list(sim_constraint(function(s) abs(s$rho) > 1, "correlation out of range",
                                       type = "invalid"))), "invalid configuration")
})

test_that("sampled designs are reproducible and record their distribution", {
  smp <- function() list(n = sample(10:100, 1), effect = runif(1))
  a <- sim_scenarios_sample(smp, n = 5, seed = 3, distribution = "n uniform on 10:100, effect U(0, 1)")
  b <- sim_scenarios_sample(smp, n = 5, seed = 3, distribution = "n uniform on 10:100, effect U(0, 1)")
  expect_identical(a$keys, b$keys)
  expect_identical(a$design, "sampled")
  expect_match(a$distribution$description, "uniform")
})

test_that("nested factors are expressed by binding designs", {
  nested <- sim_scenarios_bind(sim_scenarios(design = "parallel", k = 2L),
                               sim_scenarios(design = "crossover", k = c(2L, 4L)))
  expect_length(nested, 3L)
  df <- as.data.frame(nested)
  expect_setequal(df$design, c("parallel", "crossover"))
})

test_that("list valued factors and weights are supported", {
  sc <- sim_scenarios(R = list(diag(2), matrix(c(1, 0.5, 0.5, 1), 2)), .weights = c(1, 3))
  expect_length(sc, 2L)
  expect_equal(sc$weights, c(0.25, 0.75))
  expect_error(sim_scenarios(n = 1:2, .weights = c(-1, 1)), "weights")
})

test_that("validation catches inconsistent declarations and impossible scenarios", {
  st <- sim_study("v", scenarios = sim_scenarios(rho = c(0.5, 0.99, 1.2)),
    generate = sim_generator(function(scenario, context) rnorm(2),
      check = function(scenario) {
        sim_check_matrix(matrix(c(1, scenario$rho, scenario$rho, 1), 2))
      }),
    estimands = list(sim_estimand("m", sim_truth_exact(function(scenario) 0))),
    methods = list(a = sim_method(function(data) sim_estimate("m", 1), scenario_fields = "nope")))
  v <- sim_validate(st, error = FALSE)
  expect_true(any(grepl("unknown scenario field", v$message)))
  expect_true(any(grepl("invalid configuration", v$message)))
  expect_error(sim_validate(st), "validation failed")
  expect_error(sim_estimand("x", type = "vector"), "only scalar")
  expect_error(sim_study("dup", scenarios = sim_scenarios(n = 1L),
    generate = function(scenario, context) 1,
    estimands = list(sim_estimand("a"), sim_estimand("a")),
    methods = list(m = function(data) 1)), "unique")
})

test_that("completeness, freezing and amendments are recorded", {
  st <- normal_study(n = 10L, mu = 0)
  cmp <- sim_completeness(st)
  expect_true(any(cmp$status == "missing"))
  expect_identical(cmp$status[cmp$item == "aims"], "documented")
  fr <- sim_freeze(st)
  expect_identical(fr$status$label, "production")
  expect_identical(fr$status$frozen$hash, sim_fingerprint(st))
  am <- sim_amend(fr, alpha = 0.01, reason = "agreed stricter level")
  expect_length(am$amendments, 1L)
  expect_true(am$amendments[[1]]$after_freeze)
  expect_false(identical(am$amendments[[1]]$fingerprint_before, am$amendments[[1]]$fingerprint_after))
  expect_error(sim_amend(st, alpha = 0.01), "reason")
  # a frozen study edited directly is refused; an amended one runs
  edited <- fr
  edited$aims <- "Something else"
  expect_error(sim_run(edited, replications = 2, seed = 1, path = run_path()),
               "changed after it was frozen")
  expect_s3_class(sim_run(am, replications = 2, seed = 1, path = run_path()), "adempr_run")
})

test_that("scenario weights are part of the study identity", {
  mk <- function(w) sim_study("w", scenarios = sim_scenarios(n = c(10L, 20L), .weights = w),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(avg = function(data) sim_estimate("mean", mean(data))))
  expect_false(identical(sim_fingerprint(mk(c(1, 3))), sim_fingerprint(mk(c(3, 1)))))
  expect_identical(sim_fingerprint(mk(c(1, 3))), sim_fingerprint(mk(c(2, 6))))
  p <- run_path()
  sim_run(mk(c(1, 3)), replications = 2, seed = 1, path = p)
  r <- sim_run(mk(c(3, 1)), replications = 2, seed = 1, path = p)
  ch <- unlist(lapply(r$manifest$revisions, `[[`, "changes"))
  expect_true(any(grepl("scenario design changed", ch)))
})

test_that("the smoke test reports failing methods without stopping", {
  st <- normal_study(n = 10L, mu = 0, methods = list(bad = function(data) stop("boom")))
  v <- sim_validate(st, smoke = TRUE, error = FALSE)
  expect_true(any(v$status == "warning" & grepl("boom", v$message)))
})

test_that("a run refuses more tasks than the budget allows", {
  st <- normal_study(n = 10L, mu = 0)
  expect_error(sim_run(st, replications = 50, seed = 1, path = run_path(), max_tasks = 10),
               "max_tasks")
})

test_that("pilots project cost and keep conservative precision requirements", {
  st <- normal_study(n = 10L, mu = 0, protocol = sim_protocol(precision_targets = c(coverage = 0.01)))
  pl <- sim_pilot(st, replications = 15, seed = 2)
  expect_s3_class(pl, "adempr_pilot")
  req <- pl$precision$required_replications
  expect_true(all(req >= ceiling(0.95 * 0.05 / 0.01^2)))
  expect_equal(pl$projection$planned_replications, max(req))
  expect_identical(pl$run$manifest$label, "pilot")
})
