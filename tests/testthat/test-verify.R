test_that("sim_verify confirms a reproduction and flags a seed sensitive result", {
  st <- normal_study(n = 20L, mu = 0, measures = c("bias", "empse"))
  a <- sim_run(st, replications = 30, seed = 1, path = run_path())
  a2 <- sim_run(st, replications = 30, seed = 1, path = run_path())
  v <- sim_verify(a, a2)
  expect_equal(v$mode, "reproduce")
  expect_true(v$ok)
  expect_equal(v$table$n_identical, v$table$n_compared)
  b <- sim_run(st, replications = 30, seed = 2, path = run_path())
  s <- sim_verify(a, b)
  expect_equal(s$mode, "seed_sensitivity")
  expect_true(all(is.finite(s$table$z)))
  # a large threshold accepts, a tiny threshold flags every cell that differs
  expect_true(sim_verify(a, b, threshold = 100)$ok)
  expect_false(sim_verify(a, b, threshold = 1e-8)$ok)
  expect_output(print(s), "Seed sensitivity")
  other <- normal_study(n = 30L, mu = 0)
  o <- sim_run(other, replications = 5, seed = 1, path = run_path())
  expect_error(sim_verify(a, o), "different studies")
})

test_that("reproduction compares complete key sets, all outcomes and infinite values", {
  st <- normal_study(n = 10L, mu = 0)
  a <- sim_run(st, replications = 4, seed = 1, path = run_path())
  # 1 and 1L are the same seed
  b <- sim_run(st, replications = 4, seed = 1L, path = run_path())
  expect_equal(sim_verify(a, b)$mode, "reproduce")
  # rows in one run only fail the check
  short <- sim_run(st, replications = 2, seed = 1, path = run_path())
  v <- sim_verify(a, short)
  expect_false(v$ok)
  expect_equal(v$table$only_a, 2L)
  # equal infinite limits are equal
  inf <- normal_study(n = 10L, mu = 0, methods = list(
    m = function(data) sim_estimate("mean", mean(data$y), lower = -Inf, upper = Inf)))
  v <- sim_verify(sim_run(inf, replications = 2, seed = 1, path = run_path()),
                  sim_run(inf, replications = 2, seed = 1, path = run_path()))
  expect_true(v$ok)
  expect_equal(v$table$max_abs_diff, 0)
  # a changed test decision is a mismatch, even with the same estimate
  flag <- tempfile()
  dec <- normal_study(n = 10L, mu = 0, methods = list(
    m = function(data) sim_estimate("mean", 0, p_value = if (file.exists(flag)) 0.01 else 0.5)))
  r1 <- sim_run(dec, replications = 2, seed = 1, path = run_path())
  writeLines("x", flag)
  r2 <- sim_run(dec, replications = 2, seed = 1, path = run_path())
  v <- sim_verify(r1, r2)
  expect_false(v$ok)
  expect_equal(nrow(v$mismatches), 2L)
})

test_that("seed sensitivity treats boundary proportions and missing MCSEs carefully", {
  wide <- normal_study(n = 10L, mu = 0, measures = "coverage", methods = list(
    m = function(data) sim_estimate("mean", mean(data$y), lower = -100, upper = 100)))
  a <- sim_run(wide, replications = 20, seed = 1, path = run_path())
  b <- sim_run(wide, replications = 20, seed = 2, path = run_path())
  v <- sim_verify(a, b)
  # equal coverage of 1 in both: z = 0, not an indeterminate 0/0
  expect_equal(v$table$z[v$table$measure == "coverage"], 0)
  expect_true(v$ok)
  expect_equal(v$expected_flags, sum(!is.na(v$table$z)) * 2 * stats::pnorm(-3))
  # an estimand without any usable estimate gives no MCSE: indeterminate
  none <- normal_study(n = 10L, mu = 0, measures = "bias", methods = list(
    m = function(data) sim_estimate("mean", NA_real_)))
  v <- sim_verify(sim_run(none, replications = 3, seed = 1, path = run_path()),
                  sim_run(none, replications = 3, seed = 2, path = run_path()))
  expect_true(is.na(v$ok))
  expect_equal(v$n_indeterminate, 1L)
})

test_that("reproduction checks every scientific table, one mutated field at a time", {
  st <- sim_study("verify_schema", scenarios = sim_scenarios(n = 5L),
    generate = sim_generator(function(scenario, context) stats::rnorm(scenario$n),
                             record = function(data) list(m = mean(data))),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(a = sim_method(function(data) sim_estimate("mean", mean(data)),
                                  diagnostics = function(fit) c(d = 1))),
    measures = "bias")
  ra <- sim_results(sim_run(st, 3, 1, run_path(), progress = FALSE))
  rb <- sim_results(sim_run(st, 3, 1, run_path(), progress = FALSE))
  v <- sim_verify(ra, rb)
  expect_true(v$ok)
  expect_true(all(v$tables$identical))
  expect_true(all(c("diagnostics", "records", "attempts", "truth") %in% v$tables$table))
  expect_true("elapsed" %in% v$excluded_fields)
  mutate <- function(tab, col, val) {
    x <- rb
    x[[tab]][[col]][1] <- val
    sim_verify(ra, x)
  }
  for (m in list(c("diagnostics", "value"), c("records", "value"), c("attempts", "status"),
                 c("truth", "truth"), c("estimates", "estimate"), c("replications", "data_md5"))) {
    val <- if (is.numeric(rb[[m[1]]][[m[2]]])) 999 else "changed"
    w <- mutate(m[1], m[2], val)
    expect_false(w$ok, info = paste(m, collapse = "$"))
    expect_false(w$tables$identical[w$tables$table == m[1]], info = m[1])
  }
  # timing fields are not compared
  x <- rb
  x$attempts$elapsed <- x$attempts$elapsed + 10
  expect_true(sim_verify(ra, x)$ok)
})
