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
