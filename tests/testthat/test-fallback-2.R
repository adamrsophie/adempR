fallback_study <- function(policy = "use", fallback_warns = FALSE) {
  m <- sim_method(function(data) {
    if (data$y[1] < 0) stop("singular fit")
    sim_estimate("mean", mean(data$y))
  }, fallback = function(data) {
    if (fallback_warns) warning("fallback warned")
    Sys.sleep(0.05)
    sim_estimate("mean", median(data$y))
  })
  st <- normal_study(n = 5L, mu = 0, methods = list(m = m),
                     failure = sim_failure_policy(fallback = policy,
                       warnings = if (fallback_warns) "invalid" else "record"))
  st
}

test_that("fallback results keep the flow nonnegative and are counted", {
  r <- sim_run(fallback_study(), replications = 8, seed = 3, path = run_path())
  f <- sim_flow(r)
  e <- sim_results(r)$estimates
  n_fb <- sum(e$fallback)
  expect_gt(n_fb, 0)
  expect_true(all(f$excluded >= 0))
  expect_equal(f$fit_failed, 0L)
  expect_equal(f$fallback_used, n_fb)
  expect_equal(f$planned, f$generation_failed + f$not_generated + f$not_run +
                 f$fit_failed + f$excluded + f$usable)
  expect_equal(attr(sim_summarise(r), "failures")$fallback_used, n_fb)
})

test_that("the fallback is recorded as an attempt with its warnings and time", {
  r <- sim_run(fallback_study(fallback_warns = TRUE), replications = 8, seed = 3,
               path = run_path())
  res <- sim_results(r)
  a <- res$attempts[res$attempts$fallback %in% TRUE, ]
  expect_gt(nrow(a), 0)
  expect_true(all(a$n_warnings == 1L))
  expect_true(all(a$elapsed >= 0.04))
  e <- res$estimates[res$estimates$fallback, ]
  # the warnings policy applies to the output used: the fallback warned
  expect_true(all(e$had_warning))
  expect_false(any(e$valid))
})

test_that("an erroring fallback leaves a failed attempt and a failed fit", {
  m <- sim_method(function(data) stop("primary failed"),
                  fallback = function(data) stop("fallback failed too"))
  st <- normal_study(n = 5L, mu = 0, methods = list(m = m),
                     failure = sim_failure_policy(fallback = "use"))
  r <- sim_run(st, replications = 2, seed = 1, path = run_path())
  a <- sim_results(r)$attempts
  expect_equal(sum(a$fallback & a$status == "error"), 2L)
  expect_equal(sim_flow(r)$fit_failed, 2L)
})

test_that("failing diagnostics and malformed records are kept as errors", {
  m <- sim_method(function(data) sim_estimate("mean", mean(data$y)),
                  diagnostics = function(fit) stop("no iterations stored"))
  st <- normal_study(n = 5L, mu = 0, methods = list(m = m))
  st$generator$record <- function(data) c(a = "x")
  r <- sim_run(st, replications = 2, seed = 1, path = run_path())
  res <- sim_results(r)
  expect_equal(unique(res$diagnostics$name), ".diagnostics_error")
  expect_match(res$diagnostics$message[1], "no iterations stored")
  expect_equal(unique(res$records$name), ".record_error")
  expect_match(res$records$message[1], "must be numbers")
  d <- sim_diagnose(r)
  expect_equal(nrow(d$callback_errors), 4L)
})
