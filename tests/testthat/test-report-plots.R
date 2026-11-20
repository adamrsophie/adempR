test_that("plots leave the caller's random stream and devices alone", {
  st <- normal_study(n = 10L, mu = 0)
  run <- sim_run(st, replications = 8, seed = 3, path = run_path())
  set.seed(42)
  before <- .Random.seed
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  mine <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(mine), add = TRUE)
  sim_plot(run, "estimates")
  expect_identical(.Random.seed, before)
  out <- tempfile(fileext = ".html")
  sim_report(run, out)
  expect_identical(grDevices::dev.cur(), mine)
})

test_that("the zip plot uses coverage eligibility and has a power version", {
  st <- normal_study(n = 10L, mu = 0)
  run <- sim_run(st, replications = 8, seed = 3, path = run_path())
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(sim_plot(run, "zip"))
  expect_null(sim_plot(run, "zip", reference = "null"))
  expect_null(sim_plot(run, "precision", measure = "coverage"))
})

test_that("Bland-Altman and interval-by-truth plots draw and validate their inputs", {
  st <- normal_study(n = 20L, mu = 0, methods = c(normal_study()$methods,
    list(med = function(data) sim_estimate("mean", median(data$y), lower = median(data$y) - 0.5,
                                           upper = median(data$y) + 0.5, level = 0.95))))
  run <- sim_run(st, replications = 15, seed = 2, path = run_path())
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  dev <- grDevices::dev.cur()
  expect_null(sim_plot(run, "bland_altman"))
  expect_null(sim_plot(run, "interval_truth"))
  expect_identical(grDevices::dev.cur(), dev)
  grDevices::dev.off()
  expect_error(sim_plot(run, "bland_altman", methods = c("med", "nope")), "two different methods")
  expect_error(sim_plot(run, "bland_altman", methods = c("med", "med")), "two different methods")
})

test_that("unbounded intervals are drawn and the nominal coverage line uses eligible intervals", {
  st <- normal_study(n = 10L, mu = 0, methods = list(
    wide = function(data) sim_estimate("mean", mean(data$y), lower = -Inf, upper = Inf, level = 0.95),
    # a 0.5 level on intervals that are not eligible must not count
    half = function(data) sim_estimate("mean", mean(data$y), lower = 1, upper = 0, level = 0.5)),
    measures = "coverage")
  run <- sim_run(st, replications = 5, seed = 2, path = run_path())
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_null(sim_plot(run, "interval_truth"))
  e <- sim_results(run)$estimates
  ok <- e$interval_status == "valid" & e$estimate_status %in% adempR:::USABLE_FOR_INTERVAL
  expect_equal(unique(e$level[ok]), 0.95)
  expect_null(sim_plot(run, "precision", measure = "coverage"))
})
