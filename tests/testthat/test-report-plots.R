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
  expect_error(sim_plot(run, "bland_altman", methods = c("med", "nope")), "two methods")
})
