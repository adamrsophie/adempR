prod_run <- function() {
  st <- normal_study(n = 5L, mu = 0)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path, label = "production")
  list(st = st, path = path)
}

test_that("a production run cannot be relabeled without confirmation", {
  x <- prod_run()
  expect_error(sim_run(x$st, replications = 2, seed = 1, path = x$path, label = "draft"),
               "relabeling")
  r <- suppressWarnings(sim_run(x$st, replications = 2, seed = 1, path = x$path, label = "draft",
                                update = "production"))
  expect_equal(r$manifest$label, "draft")
  expect_match(unlist(lapply(r$manifest$revisions, `[[`, "changes")), "label changed",
               all = FALSE)
})

test_that("method names never decide whether a change alters results", {
  x <- prod_run()
  st2 <- x$st
  st2$methods$removed_outliers <- sim_method(function(data) sim_estimate("mean", mean(data$y)))
  expect_no_error(suppressWarnings(sim_run(st2, replications = 2, seed = 1, path = x$path)))
})

test_that("a retry on a production run needs confirmation", {
  flag <- tempfile(); writeLines("x", flag)
  m <- list(m = function(data) {
    if (file.exists(flag)) stop("license")
    sim_estimate("mean", 1)
  })
  st <- normal_study(n = 5L, mu = 0, methods = m)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path, label = "production")
  unlink(flag)
  expect_error(sim_run(st, replications = 2, seed = 1, path = path, retry_failed = TRUE),
               "retry supersedes")
  expect_no_error(suppressWarnings(sim_run(st, replications = 2, seed = 1, path = path,
                                           retry_failed = TRUE, update = "production")))
})

test_that("raised scenario counts and protocol text changes are revisions", {
  st <- normal_study(n = c(5L, 10L), mu = 0)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  r <- sim_run(st, replications = 2, seed = 1, path = path,
               scenario_replications = c(`n=5, mu=0, sigma=1` = 4))
  ch <- unlist(lapply(r$manifest$revisions, `[[`, "changes"))
  expect_match(ch, "replications raised", all = FALSE)
  st$aims <- "A different aim"
  r <- sim_run(st, replications = 2, seed = 1, path = path)
  ch <- unlist(lapply(r$manifest$revisions, `[[`, "changes"))
  expect_match(ch, "protocol text", all = FALSE)
})

test_that("calibration and pilot work refuse a directory on other streams", {
  st <- normal_study(n = 5L, mu = 0)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  expect_error(sim_calibrate_method(st, function(k) function(data) {
    sim_estimate("mean", mean(data$y) * k)
  }, grid = list(k = c(1, 2)), replications = 2, seed = 1, path = path),
  "needs its own directory")
  expect_error(sim_run(st, replications = 2, seed = 1, path = path, label = "pilot"),
               "needs its own directory")
})
