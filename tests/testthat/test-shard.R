test_that("sharded jobs collect into the same results as one run", {
  st <- normal_study(n = c(10L, 20L, 30L), mu = c(0, 1), methods = stochastic_methods())
  parent <- run_path()
  # each shard runs in its own R process, as an array job would
  skip_if_not_installed("callr")
  lib <- dirname(find.package("adempR"))
  # every job reads the same saved study, as recommended for array jobs
  sf <- tempfile(fileext = ".rds")
  saveRDS(st, sf)
  for (k in 1:3) {
    callr::r(function(k, parent, sf, pkg) {
      if (file.exists(file.path(pkg, "R", "run.R"))) pkgload::load_all(pkg, quiet = TRUE) else library(adempR)
      invisible(sim_run_shard(readRDS(sf), k, 3, replications = 6, seed = 11, path = parent))
    }, args = list(k = k, parent = parent, sf = sf, pkg = find.package("adempR")))
  }
  expect_length(list.dirs(parent, recursive = FALSE), 3L)
  run <- sim_collect_shards(st, 3, replications = 6, seed = 11, path = parent, into = run_path())
  one <- sim_run(st, replications = 6, seed = 11, path = run_path())
  expect_identical(est_signature(run), est_signature(one))
  expect_equal(as.data.frame(sim_summarise(run))$value, as.data.frame(sim_summarise(one))$value)
  expect_setequal(unlist(lapply(1:3, function(k) sim_shard_scenarios(st, k, 3))), st$scenarios$keys)
  expect_error(sim_collect_shards(st, 4, replications = 6, seed = 11, path = parent, into = run_path()),
               "not run yet")
  expect_error(sim_shard_scenarios(st, 4, 3), "between 1")
})

test_that("shards run as batchtools array jobs and collect", {
  skip_if_not_installed("batchtools")
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = stochastic_methods())
  parent <- run_path()
  reg <- batchtools::makeRegistry(file.dir = NA, seed = 1)
  # each job runs in its own Rscript process, as on a scheduler node
  reg$cluster.functions <- batchtools::makeClusterFunctionsInteractive(external = TRUE)
  sf <- tempfile(fileext = ".rds")
  saveRDS(st, sf)
  job <- function(k, parent, sf, pkg) {
    if (file.exists(file.path(pkg, "R", "run.R"))) pkgload::load_all(pkg, quiet = TRUE) else library(adempR)
    sim_run_shard(readRDS(sf), k, 2, replications = 5, seed = 3, path = parent)
    k
  }
  batchtools::batchMap(job, k = 1:2, more.args = list(parent = parent, sf = sf,
                                                  pkg = find.package("adempR")), reg = reg)
  suppressMessages(batchtools::submitJobs(reg = reg))
  batchtools::waitForJobs(reg = reg)
  expect_equal(nrow(batchtools::findDone(reg = reg)), 2L)
  run <- sim_collect_shards(st, 2, replications = 5, seed = 3, path = parent, into = run_path())
  one <- sim_run(st, replications = 5, seed = 3, path = run_path())
  expect_identical(est_signature(run), est_signature(one))
})
