# Operational tests 16, 17 and 18: crashes, duplicates, corruption, locks.

crash_study <- function(flag = "") {
  sim_study("crash", scenarios = sim_scenarios(n = c(10L, 20L)),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(
      avg = function(data) sim_estimate("mean", mean(data)),
      dies = sim_method(function(data, context) {
        if (nzchar(flag) && context$replication == 13 && file.exists(flag)) {
          file.remove(flag)
          tools::pskill(Sys.getpid(), tools::SIGKILL)
        }
        sim_estimate("mean", median(data))
      })))
}

crash_and_resume <- function(fault, flag = "") {
  p <- run_path()
  st <- crash_study(flag)
  res <- tryCatch(callr::r(function(st, p, fault) {
    library(adempR)
    options(adempR.fault = fault)
    sim_run(st, replications = 30, seed = 8, path = p, chunk_size = 5)
  }, list(st, p, fault), timeout = 120), error = function(e) e)
  expect_s3_class(res, "error")
  expect_true(dir.exists(file.path(p, "LOCK")))
  # resume in this process: the stale lock of the dead process is recovered
  expect_message(run <- suppressWarnings(sim_run(st, replications = 30, seed = 8, path = p)),
                 "stale lock")
  run
}

test_that("killed runs resume without double counting", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("callr")
  skip_if_dev_namespace()
  ref <- sim_run(crash_study(), replications = 30, seed = 8, path = run_path(), chunk_size = 5)
  ref_sig <- est_signature(ref)
  rownames(ref_sig) <- NULL
  for (fault in c("after_data_write", "before_sidecar", "after_commit")) {
    run <- crash_and_resume(fault)
    sig <- est_signature(run)
    rownames(sig) <- NULL
    expect_identical(sig, ref_sig, info = fault)
    expect_equal(nrow(sim_results(run)$estimates), 2 * 30 * 2)
  }
  flag <- tempfile()
  file.create(flag)
  run <- crash_and_resume("none", flag)
  sig <- est_signature(run)
  rownames(sig) <- NULL
  expect_identical(sig, ref_sig)
})

test_that("duplicate deliveries are deduplicated and conflicts are detected", {
  st <- normal_study(n = 10L, mu = 0)
  p <- run_path()
  run <- sim_run(st, replications = 6, seed = 2, path = p, chunk_size = 3)
  b <- list.files(file.path(p, "batches"), pattern = "\\.rds$", full.names = TRUE)[1]
  copy <- file.path(p, "batches", "b-duplicate.rds")
  file.copy(b, copy)
  writeLines(unname(tools::md5sum(copy)), sub("\\.rds$", ".md5", copy))
  expect_equal(nrow(sim_results(run)$estimates), 6)
  # timing differences are not conflicts
  obj <- readRDS(copy)
  obj$attempts$elapsed <- obj$attempts$elapsed + 1
  saveRDS(obj, copy)
  writeLines(unname(tools::md5sum(copy)), sub("\\.rds$", ".md5", copy))
  expect_equal(nrow(sim_results(run)$estimates), 6)
  # scientific differences are
  obj$estimates$estimate[1] <- obj$estimates$estimate[1] + 1
  saveRDS(obj, copy)
  writeLines(unname(tools::md5sum(copy)), sub("\\.rds$", ".md5", copy))
  expect_error(sim_results(run), "Conflicting duplicate")
})

test_that("truncated batches are detected and their tasks rerun", {
  st <- normal_study(n = 10L, mu = 0)
  p <- run_path()
  run <- sim_run(st, replications = 9, seed = 2, path = p, chunk_size = 3)
  ref <- est_signature(run)
  b <- list.files(file.path(p, "batches"), pattern = "\\.rds$", full.names = TRUE)[2]
  bytes <- readBin(b, "raw", file.info(b)$size)
  writeBin(bytes[seq_len(length(bytes) %/% 2)], b)
  expect_warning(res <- sim_results(run), "corrupt")
  expect_equal(nrow(res$estimates), 6)
  expect_warning(run2 <- sim_run(st, replications = 9, seed = 2, path = p), "quarantine")
  expect_identical(est_signature(run2), ref)
  expect_length(list.files(file.path(p, "quarantine")), 2L)
})

test_that("a live lock blocks a second coordinator; incompatible directories are refused", {
  st <- normal_study(n = 10L, mu = 0)
  p <- run_path()
  sim_run(st, replications = 2, seed = 1, path = p)
  dir.create(file.path(p, "LOCK"))
  jsonlite::write_json(list(pid = Sys.getpid(), host = Sys.info()[["nodename"]], time = "now"),
                       file.path(p, "LOCK", "owner.json"), auto_unbox = TRUE)
  expect_error(sim_run(st, replications = 2, seed = 1, path = p), "locked")
  if (.Platform$OS.type == "unix") {
    expect_error(sim_unlock(p), "still running")
  }
  sim_unlock(p, force = TRUE)
  expect_false(dir.exists(file.path(p, "LOCK")))
  expect_error(sim_run(st, replications = 2, seed = 99, path = p), "seed")
  other <- run_path()
  dir.create(other)
  writeLines("x", file.path(other, "notes.txt"))
  expect_error(sim_run(st, replications = 2, seed = 1, path = other), "not an adempR run")
})

test_that("an unavailable optional backend fails clearly", {
  skip_if(requireNamespace("future.apply", quietly = TRUE), "future.apply is installed")
  expect_error(sim_run(normal_study(n = 10L, mu = 0), 2, seed = 1, path = run_path(),
                       backend = "future"), "future")
})

test_that("a stale lock is taken over once, and a fresh lock is never stolen", {
  skip_on_os("windows")
  p <- adempR:::run_paths(run_path())
  dir.create(p$root)
  dead <- list(pid = 999999L, host = Sys.info()[["nodename"]], time = "t0")
  skip_if(isTRUE(adempR:::pid_alive(dead$pid)), "pid 999999 is in use")
  dir.create(p$lock)
  jsonlite::write_json(dead, file.path(p$lock, "owner.json"), auto_unbox = TRUE)
  # a second coordinator that saw the same dead owner after the first took over
  tok <- suppressMessages(adempR:::acquire_lock(p))
  expect_false(adempR:::take_over_stale_lock(p, dead))
  expect_identical(adempR:::lock_owner(p)$token, tok)
  expect_length(list.files(p$root, pattern = "stale"), 0L)
  adempR:::release_lock(p, "not-mine")
  expect_true(dir.exists(p$lock))
  adempR:::release_lock(p, tok)
  expect_false(dir.exists(p$lock))
})
