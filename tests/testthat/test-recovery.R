# Operational tests 16, 17 and 18: crashes, duplicates, corruption, locks.

crash_study <- function(flag = "", gen_flag = "") {
  sim_study("crash", scenarios = sim_scenarios(n = c(10L, 20L)),
    generate = function(scenario, context) {
      # a kill in the middle of generating a data set
      if (nzchar(gen_flag) && context$replication == 7 && file.exists(gen_flag)) {
        file.remove(gen_flag)
        tools::pskill(Sys.getpid(), tools::SIGKILL)
      }
      rnorm(scenario$n)
    },
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

crash_and_resume <- function(fault, flag = "", gen_flag = "") {
  p <- run_path()
  st <- crash_study(flag, gen_flag)
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
  gflag <- tempfile()
  file.create(gflag)
  run <- crash_and_resume("none", gen_flag = gflag)
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
  local_mocked_bindings(has_namespace = function(pkg) pkg != "future.apply")
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

test_that("a damaged manifest falls back to its backup and a bare directory restarts", {
  st <- normal_study(n = 10L, mu = 0)
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  sim_run(st, replications = 3, seed = 1, path = path)
  f <- file.path(path, "manifest.rds")
  writeBin(as.raw(1:10), f)
  expect_warning(man <- read_manifest(path), "checksum")
  expect_equal(man$replications, 3L)
  # the next write must not copy the damaged primary over the good backup
  good_bak <- readBin(paste0(f, ".bak"), "raw", 1e6)
  suppressWarnings(write_manifest(run_paths(path), man))
  expect_identical(readBin(paste0(f, ".bak"), "raw", 1e6), good_bak)
  expect_true(manifest_verifies(paste0(f, ".bak")))
  expect_equal(read_manifest(path)$replications, 3L)
  # a damaged backup is never used
  writeBin(as.raw(1:10), f)
  writeBin(as.raw(1:5), paste0(f, ".bak"))
  expect_error(read_manifest(path), "no readable backup")
  # a first call stopped after creating the folders but before the manifest
  path2 <- run_path()
  for (d in c("batches", "truth", "tmp", "quarantine")) dir.create(file.path(path2, d), recursive = TRUE)
  expect_s3_class(sim_run(st, replications = 2, seed = 1, path = path2), "adempr_run")
})

test_that("a batch deleted with its sidecar is detected and rerun", {
  st <- normal_study(n = 10L, mu = 0)
  path <- run_path()
  run <- sim_run(st, replications = 4, seed = 1, path = path, chunk_size = 1L)
  full <- est_signature(run)
  b <- list.files(file.path(path, "batches"), pattern = "\\.rds$", full.names = TRUE)[1]
  unlink(c(b, sub("\\.rds$", ".md5", b)))
  expect_warning(res <- sim_results(run), "missing from the run directory")
  expect_true(any(grepl("missing", res$problems)))
  run2 <- suppressWarnings(sim_run(st, replications = 4, seed = 1, path = path))
  expect_equal(est_signature(run2), full)
  expect_silent(sim_results(run2))
})

test_that("LPT scheduling leaves results unchanged and orders costly scenarios first", {
  st <- normal_study(n = c(10L, 40L), mu = 0, methods = stochastic_methods())
  a <- sim_run(st, replications = 6, seed = 3, path = run_path(), chunk_size = 2)
  b <- sim_run(st, replications = 6, seed = 3, path = run_path(), chunk_size = 2,
               cost = function(scenario) scenario$n^2)
  expect_identical(est_signature(a), est_signature(b))
  ctx <- list(seed = 3, domain = NULL, registry = a$manifest$rng_registry,
              cache = new.env(parent = emptyenv()), dgp_stage = "dgp", method_stage = "method")
  todo <- adempR:::plan_tasks(st, 2L, NULL)
  w <- adempR:::cost_weights(st, function(scenario) scenario$n)
  chunks <- adempR:::lpt_chunks(adempR:::build_tasks(st, ctx, todo), w, 1L)
  first <- st$scenarios$values[[adempR:::scenario_index(st$scenarios, chunks[[1]][[1]]$scenario_key)]]
  expect_equal(first$n, 40L)
  expect_error(sim_run(st, replications = 1, seed = 3, path = run_path(), cost = function(scenario) -1),
               "nonnegative")
})

test_that("LPT batches carry balanced work, and a random cost is refused before the run starts", {
  # costs 100, 1, 100, 1: contiguous batches of two would hold 101 and 101
  # only by luck; packing gives every batch one costly and one cheap task
  st <- normal_study(n = c(1L, 2L, 3L, 4L), mu = 0)
  w <- stats::setNames(c(100, 1, 100, 1), st$scenarios$keys)
  tasks <- lapply(st$scenarios$keys, function(k) list(scenario_key = k))
  ch <- adempR:::lpt_chunks(tasks, w, 2L)
  loads <- vapply(ch, function(b) sum(w[vapply(b, `[[`, "", "scenario_key")]), 0)
  expect_equal(loads, c(101, 101))
  # eight tasks into four batches of two: loads within one task of each other
  w8 <- c(9, 8, 7, 6, 5, 4, 3, 2)
  t8 <- lapply(seq_along(w8), function(i) list(scenario_key = as.character(i)))
  l8 <- vapply(adempR:::lpt_chunks(t8, stats::setNames(w8, seq_along(w8)), 2L),
               function(b) sum(w8[as.integer(vapply(b, `[[`, "", "scenario_key"))]), 0)
  expect_equal(sort(l8), c(11, 11, 11, 11))
  path <- run_path()
  set.seed(42); before <- runif(1); set.seed(42)
  expect_error(sim_run(st, replications = 1, seed = 3, path = path,
                       cost = function(scenario) runif(1)), "draws random numbers")
  expect_equal(runif(1), before)
  expect_false(dir.exists(path))
})

test_that("planning missing work scales with each scenario's own count", {
  st <- normal_study(n = c(10L, 20L), mu = 0)
  R <- stats::setNames(c(3L, 1e6L), st$scenarios$keys)
  done <- data.frame(scenario_key = st$scenarios$keys[2], replication = seq_len(1e6),
                     method = "t_interval", stringsAsFactors = FALSE)
  todo <- adempR:::plan_tasks(st, R, done)
  expect_equal(nrow(todo), 3L)
  expect_true(all(todo$scenario_key == st$scenarios$keys[1]))
  expect_equal(nrow(adempR:::plan_tasks(st, 2L, NULL)), 4L)
})

test_that("changes that alter a production run need explicit confirmation", {
  st <- normal_study(n = 10L, mu = 0, methods = stochastic_methods(), measures = "bias")
  p <- run_path()
  suppressWarnings(sim_run(st, replications = 3, seed = 1, path = p, label = "production"))
  # an extension is allowed
  # (a production run may warn about uncommitted changes in the git tree)
  expect_no_error(suppressWarnings(sim_run(st, replications = 4, seed = 1, path = p)))
  m2 <- stochastic_methods()
  m2$jitter_a <- function(data) sim_estimate("mean", mean(data$y) + rnorm(1, 0, 0.02))
  st2 <- normal_study(n = 10L, mu = 0, methods = m2, measures = "bias")
  expect_error(suppressWarnings(sim_run(st2, replications = 4, seed = 1, path = p, update = TRUE)), "production run")
  expect_s3_class(suppressWarnings(sim_run(st2, replications = 4, seed = 1, path = p, update = "production")), "adempr_run")
  expect_error(sim_run(st2, replications = 4, seed = 1, path = p, update = "yes"), "must be TRUE")
})

test_that("a takeover in progress excludes other contenders, and a dead one's is recovered", {
  skip_on_os("windows")
  p <- adempR:::run_paths(run_path())
  dir.create(p$root)
  host <- Sys.info()[["nodename"]]
  dead <- list(pid = 999999L, host = host, time = "t0")
  skip_if(isTRUE(adempR:::pid_alive(dead$pid)), "pid 999999 is in use")
  dir.create(p$lock)
  jsonlite::write_json(dead, file.path(p$lock, "owner.json"), auto_unbox = TRUE)
  # another live contender holds the takeover lock: this one must back off
  mutex <- paste0(p$lock, ".takeover")
  dir.create(mutex)
  jsonlite::write_json(list(pid = Sys.getpid(), host = host), file.path(mutex, "owner.json"),
                       auto_unbox = TRUE)
  expect_false(adempR:::take_over_stale_lock(p, dead))
  expect_identical(adempR:::lock_owner(p)$pid, dead$pid)
  # the contender died while holding it: its takeover lock is stale too
  jsonlite::write_json(list(pid = 999999L, host = host), file.path(mutex, "owner.json"),
                       auto_unbox = TRUE)
  tok <- suppressMessages(adempR:::acquire_lock(p))
  expect_identical(adempR:::lock_owner(p)$token, tok)
  expect_false(dir.exists(mutex))
  adempR:::release_lock(p, tok)
})

test_that("a batch committed before its inventory entry is added on the next resume", {
  st <- normal_study(n = 10L, mu = 0)
  path <- run_path()
  run <- sim_run(st, replications = 4, seed = 1, path = path, chunk_size = 1L)
  inv <- file.path(path, "batches", "INVENTORY")
  listed <- readLines(inv)
  # simulate a stop between the commit and the inventory append
  writeLines(listed[-1], inv)
  run2 <- sim_resume(path, replications = 4)
  expect_setequal(readLines(inv), listed)
  # the loss of that batch is now detected
  b <- file.path(path, "batches", paste0(listed[1], c(".rds", ".md5")))
  unlink(b)
  expect_warning(sim_results(run2), "missing from the run directory")
})

test_that("a truth value is published once; later writers get the committed value", {
  dir <- tempfile(); tmp <- tempfile(); dir.create(dir); dir.create(tmp)
  expect_equal(publish_first(1, dir, "t-x", tmp), 1)
  expect_equal(publish_first(2, dir, "t-x", tmp), 1)
  # a claim held by a live writer makes others wait for its value
  dir.create(file.path(dir, "t-y.claim"))
  writeLines(as.character(Sys.getpid()), file.path(dir, "t-y.claim", "pid"))
  expect_equal(publish_first(3, dir, "t-y", tmp, wait = 0.2), 3)
  expect_false(dir.exists(file.path(dir, "t-y.claim")))
  # a claim left by a dead process is recovered at once
  skip_on_os("windows")
  dir.create(file.path(dir, "t-z.claim"))
  writeLines("999999", file.path(dir, "t-z.claim", "pid"))
  expect_equal(publish_first(4, dir, "t-z", tmp, wait = 60), 4)
  expect_equal(readRDS(file.path(dir, "t-z.rds")), 4)
})

test_that("a missing value and the text NA are different content", {
  expect_false(identical(adempR:::col_text(NA_character_), adempR:::col_text("NA")))
  df <- data.frame(k = c(1, 1), note = c(NA, "NA"), stringsAsFactors = FALSE)
  expect_error(adempR:::dedup_table(df, "k", "note", "estimate"), "Conflicting duplicate")
})

test_that("a task with a missing estimand row is not counted as complete", {
  st <- sim_study("two", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
    estimands = list(sim_estimand("m1", sim_truth_exact(function(scenario) 0)),
                     sim_estimand("m2", sim_truth_exact(function(scenario) 0))),
    methods = list(a = function(data) list(sim_estimate("m1", mean(data$y)),
                                           sim_estimate("m2", median(data$y)))))
  path <- run_path()
  run <- sim_run(st, replications = 2, seed = 1, path = path, chunk_size = 1)
  full <- est_signature(run)
  # drop one estimand row from a committed batch, keeping its checksum valid
  b <- list.files(file.path(path, "batches"), pattern = "\\.rds$", full.names = TRUE)[1]
  obj <- readRDS(b)
  obj$estimates <- obj$estimates[obj$estimates$estimand == "m1", , drop = FALSE]
  saveRDS(obj, b)
  writeLines(unname(tools::md5sum(b)), sub("\\.rds$", ".md5", b))
  run2 <- sim_resume(path, replications = 2)
  expect_equal(est_signature(run2), full)
  e <- sim_results(run2)$estimates
  expect_equal(nrow(e), 4L)
})
