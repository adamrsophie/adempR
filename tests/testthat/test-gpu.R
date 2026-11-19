test_that("GPU methods see their assigned device and the variable is restored", {
  old <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA)
  Sys.setenv(CUDA_VISIBLE_DEVICES = "keep")
  on.exit(if (is.na(old)) Sys.unsetenv("CUDA_VISIBLE_DEVICES") else Sys.setenv(CUDA_VISIBLE_DEVICES = old), add = TRUE)
  m <- list(
    g = sim_method(function(data, context) {
      sim_estimate("mean", mean(data$y), note = paste(context$gpu, Sys.getenv("CUDA_VISIBLE_DEVICES")))
    }, gpu = TRUE),
    cpu = function(data, context) sim_estimate("mean", mean(data$y),
                                               note = paste(is.null(context$gpu), Sys.getenv("CUDA_VISIBLE_DEVICES"))))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = "bias")
  run <- sim_run(st, replications = 4, seed = 1, path = run_path(), chunk_size = 1, gpus = c(0, 1))
  e <- sim_results(run)$estimates
  g <- e[e$method == "g", ]
  expect_equal(g$note[order(g$replication)], c("0 0", "1 1", "0 0", "1 1"))
  expect_true(all(e$note[e$method == "cpu"] == "TRUE keep"))
  expect_equal(Sys.getenv("CUDA_VISIBLE_DEVICES"), "keep")
  # the flag is part of what the method is, so it enters its fingerprint
  expect_false(identical(adempR:::method_fingerprint(sim_method(m$cpu)),
                   adempR:::method_fingerprint(sim_method(m$cpu, gpu = TRUE))))
  expect_error(sim_run(st, replications = 1, seed = 1, path = run_path(), gpus = character()), "device")
})

test_that("devices are recorded, used in replay, restored after an error, and distinct within a wave", {
  m <- list(g = sim_method(function(data, context) {
    sim_estimate("mean", mean(data$y), note = Sys.getenv("CUDA_VISIBLE_DEVICES"))
  }, gpu = TRUE))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = "bias")
  run <- sim_run(st, replications = 2, seed = 1, path = run_path(), chunk_size = 1, gpus = c(3, 7))
  a <- sim_results(run)$attempts
  expect_equal(a$gpu[order(a$replication)], c("3", "7"))
  rp <- sim_replay(run, 1, 2)
  expect_equal(rp$estimates$note, "7")
  expect_true(rp$matches)
  expect_equal(sim_replay(run, 1, 2, gpu = "5")$estimates$note, "5")
  # an error that escapes the task still restores the variable
  old <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA)
  Sys.setenv(CUDA_VISIBLE_DEVICES = "keep")
  on.exit(if (is.na(old)) Sys.unsetenv("CUDA_VISIBLE_DEVICES") else Sys.setenv(CUDA_VISIBLE_DEVICES = old), add = TRUE)
  task <- list(scenario_key = st$scenarios$keys[1], replication = 1L, methods = "g", gpu = "9",
               seeds = adempR:::task_seeds(list(seed = 1, domain = NULL, registry = run$manifest$rng_registry,
                                                cache = new.env(), dgp_stage = "dgp", method_stage = "method"),
                                           st$scenarios$keys[1], 1L, st, "g")[[1]])
  # isolate_data() runs outside the captured fit, after the device is set
  local_mocked_bindings(isolate_data = function(...) stop("escaped"))
  expect_error(adempR:::execute_task_seeded(st, task), "escaped")
  expect_equal(Sys.getenv("CUDA_VISIBLE_DEVICES"), "keep")
})

test_that("future waves give concurrent batches different devices", {
  skip_on_cran()
  skip_if_not_installed("future.apply")
  skip_if_dev_namespace()
  m <- list(g = sim_method(function(data, context) {
    sim_estimate("mean", mean(data$y), note = paste(Sys.getpid(), Sys.getenv("CUDA_VISIBLE_DEVICES")))
  }, gpu = TRUE))
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = "bias")
  oplan <- future::plan(future::sequential)
  on.exit(future::plan(oplan), add = TRUE)
  run <- sim_run(st, replications = 4, seed = 1, path = run_path(), chunk_size = 1,
                 backend = "future", gpus = c(0, 1), workers = 2)
  a <- sim_results(run)$attempts
  # batches 1 and 2 form the first wave, 3 and 4 the second
  expect_equal(a$gpu[order(a$replication)], c("0", "1", "0", "1"))
})
