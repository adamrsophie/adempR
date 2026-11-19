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
  # without devices the method runs normally and the flag leaves fingerprints alone
  expect_identical(adempR:::method_fingerprint(sim_method(m$cpu)),
                   adempR:::method_fingerprint(sim_method(m$cpu, gpu = TRUE)))
  expect_error(sim_run(st, replications = 1, seed = 1, path = run_path(), gpus = character()), "device")
})
