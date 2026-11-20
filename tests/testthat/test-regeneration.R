test_that("a failed regeneration recorded by an older version does not block reading", {
  st <- normal_study(n = 5L, mu = 0, measures = "bias")
  p <- run_path()
  run <- sim_run(st, replications = 2, seed = 1, path = p)
  ref <- sim_results(run)$estimates
  # an older version committed a failed regeneration next to the
  # successful generation, in the same round
  rp <- adempR:::run_paths(p)
  raw <- adempR:::read_raw_results(p)$replications[1, , drop = FALSE]
  raw$generation_status <- "error"
  raw$generation_message <- "source unavailable"
  raw$data_md5 <- NA_character_
  raw$artifact <- NA_character_
  batch <- list(schema_version = adempR:::SCHEMA_VERSION, created = adempR:::now_utc(),
                replications = raw, attempts = NULL, estimates = NULL, diagnostics = NULL,
                records = NULL, dataset_truth = NULL)
  nm <- paste0("b-", adempR:::new_uuid())
  adempR:::atomic_save(batch, rp$batches, nm, rp$tmp)
  cat(nm, "\n", sep = "", file = file.path(rp$batches, "INVENTORY"), append = TRUE)
  res <- sim_results(run)
  expect_equal(res$estimates$estimate, ref$estimate)
  expect_true(all(res$replications$generation_status == "ok"))
})
