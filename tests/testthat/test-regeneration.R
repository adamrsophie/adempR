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

test_that("shards with different retry rounds collect and retry above every imported round", {
  flag <- tempfile("svc-")
  writeLines("down", flag)
  on.exit(unlink(flag), add = TRUE)
  method <- local({
    f <- flag
    function(data) {
      if (file.exists(f)) stop("service unavailable")
      sim_estimate("mean", mean(data$y))
    }
  })
  gen <- function(scenario, context) data.frame(y = stats::rnorm(scenario$n), n = scenario$n)
  st <- sim_study("shard_rounds", scenarios = sim_scenarios(n = c(5L, 6L, 7L)),
                  generate = gen,
                  estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
                  methods = list(a = method), measures = "bias")
  parent <- run_path()
  into <- run_path()
  # shard k is retried k - 1 times, so the shards reach rounds 0, 1 and 2
  for (k in 1:3) {
    sim_run_shard(st, k, 3, 2, 1, parent, progress = FALSE)
    for (i in seq_len(k - 1L)) sim_run_shard(st, k, 3, 2, 1, parent, retry_failed = TRUE,
                                             progress = FALSE)
  }
  rounds <- vapply(1:3, function(k) as.integer(adempR:::read_manifest(
    adempR:::shard_path(parent, k))$retry_round %||% 0L), integer(1))
  col <- sim_collect_shards(st, 3, 2, 1, parent, into)
  imported <- max(sim_results(col)$attempts_all$round)
  expect_equal(rounds, 0:2)
  expect_equal(col$manifest$retry_round, 2L)
  unlink(flag)
  rr <- sim_resume(into, retry_failed = TRUE, progress = FALSE)
  res <- sim_results(rr)
  expect_true(rr$manifest$retry_round > imported)
  expect_true(all(res$attempts$status == "ok"))
  expect_equal(nrow(res$estimates), 6L)
  # the failed attempts stay in the history
  expect_true(any(res$attempts_all$status != "ok"))
})
