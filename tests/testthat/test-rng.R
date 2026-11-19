# Operational tests 13 and 14: random number identity.

test_that("the caller's RNG kind and state are restored, including a missing seed", {
  old <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old)), add = TRUE)
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  set.seed(42)
  before <- .Random.seed
  run <- sim_run(normal_study(n = 10L, mu = 0), replications = 3, seed = 1, path = run_path())
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1:3], c("Mersenne-Twister", "Inversion", "Rounding"))
  rm(".Random.seed", envir = globalenv())
  sim_run(normal_study(n = 10L, mu = 0), replications = 3, seed = 1, path = run_path())
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the registry is append only and ignores display order", {
  k1 <- data.frame(stage = "dgp", scenario_key = c("b", "a"), component = "data")
  r1 <- rng_registry(k1)
  r1b <- rng_registry(k1[2:1, ])
  expect_identical(r1, r1b)
  r2 <- rng_registry(rbind(k1, data.frame(stage = "dgp", scenario_key = "0", component = "data")), r1)
  expect_identical(r2[1:2, ], r1)
  expect_identical(r2$index[r2$scenario_key == "0"], 3L)
})

test_that("results do not depend on scenario display order, method order or chunk size", {
  m <- stochastic_methods()
  s1 <- normal_study(n = c(10L, 20L), mu = c(0, 1), methods = m)
  s2 <- normal_study(n = c(20L, 10L), mu = c(1, 0), methods = rev(m))
  r1 <- sim_run(s1, replications = 12, seed = 77, path = run_path(), chunk_size = 1)
  r2 <- sim_run(s2, replications = 12, seed = 77, path = run_path(), chunk_size = 7)
  expect_identical(est_signature(r1), est_signature(r2))
})

test_that("extending replications and adding a method preserve existing streams", {
  m <- stochastic_methods()
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = m)
  p <- run_path()
  r100 <- sim_run(st, replications = 10, seed = 5, path = p)
  first <- est_signature(r100)
  reg1 <- r100$manifest$rng_registry
  m3 <- c(m, list(third = function(data) sim_estimate("mean", mean(data$y) + runif(1, -0.01, 0.01))))
  st3 <- normal_study(n = c(10L, 20L), mu = 0, methods = m3)
  r200 <- sim_run(st3, replications = 20, seed = 5, path = p)
  reg2 <- r200$manifest$rng_registry
  expect_identical(reg2[seq_len(nrow(reg1)), ], reg1)
  ext <- est_signature(r200)
  old_rows <- ext[ext$method %in% names(m) & ext$replication <= 10, ]
  rownames(old_rows) <- rownames(first) <- NULL
  expect_identical(old_rows, first)
  # the same extension in a single run of 20 replications gives the same values
  # for the original methods
  r20 <- sim_run(st, replications = 20, seed = 5, path = run_path())
  a <- est_signature(r20)
  b <- ext[ext$method %in% names(m), ]
  rownames(a) <- rownames(b) <- NULL
  expect_identical(a, b)
  # the added method is reproducible by replay
  rp <- sim_replay(r200, scenario = 1, replication = 15, methods = "third")
  expect_true(rp$matches)
})

test_that("the future backend reproduces the sequential reference", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_if_dev_namespace()
  st <- normal_study(n = c(10L, 20L), mu = c(0, 1), methods = stochastic_methods())
  seq_run <- sim_run(st, replications = 15, seed = 31, path = run_path())
  par_run <- sim_run(st, replications = 15, seed = 31, path = run_path(), backend = "future",
                     workers = 2, chunk_size = 4)
  expect_identical(est_signature(seq_run), est_signature(par_run))
  op <- future::plan(future::sequential)
  on.exit(future::plan(op), add = TRUE)
  par1 <- sim_run(st, replications = 15, seed = 31, path = run_path(), backend = "future",
                  chunk_size = 9)
  expect_identical(est_signature(seq_run), est_signature(par1))
})
