# Operational tests 19, 20 and 21.

test_that("replay restores the data set and random state", {
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = stochastic_methods(),
                     retention = sim_retention(data = TRUE, fits = TRUE))
  run <- sim_run(st, replications = 8, seed = 3, path = run_path())
  rp <- sim_replay(run, scenario = 2, replication = 6)
  expect_true(rp$matches)
  expect_true(rp$data_matches)
  expect_equal(nrow(rp$data), st$scenarios$values[[2]]$n)
  rp2 <- sim_replay(run, scenario = st$scenarios$labels[1], replication = 1, methods = "boot_b")
  expect_true(rp2$matches)
  # a missing artifact is reported, not fatal
  art <- sim_results(run)$replications$artifact[1]
  unlink(file.path(run$path, art))
  k <- sim_results(run)$replications
  rp3 <- sim_replay(run, scenario = k$scenario_key[1], replication = k$replication[1])
  expect_true(is.na(rp3$data_matches))
  expect_true(rp3$matches)
})

test_that("changes invalidate exactly the affected work", {
  m <- stochastic_methods()
  st <- normal_study(n = 10L, mu = 0, methods = m)
  p <- run_path()
  run <- sim_run(st, replications = 6, seed = 4, path = p)
  base <- est_signature(run)
  nb <- function() length(list.files(file.path(p, "batches"), pattern = "\\.md5$"))
  n0 <- nb()
  # evaluation change: no new fits, revision recorded
  st_eval <- sim_amend(st, measures = sim_measures("bias", "mae"), reason = "add MAE")
  r1 <- sim_run(st_eval, replications = 6, seed = 4, path = p)
  expect_equal(nb(), n0)
  expect_true(any(grepl("re-evaluate", unlist(lapply(r1$manifest$revisions, `[[`, "changes")))))
  expect_true("mae" %in% sim_summarise(r1)$measure)
  # truth change: no new fits, evaluation changes
  st_truth <- st
  st_truth$estimands$mean$truth <- sim_truth_exact(function(scenario) 1)
  r2 <- sim_run(st_truth, replications = 6, seed = 4, path = p)
  expect_equal(nb(), n0)
  s2 <- as.data.frame(sim_summarise(r2))
  expect_lt(s2$value[s2$measure == "bias" & s2$method == "jitter_a"], -0.5)
  # method change: refused without update, refits only that method with it
  m2 <- m
  m2$boot_b <- function(data) {
    b <- replicate(30, mean(sample(data$y, replace = TRUE)))
    sim_estimate("mean", mean(b), se = sd(b))
  }
  st_m <- normal_study(n = 10L, mu = 0, methods = m2)
  expect_error(sim_run(st_m, replications = 6, seed = 4, path = p), "update = TRUE")
  r3 <- sim_run(st_m, replications = 6, seed = 4, path = p, update = TRUE)
  sig <- est_signature(r3)
  a_old <- base[base$method == "jitter_a", ]
  a_new <- sig[sig$method == "jitter_a", ]
  rownames(a_old) <- rownames(a_new) <- NULL
  expect_identical(a_new, a_old)
  expect_false(identical(sig$estimate[sig$method == "boot_b"], base$estimate[base$method == "boot_b"]))
  raw <- adempR:::read_raw_results(p)
  expect_equal(sum(raw$attempts$method == "jitter_a"), 6)
  # data-generating mechanism change: refused
  st_g <- st_m
  st_g$generator <- sim_generator(function(scenario, context) data.frame(y = rnorm(scenario$n, 1)))
  expect_error(sim_run(st_g, replications = 6, seed = 4, path = p, update = TRUE),
               "data-generating mechanism differs")
  # declared dependency change counts as a method change
  m3 <- m2
  m3$jitter_a <- sim_method(m2$jitter_a, depends = list(table = 1:3))
  expect_error(sim_run(normal_study(n = 10L, mu = 0, methods = m3), replications = 6, seed = 4,
                       path = p), "jitter_a")
})

test_that("closure variables are part of method fingerprints", {
  mk <- function(k) function(data) sim_estimate("mean", mean(data$y) * k)
  expect_false(identical(adempR:::method_fingerprint(sim_method(mk(1))),
                         adempR:::method_fingerprint(sim_method(mk(2)))))
  expect_identical(adempR:::method_fingerprint(sim_method(mk(1))),
                   adempR:::method_fingerprint(sim_method(mk(1))))
})

test_that("retention policies determine what can be recomputed", {
  st <- normal_study(n = 10L, mu = 0)
  cap <- sim_capabilities(st)
  expect_false(cap$available[grepl("fitted model", cap$analysis)])
  expect_true(cap$available[grepl("performance measures", cap$analysis)])
  st2 <- normal_study(n = 10L, mu = 0, retention = sim_retention(fits = TRUE))
  run <- sim_run(st2, replications = 2, seed = 1, path = run_path())
  cap2 <- sim_capabilities(run)
  expect_true(cap2$available[grepl("fitted model", cap2$analysis)])
  art <- readRDS(file.path(run$path, sim_results(run)$replications$artifact[1]))
  expect_s3_class(art$fits$t_interval, "adempr_estimate")
})
