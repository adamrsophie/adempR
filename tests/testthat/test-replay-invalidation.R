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

test_that("every kind of study change is recorded as a revision", {
  changes_of <- function(run) unlist(lapply(run$manifest$revisions, `[[`, "changes"))
  m <- stochastic_methods()
  st <- normal_study(n = c(10L, 20L), mu = 0, methods = m)
  p <- run_path()
  sim_run(st, replications = 3, seed = 4, path = p)
  nb <- function() length(list.files(file.path(p, "batches"), pattern = "\\.md5$"))
  n0 <- nb()
  # estimand metadata: re-evaluate, no new fits
  st_e <- st
  st_e$estimands$mean$effect <- "conditional"
  r <- sim_run(st_e, replications = 3, seed = 4, path = p)
  expect_true(any(grepl("metadata of estimand mean changed", changes_of(r))))
  expect_equal(nb(), n0)
  # retries move from evaluation to execution
  st_r <- sim_amend(st_e, failure = sim_failure_policy(retries = 2L), reason = "retries")
  r <- sim_run(st_r, replications = 3, seed = 4, path = p)
  last <- r$manifest$revisions[[length(r$manifest$revisions)]]$changes
  expect_true(any(grepl("execution settings changed \\(retries\\)", last)))
  expect_false(any(grepl("re-evaluate", last)))
  # timeout is an execution setting
  st_t <- st_r
  st_t$methods$jitter_a <- sim_method(m$jitter_a, estimands = "mean", timeout = 30)
  r <- sim_run(st_t, replications = 3, seed = 4, path = p)
  expect_true(any(grepl("timeout of jitter_a", changes_of(r))))
  expect_equal(nb(), n0)
  # removed scenario and method are recorded and left out of summaries
  st_s <- st_t
  st_s$scenarios <- sim_scenarios(n = 10L, mu = 0, sigma = 1)
  st_s$methods$boot_b <- NULL
  r <- sim_run(st_s, replications = 3, seed = 4, path = p)
  ch <- changes_of(r)
  expect_true(any(grepl("removed 1 scenario", ch)))
  expect_true(any(grepl("removed method boot_b", ch)))
  s <- as.data.frame(sim_summarise(r))
  expect_setequal(unique(s$method), "jitter_a")
  # fewer replications: warning, nothing dropped
  expect_warning(r <- sim_run(st_s, replications = 2, seed = 4, path = p), "holds 3 replications")
  expect_equal(r$manifest$replications, 3L)
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

test_that("fingerprints distinguish formulas, calls, strings, dates and environments", {
  fp <- adempR:::fingerprint
  expect_false(fp(y ~ x) == fp(y ~ x + z))
  expect_false(fp(quote(a + b)) == fp(quote(a * b)))
  expect_false(fp("a,b") == fp(c("a", "b")))
  expect_false(fp(as.Date("2026-01-01")) == fp(as.numeric(as.Date("2026-01-01"))))
  expect_false(fp(NA_character_) == fp("NA"))
  e1 <- new.env(); e1$x <- 1
  e2 <- new.env(); e2$x <- 2
  expect_false(fp(e1) == fp(e2))
  mk <- function(fml) { force(fml); function(data) fml }
  expect_false(fp(mk(y ~ x)) == fp(mk(y ~ x + z)))
  expect_false(fp(list(y ~ x)) == fp(list(y ~ x + z)))
})

test_that("formula-valued scenario factors get distinct keys", {
  a <- sim_scenarios(n = c(10L, 20L), model = list(y ~ x))
  b <- sim_scenarios(n = c(10L, 20L), model = list(y ~ x + z))
  expect_false(any(a$keys %in% b$keys))
  expect_length(sim_scenarios(n = 10L, model = list(y ~ x, y ~ x + z))$keys, 2L)
})

test_that("a changed formula in a method factory requires an update", {
  g <- function(scenario, context) {
    x <- rnorm(20); z <- x + rnorm(20)
    data.frame(x = x, z = z, y = x + z + rnorm(20))
  }
  mk <- function(fml) { force(fml); function(data) {
    sim_estimate("b", unname(stats::coef(stats::lm(fml, data))[2])) } }
  mkst <- function(fml) sim_study("f", scenarios = sim_scenarios(n = 20L), generate = g,
    estimands = list(sim_estimand("b", sim_truth_unknown())),
    methods = list(reg = mk(fml)), measures = "empse")
  p <- tempfile("adempR-fml-")
  sim_run(mkst(y ~ x), 3, seed = 1, path = p, progress = FALSE)
  expect_error(sim_run(mkst(y ~ x + z), 3, seed = 1, path = p, progress = FALSE),
               "changed")
})

test_that("results are read by stored fingerprints, not by the current session", {
  assign("adempr_test_helper", function(x) mean(x), envir = globalenv())
  on.exit(if (exists("adempr_test_helper", envir = globalenv()))
    rm("adempr_test_helper", envir = globalenv()), add = TRUE)
  st <- sim_study("g", scenarios = sim_scenarios(n = 5L),
    generate = function(scenario, context) rnorm(scenario$n),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
    methods = list(glob = function(data) sim_estimate("mean", adempr_test_helper(data))),
    measures = "bias")
  p <- tempfile("adempR-glob-")
  sim_run(st, 4, seed = 1, path = p, progress = FALSE)
  rm("adempr_test_helper", envir = globalenv())
  res <- sim_results(sim_open(p))
  expect_equal(nrow(res$estimates), 4L)
  expect_equal(nrow(res$replications), 4L)
  expect_silent(sim_resume(p, progress = FALSE))
})

test_that("replay reports generation failures, checks artifacts and finds current fits", {
  st <- normal_study(n = 10L, mu = 0, methods = list(a = function(data) sim_estimate("mean", 1)),
                     retention = sim_retention(data = TRUE, fits = TRUE))
  path <- run_path()
  run <- sim_run(st, replications = 2, seed = 1, path = path)
  # a later method has its own artifact; replay finds its retained fit
  st2 <- normal_study(n = 10L, mu = 0,
                      methods = list(a = function(data) sim_estimate("mean", 1),
                                     b = function(data) sim_estimate("mean", 3)),
                      retention = sim_retention(data = TRUE, fits = TRUE))
  run2 <- sim_run(st2, replications = 2, seed = 1, path = path)
  rp <- sim_replay(run2, scenario = 1, replication = 1)
  expect_true(rp$matches)
  expect_equal(rp$stored_fits$b$estimate, 3)
  expect_equal(rp$stored_fits$a$estimate, 1)
  # a changed artifact stops replay
  art <- sim_results(run2)$artifacts$artifact[1]
  x <- readRDS(file.path(path, art))
  x$data$y[1] <- 99
  saveRDS(x, file.path(path, art))
  expect_error(sim_replay(run2, scenario = 1, replication = 1), "checksum")
  # a generation failure is returned, not thrown
  st3 <- sim_study("fail", aims = "x", scenarios = sim_scenarios(n = 5L),
                   generate = function(scenario, context) stop("no data"),
                   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
                   methods = list(a = function(data) sim_estimate("mean", 1)))
  run3 <- suppressWarnings(sim_run(st3, replications = 1, seed = 1, path = run_path()))
  rp3 <- sim_replay(run3, scenario = 1, replication = 1)
  expect_equal(rp3$generation$replayed, "error")
  expect_true(rp3$matches)
  expect_null(rp3$data)
  # a retention policy without retained artifacts does not advertise them
  run4 <- sim_run(normal_study(n = 10L, mu = 0), replications = 2, seed = 1, path = run_path())
  run4$manifest$study$retention$data <- TRUE
  cap <- sim_capabilities(run4)
  expect_false(cap$available[grepl("stored data", cap$analysis)])
})

test_that("ambiguous scenario labels are refused", {
  sc <- list(keys = c("k1", "k2"), labels = c("x=1", "x=1"))
  expect_error(adempR:::resolve_scenario(list(scenarios = sc), "x=1"), "several")
})
