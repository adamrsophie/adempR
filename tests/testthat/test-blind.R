test_that("blinded export, external results and a blinded report work end to end", {
  st <- normal_study(n = c(10L, 20L), mu = 0, measures = c("bias", "empse"))
  p <- run_path()
  run <- sim_run(st, replications = 6, seed = 4, path = p)
  dir <- file.path(tempdir(), paste0("blind-", basename(p)))
  key <- sim_export_blinded(run, dir)
  expect_equal(nrow(key), 12L)
  files <- list.files(file.path(dir, "data"))
  expect_setequal(sub("\\.rds$", "", files), key$code)
  # codes do not reveal the scenario or the replication
  expect_false(any(grepl("n=|^[0-9]+$", key$code)))
  expect_true(all(grepl("^[a-f][0-9a-f]{11}$", key$code)))
  expect_error(sim_export_blinded(run, dir), "not empty")
  # an external analyst computes medians on the coded files
  ext <- do.call(rbind, lapply(files, function(f) {
    d <- readRDS(file.path(dir, "data", f))
    data.frame(code = sub("\\.rds$", "", f), estimand = "mean", estimate = median(d$y))
  }))
  ext <- ext[-1, ]  # one data set was not analyzed
  m <- c(st$methods, list(external = sim_external_method(ext, paste0(dir, "-key.rds"),
                                                         description = "returned by analysts")))
  st2 <- normal_study(n = c(10L, 20L), mu = 0, measures = c("bias", "empse"), methods = m)
  run2 <- sim_run(st2, replications = 6, seed = 4, path = p)
  e <- sim_results(run2)$estimates
  ee <- e[e$method == "external", ]
  expect_equal(sum(ee$valid), 11L)
  # the imported values are the analysts' values, matched to the right data set
  k <- key[match(ext$code, key$code), ]
  idx <- match(paste(k$scenario_key, k$replication), paste(ee$scenario_key, ee$replication))
  expect_equal(ee$estimate[idx], ext$estimate)
  expect_error(sim_external_method(data.frame(code = "zz", estimand = "mean", estimate = 1), key),
               "Unknown code")
  # blinded report hides names and code; the key is written apart
  f <- tempfile(fileext = ".html")
  sim_report(run2, f, plots = FALSE, blind_methods = TRUE)
  html <- paste(readLines(f), collapse = "\n")
  # every occurrence of the name "external" (other than the protocol item
  # "external registration") would identify the method
  expect_false(grepl("t_interval|\\bexternal\\b(?! registration)|returned by analysts|t\\.test",
                     html, perl = TRUE))
  # no fingerprint of any method, at any length
  for (fp in unlist(run2$manifest$fingerprints$methods)) {
    expect_false(grepl(substr(fp, 1, 8), html, fixed = TRUE))
  }
  # a fresh order unless a seed is given; the key records it
  bk0 <- utils::read.csv(paste0(f, ".blinding-key.csv"))
  expect_true(is.numeric(bk0$seed))
  expect_true(grepl("method_A", html) && grepl("method_B", html))
  bk <- utils::read.csv(paste0(f, ".blinding-key.csv"))
  expect_setequal(bk$method, c("t_interval", "external"))
})

test_that("a blinded report hides warnings, revisions and keeps summaries unchanged", {
  m <- list(secretalpha = function(data) {
    warning("secretalpha bad warning")
    sim_estimate("mean", mean(data$y))
  })
  by_name <- sim_measure("named", contribution = function(rows) {
    as.numeric(rows$method == "secretalpha")
  })
  st <- normal_study(n = 5L, mu = 0, methods = m, measures = sim_measures("bias", by_name))
  p <- run_path()
  sim_run(st, replications = 3, seed = 1, path = p)
  run <- sim_add_methods(sim_open(p), list(secretbeta = function(data) {
    sim_estimate("mean", median(data$y))
  }), reason = "new competitor")
  f <- tempfile(fileext = ".html")
  sim_report(run, f, plots = FALSE, blind_methods = TRUE, blind_seed = 7)
  html <- paste(readLines(f), collapse = "\n")
  expect_false(grepl("secretalpha|secretbeta", html))
  s <- as.data.frame(sim_summarise(run))
  expect_equal(s$value[s$measure == "named" & s$method == "secretalpha"], 1)
  key <- utils::read.csv(paste0(f, ".blinding-key.csv"))
  lab <- key$label[key$method == "secretalpha"]
  expect_match(html, paste0(lab, "</td>"))
  expect_equal(key$seed[1], 7)
})

test_that("external results are joined by task, so identical data sets keep their own answers", {
  # a fixed design: every replication has the same data
  st <- sim_study("fixed", scenarios = sim_scenarios(n = 3L),
    generate = function(scenario, context) data.frame(y = c(1, 2, 3)),
    estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 2))),
    methods = list(avg = function(data) sim_estimate("mean", mean(data$y))), measures = "bias")
  p <- run_path()
  run <- sim_run(st, replications = 3, seed = 1, path = p)
  dir <- file.path(tempdir(), paste0("blindfix-", basename(p)))
  key <- sim_export_blinded(run, dir, seed = 7)
  expect_equal(length(unique(key$data_md5)), 1L)
  ans <- c(11, 12, 13)
  ext <- data.frame(code = key$code, estimand = "mean", estimate = ans[key$replication])
  m <- c(st$methods, list(external = sim_external_method(ext, key)))
  st2 <- sim_study("fixed", scenarios = st$scenarios, generate = st$generator,
                   estimands = st$estimands, methods = m, measures = "bias")
  e <- sim_results(sim_run(st2, replications = 3, seed = 1, path = p))$estimates
  ee <- e[e$method == "external", ]
  expect_equal(ee$estimate[order(ee$replication)], ans)
  # tables are validated when the method is built
  expect_error(sim_external_method(rbind(ext, ext[1, ]), key), "more than one row")
  bad_key <- key
  bad_key$code[2] <- bad_key$code[1]
  expect_error(sim_external_method(ext, bad_key), "duplicated codes")
  ext_chr <- ext
  ext_chr$estimate <- as.character(ext_chr$estimate)
  expect_error(sim_external_method(ext_chr, key), "must be numeric")
})

test_that("blinded export covers the committed cohort, keeps the key apart and records omissions", {
  st <- normal_study(n = c(10L, 20L), mu = 0, measures = "bias")
  p <- run_path()
  # one scenario raised to five replications
  run <- sim_run(st, replications = 3, seed = 4, path = p,
                 scenario_replications = stats::setNames(5L, st$scenarios$labels[1]))
  dir <- file.path(tempdir(), paste0("blindcoh-", basename(p)))
  k1 <- sim_export_blinded(run, dir)
  expect_equal(nrow(k1), 8L)
  expect_false(file.exists(file.path(dir, "key.rds")))
  expect_true(file.exists(paste0(dir, "-key.rds")))
  expect_false(is.null(attr(k1, "seed")))
  # a fresh seed by default: a second export has other codes
  dir2 <- paste0(dir, "b")
  k2 <- sim_export_blinded(run, dir2)
  expect_length(intersect(k1$code, k2$code), 0L)
  # asking for a replication the run lacks records it as omitted
  dir3 <- paste0(dir, "c")
  k3 <- sim_export_blinded(run, dir3, replications = 1:6, seed = 1)
  om <- attr(k3, "omitted")
  expect_true(all(om$reason == "not generated in the run"))
  expect_equal(nrow(om), 1L + 3L)
  expect_error(sim_export_blinded(run, paste0(dir, "d"), key = file.path(paste0(dir, "d"), "k.rds")),
               "outside")
})

test_that("no spelling of a path puts the key inside the export", {
  st <- normal_study(n = 5L, mu = 0, measures = "bias")
  run <- sim_run(st, replications = 2, seed = 1, path = run_path(), progress = FALSE)
  base <- tempfile("blindalias-")
  dir.create(base)
  dir <- file.path(base, "export")
  inside <- list(
    file.path(base, ".", "export", "key.rds"),
    file.path(base, "export", "sub", "..", "key.rds"),
    file.path(base, "other", "..", "export", "key.rds"),
    dir)
  for (k in inside) {
    expect_error(sim_export_blinded(run, dir, key = k, seed = 1), "outside", info = k)
    expect_false(file.exists(file.path(dir, "key.rds")))
  }
  # a symbolic link to the parent of the export
  if (.Platform$OS.type != "windows") {
    link <- file.path(tempdir(), basename(tempfile("blindlink-")))
    file.symlink(base, link)
    expect_error(sim_export_blinded(run, dir, key = file.path(link, "export", "key.rds"), seed = 1),
                 "outside")
  }
  # a genuinely separate sibling key, in a directory that does not exist yet
  sib <- file.path(base, "keys", "export-key.rds")
  dir.create(dirname(sib))
  k <- sim_export_blinded(run, dir, key = sib, seed = 1)
  expect_true(file.exists(sib))
  expect_false(any(grepl("key", list.files(dir, recursive = TRUE))))
  # a directory whose name only starts like the export is outside it
  expect_equal(adempR:::canonical_path(file.path(base, "export2", "k.rds")),
               paste0(adempR:::canonical_path(base), "/export2/k.rds"))
})

test_that("an existing key is never overwritten", {
  st <- normal_study(n = 5L, mu = 0, measures = "bias")
  run <- sim_run(st, replications = 2, seed = 1, path = run_path(), progress = FALSE)
  key <- tempfile("existing-key-", fileext = ".rds")
  saveRDS("old", key)
  dir <- tempfile("blind-exists-")
  expect_error(sim_export_blinded(run, dir, key = key, seed = 1), "exists")
  expect_identical(readRDS(key), "old")
  expect_false(dir.exists(file.path(dir, "data")))
})
