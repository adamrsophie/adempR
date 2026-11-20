# Operational test 22 and reporting.

test_that("export, import and reanalysis reproduce the summary", {
  st <- normal_study(n = 10L, mu = c(0, 1))
  run <- sim_run(st, replications = 20, seed = 6, path = run_path())
  dir <- tempfile("bundle-")
  sim_export(run, dir, report = TRUE)
  expect_true(all(file.exists(file.path(dir, c("results.rds", "manifest.json", "README.md",
                                               "report.html", "CHECKSUMS.md5",
                                               "tables/estimates.csv", "code/generator.R")))))
  txt <- paste(readLines(file.path(dir, "manifest.json")), collapse = "\n")
  expect_false(grepl(Sys.info()[["nodename"]], txt, fixed = TRUE))
  expect_false(grepl(normalizePath(tempdir()), txt, fixed = TRUE))
  res <- sim_import(dir)
  a <- as.data.frame(sim_summarise(run))
  b <- as.data.frame(sim_summarise(res))
  expect_equal(a$value, b$value)
  expect_error(sim_export(run, dir), "not empty")
  write("tampered", file.path(dir, "README.md"), append = TRUE)
  expect_error(sim_import(dir), "Checksum mismatch")
})

test_that("a clean process can reanalyse a bundle with only adempR", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_dev_namespace()
  run <- sim_run(normal_study(n = 10L, mu = 0), replications = 10, seed = 6, path = run_path())
  dir <- tempfile("bundle-")
  sim_export(run, dir, report = FALSE)
  out <- callr::r(function(dir) {
    res <- adempR::sim_import(dir)
    as.data.frame(adempR::sim_summarise(res))$value
  }, list(dir))
  expect_equal(out, as.data.frame(sim_summarise(run))$value)
  # the regeneration script also runs in a clean process
  regen <- callr::r(function(dir) {
    setwd(dir)
    source("regenerate.R", local = TRUE)
    e <- adempR::sim_results(run)$estimates
    e[order(e$scenario_key, e$replication, e$method), "estimate"]
  }, list(dir))
  e <- sim_results(run)$estimates
  expect_equal(regen, e[order(e$scenario_key, e$replication, e$method), "estimate"])
})

test_that("the report is self-contained and escapes text", {
  st <- normal_study(n = 10L, mu = 0)
  st$aims <- "Check <b>escaping</b> & ampersands"
  run <- sim_run(st, replications = 10, seed = 6, path = run_path())
  f <- tempfile(fileext = ".html")
  sim_report(run, f, plots = FALSE)
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "&lt;b&gt;escaping&lt;/b&gt; &amp; ampersands", fixed = TRUE)
  expect_match(h, "doi.org/10.1002/sim.8086", fixed = TRUE)
  # self-contained: no external scripts, stylesheets or images
  expect_false(grepl("<script[^>]*src=|<link[^>]*href=|<img[^>]*src=\"http", h))
  expect_error(sim_report(run, f), "exists")
  # with plots, every image is embedded in the file
  f2 <- tempfile(fileext = ".html")
  sim_report(run, f2)
  h2 <- paste(readLines(f2), collapse = "\n")
  imgs <- regmatches(h2, gregexpr("<img[^>]*src=\"[^\"]*", h2))[[1]]
  expect_gt(length(imgs), 0L)
  expect_true(all(grepl("src=\"data:image/", imgs)))
})

test_that("the report shows every protocol field, code and comparisons", {
  m <- stochastic_methods()
  st <- normal_study(n = 10L, mu = 0, methods = m, measures = c("bias", "empse"),
                     protocol = sim_protocol(users = c("statisticians", "trialists"),
                                             precision_targets = c(bias = 0.01),
                                             phase = "II"))
  run <- sim_run(st, replications = 6, seed = 2, path = run_path())
  f <- tempfile(fileext = ".html")
  sim_report(run, f, plots = FALSE)
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "statisticians; trialists", fixed = TRUE)
  expect_match(h, "bias = 0.01", fixed = TRUE)
  expect_match(h, ">phase<", fixed = TRUE)
  expect_match(h, "Paired comparisons with the first method", fixed = TRUE)
  expect_match(h, "Truth code of mean", fixed = TRUE)
})

test_that("bundles carry no source paths, stale files, formula cells or foreign schemas", {
  # a study built inside a function, from sourced code with source references
  make <- function() {
    secret <- "do-not-share"
    src <- tempfile(fileext = ".R")
    writeLines(c("m <- function(data) {", "  # comment", "  sim_estimate('mean', mean(data$y))", "}"), src)
    e <- new.env()
    sys.source(src, envir = e, keep.source = TRUE)
    m <- e$m
    environment(m) <- environment()
    list(study = normal_study(n = 10L, mu = 0, methods = list(m = m)), src = src)
  }
  mk <- make()
  expect_false(is.null(attr(mk$study$methods$m$fun, "srcref")))
  run <- sim_run(mk$study, replications = 3, seed = 1, path = run_path())
  expect_null(attr(run$manifest$study$methods$m$fun, "srcref"))
  dir <- tempfile("bundle-")
  expect_no_message(sim_export(run, dir, report = FALSE), message = "secret")
  res <- sim_import(dir)
  expect_false(exists("secret", envir = environment(res$manifest$study$methods$m$fun),
                      inherits = FALSE))
  expect_false(grepl("do-not-share", rawToChar(serialize(res, NULL, ascii = TRUE)), fixed = TRUE))
  expect_null(attr(res$manifest$study$methods$m$fun, "srcref"))
  expect_false(grepl(dirname(mk$src), paste(readLines(file.path(dir, "manifest.json")), collapse = "")))
  expect_false(grepl("/", res$manifest$environment$blas %||% "", fixed = TRUE))
  # overwrite replaces the earlier bundle entirely
  writeLines("old", file.path(dir, "code", "method-old.R"))
  suppressMessages(sim_export(run, dir, report = FALSE, overwrite = TRUE))
  expect_false(file.exists(file.path(dir, "code", "method-old.R")))
  expect_false(any(grepl("method-old", readLines(file.path(dir, "CHECKSUMS.md5")))))
  other <- tempfile("notbundle-")
  dir.create(other)
  writeLines("x", file.path(other, "keep.txt"))
  expect_error(suppressMessages(sim_export(run, other, overwrite = TRUE)), "not an adempR bundle")
  expect_true(file.exists(file.path(other, "keep.txt")))
  # a foreign schema is refused
  mj <- file.path(dir, "manifest.json")
  j <- jsonlite::read_json(mj)
  j$schema_version <- "9.9.9"
  jsonlite::write_json(j, mj, auto_unbox = TRUE)
  ck <- file.path(dir, "CHECKSUMS.md5")
  lines <- readLines(ck)
  lines[grepl("  manifest\\.json$", lines)] <- paste0(unname(tools::md5sum(mj)), "  manifest.json")
  writeLines(lines, ck)
  expect_error(sim_import(dir), "schema version 9.9.9")
})

test_that("text cells that spreadsheets read as formulas are neutralized", {
  x <- adempR:::csv_neutralize(c("=HYPERLINK(1)", "+a", "-1.5", "@x", "ok", NA, "-text"))
  expect_identical(x, c("'=HYPERLINK(1)", "'+a", "-1.5", "'@x", "ok", NA, "'-text"))
})

test_that("regeneration with the recorded registry reproduces an extended run", {
  m <- stochastic_methods()
  path <- run_path()
  st1 <- normal_study(n = 10L, mu = 0, methods = m["jitter_a"])
  sim_run(st1, replications = 3, seed = 4, path = path)
  st2 <- normal_study(n = 10L, mu = 0, methods = m)
  run <- sim_run(st2, replications = 3, seed = 4, path = path)
  res <- sim_results(run)
  re <- sim_run(res$manifest$study, replications = 3, seed = 4, path = run_path(),
                label = res$manifest$label, rng_registry = res$manifest$rng_registry)
  expect_equal(est_signature(re), est_signature(run))
  fresh <- sim_run(res$manifest$study, replications = 3, seed = 4, path = run_path())
  expect_false(isTRUE(all.equal(est_signature(fresh), est_signature(run))))
  dir <- tempfile()
  sim_export(run, dir, report = FALSE)
  expect_match(paste(readLines(file.path(dir, "README.md")), collapse = "\n"),
               "rng_registry = m$rng_registry", fixed = TRUE)
  expect_false(dir.exists(file.path(dir, "artifacts")))
})

test_that("exported functions keep only the bindings they use and still run", {
  make <- function(k) {
    secret <- "do-not-share"
    helper <- function(x) x * k
    function(data) sim_estimate("mean", helper(mean(data$y)))
  }
  run <- sim_run(normal_study(n = 10L, mu = 0, methods = list(m = make(2))),
                 replications = 3, seed = 1, path = run_path())
  dir <- tempfile("bundle-")
  expect_message(sim_export(run, dir, report = FALSE), "helper, k|k, helper")
  res <- sim_import(dir)
  f <- res$manifest$study$methods$m$fun
  expect_setequal(ls(environment(f)), "helper")
  expect_setequal(ls(environment(environment(f)$helper)), "k")
  expect_false(grepl("do-not-share", rawToChar(serialize(res, NULL, ascii = TRUE)), fixed = TRUE))
  # the exported study identifies and computes the same as the original
  expect_identical(method_fingerprint(res$manifest$study$methods$m),
                   method_fingerprint(run$manifest$study$methods$m))
  d <- data.frame(y = c(1, 2, 3))
  expect_equal(f(d)$estimate, 4)
  # a function that looks names up dynamically keeps its whole environment
  dyn <- local({ hidden <- 1; function(data) sim_estimate("mean", get("hidden") + mean(data$y)) })
  run2 <- sim_run(normal_study(n = 10L, mu = 0, methods = list(d = dyn)),
                  replications = 2, seed = 1, path = run_path())
  expect_message(sim_export(run2, tempfile("bundle-"), report = FALSE),
                 "whole enclosing environments are exported, including: hidden")
})

test_that("a failed overwrite keeps the earlier bundle, and import checks the bundle as a whole", {
  run <- sim_run(normal_study(n = 10L, mu = 0), replications = 3, seed = 1, path = run_path())
  dir <- tempfile("bundle-")
  suppressMessages(sim_export(run, dir, report = FALSE))
  before <- readLines(file.path(dir, "CHECKSUMS.md5"))
  local_mocked_bindings(sim_flow = function(...) stop("simulated failure"))
  expect_error(suppressMessages(sim_export(run, dir, report = FALSE, overwrite = TRUE)),
               "simulated failure")
  expect_identical(readLines(file.path(dir, "CHECKSUMS.md5")), before)
  expect_s3_class(sim_import(dir), "adempr_results")
  expect_length(list.files(dirname(dir), pattern = paste0("^\\.", basename(dir)), all.files = TRUE), 0L)
  # a file outside the checksums is refused
  writeLines("x", file.path(dir, "extra.txt"))
  expect_error(sim_import(dir), "not covered by its checksums: extra.txt")
  unlink(file.path(dir, "extra.txt"))
  # the JSON manifest must agree with the results object
  mj <- file.path(dir, "manifest.json")
  j <- jsonlite::read_json(mj)
  j$seed <- 99
  jsonlite::write_json(j, mj, auto_unbox = TRUE)
  ck <- file.path(dir, "CHECKSUMS.md5")
  lines <- readLines(ck)
  lines[grepl("  manifest\\.json$", lines)] <- paste0(unname(tools::md5sum(mj)), "  manifest.json")
  writeLines(lines, ck)
  expect_error(sim_import(dir), "disagree on: seed")
  # a checksum list naming only results.rds is not a bundle
  writeLines(lines[grepl("results\\.rds", lines)], ck)
  expect_error(sim_import(dir), "does not list manifest.json")
})

test_that("include_data copies only the artifacts of the current replications", {
  flag <- tempfile(); writeLines("x", flag)
  m <- function(data) {
    if (file.exists(flag)) stop("cannot allocate vector of size 9 Gb")
    sim_estimate("mean", mean(data$y))
  }
  st <- normal_study(n = 5L, mu = 0, methods = list(m = m),
                     retention = sim_retention(data = TRUE))
  path <- run_path()
  sim_run(st, replications = 2, seed = 1, path = path)
  unlink(flag)
  run <- sim_run(st, replications = 2, seed = 1, path = path, retry_failed = TRUE)
  on_disk <- list.files(file.path(path, "artifacts"), pattern = "\\.rds$")
  used <- basename(stats::na.omit(sim_results(run)$replications$artifact))
  dir <- tempfile("bundle-")
  suppressMessages(sim_export(run, dir, report = FALSE, include_data = TRUE))
  copied <- list.files(file.path(dir, "artifacts"), pattern = "\\.rds$")
  expect_setequal(copied, used)
  expect_true(length(on_disk) >= length(used))
  expect_s3_class(sim_import(dir), "adempr_results")
})

test_that("the bundle's regeneration script reproduces per scenario counts, and the record shows them", {
  st <- normal_study(n = 10L, mu = c(0, 1))
  path <- run_path()
  k2 <- st$scenarios$keys[2]
  run <- sim_run(st, replications = 3, seed = 4, path = path, label = "pilot",
                 scenario_replications = stats::setNames(6L, k2))
  dir <- tempfile("bundle-")
  suppressMessages(sim_export(run, dir, report = TRUE))
  j <- jsonlite::read_json(file.path(dir, "manifest.json"))
  expect_identical(j$rng_domain, "pilot")
  expect_equal(j$scenario_replications[[k2]], 6)
  expect_equal(unname(unlist(j$analyzed_replications)), c(3, 6))
  h <- paste(readLines(file.path(dir, "report.html")), collapse = "\n")
  expect_match(h, "from 3 to 6 per scenario", fixed = TRUE)
  expect_match(h, "Random number domain: pilot", fixed = TRUE)
  expect_match(h, "With 3 replications (the smallest", fixed = TRUE)
  # the bundle carries a lockfile of the creating session, covered by the
  # checksums, and the recipe checks it first
  expect_true(file.exists(file.path(dir, "renv.lock")))
  expect_match(readLines(file.path(dir, "CHECKSUMS.md5")), "renv.lock", all = FALSE)
  expect_match(readLines(file.path(dir, "regenerate.R")), "sim_check_environment", all = FALSE)
  # a plain text environment record and a Dockerfile for the recorded R
  # version, both covered by the checksums
  et <- readLines(file.path(dir, "environment.txt"))
  expect_true(any(startsWith(et, "R:") & grepl(R.version.string, et, fixed = TRUE)))
  expect_match(et, "^random number generator:\\s+L'Ecuyer-CMRG", all = FALSE)
  expect_match(et, "^adempR\\s", all = FALSE)
  df <- readLines(file.path(dir, "Dockerfile"))
  expect_true(paste0("FROM rocker/r-ver:", R.version$major, ".", R.version$minor) %in% df)
  expect_match(df, "renv::restore(lockfile = \"renv.lock\"", fixed = TRUE, all = FALSE)
  expect_match(readLines(file.path(dir, "CHECKSUMS.md5")), "  Dockerfile$", all = FALSE)
  expect_match(readLines(file.path(dir, "CHECKSUMS.md5")), "  environment.txt$", all = FALSE)
  # run the published recipe from the bundle directory
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  env <- new.env(parent = globalenv())
  suppressMessages(sys.source("regenerate.R", envir = env))
  setwd(old)
  expect_equal(est_signature(env$run), est_signature(run))
  expect_equal(table(sim_results(env$run)$replications$scenario_key)[[k2]], 6L)
})
