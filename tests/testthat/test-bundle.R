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
  expect_message(sim_export(run, dir, report = FALSE), "secret")
  res <- sim_import(dir)
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
               "rng_registry = res$manifest$rng_registry", fixed = TRUE)
  expect_false(dir.exists(file.path(dir, "artifacts")))
})
