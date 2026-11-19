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
  expect_error(sim_report(run, f), "exists")
})
