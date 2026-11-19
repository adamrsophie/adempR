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
  expect_error(sim_export_blinded(run, dir), "not empty")
  # an external analyst computes medians on the coded files
  ext <- do.call(rbind, lapply(files, function(f) {
    d <- readRDS(file.path(dir, "data", f))
    data.frame(code = sub("\\.rds$", "", f), estimand = "mean", estimate = median(d$y))
  }))
  ext <- ext[-1, ]  # one data set was not analyzed
  m <- c(st$methods, list(external = sim_external_method(ext, file.path(dir, "key.rds"),
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
