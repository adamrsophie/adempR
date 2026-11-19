test_that("git provenance is recorded without credentials", {
  expect_identical(adempR:::strip_credentials("https://user:secret@example.org/a/b.git"),
                   "https://example.org/a/b.git")
  expect_identical(adempR:::strip_credentials("git@example.org:a/b.git"), "git@example.org:a/b.git")
  skip_if(!nzchar(Sys.which("git")), "git is not installed")
  repo <- tempfile("adempR-git-")
  dir.create(repo)
  git <- function(...) system2("git", c("-C", shQuote(repo), ...), stdout = FALSE, stderr = FALSE)
  git("init", "-q")
  git("config", "user.email", "t@example.org")
  git("config", "user.name", "t")
  writeLines("x", file.path(repo, "a.txt"))
  git("add", "a.txt")
  git("commit", "-q", "-m", "init")
  old <- options(adempR.git_dir = repo)
  on.exit(options(old), add = TRUE)
  g <- adempR:::git_info()
  expect_match(g$commit, "^[0-9a-f]{40}$")
  expect_false(g$dirty)
  expect_equal(g$untracked, 0L)
  # untracked files are counted, not mistaken for changes to tracked files
  writeLines("new", file.path(repo, "helper.R"))
  g2 <- adempR:::git_info()
  expect_equal(g2$untracked, 1L)
  expect_false(g2$dirty)
  unlink(file.path(repo, "helper.R"))
  st <- sim_freeze(normal_study(n = 10L, mu = 0))
  expect_identical(st$status$frozen$git_commit, g$commit)
  run <- sim_run(st, replications = 2, seed = 1, path = run_path())
  expect_identical(run$manifest$environment$git$commit, g$commit)
  expect_length(run$manifest$sessions, 1L)
  f <- tempfile(fileext = ".html")
  sim_report(run, f, plots = FALSE)
  expect_match(paste(readLines(f), collapse = ""), "untracked files (not identified by the commit)", fixed = TRUE)
  # uncommitted changes in a production run give a warning
  writeLines("y", file.path(repo, "a.txt"))
  expect_warning(run2 <- sim_run(st, replications = 3, seed = 1, path = run$path),
                 "uncommitted changes")
  expect_length(run2$manifest$sessions, 2L)
  expect_null(adempR:::git_info(tempdir()))
})

test_that("a session with a different environment is recorded and reported", {
  st <- normal_study(n = 10L, mu = 0)
  run <- sim_run(st, replications = 2, seed = 1, path = run_path())
  man <- run$manifest
  man$environment$r_version <- "R version 0.0.0"
  expect_warning(man2 <- adempR:::record_session(man), "environment differs")
  expect_length(man2$sessions, 2L)
})

test_that("thread settings are part of the compared environment and shown in the report", {
  st <- normal_study(n = 10L, mu = 0)
  run <- sim_run(st, replications = 2, seed = 1, path = run_path())
  man <- run$manifest
  man$environment$threads$OMP_NUM_THREADS <- "64"
  expect_warning(adempR:::record_session(man), "thread settings")
  f <- tempfile(fileext = ".html")
  sim_report(run, f, plots = FALSE)
  h <- paste(readLines(f), collapse = "\n")
  expect_match(h, "OMP_NUM_THREADS", fixed = TRUE)
})

test_that("runs record package sources and write a lockfile that the check reads", {
  st <- normal_study(n = 5L, mu = 0, methods = list(
    m = function(data) sim_estimate("mean", stats::median(data$y))))
  run <- sim_run(st, replications = 2, seed = 1, path = run_path())
  env <- run$manifest$environment
  expect_true("adempR" %in% names(env$package_sources))
  expect_true("jsonlite" %in% names(env$package_sources))
  expect_false(is.null(env$locale))
  expect_type(env$renv$active, "logical")
  lock <- tempfile(fileext = ".lock")
  sim_lockfile(run, lock)
  j <- jsonlite::read_json(lock)
  expect_equal(j$R$Version, paste(R.version$major, R.version$minor, sep = "."))
  expect_equal(j$Packages$jsonlite$Version, as.character(utils::packageVersion("jsonlite")))
  chk <- sim_check_environment(lock)
  expect_true(all(chk$match[chk$package != "adempR"]))
  # a recorded version that is not installed is reported
  j$Packages$jsonlite$Version <- "0.0.1"
  jsonlite::write_json(j, lock, auto_unbox = TRUE)
  expect_warning(chk <- sim_check_environment(lock), "jsonlite 0.0.1 recorded")
  expect_false(chk$match[chk$package == "jsonlite"])
})

test_that("system requirements reach the Dockerfile and environment.txt", {
  recs <- list(xml = list(Package = "xmlish", Version = "1.0", Source = "Repository",
                          Repository = "CRAN", SystemRequirements = "libxml2 (>= 2.9)"),
               gh = list(Package = "ghpkg", Version = "0.2", Source = "GitHub",
                         RemoteUsername = "someone", RemoteRepo = "ghpkg",
                         RemoteSha = "0123456789abcdef0123"))
  txt <- adempR:::dockerfile_text(list(r_number = "4.4.1", package_sources = recs))
  expect_true("FROM rocker/r-ver:4.4.1" %in% txt)
  expect_true("#   xmlish: libxml2 (>= 2.9)" %in% txt)
  expect_false(any(grepl("ghpkg", txt)))
  et <- adempR:::environment_text(list(r_version = "R version 4.4.1", rng_kind = RNG_KIND,
                                       package_sources = recs))
  expect_match(et, "^xmlish\\s+1\\.0\\s+CRAN\\s+libxml2", all = FALSE)
  expect_match(et, "^ghpkg\\s+0\\.2\\s+GitHub someone/ghpkg@0123456789ab", all = FALSE)
  # an installed package that declares system requirements keeps them
  with_req <- Filter(function(r) !is.null(r$SystemRequirements),
                     adempR:::package_records(c("curl", "xml2", "openssl", "Rcpp")))
  skip_if(!length(with_req), "no installed package declares system requirements")
  r <- with_req[[1]]
  expect_identical(r$SystemRequirements,
                   gsub("\\s+", " ", utils::packageDescription(r$Package)$SystemRequirements))
})
