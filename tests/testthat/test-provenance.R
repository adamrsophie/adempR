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
