## Test environments

* local macOS (aarch64), R release
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1 and 4.2),
  windows-latest (release), macos-latest (release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

Tests that start separate R processes or kill them are skipped on CRAN,
and tests that run parallel workers use at most two. Examples, tests and
vignettes write only to temporary directories.
