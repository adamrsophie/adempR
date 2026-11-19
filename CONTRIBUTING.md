# Contributing to adempR

Thank you for considering a contribution to `adempR`. The package is used to
draw conclusions about statistical methods, so a change is judged on two
things: whether the code is sound, and whether the numbers it reports are
right, with their Monte Carlo error and their denominators stated.

## Before You Start

- Use GitHub issues for bug reports, feature requests, documentation gaps and
  questions.
- For larger changes, such as a new performance measure, a new backend or a
  change to the stored result schema, open an issue first so the design can be
  discussed before implementation.
- Keep pull requests focused. A new feature, a bug fix and a refactor should
  usually be separate pull requests.

## Reporting a Wrong Number

A wrong result is the most valuable bug report. Please include:

- a small study (few scenarios, few replications) that reproduces it,
- what adempR reports and what you expected, with the formula or reference
  you rely on,
- `sessionInfo()`.

## Development Setup

```r
install.packages(c("remotes", "testthat", "devtools", "roxygen2"))
remotes::install_deps(dependencies = TRUE)
```

adempR has no compiled code, so no C or C++ toolchain is needed.

## Rules

- Expected values in tests come from an independent calculation: a hand
  calculation from the formula, a closed form result, another package, or a
  bit for bit comparison between two ways of running the same study. Never
  take the expected value from the function under test.
- Every performance measure needs an eligibility rule, an MCSE (or an explicit
  statement that none is available) and a status for degenerate cases.
- Generators and methods draw random numbers only through R's generator; code
  in adempR never calls `set.seed()` on the caller's behalf without restoring
  the caller's state.
- A change to the stored result schema needs a schema version bump and a
  migration path.
- Update `vignettes/capabilities.Rmd` when a feature's status changes.
- Keep roxygen tags directly above the function they document.

## Branches and Pull Requests

- Use one branch per logical change, named after the change
  (`add-control-variates`, `fix-relbias-mcse`).
- Write commit subjects in the imperative ("Add control variate estimators").
- Keep the implementation, its tests, the roxygen source, the regenerated
  documentation and the `NEWS.md` bullet together in the same pull request.

## Testing

```r
devtools::test()
```

Tests that start child processes or parallel workers need the installed
package and are skipped under `devtools::load_all()`; run them with
`R CMD check` before a pull request:

```sh
R CMD build .
R CMD check --as-cran adempR_*.tar.gz
```

## Pull Request Checklist

- Tests relevant to the change pass, with independently derived expected
  values.
- Documentation and the capability matrix are updated.
- `NEWS.md` has an entry for user-facing changes.
- No generated build artifacts, run directories, `.Rcheck` directories or
  tarballs are committed.

## Code of Conduct

All contributors are expected to follow the project
[Code of Conduct](CODE_OF_CONDUCT.md).
