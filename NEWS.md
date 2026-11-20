# adempR 0.1.0

First release.

* Study contracts: scenarios with content keys, constraints, sampled and
  nested designs; estimands with declared population, conditioning, effect
  and scale; exact, numerical, Monte Carlo, data set specific and unknown
  truth; generators, methods with public scenario views, and estimates with
  separate estimate, SE, interval and test validity.
* Failure handling: raw attempts and estimates are stored; validity is
  classified at evaluation under a failure policy, so changing the policy
  needs no refitting. Each performance measure uses its own denominator,
  and failures are reported by default.
* Reproducibility: an append only L'Ecuyer-CMRG stream registry keyed by
  stage, scenario and component, so results do not depend on scenario
  order, method order, chunk size or backend, and adding methods or
  replications leaves existing streams unchanged. Pilots and data
  generating mechanism checks draw from their own stream families. The
  caller's RNG is restored.
* Provenance: frozen protocols with recorded amendments, study
  fingerprints that cover scenario weights and designs, the computing
  environment and git commit of every session, and a revision for every
  change to a study.
* Software environment: every session records the version, source and
  system requirements of each package in the dependency closure, the
  repositories, numerical libraries, thread settings, generator kinds,
  collation locale, time zone and renv state. `sim_lockfile()` writes an
  renv lockfile from a run and `sim_check_environment()` compares the
  installed packages with it. Export bundles carry the lockfile, a plain
  text environment record and a Dockerfile on the Rocker image for the
  recorded R version, and the regeneration script checks the environment
  first.
* Recoverable execution: atomic commits with checksums, a checksummed
  manifest with a backup, an inventory that detects lost batches, one
  coordinator lock with stale lock recovery, quarantine of corrupt files,
  duplicate detection with conflict checks, selective invalidation and
  revision records, a `future` backend that resolves declared packages and
  helper functions as the sequential backend does, a task budget checked
  before any work, replay of single replications against checksum verified
  retained data and fits, and regeneration from an export bundle with the
  recorded stream registry.
* Performance: 17 built in measures with Monte Carlo standard errors,
  exact intervals for proportions, Monte Carlo truth uncertainty, paired
  method comparisons, Pareto views, precision planning, pilots, sequential
  stopping with error spending over looks and cells (stored with the run
  and analyzed on the look prefix), and mergeable accumulators. Each
  measure states the replications it uses and a status when its value or
  MCSE is not defined, or when unattempted work leaves it incomplete
  (reported with bounds).
* Accounting: `sim_flow()` follows every planned replication to a usable
  result or the reason it is not one; `sim_acceptable()` judges a measure
  against an acceptable range with its Monte Carlo interval, or with the
  monitoring interval of a sequential run; a large
  sample sanity run in `sim_validate()`.
* Diagnostics flag outliers, data sets that repeat across replications,
  tests and intervals that disagree at the study's alpha, and degenerate
  estimates. Zip plots show the Monte Carlo band and have a null version.
* A self contained HTML report (protocol, code, flow of replications,
  results, paired comparisons, diagnostics) and reproducible export
  bundles.
* Estimands beyond scalars: categorical estimands with selection
  proportions, and vector or function valued estimands as components on a
  grid with grid averaged performance. Truths can be shared by scenarios
  that agree on given factors.
* Methods: fallback values counted only under the failure policy, per fit
  diagnostics, methods with known performance checked against Monte Carlo
  error, a substream for resampling inside a method, tuning of method
  settings on separate calibration streams, and GPU device assignment.
* Generators: per replication records of realized properties, plasmode
  resampling of real covariates, a Gaussian copula with arbitrary margins,
  and calibration of a censoring rate to a target proportion.
* Analysis: method rankings with bootstrap probabilities that state when
  they hold only for an uncertain truth (and whether that truth can
  reverse first place), a leaderboard,
  bias decomposition against a secondary reference value with MCSEs of
  every part, identity checks between estimands with tolerances in residual
  units, optional relative scales, stated truth error bounds,
  and pass, fail or unavailable statuses, a descriptive meta-model of factor influence, the
  interval score, weighted summaries within factor levels, time-uniform
  confidence sequences, simulation-based calibration of Bayesian methods,
  and batch by batch summaries for runs larger than memory.
* Variance reduction and allocation: common random numbers across
  scenarios, antithetic pairs analyzed by pair in every summary (the
  number of independent units is reported, and proportions whose pair
  means are all equal get an exact interval from the pairs), control
  variates from recorded quantities with the coefficient from the units,
  and adaptive allocation of replications to the scenarios that need them.
* Workflow: protocol export in the ADEMP-PreReg section order, blinded
  evaluation (coded export of the committed data with the key stored
  apart, external results joined by scenario and replication, a report
  with blinded method names), reproduction and seed sensitivity
  checks between runs over every stored table, sharded runs for cluster array jobs, living
  benchmarks that add methods to an existing run, longest first
  scheduling, and protection of production runs from changes that alter
  their results. Bland-Altman and interval by truth plots.
