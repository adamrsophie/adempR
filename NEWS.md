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
  and analyzed on the look prefix), and mergeable accumulators.
* Accounting: `sim_flow()` follows every planned replication to a usable
  result or the reason it is not one; `sim_acceptable()` judges a measure
  against an acceptable range with its Monte Carlo interval; a large
  sample sanity run in `sim_validate()`.
* Diagnostics flag outliers, data sets that repeat across replications,
  tests and intervals that disagree, and degenerate estimates. Zip plots
  show the Monte Carlo band and have a power version.
* A self contained HTML report (protocol, code, flow of replications,
  results, paired comparisons, diagnostics) and reproducible export
  bundles.
