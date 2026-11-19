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
  replications leaves existing streams unchanged. The caller's RNG is
  restored.
* Recoverable execution: atomic commits with checksums, one coordinator
  lock with stale lock recovery, quarantine of corrupt files, duplicate
  detection with conflict checks, selective invalidation and revision
  records, a `future` backend, replay of single replications.
* Performance: 17 built in measures with Monte Carlo standard errors,
  exact intervals for proportions, Monte Carlo truth uncertainty, paired
  method comparisons, Pareto views, precision planning, pilots, sequential
  stopping with error spending over looks and cells, and mergeable
  accumulators.
* Diagnostics, plots, a self contained HTML report and reproducible export
  bundles.
