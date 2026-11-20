# adempR: Reproducible Statistical Simulation Studies <img src="man/figures/logo.png" align="right" height="139" alt="adempR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/adempR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/adempR/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/adamrsophie/adempR/graph/badge.svg)](https://app.codecov.io/gh/adamrsophie/adempR)
[![pkgdown](https://github.com/adamrsophie/adempR/actions/workflows/pkgdown.yaml/badge.svg)](https://adamrsophie.github.io/adempR/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- DOI badge: Zenodo mints a concept DOI on the first GitHub release; replace XXXXXXXX and uncomment. -->
<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXXX) -->
[![Archived on Software Heritage](https://archive.softwareheritage.org/badge/origin/https://github.com/adamrsophie/adempR/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/adamrsophie/adempR)
<!-- badges: end -->

***adempR*** runs statistical simulation studies as reproducible scientific
experiments. A study is declared with the ADEMP structure of Morris, White and
Crowther (2019): **A**ims, **D**ata-generating mechanisms, **E**stimands,
**M**ethods and **P**erformance measures. Each part is a declared R object,
the study is fingerprinted, and the run is recorded in a directory that can be
interrupted, resumed, extended, replayed, summarized and exported.

Simulation studies are easy to get subtly wrong: a method judged against a
parameter it does not target, failed fits silently dropped, results that
change when a scenario is added or the job runs in parallel, a resumed run
that counts replications twice, performance reported without Monte Carlo
error. adempR is built to make each of these errors either impossible or
visible:

- **Estimands are first class.** Each estimand declares its population,
  conditioning set, effect type and scale; methods declare what they target,
  and declared mismatches of effect, scale, population or conditioning are
  validation errors (a method that does not declare them is not checked). Truth can be exact, numerical, Monte
  Carlo (with its own standard error carried into the summary), data set
  specific, or unknown.
- **Each measure has its own denominator.** Estimate, standard error,
  interval and test validity are classified separately, failures are reported
  by default, and composite measures over all generated data sets are
  available.
- **Random numbers are tied to content, not position.** L'Ecuyer-CMRG streams
  are keyed by stage, scenario content and component, so results do not
  change when scenarios or methods are reordered, chunks are resized, the
  backend changes, or replications and methods are added (for methods that
  draw only from the stream they are given).
- **Runs are recoverable and auditable.** Atomic commits with checksums,
  quarantine of corrupt files, duplicate and conflict detection, stale lock
  recovery, selective invalidation by fingerprint, protocol freezing and
  amendments, and replay of any single replication. Recovery is tested by
  interrupting runs at named points between commits; a crash elsewhere is
  covered by the same commit protocol but not by a dedicated test.
- **Every measure comes with its Monte Carlo error where one is defined.**
  Seventeen built in measures with MCSEs, exact intervals for proportions,
  paired method comparisons, precision planning, pilots, and sequential
  stopping with intervals that stay valid after optional stopping for
  independent replications. A cell whose MCSE is undefined says so in its
  status instead of reporting a number.

## Installation

```r
# Development version
# install.packages("remotes")
remotes::install_github("adamrsophie/adempR")
```

adempR needs R 4.2.0 or later. Its hard dependencies are `jsonlite`,
`codetools` and base R packages. It contains no compiled code. The parallel
backend uses `future` and `future.apply` when they are installed.

## Quick start

Compare the t interval with the Wald interval for a normal mean.

```r
library(adempR)

study <- sim_study(
  name = "normal_mean",
  aims = "Coverage and power of t and Wald intervals for a normal mean",
  scenarios = sim_scenarios(n = c(10L, 40L), mu = c(0, 0.5)),
  generate = function(scenario, context) rnorm(scenario$n, scenario$mu),
  estimands = list(sim_estimand("mean", population = "normal DGP",
    truth = sim_truth_exact(function(scenario) scenario$mu))),
  methods = list(
    t = function(data) {
      tt <- t.test(data)
      sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
                   lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95,
                   p_value = tt$p.value, null = 0)
    },
    wald = function(data) {
      m <- mean(data); se <- sd(data) / sqrt(length(data))
      sim_estimate("mean", m, se = se, lower = m - 1.96 * se, upper = m + 1.96 * se,
                   level = 0.95, p_value = 2 * pnorm(-abs(m / se)), null = 0)
    }),
  measures = c("bias", "empse", "coverage", "rejection"))

sim_validate(study, smoke = TRUE)
#> <adempR validation: 0 error(s), 0 warning(s)>

path <- file.path(tempdir(), "normal-mean")
run <- sim_run(study, replications = 1000, seed = 2026, path = path)
```

Performance with Monte Carlo standard errors. The rejection rate is labeled
as the type I error where the null is true and as power elsewhere.

```r
perf <- as.data.frame(sim_summarise(run))
perf[perf$measure %in% c("coverage", "rejection"),
     c("scenario_id", "method", "label", "value", "mcse")]
#>   scenario_id method        label value    mcse
#>    n=10, mu=0      t     coverage 0.946 0.00715
#>    n=10, mu=0      t type_i_error 0.054 0.00715
#>    n=10, mu=0   wald     coverage 0.909 0.00909
#>    n=10, mu=0   wald type_i_error 0.091 0.00909
#>  n=10, mu=0.5      t        power 0.285 0.01427
#>  n=10, mu=0.5   wald        power 0.382 0.01536
#>  ...
```

Both methods saw the same data sets, so their difference is estimated from
the paired replications and is far more precise than the two separate MCSEs
suggest:

```r
cmp <- sim_compare(run, c("t", "wald"), measures = "coverage")
cmp[c("scenario_id", "value", "mcse", "n_common")]
#>   scenario_id  value    mcse n_common
#>    n=10, mu=0 -0.037 0.00597     1000
#>    n=40, mu=0 -0.006 0.00244     1000
#>  n=10, mu=0.5 -0.027 0.00513     1000
#>  n=40, mu=0.5 -0.009 0.00299     1000
```

Stop R at any point and call `sim_run()` again with the same arguments: only
missing replications run. Raise `replications` to extend the study; the
first 1000 stay unchanged. Recompute any single replication from its
recorded random number streams, then write the report and a reproducible
bundle:

```r
sim_replay(run, scenario = 1, replication = 17)$matches
#> [1] TRUE

sim_report(run, file.path(tempdir(), "normal-mean.html"))
sim_export(run, file.path(tempdir(), "normal-mean-bundle"))
```

The bundle holds the results, the code, an renv lockfile of the packages
that produced them, a plain text record of the software environment and a
Dockerfile that recreates the R version, so the study can be rerun on
another machine.

## Workflow

| Step | Functions |
|------|-----------|
| Protocol and design | `sim_protocol()`, `sim_protocol_export()`, `sim_scenarios()`, `sim_constraint()`, `sim_scenarios_table()`, `sim_scenarios_sample()`, `sim_scenarios_bind()`, `sim_completeness()` |
| Estimands and truth | `sim_estimand()`, `sim_estimand_grid()`, `sim_truth_exact()`, `sim_truth_numeric()`, `sim_truth_mc()`, `sim_truth_dataset()`, `sim_truth_unknown()`, `sim_truth_value()` |
| Data and methods | `sim_generator()`, `sim_data()`, `sim_method()`, `sim_method_sbc()`, `sim_external_method()`, `sim_estimate()`, `sim_estimate_grid()`, `sim_failure_policy()`, `sim_retention()`, `sim_calibrate_method()` |
| Study | `sim_study()`, `sim_validate()`, `sim_fingerprint()`, `sim_freeze()`, `sim_amend()` |
| Planning | `sim_plan_replications()`, `sim_pilot()`, `sim_run_adaptive()` |
| Execution | `sim_run()`, `sim_resume()`, `sim_open()`, `sim_unlock()`, `sim_run_sequential()`, `sim_run_shard()`, `sim_collect_shards()`, `sim_add_methods()` |
| Results | `sim_results()`, `sim_truth()`, `sim_summarise()`, `sim_summarise_stream()`, `sim_summarise_grid()`, `sim_summarise_antithetic()`, `sim_control_variate()`, `sim_flow()`, `sim_acceptable()`, `sim_compare()`, `sim_pareto()`, `sim_rank()`, `sim_leaderboard()`, `sim_bias_decomposition()`, `sim_factor_influence()`, `sim_confidence_sequence()`, `sim_measures()`, `sim_measure()` |
| Checking | `sim_check_dgp()`, `sim_diagnose()`, `sim_check_reference()`, `sim_check_identity()`, `sim_verify()`, `sim_sbc()`, `sim_plot()`, `sim_replay()`, `sim_capabilities()` |
| Reporting | `sim_report()`, `sim_export()`, `sim_import()`, `sim_export_blinded()` |
| Generator helpers | `sim_rmvnorm()`, `sim_rcopula()`, `sim_plasmode()`, `sim_check_matrix()`, `sim_calibrate_intercept()`, `sim_calibrate_censoring()`, `sim_rweibull_ph()`, `sim_add_missing()`, `sim_rclustered()`, `sim_contaminate()` |
| Low level | `rng_registry()`, `rng_task_seed()`, `sim_accumulator()`, `sim_accumulate()`, `sim_merge()`, `sim_accumulated()` |

## Performance measures

| Measure | Eligible replications | MCSE |
|---------|-----------------------|------|
| `bias`, `mae`, `mse` | usable estimate and truth | sample SD of contributions over sqrt(n) |
| `relbias` | usable estimate, nonzero truth | MCSE(bias) / truth; delta method for varying truth |
| `rmse` | usable estimate and truth | delta method |
| `empse` | usable estimate (errors, for a data set truth) | kurtosis based variance of the sample variance |
| `mean_se`, `rmvar` | usable estimate with a valid standard error | mean, or delta method |
| `se_ratio` | usable estimate with a valid standard error | jackknife (the normal theory value of Morris et al. 2019 is in the note) |
| `coverage`, `noncoverage_lower`, `noncoverage_upper` | valid interval and truth | binomial, with Clopper-Pearson interval |
| `be_coverage` | valid interval and usable estimate | grouped jackknife over replications |
| `width` | valid interval | sample SD over sqrt(n) |
| `rejection` | valid test | binomial; labeled type I error or power only when the truth settles the null |
| `usable`, `usable_coverage` | all generated data sets | binomial; status `incomplete` while some are not analyzed, bounds when truths are missing |

Custom measures are defined with `sim_measure()`; measure factories add
selection proportions for categorical estimands
(`sim_measures_categorical()`, applied only to categorical estimands) and
the interval score (`sim_measure_interval_score()`, optionally divided by a
fixed scale or a positive column such as the truth). No MCSE is claimed for
a nonlinear aggregate unless one is supplied. Comparisons add paired
differences of bias, MAE, MSE, coverage, width, rejection and the usable
proportion, and relative precision (paired bootstrap MCSE, with the
normal theory value in the note). Printed summaries show each value to
the precision its MCSE supports.

## Documentation

The [website](https://adamrsophie.github.io/adempR/) has the reference
manual and these articles:

- [Getting started](https://adamrsophie.github.io/adempR/articles/adempR.html)
- [Protocols and scenario designs](https://adamrsophie.github.io/adempR/articles/protocol-and-scenarios.html)
- [Estimands and true values](https://adamrsophie.github.io/adempR/articles/estimands-and-truth.html)
  (including marginal and conditional odds ratios)
- [Failures and denominators](https://adamrsophie.github.io/adempR/articles/failures-and-denominators.html)
  (logistic regression with separation)
- [Random numbers and reproducibility](https://adamrsophie.github.io/adempR/articles/rng-and-reproducibility.html)
- [Execution, recovery and replay](https://adamrsophie.github.io/adempR/articles/execution-and-recovery.html)
- [Performance measures and Monte Carlo error](https://adamrsophie.github.io/adempR/articles/performance-measures.html)
- [Comparing methods](https://adamrsophie.github.io/adempR/articles/comparisons.html)
- [Planning replications and sequential stopping](https://adamrsophie.github.io/adempR/articles/planning-and-sequential.html)
- [Diagnostics and reporting](https://adamrsophie.github.io/adempR/articles/diagnostics-and-reporting.html)
- [Worked examples](https://adamrsophie.github.io/adempR/articles/examples.html):
  clustered data, survival at a horizon, missing data, meta-analysis,
  noncompliance (ITT and CACE)
- [adempR and other simulation packages](https://adamrsophie.github.io/adempR/articles/ecosystem.html)
- [Capabilities, status and evidence](https://adamrsophie.github.io/adempR/articles/capabilities.html)

## Validation

- Built in measures and their MCSEs are tested against hand calculations
  from the formulas, and against
  [rsimsum](https://CRAN.R-project.org/package=rsimsum) where the
  definitions agree (for `se_ratio`, rsimsum's point estimate and the
  normal theory MCSE given in the note; the reported MCSE is a jackknife).
- Scientific tests check known results: the t interval for a normal mean,
  informative failures and conditional versus composite denominators, the
  noncollapsible marginal odds ratio with Monte Carlo truth, paired
  comparisons, nested resampling, and the simultaneous coverage of the
  sequential stopping rule.
- Operational tests compare results bit for bit across scenario order,
  method order, chunk size, the `future` backend (for top-level and
  self-contained functions; nested closures that reach global objects are
  refused), extension and method addition; kill R with `SIGKILL` at each
  stage of a commit, inside a method and inside data generation, and check
  that the resumed run equals an uninterrupted one; and check corruption
  quarantine, lost batches, manifest damage, duplicate and conflict
  detection, stale locks, replay, selective invalidation, retention,
  regeneration from a bundle, and reanalysis of an imported bundle in a
  clean process.
- Continuous integration runs `R CMD check` on Ubuntu (R release, devel,
  the previous release and R 4.2), Windows and macOS.

The [capability matrix](https://adamrsophie.github.io/adempR/articles/capabilities.html)
lists every feature as verified, implemented or experimental, with the test
that supports each claim. Sharded runs are tested with separate R processes
and batchtools jobs on one machine, not on a real cluster; GPU support
assigns devices to methods and is not tested on GPU hardware.

## Related packages

[SimDesign](https://CRAN.R-project.org/package=SimDesign),
[simChef](https://github.com/Yu-Group/simChef),
[SimEngine](https://CRAN.R-project.org/package=SimEngine),
[simpr](https://CRAN.R-project.org/package=simpr) and
[simsalapar](https://CRAN.R-project.org/package=simsalapar) run simulation
studies; [rsimsum](https://CRAN.R-project.org/package=rsimsum),
[simhelpers](https://CRAN.R-project.org/package=simhelpers) and
[simMetric](https://CRAN.R-project.org/package=simMetric) analyze their
results; [simstudy](https://CRAN.R-project.org/package=simstudy) and
[simsurv](https://CRAN.R-project.org/package=simsurv) generate data. adempR
results are plain data frames that work with the analysis packages, and its
generators can call the data generation packages. See the
[ecosystem article](https://adamrsophie.github.io/adempR/articles/ecosystem.html)
for a comparison.

## Citing adempR

```r
citation("adempR")
```

The repository ships machine readable citation and software metadata:

- [`CITATION.cff`](CITATION.cff): Citation File Format (GitHub "Cite this
  repository", Zenodo).
- [`codemeta.json`](codemeta.json): CodeMeta 2.0 software metadata.
- [`.zenodo.json`](.zenodo.json): Zenodo deposit metadata.
- `inst/CITATION`: the R level citation.
- `inst/REFERENCES.bib`: the methodological references, checked against
  Crossref.

Please also cite the methodological work your study relies on, in
particular Morris, White and Crowther (2019),
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086).

## Contributing

Bug reports, feature requests and pull requests are welcome; see
[`CONTRIBUTING.md`](CONTRIBUTING.md). A wrong number is the most useful bug
report: include a small study that reproduces it and the value you expected.
Please follow the [Code of Conduct](CODE_OF_CONDUCT.md), and report security
issues privately as described in [`SECURITY.md`](SECURITY.md).

## Author

- Adam R Sophie (author, maintainer, copyright holder),
  <adamrsophie@gmail.com>

## License

GPL-3. See [`LICENSE.md`](LICENSE.md) or
<https://www.gnu.org/licenses/gpl-3.0> for the full license text.
