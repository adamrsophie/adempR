# adempR <img src="man/figures/logo.png" align="right" height="139" alt="adempR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/adempR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/adempR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

adempR runs statistical simulation studies as reproducible experiments.
A study is declared with the ADEMP structure of Morris, White and Crowther
(2019): **A**ims, **D**ata-generating mechanisms, **E**stimands,
**M**ethods and **P**erformance measures. adempR then

* checks that every method is evaluated against the estimand it actually
  targets, with exact, numerical, Monte Carlo or data set specific truth;
* keeps every failed fit and gives each performance measure its own
  denominator, with failures reported by default;
* assigns random number streams by content, so results do not change when
  scenarios or methods are reordered, chunks are resized, the backend
  changes, or replications and methods are added;
* records progress atomically, so interrupted runs resume without double
  counting, and any replication can be replayed;
* reports Monte Carlo standard errors for every measure, paired method
  comparisons, precision planning and sequential stopping with valid
  intervals;
* writes a self-contained HTML report and a reproducible export bundle.

## Installation

```r
# install.packages("remotes")
remotes::install_github("adamrsophie/adempR")
```

## Example

```r
library(adempR)

study <- sim_study(
  name = "normal_mean",
  aims = "Coverage of the t interval for a normal mean",
  scenarios = sim_scenarios(n = c(10L, 40L), mu = c(0, 0.5)),
  generate = function(scenario, context) rnorm(scenario$n, scenario$mu),
  estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu))),
  methods = list(t = function(data) {
    tt <- t.test(data)
    sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
                 lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95,
                 p_value = tt$p.value, null = 0)
  }),
  measures = c("bias", "empse", "coverage", "rejection"))

path <- file.path(tempdir(), "normal-mean")
run <- sim_run(study, replications = 1000, seed = 2026, path = path)
sim_summarise(run)

# interrupt at any time, then call sim_run() again to resume
sim_replay(run, scenario = 1, replication = 17)$matches
sim_report(run, file.path(tempdir(), "report.html"))
```

## Documentation

The [website](https://adamrsophie.github.io/adempR/) has a getting
started guide and articles on protocols and scenarios, estimands and
truth, failures and denominators, random numbers, recovery and replay,
performance measures and their Monte Carlo error, comparisons, planning
and sequential stopping, diagnostics and reporting, worked examples, the
relation to other simulation packages, and a
[capability matrix](https://adamrsophie.github.io/adempR/articles/capabilities.html)
that states what is verified and by which test.

## License

GPL-3.
