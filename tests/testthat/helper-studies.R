# Shared study builders for the tests.

normal_study <- function(n = c(30L, 100L), mu = c(0, 0.5), methods = NULL,
                         measures = c("bias", "mse", "rmse", "empse", "mean_se",
                                      "se_ratio", "coverage", "rejection"), ...) {
  if (is.null(methods)) {
    methods <- list(t_interval = sim_method(function(data, scenario, context) {
      fit <- t.test(data$y, mu = 0, conf.level = 0.95)
      sim_estimate("mean", estimate = mean(data$y), se = sd(data$y) / sqrt(nrow(data)),
                   lower = unname(fit$conf.int[1]), upper = unname(fit$conf.int[2]),
                   level = 0.95, interval_type = "t", df = unname(fit$parameter),
                   p_value = fit$p.value, null = 0)
    }))
  }
  sim_study(
    name = "normal_mean",
    aims = "Evaluate the t interval for a normal mean",
    scenarios = sim_scenarios(n = n, mu = mu, sigma = 1),
    generate = function(scenario, context) {
      data.frame(y = rnorm(scenario$n, scenario$mu, scenario$sigma))
    },
    estimands = list(sim_estimand("mean", population = "normal DGP", effect = "marginal",
      truth = sim_truth_exact(function(scenario) scenario$mu))),
    methods = methods, measures = measures, ...)
}

# Two stochastic methods: their output depends on their own random stream.
stochastic_methods <- function() {
  list(
    jitter_a = function(data) sim_estimate("mean", mean(data$y) + rnorm(1, 0, 0.01)),
    boot_b = function(data) {
      b <- replicate(20, mean(sample(data$y, replace = TRUE)))
      sim_estimate("mean", mean(b), se = sd(b))
    })
}

est_signature <- function(run) {
  e <- sim_results(run)$estimates
  e <- e[order(e$scenario_key, e$method, e$estimand, e$replication), ]
  out <- e[c("scenario_key", "replication", "method", "estimand", "estimate", "se", "lower",
             "upper", "p_value", "estimate_status")]
  rownames(out) <- NULL
  out
}

run_path <- function() tempfile("adempR-test-")

# Parallel workers and child processes load the installed package, so tests
# that use them are skipped under devtools::load_all().
skip_if_dev_namespace <- function() {
  testthat::skip_if(exists(".__DEVTOOLS__", envir = asNamespace("adempR"), inherits = FALSE),
                    "needs the installed package")
}
