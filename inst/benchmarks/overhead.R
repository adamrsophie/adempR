# Overhead of adempR compared with a plain R loop and two other frameworks.
#
# The same study in every implementation: normal data with n in {20, 100},
# two methods (mean and median), R replications per scenario. The work per
# replication is tiny on purpose, so the timings measure framework overhead
# (bookkeeping, storage, random number streams), not statistics. Run with
#   Rscript inst/benchmarks/overhead.R [replications]
# Frameworks that are not installed are skipped.

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args)) as.integer(args[1]) else 1000L
ns <- c(20L, 100L)
library(adempR)

time_it <- function(expr) unname(system.time(expr)[["elapsed"]])
results <- list()
add <- function(name, seconds, note = "") {
  results[[length(results) + 1L]] <<- data.frame(implementation = name, seconds = seconds,
                                                 per_replication_ms = 1000 * seconds / (R * length(ns)),
                                                 note = note)
}

# 1. Plain R: a loop with set.seed, results kept in memory
add("plain R loop", time_it({
  set.seed(1)
  out <- vector("list", length(ns))
  for (i in seq_along(ns)) {
    m <- matrix(NA_real_, R, 2)
    for (r in seq_len(R)) {
      y <- rnorm(ns[i])
      m[r, ] <- c(mean(y), median(y))
    }
    out[[i]] <- colMeans(m)
  }
}), "no storage, one stream")

# 2. adempR: committed batches on disk, one stream per task and method
study <- sim_study("overhead", scenarios = sim_scenarios(n = ns),
  generate = function(scenario, context) data.frame(y = rnorm(scenario$n)),
  estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
  methods = list(mean = function(data) sim_estimate("mean", mean(data$y)),
                 median = function(data) sim_estimate("mean", median(data$y))),
  measures = c("bias", "empse"))
for (cs in c(25L, 250L)) {
  path <- tempfile("adempR-bench-")
  add(sprintf("adempR (chunk_size = %d)", cs),
      time_it(sim_summarise(sim_run(study, replications = R, seed = 1, path = path,
                                    chunk_size = cs, progress = FALSE))),
      "atomic commits, per-task streams, summary included")
  unlink(path, recursive = TRUE)
}

# 3. SimDesign
if (requireNamespace("SimDesign", quietly = TRUE)) {
  design <- SimDesign::createDesign(n = ns)
  gen <- function(condition, fixed_objects) rnorm(condition$n)
  ana <- function(condition, dat, fixed_objects) c(mean = mean(dat), median = median(dat))
  summ <- function(condition, results, fixed_objects) colMeans(results)
  add("SimDesign", time_it(suppressMessages(SimDesign::runSimulation(design, replications = R,
      generate = gen, analyse = ana, summarise = summ, seed = rep(1L, length(ns)),
      verbose = FALSE, progress = FALSE))), "in memory")
}

# 4. simChef
if (requireNamespace("simChef", quietly = TRUE)) {
  dgp <- simChef::create_dgp(.dgp_fun = function(n) list(y = rnorm(n)), .name = "normal")
  meth <- simChef::create_method(.method_fun = function(y) list(mean = mean(y), median = median(y)),
                                 .name = "both")
  exper <- simChef::create_experiment(name = "overhead", save_dir = tempfile("simchef-")) |>
    simChef::add_dgp(dgp) |> simChef::add_method(meth) |>
    simChef::add_vary_across(.dgp = "normal", n = ns)
  add("simChef", time_it(suppressMessages(simChef::run_experiment(exper, n_reps = R, save = FALSE,
                                                                 verbose = 0))), "in memory")
}

out <- do.call(rbind, results)
out$relative_to_plain <- out$seconds / out$seconds[1]
print(out, row.names = FALSE, digits = 3)
cat(sprintf("\nR %s, %s, %d replications x %d scenarios x 2 methods\n",
            getRversion(), R.version$platform, R, length(ns)))
