# The one-pass task builder must give exactly the tasks, seeds and order of
# the straightforward per-replication scan it replaced.
old_build_tasks <- function(study, ctx, todo, data_md5 = NULL) {
  out <- list()
  for (key in unique(todo$scenario_key)) {
    sub <- todo[todo$scenario_key == key, , drop = FALSE]
    reps <- sort(unique(sub$replication))
    seeds <- task_seeds(ctx, key, reps, study, names(study$methods))
    for (j in seq_along(reps)) {
      ms <- sub$method[sub$replication == reps[j]]
      ms <- names(study$methods)[names(study$methods) %in% ms]
      out[[length(out) + 1L]] <- list(scenario_key = key, replication = reps[j],
                                      seeds = seeds[[j]], methods = ms,
                                      expect_md5 = unname(data_md5[paste(key, reps[j], sep = "\r")]))
    }
  }
  out
}

test_that("tasks are identical to the per-replication scan", {
  mk <- function(methods, n = c(3L, 5L, 8L)) {
    sim_study("plan", scenarios = sim_scenarios(n = n),
      generate = function(scenario, context) rnorm(scenario$n),
      estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
      methods = methods)
  }
  m <- list(a = function(data) sim_estimate("mean", mean(data)),
            b = function(data) sim_estimate("mean", median(data)),
            c = function(data) sim_estimate("mean", data[1]))
  ctx <- function(s) list(seed = 7, domain = NULL, registry = rng_registry(rng_keys(s)),
                          cache = new.env(parent = emptyenv()), dgp_stage = "dgp",
                          method_stage = "method")
  s <- mk(m)
  # fresh run
  todo <- plan_tasks(s, 12L, NULL)
  expect_identical(build_tasks(s, ctx(s), todo), old_build_tasks(s, ctx(s), todo))
  # resume: some cells done, in scattered order, with data hashes
  set.seed(3)
  all <- plan_tasks(s, 12L, NULL)
  done <- all[sample(nrow(all), 40), ]
  todo <- plan_tasks(s, 12L, done)
  md5 <- stats::setNames(sprintf("h%d", seq_len(nrow(done))),
                         paste(done$scenario_key, done$replication, sep = "\r"))
  expect_identical(build_tasks(s, ctx(s), todo, md5), old_build_tasks(s, ctx(s), todo, md5))
  # a method added to a complete run
  s2 <- mk(c(m, list(d = function(data) sim_estimate("mean", max(data)))))
  done2 <- all
  todo <- plan_tasks(s2, 12L, done2)
  expect_identical(build_tasks(s2, ctx(s2), todo), old_build_tasks(s2, ctx(s2), todo))
  # reordered scenarios and unsorted pending rows
  s3 <- mk(m, n = c(8L, 3L, 5L))
  todo <- plan_tasks(s3, 6L, NULL)
  todo <- todo[sample(nrow(todo)), ]
  expect_identical(build_tasks(s3, ctx(s3), todo), old_build_tasks(s3, ctx(s3), todo))
  expect_identical(build_tasks(s, ctx(s), todo[0, ]), list())
})
