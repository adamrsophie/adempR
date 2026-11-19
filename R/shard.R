# Cluster schedulers: shards and collection -------------------------------

#' Run a study as independent jobs on a cluster scheduler
#'
#' For a scheduler with array jobs (Slurm, PBS, SGE, LSF, through
#' \pkg{batchtools}, \pkg{targets} or a plain job script), the scenarios are
#' split into `n_shards` shards. Job `shard` runs `sim_run_shard()`, which
#' writes a run directory of its own (`<path>/shard-<k>`), with the stream
#' registry of the full study, so each scenario draws exactly the streams it
#' would draw in one run. No job writes into another job's directory, so no
#' locking across nodes is needed. After all jobs finish,
#' `sim_collect_shards()` assembles a run directory for the full study from
#' the committed batches of every shard, verifying each against its
#' checksum. The result is identical to a single [sim_run()] with the same
#' seed. A failed or unfinished job is rerun by calling `sim_run_shard()`
#' again for that shard; it resumes.
#'
#' With a scheduler that \pkg{future} supports (for example through
#' \pkg{future.batchtools}), `sim_run(backend = "future")` is an
#' alternative in which one coordinating R session commits all results.
#'
#' @param study A study.
#' @param shard Shard number, from 1 to `n_shards` (for example
#'   `as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))`).
#' @param n_shards Number of shards.
#' @param replications Replications per scenario.
#' @param seed Seed.
#' @param path Parent directory of the shard directories.
#' @param ... Passed to [sim_run()] (not `rng_registry`).
#' @return `sim_run_shard()`: the shard's run. `sim_shard_scenarios()`: the
#'   scenario keys of a shard. `sim_collect_shards()`: the collected run.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 20L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data))))
#' parent <- tempfile("adempR-shards-")
#' for (k in 1:2) sim_run_shard(study, k, 2, replications = 10, seed = 1, path = parent)
#' run <- sim_collect_shards(study, 2, replications = 10, seed = 1, path = parent,
#'                           into = tempfile("adempR-run-"))
#' run
#' @export
sim_run_shard <- function(study, shard, n_shards, replications, seed, path, ...) {
  stopifnot(inherits(study, "adempr_study"))
  keys <- sim_shard_scenarios(study, shard, n_shards)
  sub <- study
  idx <- match(keys, study$scenarios$keys)
  sc <- study$scenarios
  sc$values <- sc$values[idx]
  sc$keys <- sc$keys[idx]
  sc$labels <- sc$labels[idx]
  if (!is.null(sc$weights)) sc$weights <- sc$weights[idx] / sum(sc$weights[idx])
  sub$scenarios <- sc
  sim_run(sub, replications = replications, seed = seed, path = shard_path(path, shard),
          rng_registry = rng_registry(rng_keys(study)), ...)
}

#' @rdname sim_run_shard
#' @export
sim_shard_scenarios <- function(study, shard, n_shards) {
  if (!is_count(n_shards, 1)) stopf("`n_shards` must be a positive whole number.")
  if (!is_count(shard, 1) || shard > n_shards) stopf("`shard` must be between 1 and `n_shards`.")
  keys <- study$scenarios$keys
  if (n_shards > length(keys)) stopf("There are only %d scenarios; use at most that many shards.", length(keys))
  keys[(seq_along(keys) - 1L) %% n_shards + 1L == shard]
}

#' @rdname sim_run_shard
#' @param into Directory of the collected run (new).
#' @export
sim_collect_shards <- function(study, n_shards, replications, seed, path, into) {
  stopifnot(inherits(study, "adempr_study"))
  if (dir.exists(into) && length(list.files(into, all.files = TRUE, no.. = TRUE))) {
    stopf("%s is not empty.", into)
  }
  dirs <- vapply(seq_len(n_shards), function(k) shard_path(path, k), character(1))
  missing <- dirs[!file.exists(file.path(dirs, "manifest.rds"))]
  if (length(missing)) stopf("Shard(s) not run yet: %s.", paste(basename(missing), collapse = ", "))
  registry <- rng_registry(rng_keys(study))
  for (d in dirs) {
    m <- read_manifest(d)
    if (!identical(m$seed, seed)) stopf("Shard %s used seed %s, not %s.", basename(d), m$seed, seed)
    if (!identical(m$fingerprints$generator, generator_fingerprint(study$generator)) ||
        !identical(unlist(m$fingerprints$methods), unlist(lapply(study$methods, method_fingerprint)))) {
      stopf("Shard %s was run with a different generator or methods.", basename(d))
    }
    if (m$replications < replications) {
      stopf("Shard %s holds %d replications, fewer than %d.", basename(d), m$replications, replications)
    }
  }
  p <- run_paths(into)
  ensure_dirs(p)
  man <- new_manifest(study, seed, replications, "collected", registry)
  man$collected_from <- basename(dirs)
  write_manifest(p, man)
  inv <- character()
  for (d in dirs) {
    rc <- read_committed(file.path(d, "batches"), file.path(d, "quarantine"))
    for (nm in names(rc$objects)) {
      b <- sub("\\.rds$", "", nm)
      file.copy(file.path(d, "batches", c(nm, paste0(b, ".md5"))), p$batches)
      inv <- c(inv, b)
    }
    art <- file.path(d, "artifacts")
    if (dir.exists(art)) {
      dir.create(file.path(into, "artifacts"), showWarnings = FALSE)
      file.copy(list.files(art, full.names = TRUE), file.path(into, "artifacts"))
    }
  }
  writeLines(inv, file.path(p$batches, "INVENTORY"))
  # a final call computes truths, records the session and finds nothing to run
  sim_run(study, replications = replications, seed = seed, path = into, progress = FALSE)
}

shard_path <- function(path, shard) file.path(path, sprintf("shard-%03d", as.integer(shard)))
