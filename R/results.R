# Reduction of committed results -------------------------------------------

col_text <- function(x) {
  if (is.double(x)) {
    out <- formatC(x, digits = 17, format = "g")
    out[is.na(x)] <- ifelse(is.nan(x[is.na(x)]), "NaN", "NA")
    return(out)
  }
  out <- as.character(x)
  out[is.na(out)] <- "NA"
  out
}

# Remove exact duplicates (repeated delivery of the same task). Duplicates
# whose scientific content differs are a conflict and stop with an error;
# volatile columns (timings, artifact paths) are not compared.
dedup_table <- function(df, keys, scientific, what) {
  if (is.null(df) || !nrow(df)) return(df)
  id <- do.call(paste, c(lapply(df[keys], col_text), sep = "\r"))
  dup <- duplicated(id)
  if (!any(dup)) return(df)
  sci <- intersect(scientific, names(df))
  sig <- do.call(paste, c(lapply(df[sci], col_text), sep = "\r"))
  first_sig <- sig[match(id, id)]
  conflict <- dup & sig != first_sig
  if (any(conflict)) {
    bad <- df[conflict, keys, drop = FALSE][1, , drop = FALSE]
    stopf("Conflicting duplicate %s results: the same task (%s) was committed twice with different scientific content. The run directory is inconsistent; inspect the batches before continuing.",
          what, paste(names(bad), vapply(bad, col_text, character(1)), sep = "=", collapse = ", "))
  }
  df[!dup, , drop = FALSE]
}

# Keep, for every task identified by `keys`, only the rows of its latest
# retry round.
latest_round <- function(df, keys) {
  if (is.null(df) || !nrow(df)) return(df)
  id <- do.call(paste, c(lapply(df[keys], col_text), sep = "\r"))
  top <- tapply(df$round, id, max)
  df[df$round == top[id], , drop = FALSE]
}

read_raw_results <- function(path, repair = FALSE) {
  p <- run_paths(path)
  if (!dir.exists(p$batches)) {
    return(list(replications = NULL, attempts = NULL, estimates = NULL,
                dataset_truth = NULL, problems = character()))
  }
  rc <- read_committed(p$batches, p$quarantine, repair)
  objs <- rc$objects
  for (o in objs) {
    if (!identical(o$schema_version, SCHEMA_VERSION)) {
      stopf("A committed batch uses result schema %s; expected %s.", o$schema_version,
            SCHEMA_VERSION)
    }
  }
  # Every row carries the retry round in which it was computed (0 for the
  # first execution). Within a round, a repeated delivery must be identical;
  # across rounds, the latest round supersedes earlier ones.
  get <- function(nm) {
    df <- bind_rows(lapply(objs, `[[`, nm))
    if (!is.null(df) && nrow(df) && is.null(df$round)) df$round <- 0L
    df
  }
  # every artifact ever committed, before duplicates of a generation are
  # removed: a later batch (a new or updated method) has its own artifact
  all_reps <- get("replications")
  arts <- NULL
  if (!is.null(all_reps) && "artifact" %in% names(all_reps)) {
    arts <- all_reps[!is.na(all_reps$artifact), intersect(
      c("scenario_key", "replication", "generator_fp", "round", "artifact", "artifact_fits"),
      names(all_reps)), drop = FALSE]
    rownames(arts) <- NULL
  }
  reps <- latest_round(dedup_table(get("replications"),
                                   c("scenario_key", "replication", "generator_fp", "round"),
                                   c("generation_status", "generation_message", "data_md5"),
                                   "generation"),
                       c("scenario_key", "replication", "generator_fp"))
  task_keys <- c("scenario_key", "replication", "method", "method_fp")
  att <- latest_round(dedup_table(get("attempts"), c(task_keys, "attempt", "round"),
                                  c("status", "failed_stage", "message", "n_warnings", "warnings"),
                                  "method attempt"), task_keys)
  est <- latest_round(dedup_table(get("estimates"), c(task_keys, "estimand", "round"),
                                  c(RAW_ESTIMATE_COLUMNS, "attempt_status", "had_warning",
                                    "fallback"),
                                  "estimate"), task_keys)
  dtr <- latest_round(dedup_table(get("dataset_truth"),
                                  c("scenario_key", "replication", "estimand", "round"),
                                  "truth", "data set truth"),
                      c("scenario_key", "replication", "estimand"))
  dgn <- latest_round(dedup_table(get("diagnostics"), c(task_keys, "name", "round"), "value",
                                  "fit diagnostic"), task_keys)
  rcd <- latest_round(dedup_table(get("records"), c("scenario_key", "replication", "name", "round"),
                                  "value", "data set record"),
                      c("scenario_key", "replication", "name"))
  list(replications = reps, attempts = att, estimates = est, dataset_truth = dtr,
       diagnostics = dgn, records = rcd, artifacts = arts, problems = rc$problems)
}

# Read a retained artifact after checking it against its checksum sidecar;
# `NULL` when it is missing (a missing copy is reported, a changed one stops).
read_artifact <- function(root, rel) {
  f <- file.path(root, rel)
  side <- sub("\\.rds$", ".md5", f)
  if (!file.exists(f) || !file.exists(side)) return(NULL)
  if (!identical(readLines(side, warn = FALSE)[1], file_md5(f))) {
    stopf("Retained artifact %s does not match its checksum; it was changed or damaged.", rel)
  }
  readRDS(f)
}

#' Compute or load the truth of every estimand
#'
#' Exact and numerical truths are computed directly. Monte Carlo truths are
#' computed in their own stage with a registered random number stream and
#' cached in the run directory under their fingerprint, so they are
#' recomputed only when the truth definition or the scenario changes.
#' Data set specific truths come from the generator.
#'
#' @param run A run from [sim_run()] or [sim_open()].
#' @return A data frame with one row per scenario and estimand.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L, mu = c(0, 1)),
#'   generate = function(scenario, context) rnorm(scenario$n, scenario$mu),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) scenario$mu))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 5, seed = 1, path = tempfile("adempR-run-"))
#' sim_truth(run)
#' @export
sim_truth <- function(run) {
  stopifnot(inherits(run, "adempr_run"))
  man <- run$manifest
  study <- man$study
  p <- run_paths(run$path)
  ensure_dirs(p)
  rows <- list()
  base <- NULL
  cache <- new.env(parent = emptyenv())
  for (en in names(study$estimands)) {
    tr <- study$estimands[[en]]$truth
    tfp <- truth_fingerprint(tr)
    for (i in seq_along(study$scenarios$values)) {
      key <- study$scenarios$keys[i]
      sc <- study$scenarios$values[[i]]
      row <- NULL
      if (tr$kind == "mc") {
        tk <- truth_key(tr, sc, key)
        nm <- paste0("t-", md5_string(c(tfp, tk, en, man$seed)))
        f <- file.path(p$truth, paste0(nm, ".rds"))
        sidecar <- file.path(p$truth, paste0(nm, ".md5"))
        if (file.exists(f) && file.exists(sidecar) &&
            identical(file_md5(f), readLines(sidecar, warn = FALSE)[1])) {
          row <- readRDS(f)
        } else {
          if (is.null(base)) base <- rng_base_seed(man$seed, man$rng_domain)
          k <- rng_lookup(man$rng_registry, "truth", tk, en)
          seed <- rng_substreams(rng_stream_seeds(base, k, cache)[[1]], 1L)[[1]]
          t0 <- proc.time()[["elapsed"]]
          row <- compute_truth_one(tr, sc, seed)
          row$truth_elapsed <- proc.time()[["elapsed"]] - t0
          # another reader may have published the value meanwhile: the first
          # committed value wins, so every reader reports the same row
          if (file.exists(f) && file.exists(sidecar) &&
              identical(file_md5(f), readLines(sidecar, warn = FALSE)[1])) {
            row <- readRDS(f)
          } else {
            atomic_save(row, p$truth, nm, p$tmp)
          }
        }
      } else {
        row <- compute_truth_one(tr, sc)
        row$truth_elapsed <- NA_real_
      }
      row$scenario_key <- key
      row$estimand <- en
      row$truth_fp <- tfp
      rows[[length(rows) + 1L]] <- row
    }
  }
  out <- bind_rows(rows)
  out[c("scenario_key", "estimand", setdiff(names(out), c("scenario_key", "estimand")))]
}

#' Results of a run as validated tables
#'
#' Reads the committed batches, removes duplicate deliveries, detects
#' conflicting duplicates, keeps the results that belong to the current
#' definition of the study (current method fingerprints, current scenarios,
#' replications up to the planned number), joins the truth and classifies
#' the validity of every output under the study's failure policy.
#'
#' @param run A run from [sim_run()] or [sim_open()].
#' @return An object of class `adempr_results`: a list of data frames
#'   `scenarios`, `replications`, `attempts`, `estimates` and `truth`,
#'   plus the `manifest`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 5, seed = 1, path = tempfile("adempR-run-"))
#' res <- sim_results(run)
#' head(res$estimates)
#' @export
sim_results <- function(run) {
  if (inherits(run, "adempr_results")) return(run)
  stopifnot(inherits(run, "adempr_run"))
  man <- run$manifest
  study <- man$study
  raw <- read_raw_results(run$path)
  keys <- study$scenarios$keys
  R <- man$replications
  # A sequential run is analyzed at the replication count of its stopping
  # look, even when later replications are on disk.
  if (!is.null(man$sequential$analyzed)) R <- min(R, man$sequential$analyzed)
  # Results are selected by the fingerprints stored when they were run, never
  # by fingerprints recomputed in the current session.
  fps <- unlist(man$fingerprints$methods)
  gen_fp <- man$fingerprints$generator
  keep_rep <- function(df) {
    if (is.null(df)) return(df)
    df[df$scenario_key %in% keys & df$replication <= R, , drop = FALSE]
  }
  keep_meth <- function(df) {
    if (is.null(df)) return(df)
    df <- keep_rep(df)
    df[df$method %in% names(fps) & df$method_fp == fps[df$method], , drop = FALSE]
  }
  reps <- keep_rep(raw$replications)
  if (!is.null(reps)) reps <- reps[reps$generator_fp == gen_fp, , drop = FALSE]
  att <- keep_meth(raw$attempts)
  est <- keep_meth(raw$estimates)
  dgn <- keep_meth(raw$diagnostics)
  truth <- sim_truth(run)
  if (!is.null(est) && nrow(est)) {
    common <- truth[truth$truth_status == "ok", c("scenario_key", "estimand", "truth",
                                                    "truth_se", "truth_kind")]
    est <- merge(est, common, by = c("scenario_key", "estimand"), all.x = TRUE, sort = FALSE)
    dtr <- keep_rep(raw$dataset_truth)
    if (!is.null(dtr) && nrow(dtr)) {
      ds_est <- names(study$estimands)[vapply(study$estimands, function(e) {
        e$truth$kind == "dataset"
      }, logical(1))]
      idx <- match(paste(est$scenario_key, est$replication, est$estimand),
                   paste(dtr$scenario_key, dtr$replication, dtr$estimand))
      use <- est$estimand %in% ds_est & !is.na(idx)
      est$truth[use] <- dtr$truth[idx[use]]
      est$truth_se[use] <- 0
      est$truth_kind[use] <- "dataset"
    }
    est$truth_kind[is.na(est$truth_kind)] <- vapply(est$estimand[is.na(est$truth_kind)],
      function(id) study$estimands[[id]]$truth$kind, character(1))
    est <- classify_rows(est, study$estimands, study$failure)
    est <- est[order(match(est$scenario_key, keys), est$method, est$estimand,
                     est$replication), , drop = FALSE]
    rownames(est) <- NULL
  } else {
    est <- empty_estimates(study)
  }
  sc <- as.data.frame(study$scenarios)
  structure(list(scenarios = sc, replications = reps, attempts = att, estimates = est,
                 truth = truth, manifest = man, path = run$path,
                 diagnostics = dgn, records = keep_rep(raw$records),
                 artifacts = raw$artifacts, problems = raw$problems),
            class = "adempr_results")
}

# A typed estimates table with no rows, so that summaries, comparisons,
# pilots and reports of a run without any estimate (for example when every
# generation failed) report unavailable cells instead of failing.
empty_estimates <- function(study) {
  en <- names(study$estimands)[1]
  row <- cbind(data.frame(scenario_key = NA_character_, replication = NA_integer_,
                          method = NA_character_, stringsAsFactors = FALSE),
               missing_row(en))
  row$attempt_status <- "ok"
  row$had_warning <- FALSE
  row$fallback <- FALSE
  row$method_fp <- NA_character_
  row$round <- 0L
  row$truth <- NA_real_
  row$truth_se <- NA_real_
  row$truth_kind <- NA_character_
  classify_rows(row, study$estimands, study$failure)[0, , drop = FALSE]
}

#' @export
print.adempr_results <- function(x, ...) {
  cat(sprintf("<adempR results: %d data set(s), %d estimate row(s)>\n",
              NROW(x$replications), NROW(x$estimates)))
  invisible(x)
}
