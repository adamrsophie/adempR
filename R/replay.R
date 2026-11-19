# Replay of individual tasks ----------------------------------------------

resolve_scenario <- function(study, scenario) {
  sc <- study$scenarios
  if (is.numeric(scenario)) {
    if (!is_count(scenario, 1) || scenario > length(sc)) stopf("Scenario index out of range.")
    return(sc$keys[scenario])
  }
  if (!is_scalar_character(scenario)) stopf("`scenario` must be a key, a label or an index.")
  if (scenario %in% sc$keys) return(scenario)
  if (scenario %in% sc$labels) return(sc$keys[match(scenario, sc$labels)])
  hit <- which(startsWith(sc$keys, scenario))
  if (length(hit) == 1L) return(sc$keys[hit])
  stopf("Unknown scenario `%s`.", scenario)
}

#' Replay one replication
#'
#' Regenerates the data set of one replication from its recorded random
#' number state and reruns the requested methods with their own recorded
#' states, without writing anything to the run directory. The replayed
#' estimates are compared with the committed ones (timings are excluded,
#' since they are inherently variable). If the generated data were retained
#' ([sim_retention()]), the regenerated data are also compared with the
#' stored copy.
#'
#' Replay is exact only under the same study definition, R version,
#' package versions and numerical libraries; external random number
#' generators and multithreaded native code may need additional control.
#'
#' @param run A run from [sim_run()] or [sim_open()].
#' @param scenario Scenario key, label or index.
#' @param replication Replication number.
#' @param methods Methods to rerun; defaults to all methods.
#' @return A list with `scenario` (the full scenario), `data` (the analysis
#'   data), `meta` (protected metadata), `estimates` (replayed), `committed`
#'   (stored rows), `attempts`, and logical `matches` and `data_matches`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 5, seed = 1, path = tempfile("adempR-run-"))
#' rp <- sim_replay(run, scenario = 1, replication = 3)
#' rp$matches
#' @export
sim_replay <- function(run, scenario, replication, methods = NULL) {
  stopifnot(inherits(run, "adempr_run"))
  man <- run$manifest
  study <- man$study
  key <- resolve_scenario(study, scenario)
  if (!is_count(replication, 1)) stopf("`replication` must be a positive whole number.")
  methods <- methods %||% names(study$methods)
  if (!all(methods %in% names(study$methods))) stopf("Unknown method(s).")
  ctx <- list(seed = man$seed, registry = man$rng_registry,
              cache = new.env(parent = emptyenv()), dgp_stage = "dgp",
              method_stage = "method")
  seeds <- task_seeds(ctx, key, replication, study, methods)[[1]]
  idx <- scenario_index(study$scenarios, key)
  scenario_full <- study$scenarios$values[[idx]]
  gen <- with_rng_state(seeds$data, study$generator$fun(
    scenario_full, list(replication = as.integer(replication), scenario_key = key, stage = "dgp")))
  data <- if (inherits(gen, "adempr_data")) gen$data else gen
  meta <- if (inherits(gen, "adempr_data")) gen[c("truth", "meta")] else list()
  out <- execute_task_seeded(study, list(scenario_key = key, replication = replication,
                                         seeds = seeds, methods = methods))
  raw <- read_raw_results(run$path)
  fps <- vapply(study$methods, method_fingerprint, character(1))
  committed <- raw$estimates
  if (!is.null(committed)) {
    committed <- committed[committed$scenario_key == key &
                             committed$replication == replication &
                             committed$method %in% methods &
                             committed$method_fp == fps[committed$method], , drop = FALSE]
  }
  matches <- NA
  if (!is.null(committed) && nrow(committed) && !is.null(out$estimates)) {
    cols <- c("method", "estimand", RAW_ESTIMATE_COLUMNS, "attempt_status")
    sig <- function(df) {
      df <- df[order(df$method, df$estimand), cols, drop = FALSE]
      do.call(paste, c(lapply(df, col_text), sep = "\r"))
    }
    matches <- identical(sig(committed), sig(out$estimates))
  }
  data_matches <- NA
  rep_row <- raw$replications
  if (!is.null(rep_row)) {
    rep_row <- rep_row[rep_row$scenario_key == key & rep_row$replication == replication &
                         !is.na(rep_row$artifact), , drop = FALSE]
    if (nrow(rep_row)) {
      f <- file.path(run$path, rep_row$artifact[1])
      if (file.exists(f)) data_matches <- identical(readRDS(f)$data, data)
    }
  }
  list(scenario = scenario_full, data = data, meta = meta, estimates = out$estimates,
       attempts = out$attempts, committed = committed, matches = matches,
       data_matches = data_matches)
}

#' Analyses available under the retention policy
#'
#' Reports which analyses remain possible for a run, given what was
#' retained. Estimates, intervals, tests, statuses and timings are always
#' kept, so all performance measures, new truths, new measures, paired
#' comparisons and bootstrap comparisons can be recomputed. Inspecting
#' fitted objects requires `fits = TRUE`; replay against stored data
#' requires `data = TRUE` (otherwise replay regenerates the data, which is
#' exact only under the same software environment).
#'
#' @param run A run or a study.
#' @return A data frame with columns `analysis`, `available` and `reason`.
#' @examples
#' sim_capabilities(sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data)))))
#' @export
sim_capabilities <- function(run) {
  study <- if (inherits(run, "adempr_run")) run$manifest$study else run
  stopifnot(inherits(study, "adempr_study"))
  r <- study$retention
  rows <- list(
    c("performance measures (built in and custom per replication)", TRUE, "estimate rows retained"),
    c("re-evaluation against a changed truth", TRUE, "estimate rows retained"),
    c("paired comparisons and paired bootstrap", TRUE, "replication identity retained"),
    c("failure and validity reclassification under a new policy", TRUE,
      "raw outputs, statuses and warning flags retained"),
    c("runtime distributions", TRUE, "attempt timings retained"),
    c("replay by regeneration", TRUE, "random number states recorded in the registry"),
    c("replay checked against stored data", r$data,
      if (r$data) "generated data retained" else "set sim_retention(data = TRUE)"),
    c("DGP diagnostics on the analyzed data sets", r$data,
      if (r$data) "generated data retained" else "regenerate with sim_replay() or sim_check_dgp()"),
    c("inspection of fitted model objects", r$fits,
      if (r$fits) "raw method outputs retained" else "set sim_retention(fits = TRUE)"),
    c("new estimands from existing fits", r$fits,
      if (r$fits) "refit extraction on retained outputs" else "requires refitting")
  )
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- c("analysis", "available", "reason")
  df$available <- as.logical(df$available)
  df
}
