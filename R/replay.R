# Replay of individual tasks ----------------------------------------------

resolve_scenario <- function(study, scenario) {
  sc <- study$scenarios
  if (is.numeric(scenario)) {
    if (!is_count(scenario, 1) || scenario > length(sc)) stopf("Scenario index out of range.")
    return(sc$keys[scenario])
  }
  if (!is_scalar_character(scenario)) stopf("`scenario` must be a key, a label or an index.")
  if (scenario %in% sc$keys) return(scenario)
  if (scenario %in% sc$labels) {
    if (sum(sc$labels == scenario) > 1L) {
      stopf("The label `%s` names several scenarios; use the scenario key or index.", scenario)
    }
    return(sc$keys[match(scenario, sc$labels)])
  }
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
#' @param gpu Device for methods declared with `gpu = TRUE`; defaults to the
#'   device recorded for this replication, if any.
#' @return A list with `scenario` (the full scenario), `data` (the analysis
#'   data, `NULL` when generation failed), `meta` (protected metadata),
#'   `generation` (recorded and replayed generation status and message),
#'   `estimates` (replayed), `committed` (stored rows), `attempts`,
#'   `stored_fits` (retained fits of the current method versions, when
#'   kept, from the round whose estimates are committed), and logical
#'   `matches` (same generation outcome, same data checksum as recorded,
#'   same data as the retained copy when there is one, and the same raw
#'   outputs, attempt statuses, fallback flags and warning flags),
#'   `data_checksum_matches` (replayed data have the checksum recorded with
#'   the replication) and `data_matches` (replayed data identical to the
#'   retained copy). A retained artifact whose checksum does not match
#'   stops with an error.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 5, seed = 1, path = tempfile("adempR-run-"))
#' rp <- sim_replay(run, scenario = 1, replication = 3)
#' rp$matches
#' @export
sim_replay <- function(run, scenario, replication, methods = NULL, gpu = NULL) {
  stopifnot(inherits(run, "adempr_run"))
  man <- run$manifest
  study <- man$study
  key <- resolve_scenario(study, scenario)
  if (!is_count(replication, 1)) stopf("`replication` must be a positive whole number.")
  methods <- methods %||% names(study$methods)
  if (!all(methods %in% names(study$methods))) stopf("Unknown method(s).")
  ctx <- list(seed = man$seed, domain = man$rng_domain, registry = man$rng_registry,
              cache = new.env(parent = emptyenv()), dgp_stage = "dgp",
              method_stage = "method")
  seeds <- task_seeds(ctx, key, replication, study, methods)[[1]]
  idx <- scenario_index(study$scenarios, key)
  scenario_full <- study$scenarios$values[[idx]]
  raw <- read_raw_results(run$path)
  # the device recorded for this replication is used unless another is given
  if (is.null(gpu) && !is.null(raw$attempts$gpu)) {
    a <- raw$attempts
    d <- a$gpu[a$scenario_key == key & a$replication == replication & !is.na(a$gpu)]
    if (length(d)) gpu <- d[1]
  }
  # generation runs once, through the same captured path as the run
  out <- execute_task_seeded(study,
    list(scenario_key = key, replication = replication, seeds = seeds, methods = methods,
         gpu = gpu),
    keep_data = TRUE)
  data <- out$data
  meta <- out$meta
  gen_fp <- man$fingerprints$generator
  recorded <- raw$replications
  if (!is.null(recorded)) {
    recorded <- recorded[recorded$scenario_key == key & recorded$replication == replication &
                           recorded$generator_fp == gen_fp, , drop = FALSE]
  }
  generation <- data.frame(
    recorded = if (!is.null(recorded) && nrow(recorded)) recorded$generation_status[1] else NA_character_,
    replayed = out$replication$generation_status,
    recorded_message = if (!is.null(recorded) && nrow(recorded)) recorded$generation_message[1] else NA_character_,
    replayed_message = out$replication$generation_message, stringsAsFactors = FALSE)
  fps <- unlist(man$fingerprints$methods)
  committed <- raw$estimates
  if (!is.null(committed)) {
    committed <- committed[committed$scenario_key == key &
                             committed$replication == replication &
                             committed$method %in% methods &
                             committed$method_fp == fps[committed$method], , drop = FALSE]
  }
  # retained artifacts are verified against their checksums before use;
  # data come from the newest artifact holding data of this generation, and
  # each method's fit only from the artifact of the round whose estimates
  # are committed, never from an earlier attempt
  data_matches <- NA
  fits <- list()
  arts <- raw$artifacts
  if (!is.null(arts) && nrow(arts)) {
    arts <- arts[arts$scenario_key == key & arts$replication == replication &
                   arts$generator_fp == gen_fp, , drop = FALSE]
    arts <- arts[order(arts$round, decreasing = TRUE), , drop = FALSE]
    fit_round <- if (!is.null(committed) && nrow(committed))
      tapply(committed$round, committed$method, max) else integer()
    for (i in seq_len(nrow(arts))) {
      a <- read_artifact(run$path, arts$artifact[i])
      if (is.null(a)) next
      if (is.na(data_matches) && !is.null(a$data)) data_matches <- identical(a$data, data)
      held <- if (is.null(arts$artifact_fits[i]) || is.na(arts$artifact_fits[i])) character() else
        strsplit(arts$artifact_fits[i], ";", fixed = TRUE)[[1]]
      for (mn in methods) {
        if (is.null(fits[[mn]]) && paste0(mn, "=", fps[[mn]]) %in% held &&
            identical(as.integer(fit_round[mn]), as.integer(arts$round[i]))) {
          fits[[mn]] <- a$fits[[mn]]
        }
      }
    }
  }
  # a match covers the generation outcome, the checksum of the generated
  # data (recorded with every replication, retained or not), the retained
  # copy when there is one, every raw output, the attempt status, the
  # fallback flag and the warning flag (which decides validity under a
  # warnings policy); messages and timings are compared only through the
  # status
  rec_md5 <- if (!is.null(recorded) && nrow(recorded) && !is.null(recorded$data_md5))
    recorded$data_md5[1] else NA_character_
  data_checksum_matches <- if (is.na(rec_md5)) NA else
    identical(rec_md5, out$replication$data_md5)
  matches <- NA
  if (!is.na(generation$recorded)) {
    matches <- identical(generation$recorded, generation$replayed) &&
      !isFALSE(data_checksum_matches) && !isFALSE(data_matches)
    if (matches && !is.null(committed) && nrow(committed)) {
      cols <- c("method", "estimand", RAW_ESTIMATE_COLUMNS, "attempt_status", "had_warning",
                "fallback")
      sig <- function(df) {
        if (is.null(df)) return(character())
        for (cl in setdiff(cols, names(df))) df[[cl]] <- NA
        df <- df[order(df$method, df$estimand), cols, drop = FALSE]
        if (!is.null(df$fallback)) df$fallback <- df$fallback %in% TRUE
        do.call(paste, c(lapply(df, col_text), sep = "\r"))
      }
      matches <- identical(sig(committed), sig(out$estimates))
    }
  }
  list(scenario = scenario_full, data = data, meta = meta, generation = generation,
       estimates = out$estimates, attempts = out$attempts, committed = committed,
       stored_fits = fits, matches = matches, data_matches = data_matches,
       data_checksum_matches = data_checksum_matches)
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
  if (inherits(run, "adempr_run")) {
    # for a run, what was actually retained decides, not the current policy
    # only artifacts of the current generator and method versions whose
    # files are present count, and data only where an artifact holds data
    raw <- read_raw_results(run$path)
    man <- run$manifest
    gen_fp <- man$fingerprints$generator
    reps <- raw$replications
    if (!is.null(reps)) reps <- reps[reps$generator_fp == gen_fp & reps$generation_status == "ok", , drop = FALSE]
    ok <- if (is.null(reps)) character() else unique(paste(reps$scenario_key, reps$replication))
    arts <- raw$artifacts
    if (!is.null(arts)) {
      arts <- arts[arts$generator_fp == gen_fp & file.exists(file.path(run$path, arts$artifact)), , drop = FALSE]
    }
    has_data <- if (is.null(arts) || !nrow(arts)) logical() else if (!is.null(arts$artifact_data) &&
                                                                       !anyNA(arts$artifact_data)) {
      arts$artifact_data
    } else {
      vapply(arts$artifact, function(a) !is.null(read_artifact(run$path, a)$data), logical(1))
    }
    have <- if (!length(has_data)) character() else
      unique(paste(arts$scenario_key, arts$replication)[has_data])
    current <- paste0(names(man$fingerprints$methods), "=", unlist(man$fingerprints$methods))
    held <- if (is.null(arts) || is.null(arts$artifact_fits)) character() else
      unlist(strsplit(arts$artifact_fits[!is.na(arts$artifact_fits)], ";", fixed = TRUE))
    r$data <- length(ok) > 0L && all(ok %in% have)
    r$fits <- any(held %in% current)
  }
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
