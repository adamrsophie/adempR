# Living benchmarks -------------------------------------------------------

#' Add methods to an existing benchmark
#'
#' A living benchmark keeps one run directory and adds methods as they are
#' proposed. `sim_add_methods()` takes the study stored in the run, adds the
#' new methods as a recorded amendment, and extends the run: the new
#' methods are fitted to exactly the data sets the earlier methods saw
#' (each method has its own random number stream, so nothing earlier
#' changes), and the methods already present are not refitted. The date of
#' each addition is kept in the manifest and shown by [sim_leaderboard()].
#'
#' @param run A run or a run directory. The current manifest in the
#'   directory is used, so a handle from before an earlier addition is safe.
#' @param methods Named list of new methods (functions or [sim_method()]).
#' @param reason Reason recorded with the amendment.
#' @param ... Passed to [sim_run()].
#' @return The extended run.
#' @export
sim_add_methods <- function(run, methods, reason = NULL, ...) {
  if (inherits(run, "adempr_run")) run <- run$path
  if (!is_scalar_character(run)) stopf("`run` must be a run or a run directory.")
  # the current manifest, never the snapshot an earlier handle holds, so an
  # addition made since is kept
  run <- sim_open(run)
  man <- run$manifest
  study <- man$study
  attr(study, "adempr_stored") <- NULL
  if (!is.list(methods) || !length(methods) || is.null(names(methods))) {
    stopf("`methods` must be a named list.")
  }
  clash <- intersect(names(methods), names(study$methods))
  if (length(clash)) stopf("Method(s) already in the benchmark: %s.", paste(clash, collapse = ", "))
  reason <- reason %||% sprintf("added method(s) %s to the benchmark", paste(names(methods), collapse = ", "))
  new <- sim_amend(study, methods = c(study$methods, methods), reason = reason)
  # sim_run() refuses the change if another addition was committed between
  # reading the manifest and taking the run lock
  attr(new, "adempr_base") <- man$study_fingerprint
  sim_run(new, replications = man$replications, seed = man$seed, path = run$path,
          scenario_replications = man$scenario_replications, ...)
}

#' Leaderboard of a benchmark
#'
#' Ranks the methods in every scenario with [sim_rank()] and summarizes per
#' method: the mean rank over the ranked scenarios, the number of scenarios
#' ranked first (a tie for first credits each tied method with its share),
#' the mean share of first place over bootstrap resamples, the number of
#' scenarios that could not be ranked, and when the method entered the
#' benchmark. A leaderboard is a summary; the scenario level
#' results and their Monte Carlo errors remain the evidence.
#'
#' @param run A run or results.
#' @param measure Ranking measure (see [sim_rank()]).
#' @param estimand Estimand; defaults to the first.
#' @param B Bootstrap resamples for the probability of ranking first.
#' @return A data frame, best mean rank first.
#' @export
sim_leaderboard <- function(run, measure = "mse", estimand = NULL, B = 200L) {
  res <- sim_results(run)
  man <- res$manifest
  r <- sim_rank(res, measure = measure, estimand = estimand, B = B)
  agg <- do.call(rbind, lapply(split(r, r$method), function(d) {
    ok <- d$status == "ok"
    data.frame(method = d$method[1], scenarios = sum(ok), not_ranked = sum(!ok),
               mean_rank = if (any(ok)) mean(d$rank[ok]) else NA_real_,
               first = sum(d$first[ok]),
               mean_p_first = if (any(ok)) mean(d$p_first[ok]) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  added <- vapply(agg$method, function(mn) {
    hit <- Filter(function(rv) {
      if (!is.null(rv$added_methods)) mn %in% rv$added_methods else
        paste("added method", mn) %in% rv$changes
    }, man$revisions)
    if (length(hit)) hit[[1]]$time else man$created
  }, character(1))
  agg$added <- unname(added)
  agg$fingerprint <- vapply(agg$method, function(mn) short_hash(man$fingerprints$methods[[mn]] %||% ""), character(1))
  agg <- agg[order(agg$mean_rank, -agg$first), , drop = FALSE]
  rownames(agg) <- NULL
  agg
}
