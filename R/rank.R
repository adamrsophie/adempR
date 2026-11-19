# Ranking of methods ------------------------------------------------------

RANK_MEASURES <- c("bias", "mae", "mse", "rmse", "empse", "coverage", "width")

rank_value <- function(measure, e, t, lo, up) {
  switch(measure,
         bias = mean(e - t), mae = mean(abs(e - t)), mse = mean((e - t)^2),
         rmse = sqrt(mean((e - t)^2)), empse = stats::sd(e),
         coverage = mean(lo <= t & t <= up), width = mean(up - lo))
}

# Loss used for ranking: smaller is better.
rank_loss <- function(measure, value, nominal) {
  switch(measure, bias = abs(value), coverage = abs(value - nominal), value)
}

#' Rank methods with their uncertainty
#'
#' Ranks the methods in each scenario and estimand by a performance
#' measure, on the replications where every ranked method has a usable
#' result (paired, so the ranking compares methods on the same data sets).
#' The uncertainty of the ranking comes from a paired bootstrap over
#' replications: `p_first` is the proportion of resamples in which the
#' method ranks first, and `mean_rank` its average rank. A ranking whose
#' first place changes often between resamples is not established by the
#' study.
#'
#' Smaller is better for MAE, MSE, RMSE, empirical SE and width; for bias
#' the absolute bias is ranked, and for coverage the distance from the
#' nominal level of the intervals.
#'
#' @param run A run or results.
#' @param measure One of `"bias"`, `"mae"`, `"mse"`, `"rmse"`, `"empse"`,
#'   `"coverage"` or `"width"`.
#' @param methods Methods to rank; defaults to all.
#' @param estimand Estimand; defaults to the first.
#' @param B Number of bootstrap resamples.
#' @param seed Seed of the bootstrap; the caller's random state is left
#'   untouched.
#' @param reference Optional character vector giving a reference order of
#'   the methods (best first). The Kendall distance between the observed
#'   ranking and the reference (proportion of discordant pairs, 0 = same
#'   order, 1 = reversed) is returned per cell, with its bootstrap mean.
#' @return A data frame with one row per scenario and method: `value`,
#'   `loss`, `rank`, `mean_rank`, `p_first`, `n_common`, and, with a
#'   reference, `kendall` and `kendall_boot`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data)),
#'                  median = function(data) sim_estimate("mean", median(data))))
#' run <- sim_run(study, replications = 50, seed = 1, path = tempfile("adempR-run-"))
#' sim_rank(run, "mse", B = 100)
#' @export
sim_rank <- function(run, measure = "mse", methods = NULL, estimand = NULL, B = 500L,
                     seed = 1L, reference = NULL) {
  measure <- match.arg(measure, RANK_MEASURES)
  if (!is_count(B, 10)) stopf("`B` must be a whole number of at least 10.")
  res <- sim_results(run)
  study <- res$manifest$study
  estimand <- estimand %||% names(study$estimands)[1]
  methods <- methods %||% names(study$methods)[vapply(study$methods, function(m) {
    estimand %in% m$estimands
  }, logical(1))]
  if (length(methods) < 2L) stopf("Ranking needs at least two methods.")
  if (!is.null(reference) && !setequal(reference, methods)) {
    stopf("`reference` must order exactly the ranked methods.")
  }
  est <- res$estimates
  est <- est[est$estimand == estimand & est$method %in% methods, , drop = FALSE]
  need_int <- measure %in% c("coverage", "width")
  ok <- est$valid & is.finite(est$estimate)
  if (need_int) ok <- est$interval_status == "valid" & est$estimate_status %in% USABLE_FOR_INTERVAL
  if (measure %in% c("bias", "mae", "mse", "rmse", "coverage")) ok <- ok & !is.na(est$truth)
  est <- est[ok, , drop = FALSE]
  out <- list()
  with_rng_state(rng_base_seed(seed, "rank"), {
    for (key in study$scenarios$keys) {
      sub <- est[est$scenario_key == key, , drop = FALSE]
      common <- Reduce(intersect, lapply(methods, function(m) sub$replication[sub$method == m]))
      common <- sort(common)
      n <- length(common)
      if (n < 2L) next
      get <- function(m, col) {
        r <- sub[sub$method == m, , drop = FALSE]
        r[[col]][match(common, r$replication)]
      }
      cols <- lapply(methods, function(m) list(e = get(m, "estimate"), t = get(m, "truth"),
                                               lo = get(m, "lower"), up = get(m, "upper"),
                                               lv = get(m, "level")))
      names(cols) <- methods
      nominal <- vapply(cols, function(cc) {
        l <- unique(cc$lv[!is.na(cc$lv)])
        if (length(l) == 1L) l else NA_real_
      }, numeric(1))
      if (measure == "coverage" && anyNA(nominal)) next
      losses <- function(idx) {
        vapply(methods, function(m) {
          cc <- cols[[m]]
          rank_loss(measure, rank_value(measure, cc$e[idx], cc$t[idx], cc$lo[idx], cc$up[idx]),
                    nominal[[m]])
        }, numeric(1))
      }
      full <- seq_len(n)
      loss <- losses(full)
      value <- vapply(methods, function(m) {
        cc <- cols[[m]]
        rank_value(measure, cc$e, cc$t, cc$lo, cc$up)
      }, numeric(1))
      rk <- rank(loss, ties.method = "average")
      units <- rep_units(study, common)
      boot <- vapply(seq_len(B), function(b) rank(losses(boot_index(units)),
                                                  ties.method = "average"), numeric(length(methods)))
      boot <- matrix(boot, nrow = length(methods))
      df <- data.frame(scenario_key = key,
                       scenario_id = study$scenarios$labels[match(key, study$scenarios$keys)],
                       estimand = estimand, method = methods, measure = measure,
                       value = unname(value), loss = unname(loss), rank = unname(rk),
                       mean_rank = rowMeans(boot), p_first = rowMeans(boot == 1),
                       n_common = n, stringsAsFactors = FALSE)
      if (!is.null(reference)) {
        df$kendall <- kendall_distance(rk, match(methods, reference))
        df$kendall_boot <- mean(apply(boot, 2, kendall_distance, ref = match(methods, reference)))
      }
      out[[length(out) + 1L]] <- df
    }
  })
  out <- bind_rows(out)
  if (is.null(out)) stopf("No scenario has at least two replications usable by every method.")
  rownames(out) <- NULL
  out
}

# Proportion of discordant pairs between two rankings (ties count half).
kendall_distance <- function(r, ref) {
  k <- length(r)
  d <- 0
  for (i in seq_len(k - 1L)) for (j in (i + 1L):k) {
    s <- sign(r[i] - r[j]) * sign(ref[i] - ref[j])
    d <- d + if (s < 0) 1 else if (s == 0) 0.5 else 0
  }
  d / (k * (k - 1) / 2)
}
