# Performance summaries ---------------------------------------------------

USABLE_FOR_INTERVAL <- c("valid", "infinite", "missing")

# Reject when p <= alpha (Morris, White and Crowther 2019, Table 6); the
# equality matters for discrete tests.
reject_indicator <- function(rows, alpha) {
  a <- ifelse(is.na(rows$alpha), alpha, rows$alpha)
  ifelse(!is.na(rows$reject), rows$reject, rows$p_value <= a)
}

rejection_label <- function(rows, estimand) {
  h0 <- ifelse(is.na(rows$null), estimand$null %||% NA_real_, rows$null)
  if (all(is.na(rows$truth)) || all(is.na(h0))) return("rejection")
  alt <- ifelse(is.na(rows$alternative), "two.sided", rows$alternative)
  # A Monte Carlo truth within 3 standard errors of the null cannot tell
  # whether the null holds: use the neutral label.
  tse <- rows$truth_se
  if (any(!is.na(tse) & tse > 0 & abs(rows$truth - h0) <= 3 * tse)) return("rejection")
  # likewise a numerical truth whose error bound reaches the null
  ter <- rows$truth_error
  if (!is.null(ter) && any(!is.na(ter) & ter > 0 & abs(rows$truth - h0) <= ter)) return("rejection")
  nt <- null_true(rows$truth, h0, alt)
  if (anyNA(nt)) return("rejection")
  if (all(nt)) "type_i_error" else if (!any(nt)) "power" else "rejection"
}

# One measure on the rows of one cell. `shift` moves the truth (used for
# the truth uncertainty sensitivity analysis).
measure_cell <- function(m, rows, study, estimand, shift = 0) {
  level <- study$level
  truth <- rows$truth + shift
  usable <- rows$valid
  has_truth <- !is.na(truth)
  se_ok <- rows$se_status == "valid"
  int_ok <- rows$estimate_status %in% USABLE_FOR_INTERVAL & rows$interval_status == "valid"
  test_ok <- rows$estimate_status %in% USABLE_FOR_INTERVAL & rows$test_status == "valid"
  e <- rows$estimate - truth
  # A data set truth (declared, not inferred from the values seen) makes
  # the precision measures use the errors, whose spread excludes the
  # spread of the truth.
  varying_truth <- identical(estimand$truth$kind, "dataset")
  if (inherits(m, "adempr_measure")) {
    elig <- rep(TRUE, nrow(rows))
    if ("estimate" %in% m$requires) elig <- elig & usable
    if ("truth" %in% m$requires) elig <- elig & has_truth
    if ("se" %in% m$requires) elig <- elig & se_ok
    if ("interval" %in% m$requires) elig <- elig & int_ok
    if ("test" %in% m$requires) elig <- elig & test_ok
    sub <- rows[elig, , drop = FALSE]
    sub$truth <- truth[elig]
    return(k_custom(m, sub, level))
  }
  switch(m,
    bias = k_mean(e[usable & has_truth], level),
    mae = k_mean(abs(e[usable & has_truth]), level),
    mse = k_mean(e[usable & has_truth]^2, level),
    rmse = k_rmse(e[usable & has_truth], level),
    relbias = {
      sel <- usable & has_truth
      t <- truth[sel]
      if (!length(t)) return(k_mean(numeric(), level))
      if (!all(is.finite(e[sel]))) {
        return(kernel_result(n = length(t), status = "infinite",
                             note = "nonfinite errors"))
      }
      # Work on a common scale so that sums and spreads of values near the
      # largest double cannot overflow.
      sc <- max(abs(c(t, e[sel], rows$estimate[sel])))
      if (!is.finite(sc) || sc == 0) sc <- 1
      ts <- t / sc
      es <- e[sel] / sc
      # A truth that is small next to the spread of the estimates makes the
      # ratio unstable: it is reported, but flagged.
      spread <- if (length(t) > 1L) stats::sd(rows$estimate[sel] / sc) else 0
      unstable <- function(r, ref) {
        if (abs(ref) < spread) {
          r$status <- "unstable"
          r$note <- paste0(r$note, "; the truth is smaller than the empirical SE, so the ratio is unstable")
        }
        r
      }
      # only an exact zero makes the ratio undefined; a tiny truth gives a
      # value flagged as unstable
      if (all(t == t[1])) {
        if (t[1] == 0) {
          kernel_result(n = length(t), status = "undefined",
                        note = "relative bias is undefined at zero truth")
        } else {
          r <- k_mean(es / ts[1], level)
          r$note <- "bias divided by the (constant) truth"
          unstable(r, ts[1])
        }
      } else {
        mt <- mean(ts)
        if (mt == 0) {
          return(kernel_result(n = length(t), status = "undefined",
                               note = "mean truth is zero"))
        }
        val <- sum(es) / sum(ts)
        lin <- (es - val * ts) / mt
        mc <- if (length(t) > 1L) stats::sd(lin) / sqrt(length(t)) else NA_real_
        z <- stats::qnorm(1 - (1 - level) / 2)
        unstable(kernel_result(value = val, mcse = mc, lower = val - z * mc, upper = val + z * mc,
                               n = length(t),
                               note = "varying truth: ratio of averages, sum(error) / sum(truth), delta method MCSE"),
                 mt)
      }
    },
    # With a truth that varies between data sets, the spread of the
    # estimates includes the spread of the truth; the errors are used.
    empse = if (varying_truth) {
      r <- k_empse(e[usable & has_truth], level)
      r$note <- paste(c("SD of the errors (data set truth)",
                        if (!is.na(r$note)) r$note), collapse = "; ")
      r
    } else k_empse(rows$estimate[usable], level),
    mean_se = k_mean(rows$se[usable & se_ok], level),
    rmvar = k_rmvar(rows$se[usable & se_ok], level),
    se_ratio = if (varying_truth) {
      k_se_ratio(e[usable & se_ok & has_truth], rows$se[usable & se_ok & has_truth], level)
    } else k_se_ratio(rows$estimate[usable & se_ok], rows$se[usable & se_ok], level),
    coverage = {
      sel <- int_ok & has_truth
      r <- k_prop(rows$lower[sel] <= truth[sel] & truth[sel] <= rows$upper[sel], level)
      lv <- sort(unique(rows$level[sel][!is.na(rows$level[sel])]))
      r$note <- paste(c(if (length(lv) == 1L) paste("nominal level", format(lv)),
                        if (length(lv) > 1L) paste("intervals at several nominal levels:",
                                                   paste(format(lv), collapse = ", ")),
                        if (!is.na(r$note)) r$note), collapse = "; ")
      if (!nzchar(r$note)) r$note <- NA_character_
      r
    },
    noncoverage_lower = {
      sel <- int_ok & has_truth
      k_prop(truth[sel] < rows$lower[sel], level)
    },
    noncoverage_upper = {
      sel <- int_ok & has_truth
      k_prop(truth[sel] > rows$upper[sel], level)
    },
    be_coverage = if (identical(estimand$truth$kind, "unknown")) {
      # a fixed truth cancels: truth plus the mean error is the mean estimate
      r <- k_be_coverage(rows, int_ok & usable, usable, rep(0, nrow(rows)), level)
      r$note <- paste(r$note, "(truth unknown but fixed, so the target is the mean estimate)")
      r
    } else k_be_coverage(rows, int_ok & usable & has_truth, usable & has_truth, truth, level),
    width = k_mean((rows$upper - rows$lower)[int_ok], level),
    rejection = k_prop(reject_indicator(rows[test_ok, , drop = FALSE], study$alpha), level),
    usable = k_prop(usable, level),
    usable_coverage = {
      hit <- int_ok & has_truth & rows$lower <= truth & truth <= rows$upper
      unknown <- int_ok & !has_truth
      if (!any(unknown)) k_prop(hit, level) else {
        # A usable interval whose truth is missing is neither a hit nor a
        # miss: report bounds that count it both ways, and no point value.
        lo <- k_prop(hit, level)
        hi <- k_prop(hit | unknown, level)
        kernel_result(lower = lo$lower, upper = hi$upper, n = length(hit),
                      status = "indeterminate",
                      note = sprintf("truth missing for %d usable interval(s); bounds count them as misses and as hits",
                                     sum(unknown)))
      }
    },
    stopf("Unknown measure %s.", m)
  )
}

# Smallest and largest value of a measure while the common truth moves over
# [truth - 2 SE, truth + 2 SE]. Coverage changes only where the truth
# crosses an interval endpoint, so every crossing is evaluated (by counting
# with sorted endpoints). For other measures the range is evaluated on a
# grid of 41 points plus the shifts where bias is zero and where the
# absolute error is smallest; this is exact for bias, MSE, RMSE and MAE.
truth_sensitivity <- function(m, rows, study, estimand, lim) {
  if (identical(m, "coverage")) {
    usable_int <- rows$estimate_status %in% USABLE_FOR_INTERVAL &
      rows$interval_status == "valid" & !is.na(rows$truth)
    t <- rows$truth[usable_int]
    if (!length(t)) return(c(NA_real_, NA_real_))
    # Interval i covers the shifted truth t + s when l_i <= s <= u_i. The
    # count is constant between endpoints, so the endpoints themselves and
    # one point inside every gap between them (midpoints, no fixed epsilon
    # that could step over a narrow gap) give its exact range.
    l <- sort(rows$lower[usable_int] - t)
    u <- sort(rows$upper[usable_int] - t)
    b <- unique(sort(c(-lim, lim, l, u)))
    b <- b[is.finite(b) & b >= -lim & b <= lim]
    s <- c(b, (b[-1] + b[-length(b)]) / 2)
    covered <- findInterval(s, l) - findInterval(s, u, left.open = TRUE)
    return(range(covered) / length(t))
  }
  e <- rows$estimate - rows$truth
  ok <- rows$valid & is.finite(e)
  # shifting the truth by s turns the errors into e - s, so a squared or
  # absolute error is smallest at s = mean(e) or s = median(e)
  extra <- if (any(ok)) c(mean(e[ok]), stats::median(e[ok])) else numeric()
  s <- c(seq(-lim, lim, length.out = 41L), extra[abs(extra) <= lim])
  v <- vapply(s, function(sh) measure_cell(m, rows, study, estimand, shift = sh)$value,
              numeric(1))
  if (all(is.na(v))) c(NA_real_, NA_real_) else range(v, na.rm = TRUE)
}

measure_needs_truth <- function(m) {
  if (inherits(m, "adempr_measure")) return("truth" %in% m$requires)
  MEASURE_INFO$needs_truth[MEASURE_INFO$measure == m]
}

#' Summarize performance
#'
#' Computes every selected performance measure for every scenario, method
#' and estimand, with its Monte Carlo standard error, a confidence interval
#' (exact Clopper-Pearson for proportions), and the counts that define its
#' denominator. The probability of a usable result is always included.
#'
#' Each measure uses its own eligible replications: bias and MSE use usable
#' estimates; SE based measures also need a valid SE; coverage needs a
#' valid interval (a point estimate is not required); rejection needs a
#' valid test result. These are *conditional* summaries, conditional on the
#' method producing the required output, and they are shown next to the
#' usable probability and the failure counts. The composite measure
#' `usable_coverage` (opt in) is unconditional over generated data sets and
#' is a different quantity from coverage.
#'
#' When a common truth is a Monte Carlo approximation with standard error
#' `truth_se`, the outer Monte Carlo standard error is conditional on that
#' approximation. The output then also gives `mcse_total` for bias (outer
#' and truth uncertainty combined, which is valid because the truth stage
#' uses an independent random number stream) and, for every measure that
#' depends on the truth, a sensitivity range: the smallest and largest value
#' of the measure while the truth moves over plus and minus two standard
#' errors. A numerical truth with an error bound gets the same range over
#' plus and minus the bound. For coverage every point where the truth crosses an interval
#' endpoint is evaluated; for bias, MSE, RMSE and MAE the range is exact;
#' for other measures it is evaluated on a grid of 41 points. The range is
#' a sensitivity analysis, not a confidence interval, and the point value
#' and MCSE are conditional on the approximate truth. Nonlinear measures do
#' not get a combined standard error.
#'
#' @param run A run, or results from [sim_results()].
#' @param measures Measures to compute; defaults to the study's measures.
#' @param weighted If `TRUE` and the scenario design has weights, adds
#'   prespecified weighted averages over scenarios for each method,
#'   estimand and measure (scenarios treated as independent; scenarios of
#'   weight 0 are left out). When a scenario has a boundary proportion, the
#'   weighted interval comes from exact per scenario intervals with a
#'   Bonferroni split. An average over data generating mechanisms hides
#'   where a method works and where it fails; report it next to the
#'   scenario results, never instead of them. Each weighted row's note
#'   gives the number of scenarios and the range of their values.
#' @param by Optional scenario fields within whose levels the weighted
#'   averages are formed (weights renormalized within each level), for
#'   example a treatment arm or design factor that must not be averaged
#'   over.
#' @return An object of class `adempr_summary`: a data frame with one row
#'   per scenario, method, estimand and measure, and attributes `failures`
#'   and `runtime`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) {
#'     tt <- t.test(data)
#'     sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
#'                  lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95)
#'   }))
#' run <- sim_run(study, replications = 50, seed = 1, path = tempfile("adempR-run-"))
#' sim_summarise(run)
#' @export
sim_summarise <- function(run, measures = NULL, weighted = FALSE, by = NULL) {
  res <- sim_results(run)
  man <- res$manifest
  study <- man$study
  if (!is.null(by)) {
    # `by` only shapes the weighted rows; asking for it otherwise is a mistake
    if (!isTRUE(weighted)) stopf("`by` needs `weighted = TRUE`.")
    if (is.null(study$scenarios$weights)) stopf("`by` needs a scenario design with weights.")
    bad <- setdiff(by, names(as.data.frame(study$scenarios)))
    if (length(bad)) stopf("Unknown scenario field(s) in `by`: %s.", paste(bad, collapse = ", "))
  }
  if (is.null(measures)) measures <- study$measures
  if (is.character(measures)) measures <- do.call(sim_measures, as.list(measures))
  mlist <- unclass(measures)
  if (!("usable" %in% names(mlist))) mlist <- c(mlist, list(usable = "usable"))
  est <- res$estimates
  reps <- res$replications
  out <- list()
  keys <- study$scenarios$keys
  labels <- study$scenarios$labels
  # measures whose contributions are events (built in proportions, and
  # custom measures returning logical contributions)
  prop_ids <- names(mlist)[vapply(mlist, function(m) is.character(m) &&
                                    m %in% PROPORTION_MEASURES, logical(1))]
  for (si in seq_along(keys)) {
    key <- keys[si]
    gen_ok <- if (is.null(reps)) integer() else
      reps$replication[reps$scenario_key == key & reps$generation_status == "ok"]
    n_gen_attempted <- if (is.null(reps)) 0L else sum(reps$scenario_key == key)
    for (mn in names(study$methods)) {
      for (en in study$methods[[mn]]$estimands) {
        rows <- if (is.null(est)) NULL else
          est[est$scenario_key == key & est$method == mn & est$estimand == en, , drop = FALSE]
        if (is.null(rows)) rows <- est[0, , drop = FALSE]
        rows <- rows[rows$replication %in% gen_ok, , drop = FALSE]
        estimand <- study$estimands[[en]]
        tse <- if (nrow(rows)) unique(rows$truth_se[!is.na(rows$truth_se)]) else numeric()
        tse <- if (length(tse) == 1L) tse else NA_real_
        tkind <- estimand$truth$kind
        for (mid in names(mlist)) {
          m <- mlist[[mid]]
          needs_truth <- measure_needs_truth(m)
          r <- if (!measure_applies(m, estimand$type)) {
            kernel_result(status = "not_applicable",
                          note = sprintf("measure does not apply to a %s estimand", estimand$type %||% "scalar"))
          } else if (is.null(rows) || !nrow(rows)) {
            kernel_result(status = "no_eligible", note = "no results")
          } else if (needs_truth && all(is.na(rows$truth)) && !identical(m, "usable_coverage")) {
            # (usable coverage still knows its misses: a data set without a
            # usable interval fails whatever the truth)
            kernel_result(status = "unavailable", note = sprintf("truth is %s", tkind))
          } else {
            measure_cell(m, rows, study, estimand)
          }
          if (isTRUE(r$proportion)) prop_ids <- union(prop_ids, mid)
          # The usable proportions are over all generated data sets. While
          # some have not been analyzed yet, the value covers only the
          # attempted ones; the bounds over all generated data sets count
          # the outstanding ones as all failing or all succeeding.
          n_out <- length(gen_ok) - nrow(rows)
          if (mid %in% c("usable", "usable_coverage") && n_out > 0L &&
              r$status %in% c("ok", "boundary", "no_eligible")) {
            k <- if (is.na(r$value)) 0 else round(r$value * r$n)
            r$status <- "incomplete"
            r$note <- paste(c(sprintf("%d of %d generated data sets not yet analyzed; the value covers the analyzed ones; over all generated data sets it lies in [%s, %s]",
                                      n_out, length(gen_ok), format(signif(k / length(gen_ok), 4)),
                                      format(signif((k + n_out) / length(gen_ok), 4))),
                              if (!is.na(r$note)) r$note), collapse = "; ")
          }
          # antithetic replications are dependent in pairs: the uncertainty
          # comes from the pairs, the independent units
          if (isTRUE(study$generator$antithetic) && nrow(rows)) {
            r <- pair_jackknife(m, rows, study, estimand, r)
          }
          label <- if (identical(m, "rejection") && nrow(rows)) {
            # the label uses exactly the rows the rejection rate uses
            rejection_label(rows[rows$estimate_status %in% USABLE_FOR_INTERVAL &
                                   rows$test_status == "valid", , drop = FALSE], estimand)
          } else mid
          sens_lo <- sens_hi <- mcse_total <- NA_real_
          # the truth moves over two standard errors of a Monte Carlo truth,
          # or over the error bound of a numerical one
          ter <- if (tkind == "numeric" && nrow(rows)) unique(rows$truth_error[!is.na(rows$truth_error)]) else numeric()
          lim <- if (!is.na(tse) && tse > 0 && tkind == "mc") 2 * tse else
            if (length(ter) == 1L && ter > 0) ter else NA_real_
          if (needs_truth && !is.na(lim) && r$status %in% c("ok", "boundary")) {
            rg <- truth_sensitivity(m, rows, study, estimand, lim)
            sens_lo <- rg[1]
            sens_hi <- rg[2]
            if (identical(m, "bias") && tkind == "mc") mcse_total <- sqrt(r$mcse^2 + tse^2)
          }
          out[[length(out) + 1L]] <- data.frame(
            scenario_key = key, scenario_id = labels[si], method = mn, estimand = en,
            measure = mid, label = label, value = r$value, mcse = r$mcse,
            lower = r$lower, upper = r$upper, n_eligible = r$n,
            n_units = as.integer(r$n_units %||% if (isTRUE(study$generator$antithetic)) NA_integer_ else r$n),
            n_planned = unname(planned_replications(man, key)), n_generated = length(gen_ok),
            n_generation_failed = n_gen_attempted - length(gen_ok),
            n_attempted = nrow(rows), n_completed = sum(rows$attempt_status == "ok"),
            n_valid = sum(rows$valid), truth_kind = tkind, truth_se = tse,
            mcse_total = mcse_total, truth_sens_lower = sens_lo,
            truth_sens_upper = sens_hi, status = r$status, note = r$note,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  summ <- bind_rows(out)
  if (weighted && !is.null(study$scenarios$weights)) {
    summ <- rbind(summ, weighted_rows(summ, study, by, prop_ids))
  }
  attr(summ, "failures") <- failure_table(res)
  attr(summ, "runtime") <- runtime_table(res)
  attr(summ, "study") <- study$name
  attr(summ, "sequential") <- man$sequential
  class(summ) <- c("adempr_summary", "data.frame")
  summ
}

#' @rdname sim_summarise
#' @export
sim_summarize <- sim_summarise

# The one registry of built in measures whose contributions are events (0
# or 1 per replication). Planning, adaptive allocation, sequential
# stopping, weighted averages and exact intervals all read it.
PROPORTION_MEASURES <- c("coverage", "noncoverage_lower", "noncoverage_upper", "rejection",
                         "usable", "usable_coverage")

weighted_rows <- function(summ, study, by = NULL, prop_ids = PROPORTION_MEASURES) {
  if (!is.null(by)) {
    sc <- as.data.frame(study$scenarios)
    miss <- setdiff(by, names(sc))
    if (length(miss)) stopf("Unknown scenario field(s) in `by`: %s.", paste(miss, collapse = ", "))
    # groups are typed tuples of the factor values, so ("a, b", "c") and
    # ("a", "b, c") stay apart
    lev <- do.call(paste, c(lapply(sc[by], col_text), sep = "\r"))
    lab <- do.call(paste, c(Map(function(n, v) paste0(n, "=", v), by, lapply(sc[by], as.character)),
                            sep = ", "))
    out <- lapply(unique(lev), function(l) {
      keys <- study$scenarios$keys[lev == l]
      sub <- summ[summ$scenario_key %in% keys, , drop = FALSE]
      st <- study
      w <- study$scenarios$weights[lev == l]
      if (!sum(w)) return(NULL)
      st$scenarios$keys <- keys
      st$scenarios$weights <- w / sum(w)
      r <- weighted_rows(sub, st, prop_ids = prop_ids)
      r$scenario_key <- paste0("weighted: ", lab[lev == l][1])
      r$scenario_id <- paste("prespecified weighted average over scenarios with", lab[lev == l][1])
      r
    })
    return(bind_rows(out))
  }
  w <- stats::setNames(study$scenarios$weights, study$scenarios$keys)
  grp <- unique(summ[c("method", "estimand", "measure")])
  z <- stats::qnorm(1 - (1 - study$level) / 2)
  # Scenarios that share random numbers (common = ) are dependent. Their
  # covariances are not estimated here, so the MCSE is the bound for
  # perfectly correlated scenarios, sum(w * MCSE), which holds for any
  # correlation.
  shared <- length(study$generator$common) > 0L
  rows <- lapply(seq_len(nrow(grp)), function(i) {
    s <- summ[summ$method == grp$method[i] & summ$estimand == grp$estimand[i] &
                summ$measure == grp$measure[i], , drop = FALSE]
    ww <- w[s$scenario_key]
    # scenarios with weight 0 do not enter the average
    s <- s[ww > 0, , drop = FALSE]
    ww <- ww[ww > 0]
    # the point value needs a value in every weighted scenario, whatever
    # its uncertainty status
    have <- !is.na(s$value) & is.finite(s$value)
    val <- if (all(have)) sum(ww * s$value) else NA_real_
    mc <- if (all(have) && !anyNA(s$mcse)) {
      if (shared) sum(ww * s$mcse) else sqrt(sum(ww^2 * s$mcse^2))
    } else NA_real_
    lo <- val - z * mc
    hi <- val + z * mc
    status <- if (all(have)) "ok" else "unavailable"
    note <- if (!all(have)) "some weighted scenarios have no value" else if (shared)
      "weighted average; scenarios share random numbers, so the MCSE is the bound for perfectly correlated scenarios" else
      "weighted average; scenarios treated as independent"
    boundary <- s$status %in% "boundary"
    if (all(have) && any(boundary)) {
      if (grp$measure[i] %in% prop_ids) {
        # A plug-in MCSE of 0 is not certainty. Per scenario intervals at
        # level 1 - (1 - level) / S all hold with probability at least
        # `level`, so their weighted sum bounds the weighted proportion.
        S <- nrow(s)
        conf <- 1 - (1 - study$level) / S
        if (isTRUE(study$generator$antithetic)) {
          # dependent pairs: exact limits are taken over the independent
          # pairs, for the shares of pairs at 1 and above 0 (scenarios
          # with all pair means equal), and Wald limits from the pair
          # jackknife MCSE elsewhere
          G <- s$n_units
          zs <- stats::qnorm(1 - (1 - conf) / 2)
          bnd <- s$status %in% "boundary"
          lo_s <- pmax(0, s$value - zs * s$mcse)
          hi_s <- pmin(1, s$value + zs * s$mcse)
          kone <- ifelse(s$value >= 1, G, 0)
          kpos <- ifelse(s$value > 0, G, 0)
          lo_s[bnd] <- cp_interval(kone[bnd], G[bnd], conf)[, "lower"]
          hi_s[bnd] <- cp_interval(kpos[bnd], G[bnd], conf)[, "upper"]
          lo <- if (anyNA(G[bnd])) NA_real_ else sum(ww * lo_s)
          hi <- if (anyNA(G[bnd])) NA_real_ else sum(ww * hi_s)
          note <- "weighted average; interval from per scenario intervals over independent antithetic pairs with a Bonferroni split, because some scenarios have all pair means equal"
        } else {
          k <- round(s$value * s$n_eligible)
          ci <- cp_interval(k, s$n_eligible, conf)
          lo <- sum(ww * ci[, "lower"])
          hi <- sum(ww * ci[, "upper"])
          note <- "weighted average; interval from exact per scenario intervals with a Bonferroni split, because some scenarios observed no or only events"
        }
      } else {
        lo <- hi <- mc <- NA_real_
        note <- "weighted average; no interval, because some scenarios have a degenerate MCSE"
      }
      status <- "boundary"
    }
    s1 <- s[1, , drop = FALSE]
    s1$scenario_key <- "weighted"
    s1$scenario_id <- "prespecified weighted average over scenarios"
    s1$label <- s1$measure
    s1$value <- val
    s1$mcse <- mc
    s1$lower <- lo
    s1$upper <- hi
    # counts are totals over the weighted scenarios
    for (cn in c("n_eligible", "n_units", "n_planned", "n_generated", "n_attempted", "n_completed",
                 "n_valid", "n_generation_failed")) {
      s1[[cn]] <- sum(s[[cn]])
    }
    s1$truth_se <- s1$mcse_total <- s1$truth_sens_lower <- s1$truth_sens_upper <- NA_real_
    # a Monte Carlo truth adds its own uncertainty to the bias; scenarios
    # may share one truth, so the truth terms are added as perfectly
    # correlated
    tse <- s$truth_se
    if (identical(grp$measure[i], "bias") && is.finite(mc) && any(is.finite(tse) & tse > 0)) {
      tse[!is.finite(tse)] <- 0
      s1$mcse_total <- sqrt(mc^2 + sum(ww * tse)^2)
    }
    s1$status <- status
    rngv <- range(s$value[have])
    s1$note <- if (any(have)) sprintf("%s; %d scenarios, values from %s to %s", note, nrow(s),
                                      format(rngv[1], digits = 4), format(rngv[2], digits = 4)) else note
    s1
  })
  bind_rows(rows)
}

failure_table <- function(res) {
  est <- res$estimates
  study <- res$manifest$study
  rows <- list()
  for (key in study$scenarios$keys) {
    reps <- res$replications
    n_gen <- if (is.null(reps)) 0L else sum(reps$scenario_key == key & reps$generation_status == "ok")
    n_genfail <- if (is.null(reps)) 0L else sum(reps$scenario_key == key & reps$generation_status != "ok")
    for (mn in names(study$methods)) {
      for (en in study$methods[[mn]]$estimands) {
        r <- if (is.null(est)) NULL else
          est[est$scenario_key == key & est$method == mn & est$estimand == en, , drop = FALSE]
        tab <- if (is.null(r) || !nrow(r)) character() else table(r$estimate_status)
        rows[[length(rows) + 1L]] <- data.frame(
          scenario_key = key, scenario_id = study$scenarios$labels[match(key, study$scenarios$keys)],
          method = mn, estimand = en, generated = n_gen, generation_failed = n_genfail,
          attempted = NROW(r), usable = sum(r$valid %in% TRUE),
          statuses = if (length(tab)) paste(names(tab), tab, sep = ":", collapse = ", ") else "",
          # counted among outputs that pass the failure policy, so that a
          # structurally valid SE of an excluded fit is not counted
          se_valid = sum(r$valid %in% TRUE & r$se_status %in% "valid"),
          interval_valid = sum(r$valid %in% TRUE & r$interval_status %in% "valid"),
          test_valid = sum(r$valid %in% TRUE & r$test_status %in% "valid"),
          with_warnings = sum(r$had_warning %in% TRUE),
          fallback_used = sum(r$fallback %in% TRUE),
          stringsAsFactors = FALSE)
      }
    }
  }
  bind_rows(rows)
}

runtime_table <- function(res) {
  # all attempts, including those of earlier retry rounds: the cost of the
  # run is what was spent, not only the attempts whose results are used
  att <- res$attempts_all %||% res$attempts
  if (is.null(att) || !nrow(att)) return(NULL)
  grp <- unique(att[c("scenario_key", "method")])
  rows <- lapply(seq_len(nrow(grp)), function(i) {
    a <- att[att$scenario_key == grp$scenario_key[i] & att$method == grp$method[i], ,
             drop = FALSE]
    ok <- a$status == "ok"
    data.frame(scenario_key = grp$scenario_key[i], method = grp$method[i],
               attempts = nrow(a), failed = sum(!ok), timeouts = sum(a$status == "timeout"),
               overruns = sum(a$overrun %in% TRUE),
               timeout_limit = a$timeout[1], total_seconds = sum(a$elapsed),
               mean_seconds = mean(a$elapsed), median_seconds = stats::median(a$elapsed),
               p90_seconds = unname(stats::quantile(a$elapsed, 0.9)),
               max_seconds = max(a$elapsed),
               mean_seconds_successful = if (any(ok)) mean(a$elapsed[ok]) else NA_real_,
               total_cpu_seconds = sum(a$cpu), stringsAsFactors = FALSE)
  })
  bind_rows(rows)
}

#' @export
print.adempr_summary <- function(x, digits = "mcse", ...) {
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  need <- c("scenario_id", "method", "estimand", "label", "value", "mcse", "n_eligible",
            "n_generated", "status")
  # a column subset is no longer a full summary: print it as a data frame
  if (!all(need %in% names(df))) return(print(df, ...))
  cat(sprintf("<adempR performance summary of `%s`>\n", attr(x, "study")))
  show <- data.frame(scenario = df$scenario_id, method = df$method, estimand = df$estimand,
                     measure = df$label,
                     value = if (identical(digits, "mcse")) fmt_mcse(df$value, df$mcse) else
                       fmt_num(df$value, digits),
                     mcse = fmt_num(df$mcse, 2), n = df$n_eligible,
                     of = df$n_generated, status = df$status,
                     stringsAsFactors = FALSE)
  print(show, row.names = FALSE, right = FALSE)
  sq <- attr(x, "sequential")
  if (!is.null(sq$analyzed)) {
    cat(sequential_note(sq), "\n", sep = "")
  }
  f <- attr(x, "failures")
  if (!is.null(f) && any(f$usable < f$generated | f$generation_failed > 0)) {
    cat("\nFailures and invalid outputs (see attr(x, \"failures\")):\n")
    ff <- f[f$usable < f$generated | f$generation_failed > 0, , drop = FALSE]
    print(ff[c("scenario_id", "method", "estimand", "generated", "usable", "statuses")],
          row.names = FALSE, right = FALSE)
  }
  invisible(x)
}

#' @export
as.data.frame.adempr_summary <- function(x, ...) {
  attr(x, "failures") <- NULL
  attr(x, "runtime") <- NULL
  attr(x, "study") <- NULL
  class(x) <- "data.frame"
  x
}

#' Flow of replications through a run
#'
#' A CONSORT style account of every planned replication, per scenario,
#' method and estimand: how many data sets were generated or failed to
#' generate, how many were not yet analyzed, how many method fits failed,
#' how many results were excluded by the validity rules, and how many were
#' usable. The counts add up: `planned = generation_failed + not_generated +
#' not_run + fit_failed + excluded + usable`. A primary fit that failed and
#' was replaced by its fallback ([sim_failure_policy()]) is counted with the
#' fallback output (usable or excluded), and `fallback_used` says how many
#' results came from a fallback. The report always prints this table, also
#' when nothing failed.
#'
#' @param run A run or results.
#' @return A data frame with one row per scenario, method and estimand.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 20, seed = 1, path = tempfile("adempR-run-"))
#' sim_flow(run)
#' @export
sim_flow <- function(run) {
  res <- sim_results(run)
  man <- res$manifest
  f <- failure_table(res)
  planned <- unname(planned_replications(man, man$study$scenarios$keys)[f$scenario_key])
  est <- res$estimates
  count <- function(cond) vapply(seq_len(nrow(f)), function(i) {
    if (is.null(est) || !nrow(est)) return(0L)
    sum(est$scenario_key == f$scenario_key[i] & est$method == f$method[i] &
          est$estimand == f$estimand[i] & cond(est))
  }, integer(1))
  # a failed primary fit replaced by a fallback output is not a failed fit:
  # the fallback result is usable or excluded like any other output
  fallback <- count(function(e) e$fallback %in% TRUE)
  fit_failed <- count(function(e) e$attempt_status != "ok" & !(e$fallback %in% TRUE))
  out <- data.frame(
    scenario_id = f$scenario_id, method = f$method, estimand = f$estimand,
    planned = planned, generation_failed = f$generation_failed,
    not_generated = pmax(planned - f$generated - f$generation_failed, 0L),
    not_run = f$generated - f$attempted, fit_failed = fit_failed,
    excluded = f$attempted - fit_failed - f$usable, usable = f$usable,
    fallback_used = fallback, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

#' Acceptability of a performance measure
#'
#' Classifies each cell of a measure against an acceptable range, taking
#' the Monte Carlo interval into account: `"inside"` when the whole
#' interval lies in the range, `"outside"` when it lies entirely outside,
#' and `"uncertain"` otherwise (more replications are needed to decide).
#' For the rejection rate the default range is Bradley's (1978) liberal
#' criterion, `alpha / 2` to `3 alpha / 2`, with the study's `alpha` by
#' default. It applies to type I error, so it is used only on cells whose
#' label is `type_i_error` (the null hypothesis holds); other cells get no
#' verdict and a note. Other measures need an explicit range. For a run
#' made by [sim_run_sequential()], a monitored measure is judged with the
#' design's monitoring intervals, which stay valid after the stopping.
#'
#' @param run A run or results.
#' @param measure A built in measure.
#' @param range Two numbers, the acceptable range.
#' @param alpha Nominal level for the default rejection range; defaults to
#'   the study's `alpha`.
#' @return The rows of [sim_summarise()] for the measure with columns
#'   `range_lower`, `range_upper`, `verdict` and `verdict_note`.
#' @references Bradley JV (1978). Robustness? British Journal of
#'   Mathematical and Statistical Psychology 31, 144-152.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(t = function(data) {
#'     sim_estimate("mean", mean(data), p_value = t.test(data)$p.value, null = 0)
#'   }), measures = "rejection")
#' run <- sim_run(study, replications = 100, seed = 1, path = tempfile("adempR-run-"))
#' sim_acceptable(run)
#' @export
sim_acceptable <- function(run, measure = "rejection", range = NULL, alpha = NULL) {
  if (!is_scalar_character(measure) || !(measure %in% MEASURE_INFO$measure)) {
    stopf("`measure` must be one built in measure; see sim_measure_info().")
  }
  res <- sim_results(run)
  null_only <- is.null(range)
  if (null_only) {
    if (measure != "rejection") stopf("Give an acceptable `range` for %s.", measure)
    alpha <- alpha %||% res$manifest$study$alpha
    if (!is_scalar_number(alpha) || alpha <= 0 || alpha >= 1) stopf("`alpha` must be in (0, 1).")
    range <- c(alpha / 2, 3 * alpha / 2)
  }
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) || range[1] > range[2]) {
    stopf("`range` must be two increasing numbers.")
  }
  s <- as.data.frame(sim_summarise(res, measures = measure))
  s <- s[s$measure == measure, , drop = FALSE]
  # A sequential run was stopped by its monitoring intervals; only those
  # (Clopper-Pearson at the per cell level of the design) stay valid after
  # the stopping, so the verdict uses them.
  design <- res$manifest$sequential$design
  if (!is.null(design) && measure %in% design$measures && nrow(s)) {
    ok <- s$n_eligible > 0 & !is.na(s$value)
    ci <- cp_interval(round(s$value[ok] * s$n_eligible[ok]), s$n_eligible[ok],
                      design$conf_level_per_interval)
    s$lower[ok] <- ci[, "lower"]
    s$upper[ok] <- ci[, "upper"]
    s$interval <- "sequential monitoring"
  }
  s$range_lower <- range[1]
  s$range_upper <- range[2]
  s$verdict <- ifelse(is.na(s$lower) | is.na(s$upper), NA_character_,
                      ifelse(s$lower >= range[1] & s$upper <= range[2], "inside",
                             ifelse(s$upper < range[1] | s$lower > range[2], "outside",
                                    "uncertain")))
  s$verdict_note <- NA_character_
  if (null_only) {
    off <- s$label != "type_i_error"
    s$verdict[off] <- NA_character_
    s$verdict_note[off] <- "the type I error range applies only where the null hypothesis holds"
  }
  rownames(s) <- NULL
  s
}
