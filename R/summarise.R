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
      if (all(t == t[1])) {
        if (abs(t[1]) < .Machine$double.eps^0.5) {
          kernel_result(n = length(t), status = "undefined",
                        note = "relative bias is undefined at zero truth")
        } else {
          r <- k_mean(e[sel] / t[1], level)
          r$note <- "bias divided by the (constant) truth"
          r
        }
      } else {
        mt <- mean(t)
        if (abs(mt) < .Machine$double.eps^0.5) {
          return(kernel_result(n = length(t), status = "undefined",
                               note = "mean truth is zero"))
        }
        val <- sum(e[sel]) / sum(t)
        lin <- (e[sel] - val * t) / mt
        mc <- if (length(t) > 1L) stats::sd(lin) / sqrt(length(t)) else NA_real_
        kernel_result(value = val, mcse = mc, n = length(t),
                      note = "varying truth: ratio of averages, sum(error) / sum(truth), delta method MCSE")
      }
    },
    empse = k_empse(rows$estimate[usable], level),
    mean_se = k_mean(rows$se[usable & se_ok], level),
    rmvar = k_rmvar(rows$se[usable & se_ok], level),
    se_ratio = k_se_ratio(rows$estimate[usable & se_ok], rows$se[usable & se_ok], level),
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
    be_coverage = {
      sel <- int_ok & usable
      if (!any(sel)) return(k_prop(logical(), level))
      mbar <- mean(rows$estimate[usable])
      r <- k_prop(rows$lower[sel] <= mbar & mbar <= rows$upper[sel], level)
      r$note <- paste(c("relative to the mean estimate; MCSE ignores its uncertainty",
                        if (!is.na(r$note)) r$note), collapse = "; ")
      r
    },
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
truth_sensitivity <- function(m, rows, study, estimand, tse) {
  lim <- 2 * tse
  if (identical(m, "coverage")) {
    usable_int <- rows$estimate_status %in% USABLE_FOR_INTERVAL &
      rows$interval_status == "valid" & !is.na(rows$truth)
    t <- rows$truth[usable_int]
    if (!length(t)) return(c(NA_real_, NA_real_))
    l <- sort(rows$lower[usable_int] - t)
    u <- sort(rows$upper[usable_int] - t)
    eps <- 1e-9 * max(1, abs(c(l, u)[is.finite(c(l, u))]))
    s <- c(-lim, lim, l, l - eps, u, u + eps)
    s <- s[is.finite(s) & s >= -lim & s <= lim]
    covered <- findInterval(s, l) - findInterval(s, u, left.open = TRUE)
    return(range(covered) / length(t))
  }
  e <- rows$estimate - rows$truth
  ok <- rows$valid & is.finite(e)
  extra <- if (any(ok)) c(-mean(e[ok]), -stats::median(e[ok])) else numeric()
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
#' errors. For coverage every point where the truth crosses an interval
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
#'   estimand and measure (scenarios treated as independent).
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
sim_summarise <- function(run, measures = NULL, weighted = FALSE) {
  res <- sim_results(run)
  man <- res$manifest
  study <- man$study
  if (is.null(measures)) measures <- study$measures
  if (is.character(measures)) measures <- do.call(sim_measures, as.list(measures))
  mlist <- unclass(measures)
  if (!("usable" %in% names(mlist))) mlist <- c(mlist, list(usable = "usable"))
  est <- res$estimates
  reps <- res$replications
  out <- list()
  keys <- study$scenarios$keys
  labels <- study$scenarios$labels
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
          r <- if (is.null(rows) || !nrow(rows)) {
            kernel_result(status = "no_eligible", note = "no results")
          } else if (needs_truth && all(is.na(rows$truth))) {
            kernel_result(status = "unavailable", note = sprintf("truth is %s", tkind))
          } else {
            measure_cell(m, rows, study, estimand)
          }
          label <- if (identical(m, "rejection") && nrow(rows)) {
            # the label uses exactly the rows the rejection rate uses
            rejection_label(rows[rows$estimate_status %in% USABLE_FOR_INTERVAL &
                                   rows$test_status == "valid", , drop = FALSE], estimand)
          } else mid
          sens_lo <- sens_hi <- mcse_total <- NA_real_
          if (needs_truth && !is.na(tse) && tse > 0 && tkind == "mc" &&
              r$status %in% c("ok", "boundary")) {
            rg <- truth_sensitivity(m, rows, study, estimand, tse)
            sens_lo <- rg[1]
            sens_hi <- rg[2]
            if (identical(m, "bias")) mcse_total <- sqrt(r$mcse^2 + tse^2)
          }
          out[[length(out) + 1L]] <- data.frame(
            scenario_key = key, scenario_id = labels[si], method = mn, estimand = en,
            measure = mid, label = label, value = r$value, mcse = r$mcse,
            lower = r$lower, upper = r$upper, n_eligible = r$n,
            n_planned = man$replications, n_generated = length(gen_ok),
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
    summ <- rbind(summ, weighted_rows(summ, study))
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

weighted_rows <- function(summ, study) {
  w <- stats::setNames(study$scenarios$weights, study$scenarios$keys)
  grp <- unique(summ[c("method", "estimand", "measure")])
  rows <- lapply(seq_len(nrow(grp)), function(i) {
    s <- summ[summ$method == grp$method[i] & summ$estimand == grp$estimand[i] &
                summ$measure == grp$measure[i], , drop = FALSE]
    ok <- s$status %in% c("ok", "boundary") & !is.na(s$value)
    allok <- all(ok)
    ww <- w[s$scenario_key]
    val <- if (allok) sum(ww * s$value) else NA_real_
    mc <- if (allok && !anyNA(s$mcse)) sqrt(sum(ww^2 * s$mcse^2)) else NA_real_
    s1 <- s[1, , drop = FALSE]
    s1$scenario_key <- "weighted"
    s1$scenario_id <- "prespecified weighted average over scenarios"
    s1$label <- s1$measure
    s1$value <- val
    s1$mcse <- mc
    s1$lower <- val - stats::qnorm(1 - (1 - study$level) / 2) * mc
    s1$upper <- val + stats::qnorm(1 - (1 - study$level) / 2) * mc
    s1[c("n_eligible", "n_generated", "n_attempted", "n_completed", "n_valid",
         "n_generation_failed")] <- NA
    s1$truth_se <- s1$mcse_total <- s1$truth_sens_lower <- s1$truth_sens_upper <- NA_real_
    s1$status <- if (allok) "ok" else "unavailable"
    s1$note <- if (allok) "weighted average; scenarios treated as independent" else
      "some scenarios have no usable value"
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
          se_valid = sum(r$se_status %in% "valid"),
          interval_valid = sum(r$interval_status %in% "valid"),
          test_valid = sum(r$test_status %in% "valid"),
          with_warnings = sum(r$had_warning %in% TRUE),
          stringsAsFactors = FALSE)
      }
    }
  }
  bind_rows(rows)
}

runtime_table <- function(res) {
  att <- res$attempts
  if (is.null(att) || !nrow(att)) return(NULL)
  grp <- unique(att[c("scenario_key", "method")])
  rows <- lapply(seq_len(nrow(grp)), function(i) {
    a <- att[att$scenario_key == grp$scenario_key[i] & att$method == grp$method[i], ,
             drop = FALSE]
    ok <- a$status == "ok"
    data.frame(scenario_key = grp$scenario_key[i], method = grp$method[i],
               attempts = nrow(a), failed = sum(!ok), timeouts = sum(a$status == "timeout"),
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
print.adempr_summary <- function(x, digits = 4, ...) {
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  need <- c("scenario_id", "method", "estimand", "label", "value", "mcse", "n_eligible",
            "n_generated", "status")
  # a column subset is no longer a full summary: print it as a data frame
  if (!all(need %in% names(df))) return(print(df, ...))
  cat(sprintf("<adempR performance summary of `%s`>\n", attr(x, "study")))
  show <- data.frame(scenario = df$scenario_id, method = df$method, estimand = df$estimand,
                     measure = df$label, value = fmt_num(df$value, digits),
                     mcse = fmt_num(df$mcse, 2), n = df$n_eligible,
                     of = df$n_generated, status = df$status,
                     stringsAsFactors = FALSE)
  print(show, row.names = FALSE, right = FALSE)
  sq <- attr(x, "sequential")
  if (!is.null(sq)) {
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
