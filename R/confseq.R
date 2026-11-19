# Time-uniform confidence sequences ---------------------------------------

#' Confidence sequences for performance measures
#'
#' A confidence interval for coverage is valid at one planned number of
#' replications; looking at it repeatedly while replications accumulate and
#' stopping when it looks good inflates the error rate. A confidence
#' sequence is valid at every number of replications simultaneously, so a
#' run can be monitored and stopped at any time.
#'
#' For proportions (`"coverage"`, `"rejection"`) the contributions are
#' bounded between 0 and 1, hence 1/2 sub-Gaussian, and the two-sided normal
#' mixture boundary of Robbins (1970) gives a sequence that is valid for
#' every sample size without approximation:
#' `|S_t - t p| < sqrt((t/4 + rho) log((t/4 + rho) / (rho alpha^2)))`.
#' For unbounded measures (`"bias"`, `"mse"`) the asymptotic confidence
#' sequence of Waudby-Smith et al. (2024) is used, with the running
#' standard deviation in place of the unknown one. The mixing parameter
#' `rho` is tuned to be tightest at `t_opt` replications.
#'
#' The guarantee holds per cell (scenario, method and estimand), with the
#' replication as the independent unit, and only for tuning fixed before
#' the first look. The tuning (`t_opt` and `alpha`) is therefore recorded
#' in the run at the first call for a measure and reused by later calls,
#' also after the run is extended; a call with different tuning is an
#' error. The asymptotic sequence assumes a finite positive variance of the
#' contributions (for `"mse"`, a finite fourth moment of the error) and
#' becomes valid as the number of replications grows; cells whose running
#' standard deviation is 0 get no interval (`NA`) because the asymptotic
#' boundary is not defined there. Rejection uses the decision and the
#' eligibility of [sim_summarise()]: an explicit `reject` wins over the
#' p-value, and tests from fits excluded by the failure policy do not
#' count.
#'
#' @param run A run or results.
#' @param measure `"coverage"`, `"rejection"`, `"bias"` or `"mse"`.
#' @param alpha Error probability over the whole sequence of one cell.
#' @param t_opt Replication count at which the sequence is tightest. At
#'   the first call for a measure it defaults to the number of
#'   replications of the run; later calls reuse the recorded value.
#' @param path Return the sequence at every replication count instead of
#'   only the last.
#' @return A data frame per scenario, method and estimand (and replication
#'   count with `path = TRUE`): `n`, `estimate`, `lower`, `upper` and
#'   `type` (`"exact"` or `"asymptotic"`).
#' @references Robbins H (1970). Statistical methods related to the law of
#'   the iterated logarithm. Annals of Mathematical Statistics 41:1397-1409.
#'   Howard SR, Ramdas A, McAuliffe J, Sekhon J (2021). Time-uniform,
#'   nonparametric, nonasymptotic confidence sequences. Annals of
#'   Statistics 49:1055-1080. Waudby-Smith I, Arbour D, Sinha R, Kennedy
#'   EH, Ramdas A (2024). Time-uniform central limit theory and asymptotic
#'   confidence sequences. Annals of Statistics 52:2613-2640.
#' @export
sim_confidence_sequence <- function(run, measure = c("coverage", "rejection", "bias", "mse"),
                                    alpha = 0.05, t_opt = NULL, path = FALSE) {
  measure <- match.arg(measure)
  if (!is_scalar_number(alpha) || alpha <= 0 || alpha >= 1) stopf("`alpha` must be in (0, 1).")
  res <- sim_results(run)
  study <- res$manifest$study
  if (isTRUE(study$generator$antithetic)) {
    stopf("Confidence sequences need independent replications; antithetic pairs are dependent. Use sim_summarise(), whose MCSEs come from the pairs.")
  }
  tuning <- confseq_tuning(res, measure, alpha, t_opt)
  t_opt <- tuning$t_opt
  e <- res$estimates
  if (measure == "coverage") {
    e <- e[e$interval_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(e$truth), , drop = FALSE]
    e$x <- as.numeric(e$lower <= e$truth & e$truth <= e$upper)
  } else if (measure == "rejection") {
    e <- e[e$test_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL, , drop = FALSE]
    e$x <- as.numeric(reject_indicator(e, study$alpha))
    e <- e[!is.na(e$x), , drop = FALSE]
  } else {
    e <- e[e$valid & is.finite(e$estimate) & !is.na(e$truth), , drop = FALSE]
    e$x <- if (measure == "bias") e$estimate - e$truth else (e$estimate - e$truth)^2
  }
  bounded <- measure %in% c("coverage", "rejection")
  out <- list()
  for (grp in split(e, list(e$scenario_key, e$method, e$estimand), drop = TRUE)) {
    grp <- grp[order(grp$replication), , drop = FALSE]
    x <- grp$x
    t <- seq_along(x)
    m <- cumsum(x) / t
    if (bounded) {
      # Robbins normal mixture for 1/2 sub-Gaussian sums; rho minimizing the
      # boundary near t_opt (Howard et al. 2021, section 3.5)
      half <- robbins_halfwidth(t, alpha, t_opt)
      lo <- pmax(0, m - half)
      up <- pmin(1, m + half)
    } else {
      s <- running_sd(x)
      r2 <- (-2 * log(alpha) + log(-2 * log(alpha) + 1)) / t_opt
      half <- s * sqrt(2 * (t * r2 + 1) / (t^2 * r2) * log(sqrt(t * r2 + 1) / alpha))
      # a zero running SD gives no asymptotic interval, not a zero width one
      half[!(s > 0)] <- NA_real_
      lo <- m - half
      up <- m + half
    }
    df <- data.frame(scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
                     method = grp$method[1], estimand = grp$estimand[1], measure = measure,
                     n = t, estimate = m, lower = lo, upper = up,
                     type = if (bounded) "exact" else "asymptotic", t_opt = t_opt, alpha = alpha,
                     stringsAsFactors = FALSE)
    out[[length(out) + 1L]] <- if (path) df else df[nrow(df), , drop = FALSE]
  }
  out <- bind_rows(out)
  if (!is.null(out)) rownames(out) <- NULL
  out
}

# Half width of the two-sided normal mixture confidence sequence for the
# mean of [0, 1] bounded (1/2 sub-Gaussian) contributions after t of them.
robbins_halfwidth <- function(t, alpha, t_opt) {
  v <- t / 4
  rho <- (t_opt / 4) / (2 * log(1 / alpha) + log(1 + 2 * log(1 / alpha)))
  sqrt((v + rho) * log((v + rho) / (rho * alpha^2))) / t
}

# Running standard deviation of x[1:k] for every k in one pass (Welford's
# algorithm).
running_sd <- function(x) {
  n <- length(x)
  out <- rep(NA_real_, n)
  mu <- 0
  m2 <- 0
  for (k in seq_len(n)) {
    d <- x[k] - mu
    mu <- mu + d / k
    m2 <- m2 + d * (x[k] - mu)
    if (k > 1) out[k] <- sqrt(max(m2, 0) / (k - 1))
  }
  out
}

# The tuning of a measure's confidence sequence is fixed at its first use
# and recorded in the run, so later looks (after extensions) belong to the
# same sequence.
confseq_tuning <- function(res, measure, alpha, t_opt) {
  if (!is.null(t_opt) && (!is_scalar_number(t_opt) || t_opt < 1)) stopf("`t_opt` must be a positive number.")
  rec <- res$manifest$monitoring[[measure]]
  if (!is.null(rec)) {
    if (!is.null(t_opt) && !isTRUE(all.equal(t_opt, rec$t_opt))) {
      stopf("The confidence sequence for `%s` was tuned with t_opt = %s at its first use; a sequence is valid only with tuning fixed before the first look.",
            measure, format(rec$t_opt))
    }
    if (!isTRUE(all.equal(alpha, rec$alpha))) {
      stopf("The confidence sequence for `%s` was started with alpha = %s; use the same alpha at every look.",
            measure, format(rec$alpha))
    }
    return(rec)
  }
  rec <- list(t_opt = t_opt %||% res$manifest$replications, alpha = alpha,
              unit = "replication", first_look = as.character(Sys.time()))
  if (!is.null(res$path) && dir.exists(res$path)) {
    update_manifest(res$path, function(man) {
      man$monitoring[[measure]] <- man$monitoring[[measure]] %||% rec
      man
    })
  }
  rec
}
