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
#' @param run A run or results.
#' @param measure `"coverage"`, `"rejection"`, `"bias"` or `"mse"`.
#' @param alpha Error probability over the whole sequence.
#' @param t_opt Replication count at which the sequence is tightest;
#'   defaults to the number of replications of the run.
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
  t_opt <- t_opt %||% res$manifest$replications
  if (!is_scalar_number(t_opt) || t_opt < 1) stopf("`t_opt` must be a positive number.")
  e <- res$estimates
  if (measure == "coverage") {
    e <- e[e$interval_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(e$truth), , drop = FALSE]
    e$x <- as.numeric(e$lower <= e$truth & e$truth <= e$upper)
  } else if (measure == "rejection") {
    e <- e[e$test_status == "valid" & !is.na(e$p_value), , drop = FALSE]
    a <- ifelse(is.na(e$alpha), study$alpha, e$alpha)
    e$x <- as.numeric(e$p_value <= a)
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
      s <- vapply(t, function(k) if (k > 1) stats::sd(x[seq_len(k)]) else NA_real_, numeric(1))
      r2 <- (-2 * log(alpha) + log(-2 * log(alpha) + 1)) / t_opt
      half <- s * sqrt(2 * (t * r2 + 1) / (t^2 * r2) * log(sqrt(t * r2 + 1) / alpha))
      lo <- m - half
      up <- m + half
    }
    df <- data.frame(scenario_id = study$scenarios$labels[match(grp$scenario_key[1], study$scenarios$keys)],
                     method = grp$method[1], estimand = grp$estimand[1], measure = measure,
                     n = t, estimate = m, lower = lo, upper = up,
                     type = if (bounded) "exact" else "asymptotic", stringsAsFactors = FALSE)
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
