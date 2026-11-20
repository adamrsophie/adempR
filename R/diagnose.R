# Diagnostics and plots ----------------------------------------------------

#' Diagnose a simulation study
#'
#' Collects the checks recommended by White et al. (2024): failure and
#' validity counts by scenario and method, the most frequent errors and
#' warnings, outlying estimates and standard errors (flagged, never
#' removed, with the scenario and replication needed to replay them), the
#' association between estimates and their standard errors, nonfinite
#' outputs, storage problems and recorded revisions. It also flags data
#' sets that repeat across the replications of a scenario (a generator that
#' fixes the complete data removes a source of variation; Morris et al.),
#' methods whose two-sided test at `alpha` and interval at level
#' `1 - alpha` disagree about the null (a calibration problem), and
#' degenerate estimates: constant within a cell, or identical for two
#' estimands of one method in every replication.
#'
#' @param run A run or results.
#' @param outlier_z Robust z threshold, `|x - median| / mad(x)` (`mad()`
#'   already includes the factor 1.4826), above
#'   which an estimate or standard error is flagged.
#' @return An object of class `adempr_diagnosis`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(5L, 50L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) {
#'     sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)))
#'   }))
#' run <- sim_run(study, replications = 50, seed = 1, path = tempfile("adempR-run-"))
#' sim_diagnose(run)
#' @export
sim_diagnose <- function(run, outlier_z = 5) {
  res <- sim_results(run)
  study <- res$manifest$study
  est <- res$estimates
  lab <- function(k) study$scenarios$labels[match(k, study$scenarios$keys)]
  outliers <- list()
  assoc <- list()
  if (!is.null(est) && nrow(est)) {
    cells <- unique(est[c("scenario_key", "method", "estimand")])
    for (i in seq_len(nrow(cells))) {
      r <- est[est$scenario_key == cells$scenario_key[i] & est$method == cells$method[i] &
                 est$estimand == cells$estimand[i] & est$valid, , drop = FALSE]
      if (nrow(r) < 5L) next
      for (what in c("estimate", "se")) {
        x <- r[[what]]
        ok <- is.finite(x)
        if (sum(ok) < 5L) next
        md <- stats::median(x[ok])
        s <- stats::mad(x[ok], center = md)  # mad() already scales by 1.4826
        z <- if (s > 0) (x - md) / s else ifelse(x == md, 0, Inf)
        flag <- !is.finite(x) | abs(z) > outlier_z
        flag[is.na(flag)] <- FALSE
        if (any(flag)) {
          outliers[[length(outliers) + 1L]] <- data.frame(
            scenario_id = lab(cells$scenario_key[i]), scenario_key = cells$scenario_key[i],
            method = cells$method[i], estimand = cells$estimand[i], quantity = what,
            replication = r$replication[flag], value = x[flag], robust_z = z[flag],
            stringsAsFactors = FALSE)
        }
      }
      ok <- is.finite(r$estimate) & is.finite(r$se)
      if (sum(ok) >= 5L) {
        dev <- abs(r$estimate[ok] - (if (all(is.na(r$truth[ok]))) mean(r$estimate[ok]) else r$truth[ok]))
        assoc[[length(assoc) + 1L]] <- data.frame(
          scenario_id = lab(cells$scenario_key[i]), method = cells$method[i],
          estimand = cells$estimand[i],
          cor_estimate_se = suppressWarnings(stats::cor(r$estimate[ok], r$se[ok], method = "spearman")),
          cor_abs_error_se = suppressWarnings(stats::cor(dev, r$se[ok], method = "spearman")),
          stringsAsFactors = FALSE)
      }
    }
  }
  msgs <- NULL
  att <- res$attempts
  if (!is.null(att) && nrow(att)) {
    m <- att[att$status != "ok" | att$n_warnings > 0, , drop = FALSE]
    if (nrow(m)) {
      txt <- ifelse(m$status != "ok", paste0(m$status, ": ", m$message),
                    paste0("warning: ", m$warnings))
      tab <- sort(table(paste(m$method, txt, sep = " | ")), decreasing = TRUE)
      msgs <- data.frame(method_message = names(tab), count = as.integer(tab),
                         stringsAsFactors = FALSE)
      msgs <- utils::head(msgs, 20L)
    }
  }
  gen_fail <- NULL
  if (!is.null(res$replications)) {
    g <- res$replications[res$replications$generation_status != "ok", , drop = FALSE]
    if (nrow(g)) {
      gen_fail <- data.frame(scenario_id = lab(g$scenario_key), replication = g$replication,
                             message = g$generation_message, stringsAsFactors = FALSE)
    }
  }
  # data sets that repeat across replications of a scenario: the outer
  # replication is then not a fresh draw from the mechanism
  repeated <- NULL
  rp <- res$replications
  if (!is.null(rp) && "data_md5" %in% names(rp)) {
    rp <- rp[rp$generation_status == "ok" & !is.na(rp$data_md5), , drop = FALSE]
    rows <- lapply(unique(rp$scenario_key), function(k) {
      d <- rp$data_md5[rp$scenario_key == k]
      if (length(d) < 2L || !anyDuplicated(d)) return(NULL)
      data.frame(scenario_id = lab(k), replications = length(d), distinct = length(unique(d)),
                 stringsAsFactors = FALSE)
    })
    repeated <- bind_rows(rows)
  }
  # an interval and a test from the same method that disagree (Morris):
  # two-sided test at alpha against an interval at level 1 - alpha
  disagree <- NULL
  degenerate <- NULL
  if (!is.null(est) && nrow(est)) {
    a <- ifelse(is.na(est$alpha), study$alpha, est$alpha)
    both <- est$test_status == "valid" & est$interval_status == "valid" &
      est$alternative == "two.sided" & !is.na(est$null) & !is.na(est$level) &
      !is.na(est$p_value) & abs(a - (1 - est$level)) < 1e-12
    both[is.na(both)] <- FALSE
    if (any(both)) {
      b <- est[both, , drop = FALSE]
      excl <- b$null < b$lower | b$null > b$upper
      b$disagree <- reject_indicator(b, study$alpha) != excl
      cells <- unique(b[c("scenario_key", "method", "estimand")])
      disagree <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
        r <- b[b$scenario_key == cells$scenario_key[i] & b$method == cells$method[i] &
                 b$estimand == cells$estimand[i], , drop = FALSE]
        if (!any(r$disagree)) return(NULL)
        data.frame(scenario_id = lab(cells$scenario_key[i]), method = cells$method[i],
                   estimand = cells$estimand[i], n_both = nrow(r), n_disagree = sum(r$disagree),
                   replications = paste(utils::head(r$replication[r$disagree], 5), collapse = ", "),
                   stringsAsFactors = FALSE)
      }))
    }
    # degenerate estimates: constant within a cell, or equal for two
    # estimands of one method in every replication
    u <- est[est$valid & is.finite(est$estimate), , drop = FALSE]
    cells <- unique(u[c("scenario_key", "method")])
    degenerate <- bind_rows(lapply(seq_len(nrow(cells)), function(i) {
      r <- u[u$scenario_key == cells$scenario_key[i] & u$method == cells$method[i], , drop = FALSE]
      out <- list()
      for (en in unique(r$estimand)) {
        x <- r$estimate[r$estimand == en]
        if (length(x) >= 2L && all(x == x[1])) {
          out[[length(out) + 1L]] <- data.frame(scenario_id = lab(cells$scenario_key[i]),
            method = cells$method[i], estimand = en,
            problem = sprintf("the same estimate (%g) in all %d replications", x[1], length(x)),
            stringsAsFactors = FALSE)
        }
      }
      ens <- unique(r$estimand)
      if (length(ens) >= 2L) {
        for (a in seq_along(ens)) for (bb in seq_along(ens)) if (a < bb) {
          xa <- r[r$estimand == ens[a], c("replication", "estimate")]
          xb <- r[r$estimand == ens[bb], c("replication", "estimate")]
          mm <- merge(xa, xb, by = "replication")
          if (nrow(mm) >= 2L && all(mm$estimate.x == mm$estimate.y)) {
            out[[length(out) + 1L]] <- data.frame(scenario_id = lab(cells$scenario_key[i]),
              method = cells$method[i], estimand = paste(ens[a], ens[bb], sep = " = "),
              problem = sprintf("identical estimates for two estimands in all %d replications",
                                nrow(mm)), stringsAsFactors = FALSE)
          }
        }
      }
      bind_rows(out)
    }))
  }
  # methods with known performance: a miss beyond Monte Carlo error points
  # to noise or a bug in the generator or the evaluation
  ref_check <- if (!is.null(est) && nrow(est)) sim_check_reference(res) else NULL
  # record and diagnostics callbacks that failed
  cb <- bind_rows(list(
    if (!is.null(res$records) && nrow(res$records)) {
      r <- res$records[res$records$name == ".record_error", , drop = FALSE]
      if (nrow(r)) data.frame(callback = "record", scenario_id = lab(r$scenario_key),
                              method = NA_character_, replication = r$replication,
                              message = r$message %||% NA_character_, stringsAsFactors = FALSE)
    },
    if (!is.null(res$diagnostics) && nrow(res$diagnostics)) {
      d <- res$diagnostics[res$diagnostics$name == ".diagnostics_error", , drop = FALSE]
      if (nrow(d)) data.frame(callback = "diagnostics", scenario_id = lab(d$scenario_key),
                              method = d$method, replication = d$replication,
                              message = d$message %||% NA_character_, stringsAsFactors = FALSE)
    }))
  structure(list(
    failures = failure_table(res), messages = msgs, generation_failures = gen_fail,
    callback_errors = cb,
    reference_check = ref_check, records = record_summary(res, lab),
    repeated_data = repeated, test_interval_disagreement = disagree,
    degenerate = degenerate,
    outliers = bind_rows(outliers), association = bind_rows(assoc),
    storage_problems = res$problems, revisions = res$manifest$revisions,
    amendments = study$amendments, outlier_z = outlier_z
  ), class = "adempr_diagnosis")
}

#' @export
print.adempr_diagnosis <- function(x, ...) {
  cat("<adempR diagnosis>\n")
  f <- x$failures
  cat("\nUsable results by scenario and method:\n")
  print(f[c("scenario_id", "method", "estimand", "generated", "usable", "statuses")],
        row.names = FALSE, right = FALSE)
  if (!is.null(x$generation_failures)) {
    cat(sprintf("\nData generation failed in %d replication(s).\n", nrow(x$generation_failures)))
  }
  if (!is.null(x$callback_errors) && nrow(x$callback_errors)) {
    cat(sprintf("\nRecord or diagnostics callbacks failed in %d case(s) (see $callback_errors).\n",
                nrow(x$callback_errors)))
  }
  if (!is.null(x$messages)) {
    cat("\nMost frequent errors and warnings:\n")
    print(x$messages, row.names = FALSE, right = FALSE)
  }
  if (!is.null(x$outliers) && nrow(x$outliers)) {
    cat(sprintf("\n%d outlying value(s) (robust |z| > %g); inspect with sim_replay():\n",
                nrow(x$outliers), x$outlier_z))
    print(utils::head(x$outliers[c("scenario_id", "method", "quantity", "replication",
                                   "value", "robust_z")], 10L), row.names = FALSE)
  } else {
    cat("\nNo outlying estimates or standard errors.\n")
  }
  if (NROW(x$records)) {
    cat("\nRealized properties of the generated data sets:\n")
    print(x$records, row.names = FALSE)
  }
  if (NROW(x$reference_check) && !all(x$reference_check$flag %in% FALSE)) {
    cat("\nMethods with known performance that miss it by more than 3 MCSE (Monte Carlo noise or a bug):\n")
    print(x$reference_check[!x$reference_check$flag %in% FALSE, ], row.names = FALSE)
  }
  if (NROW(x$repeated_data)) {
    cat("\nScenarios where data sets repeat across replications (fixed rather than redrawn data?):\n")
    print(x$repeated_data, row.names = FALSE)
  }
  if (NROW(x$test_interval_disagreement)) {
    cat("\nTests and intervals of one method that disagree (p <= alpha but the interval contains the null, or the reverse):\n")
    print(x$test_interval_disagreement, row.names = FALSE)
  }
  if (NROW(x$degenerate)) {
    cat("\nDegenerate estimates:\n")
    print(x$degenerate, row.names = FALSE)
  }
  if (length(x$storage_problems)) {
    cat("\nStorage problems:\n", paste(" ", x$storage_problems, collapse = "\n"), "\n")
  }
  if (length(x$revisions)) cat(sprintf("\nRun revisions recorded: %d\n", length(x$revisions)))
  invisible(x)
}

#' Check a data-generating mechanism against its intended targets
#'
#' Generates `n` data sets for each scenario from a dedicated random number
#' stream and compares the mean of user supplied summaries (means,
#' prevalences, correlations, censoring or missingness proportions,
#' treatment allocation) with their intended targets. The tolerance is
#' uncertainty aware: a summary is flagged when its z statistic,
#' `(mean - target) / (sd / sqrt(n))`, exceeds a Bonferroni corrected
#' normal quantile, so ordinary sampling fluctuation is not reported as an
#' error. A flagged summary means the mechanism does not do what was
#' intended (or the target is wrong); deterministic configuration errors
#' are caught earlier by [sim_validate()].
#'
#' @param study A [sim_study()].
#' @param summaries Function of one generated data set returning a named
#'   numeric vector.
#' @param targets Function of `scenario` returning a named numeric vector of
#'   intended values (a subset of the summaries' names).
#' @param n Data sets per scenario.
#' @param seed Seed for the check.
#' @param alpha Familywise error rate for flagging.
#' @return A data frame with one row per scenario and summary.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 50L, p = c(0.1, 0.5)),
#'   generate = function(scenario, context) rbinom(scenario$n, 1, scenario$p),
#'   estimands = list(sim_estimand("p", sim_truth_exact(function(scenario) scenario$p))),
#'   methods = list(prop = function(data) sim_estimate("p", mean(data))))
#' sim_check_dgp(study, summaries = function(data) c(prevalence = mean(data)),
#'               targets = function(scenario) c(prevalence = scenario$p), n = 100, seed = 1)
#' @export
sim_check_dgp <- function(study, summaries, targets, n = 200L, seed = 1L, alpha = 0.01) {
  stopifnot(inherits(study, "adempr_study"))
  check_function(summaries, "summaries")
  check_function(targets, "targets", "scenario")
  # the same stream layout as a run (common random number groups and
  # antithetic pairs included), in a domain of its own
  ctx <- list(seed = seed, domain = "dgp_check",
              registry = rng_registry(rng_keys(study, "dgp_check", "dgp_check_method")),
              cache = new.env(parent = emptyenv()), dgp_stage = "dgp_check",
              method_stage = "dgp_check_method")
  anti <- isTRUE(study$generator$antithetic)
  out <- list()
  for (i in seq_along(study$scenarios$values)) {
    sc <- study$scenarios$values[[i]]
    key <- study$scenarios$keys[i]
    seeds <- task_seeds(ctx, key, seq_len(n), study, character())
    vals <- lapply(seq_len(n), function(r) {
      gctx <- list(replication = r, scenario_key = key, stage = "dgp_check")
      if (anti) gctx$antithetic <- r %% 2L == 0L
      g <- with_rng_state(seeds[[r]]$data, with_packages(study$generator$packages,
                                                           study$generator$fun(sc, gctx)))
      d <- if (inherits(g, "adempr_data")) g$data else g
      # user summaries never touch the caller's random numbers
      with_preserved_rng(summaries(d))
    })
    m <- do.call(rbind, vals)
    tg <- with_preserved_rng(targets(sc))
    for (nm in names(tg)) {
      if (!(nm %in% colnames(m))) stopf("Summary `%s` is not returned by `summaries`.", nm)
      x <- m[, nm]
      # antithetic pairs are dependent: the independent unit is the pair
      u <- if (anti) as.numeric(tapply(x, (seq_along(x) + 1L) %/% 2L, mean)) else x
      se <- stats::sd(u) / sqrt(length(u))
      zval <- if (se > 0) (mean(x) - tg[[nm]]) / se else if (mean(x) == tg[[nm]]) 0 else Inf
      out[[length(out) + 1L]] <- data.frame(
        scenario_id = study$scenarios$labels[i], summary = nm, target = tg[[nm]],
        mean = mean(x), mcse = se, z = zval, stringsAsFactors = FALSE)
    }
  }
  df <- bind_rows(out)
  crit <- stats::qnorm(1 - alpha / (2 * nrow(df)))
  df$flag <- abs(df$z) > crit
  attr(df, "critical_z") <- crit
  df
}

#' Plot simulation results
#'
#' Base graphics displays of a run:
#' * `"estimates"`: distribution of estimates by method, with the truth;
#' * `"est_se"`: estimates against their standard errors;
#' * `"zip"`: zip plot of confidence intervals ranked by their z statistic
#'   (Morris et al. 2019), showing which intervals miss the truth, on the
#'   same rows as coverage in [sim_summarise()];
#' * `"failures"`: proportion of usable results by method and scenario,
#'   with failures shown explicitly;
#' * `"precision"`: a performance measure with its Monte Carlo interval
#'   across scenarios and methods;
#' * `"runtime"`: distribution of elapsed time per attempt by method;
#' * `"bland_altman"`: agreement of two methods on the same data sets:
#'   difference of their estimates against their mean, with the mean
#'   difference and limits of agreement (mean +/- 1.96 SD);
#' * `"interval_truth"`: each interval drawn at its true value (useful when
#'   the truth varies between data sets), with the line of equality;
#'   intervals that miss the truth are highlighted.
#'
#' @param run A run or results.
#' @param type Plot type.
#' @param scenario Scenario (key, label or index) for `"estimates"`,
#'   `"est_se"` and `"zip"`; defaults to the first.
#' @param estimand Estimand; defaults to the first.
#' @param measure Measure for `"precision"`.
#' @param methods For `"bland_altman"`, the two methods compared (defaults
#'   to the first two that target the estimand).
#' @param reference For `"zip"`: rank intervals against the true value
#'   (`"truth"`, a coverage zip plot) or against the declared null value
#'   (`"null"`, a power zip plot of the intervals that include the null;
#'   they equal the nonrejections only when the method's test and interval
#'   agree). The dotted lines give the Monte Carlo confidence interval of
#'   the proportion of intervals that contain the reference.
#' @param ... Graphical parameters passed to the main plotting call.
#' @return `NULL`, invisibly.
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
#' sim_plot(run, "zip")
#' sim_plot(run, "precision", measure = "coverage")
#' @export
sim_plot <- function(run, type = c("estimates", "est_se", "zip", "failures", "precision",
                                   "runtime", "bland_altman", "interval_truth"),
                     scenario = 1L, estimand = NULL, measure = "bias",
                     reference = c("truth", "null"), methods = NULL, ...) {
  type <- match.arg(type)
  res <- sim_results(run)
  study <- res$manifest$study
  estimand <- estimand %||% names(study$estimands)[1]
  est <- res$estimates
  cols <- grDevices::hcl.colors(max(3L, length(study$methods)), "Dark 3")
  if (type %in% c("estimates", "est_se", "zip", "bland_altman", "interval_truth")) {
    key <- resolve_scenario(study, scenario)
    lab <- study$scenarios$labels[match(key, study$scenarios$keys)]
    sub <- est[est$scenario_key == key & est$estimand == estimand, , drop = FALSE]
    truth <- unique(sub$truth[!is.na(sub$truth)])
  }
  if (type == "estimates") {
    ok <- sub$valid & is.finite(sub$estimate)
    graphics::boxplot(estimate ~ method, data = sub[ok, , drop = FALSE], col = "grey92",
                      main = lab, ylab = paste("estimate of", estimand), xlab = "", ...)
    # fixed jitter that leaves the caller's random number stream untouched
    state <- rng_save()
    on.exit(rng_restore(state), add = TRUE)
    set.seed(1L)
    graphics::stripchart(estimate ~ method, data = sub[ok, , drop = FALSE], vertical = TRUE,
                         method = "jitter", add = TRUE, pch = 16, cex = 0.4,
                         col = grDevices::adjustcolor("black", 0.35))
    if (length(truth) == 1L) graphics::abline(h = truth, col = "firebrick", lwd = 2)
    nbad <- tapply(!ok, sub$method, sum)
    if (any(nbad > 0)) {
      graphics::mtext(paste("not usable:", paste(names(nbad), nbad, sep = " ", collapse = ", ")),
                      side = 3, line = 0.2, cex = 0.8)
    }
  } else if (type == "est_se") {
    ok <- sub$valid & is.finite(sub$estimate) & is.finite(sub$se)
    ms <- unique(sub$method)
    graphics::plot(sub$estimate[ok], sub$se[ok], col = cols[match(sub$method[ok], ms)],
                   pch = 16, cex = 0.6, xlab = "estimate", ylab = "standard error",
                   main = lab, ...)
    if (length(truth) == 1L) graphics::abline(v = truth, col = "firebrick", lty = 2)
    graphics::legend("topright", legend = ms, col = cols[seq_along(ms)], pch = 16, bty = "n")
  } else if (type == "bland_altman") {
    ok <- sub$valid & is.finite(sub$estimate)
    avail <- unique(sub$method[ok])
    methods <- methods %||% utils::head(avail, 2L)
    if (length(methods) != 2L || methods[1] == methods[2] || !all(methods %in% avail)) {
      stopf("`methods` must name two different methods with usable estimates for this scenario.")
    }
    a <- sub[ok & sub$method == methods[1], c("replication", "estimate")]
    b <- sub[ok & sub$method == methods[2], c("replication", "estimate")]
    m <- merge(a, b, by = "replication")
    if (nrow(m) < 2L) stopf("Fewer than two data sets have usable estimates from both methods.")
    avg <- (m$estimate.x + m$estimate.y) / 2
    dif <- m$estimate.x - m$estimate.y
    md <- mean(dif)
    lim <- md + c(-1.96, 1.96) * stats::sd(dif)
    graphics::plot(avg, dif, pch = 16, cex = 0.6, col = grDevices::adjustcolor("black", 0.5),
                   xlab = sprintf("mean of %s and %s", methods[1], methods[2]),
                   ylab = sprintf("%s minus %s", methods[1], methods[2]), main = lab,
                   ylim = range(c(dif, lim)), ...)
    graphics::abline(h = 0, col = "grey60")
    graphics::abline(h = md, col = "firebrick", lwd = 2)
    graphics::abline(h = lim, col = "firebrick", lty = 2)
    graphics::mtext(sprintf("%d paired data sets; mean difference %.3g, limits %.3g to %.3g",
                            nrow(m), md, lim[1], lim[2]), side = 3, line = 0.2, cex = 0.8)
  } else if (type == "interval_truth") {
    ok <- sub$interval_status == "valid" & sub$estimate_status %in% USABLE_FOR_INTERVAL &
      !is.na(sub$truth)
    r <- sub[ok, , drop = FALSE]
    if (!nrow(r)) stopf("No valid intervals with a truth to plot.")
    ms <- unique(r$method)
    op <- graphics::par(mfrow = c(1L, length(ms)), mar = c(4, 4, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    for (mm in ms) {
      q <- r[r$method == mm, , drop = FALSE]
      miss <- q$truth < q$lower | q$truth > q$upper
      # unbounded limits are drawn to the edge of the plot and marked
      fin <- c(q$lower[is.finite(q$lower)], q$upper[is.finite(q$upper)], q$truth)
      yr <- range(fin)
      pad <- max(diff(yr), abs(yr), 1) * 0.1
      lo_edge <- yr[1] - pad
      hi_edge <- yr[2] + pad
      lo <- ifelse(is.finite(q$lower), q$lower, lo_edge)
      hi <- ifelse(is.finite(q$upper), q$upper, hi_edge)
      graphics::plot(range(q$truth), c(lo_edge, hi_edge), type = "n",
                     xlab = "true value", ylab = "interval", main = paste(mm, lab, sep = "\n"), ...)
      graphics::segments(q$truth, lo, q$truth, hi,
                         col = ifelse(miss, "firebrick", grDevices::adjustcolor("grey30", 0.4)))
      unb <- !is.finite(q$lower) | !is.finite(q$upper)
      if (any(!is.finite(q$lower))) graphics::points(q$truth[!is.finite(q$lower)], lo[!is.finite(q$lower)], pch = 6, cex = 0.6)
      if (any(!is.finite(q$upper))) graphics::points(q$truth[!is.finite(q$upper)], hi[!is.finite(q$upper)], pch = 2, cex = 0.6)
      graphics::abline(0, 1, col = "steelblue", lwd = 2)
      graphics::mtext(sprintf("%d of %d miss%s", sum(miss), nrow(q),
                              if (any(unb)) sprintf("; %d unbounded (triangles)", sum(unb)) else ""),
                      side = 3, line = 0.2, cex = 0.8)
    }
  } else if (type == "zip") {
    reference <- match.arg(reference)
    # the rows coverage uses in sim_summarise(): generated data sets, an
    # estimate status that allows an interval, a valid interval
    reps <- res$replications
    gen_ok <- if (is.null(reps)) integer() else
      reps$replication[reps$scenario_key == key & reps$generation_status == "ok"]
    sub$ref <- if (reference == "truth") sub$truth else sub$null
    elig <- sub$replication %in% gen_ok & sub$estimate_status %in% USABLE_FOR_INTERVAL &
      sub$interval_status == "valid" & !is.na(sub$ref)
    ms <- unique(sub$method[elig])
    if (!length(ms)) stopf("No eligible intervals to plot.")
    op <- graphics::par(mfrow = c(1L, length(ms)), mar = c(4, 4, 3, 1))
    on.exit(graphics::par(op), add = TRUE)
    for (m in ms) {
      r <- sub[elig & sub$method == m, , drop = FALSE]
      # rank by the distance of the reference from the interval center in
      # half widths: the folded p value order for intervals of one level
      half <- (r$upper - r$lower) / 2
      zz <- abs((r$lower + r$upper) / 2 - r$ref) / ifelse(half > 0, half, NA)
      r <- r[order(zz), , drop = FALSE]
      n <- nrow(r)
      rank <- 100 * seq_len(n) / n
      cover <- r$lower <= r$ref & r$ref <= r$upper
      lo <- pmax(r$lower, -1e300)
      hi <- pmin(r$upper, 1e300)
      rng <- range(c(lo[is.finite(lo)], hi[is.finite(hi)], r$ref))
      graphics::plot(NA, xlim = rng, ylim = c(0, 100), xlab = "interval",
                     ylab = "centile of |z|", main = m, ...)
      graphics::segments(lo, rank, hi, rank, col = ifelse(cover, "grey60", "firebrick"))
      graphics::abline(v = unique(r$ref), lwd = 2)
      lv <- unique(r$level[!is.na(r$level)])
      if (length(lv) == 1L) graphics::abline(h = 100 * lv, lty = 2)
      ci <- cp_interval(sum(cover), n, study$level)
      graphics::abline(h = 100 * ci, lty = 3, col = "steelblue")
      graphics::mtext(sprintf("%s %.3f (MC %g%% CI %.3f to %.3f)%s",
                              if (reference == "truth") "coverage" else "intervals including the null",
                              mean(cover), 100 * study$level, ci[1], ci[2],
                              if (length(lv) > 1L) ", mixed levels" else ""),
                      side = 3, line = 0.2, cex = 0.7)
    }
  } else if (type == "failures") {
    f <- failure_table(res)
    f <- f[f$estimand == estimand, , drop = FALSE]
    prop <- ifelse(f$generated > 0, f$usable / f$generated, NA)
    tab <- tapply(prop, list(f$method, f$scenario_id), mean)
    graphics::barplot(tab, beside = TRUE, ylim = c(0, 1), col = cols[seq_len(nrow(tab))],
                      ylab = "proportion of usable results", las = 2, cex.names = 0.7, ...)
    graphics::legend("bottomright", legend = rownames(tab), fill = cols[seq_len(nrow(tab))],
                     bty = "n")
    graphics::abline(h = 1, lty = 3)
  } else if (type == "precision") {
    s <- as.data.frame(sim_summarise(res, measures = measure))
    s <- s[s$measure == measure & s$estimand == estimand, , drop = FALSE]
    y <- seq_len(nrow(s))
    xr <- range(c(s$lower, s$upper, s$value), finite = TRUE)
    graphics::plot(s$value, y, xlim = xr, yaxt = "n", pch = 16,
                   col = cols[match(s$method, names(study$methods))],
                   xlab = paste(measure, "with Monte Carlo interval"), ylab = "", ...)
    graphics::segments(s$lower, y, s$upper, y, col = cols[match(s$method, names(study$methods))])
    graphics::axis(2, at = y, labels = paste(s$scenario_id, s$method, sep = " | "), las = 1,
                   cex.axis = 0.6)
    if (measure == "coverage") {
      # the nominal level is the level the methods declared, not the Monte
      # Carlo confidence of the study
      # the levels of the intervals coverage counts, as in sim_summarise()
      ok_int <- est$estimand == estimand & est$interval_status == "valid" &
        est$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(est$truth)
      lv <- unique(est$level[ok_int & !is.na(est$level)])
      if (length(lv) == 1L) graphics::abline(v = lv, lty = 2) else
        graphics::mtext("methods declare different nominal levels", side = 3, cex = 0.8)
    }
    if (measure %in% c("bias", "relbias", "se_ratio")) graphics::abline(v = 0, lty = 2)
  } else if (type == "runtime") {
    att <- res$attempts
    att$elapsed <- pmax(att$elapsed, 1e-4)
    graphics::boxplot(elapsed ~ method, data = att, log = "y",
                      ylab = "seconds per attempt (floored at 0.1 ms)",
                      xlab = "", col = "grey92", ...)
  }
  invisible(NULL)
}

# Summary of the realized properties recorded for each generated data set.
record_summary <- function(res, lab) {
  r <- res$records
  if (is.null(r) || !nrow(r)) return(NULL)
  grp <- unique(r[c("scenario_key", "name")])
  out <- lapply(seq_len(nrow(grp)), function(i) {
    v <- r$value[r$scenario_key == grp$scenario_key[i] & r$name == grp$name[i]]
    ok <- is.finite(v)
    data.frame(scenario_id = lab(grp$scenario_key[i]), name = grp$name[i], n = sum(ok),
               mean = if (any(ok)) mean(v[ok]) else NA_real_,
               sd = if (sum(ok) > 1) stats::sd(v[ok]) else NA_real_,
               min = if (any(ok)) min(v[ok]) else NA_real_,
               max = if (any(ok)) max(v[ok]) else NA_real_,
               missing = sum(!ok), stringsAsFactors = FALSE)
  })
  bind_rows(out)
}

#' Check methods with known performance against Monte Carlo error
#'
#' A method declared with `sim_method(reference = c(bias = 0, coverage =
#' 0.95))` has performance known from theory (for example the sample mean
#' under a normal model). Its estimated performance must agree with the
#' declared value within Monte Carlo error; a miss beyond `threshold` MCSE
#' points to Monte Carlo noise larger than reported, or to a bug in the
#' generator, the truth or the evaluation. `sim_diagnose()` runs this
#' check.
#'
#' A proportion (coverage, rejection and the like) is checked against its
#' exact (Clopper and Pearson) interval at the confidence level that
#' corresponds to `threshold`, so a proportion of 0 or 1, whose plug-in
#' MCSE is 0, is still checked. Other measures use `mcse_total` when the
#' truth is uncertain and `mcse` otherwise. An MCSE of 0 is a pass when the
#' value equals the reference and a flag otherwise. A cell without an MCSE
#' gives `flag = NA` (insufficient evidence). A reference given as a data
#' frame with columns `measure`, `value` and optionally `estimand` and
#' `scenario_id` applies only to the rows it names.
#'
#' @param run A run or results.
#' @param threshold Absolute z above which a cell is flagged.
#' @return A data frame with one row per scenario, reference method,
#'   estimand and declared measure: `value`, `mcse` (the one used),
#'   `reference`, `z`, the exact interval `lower` and `upper` for
#'   proportions, `check` and `flag`; `NULL` when no method declares a
#'   reference.
#' @export
sim_check_reference <- function(run, threshold = 3) {
  if (!is_scalar_number(threshold) || threshold <= 0) stopf("`threshold` must be a positive number.")
  res <- sim_results(run)
  study <- res$manifest$study
  refm <- Filter(function(m) !is.null(m$reference), study$methods)
  if (!length(refm)) return(NULL)
  refs <- lapply(refm, function(m) reference_table(m$reference))
  ms <- unique(unlist(lapply(refs, `[[`, "measure")))
  s <- as.data.frame(sim_summarise(res, measures = ms))
  alpha <- 2 * stats::pnorm(-threshold)
  rows <- list()
  for (mn in names(refm)) {
    ref <- refs[[mn]]
    bad_e <- setdiff(stats::na.omit(ref$estimand), names(study$estimands))
    bad_s <- setdiff(stats::na.omit(ref$scenario_id), study$scenarios$labels)
    if (length(bad_e) || length(bad_s)) {
      stopf("The reference of method `%s` names unknown estimand(s) or scenario(s): %s.", mn,
            paste(c(bad_e, bad_s), collapse = ", "))
    }
    for (k in seq_len(nrow(ref))) {
      r <- s[s$method == mn & s$measure == ref$measure[k], , drop = FALSE]
      if (!is.na(ref$estimand[k])) r <- r[r$estimand == ref$estimand[k], , drop = FALSE]
      if (!is.na(ref$scenario_id[k])) r <- r[r$scenario_id == ref$scenario_id[k], , drop = FALSE]
      if (!nrow(r)) next
      r$reference <- ref$value[k]
      rows[[length(rows) + 1L]] <- r
    }
  }
  out <- bind_rows(rows)
  if (is.null(out)) return(NULL)
  se <- if (is.null(out$mcse_total)) out$mcse else ifelse(is.na(out$mcse_total), out$mcse, out$mcse_total)
  d <- out$value - out$reference
  z <- ifelse(se > 0, d / se, ifelse(abs(d) <= 1e-12 * pmax(1, abs(out$reference)), 0, sign(d) * Inf))
  flag <- !is.na(z) & abs(z) > threshold
  flag[is.na(z) | is.na(out$value)] <- NA
  check <- ifelse(is.na(flag), "insufficient evidence", ifelse(se %in% 0, "zero MCSE", "normal"))
  lo <- up <- rep(NA_real_, nrow(out))
  prop <- out$measure %in% PROPORTION_MEASURES & !is.na(out$value) & out$n_eligible > 0
  if (isTRUE(study$generator$antithetic)) {
    # dependent pairs: a binomial interval over the replications does not
    # hold. A cell with all pair means equal uses the bounded pair interval
    # of the summary (at the study level); the others are judged by the
    # pair based MCSE above.
    bnd <- prop & out$status %in% "boundary" & !is.na(out$lower) & !is.na(out$upper)
    lo[bnd] <- out$lower[bnd]
    up[bnd] <- out$upper[bnd]
    flag[bnd] <- out$reference[bnd] < lo[bnd] | out$reference[bnd] > up[bnd]
    check[bnd] <- "bounded pair interval"
    prop <- rep(FALSE, length(prop))
  }
  if (any(prop)) {
    n <- out$n_eligible[prop]
    x <- round(out$value[prop] * n)
    lo[prop] <- ifelse(x == 0, 0, stats::qbeta(alpha / 2, x, n - x + 1))
    up[prop] <- ifelse(x == n, 1, stats::qbeta(1 - alpha / 2, x + 1, n - x))
    flag[prop] <- out$reference[prop] < lo[prop] | out$reference[prop] > up[prop]
    check[prop] <- "exact interval"
  }
  out <- data.frame(out[c("scenario_id", "method", "estimand", "measure", "value")],
                    mcse = se, reference = out$reference, z = z, lower = lo, upper = up,
                    check = check, flag = flag, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# A method's reference as a table: a named vector applies to every
# scenario and estimand; a data frame may restrict rows to an estimand
# and a scenario label.
reference_table <- function(ref) {
  if (is.list(ref) && all(c("measure", "value") %in% names(ref)) && !is.list(ref$measure)) {
    data.frame(measure = as.character(ref$measure), value = as.numeric(ref$value),
               estimand = if (is.null(ref$estimand)) NA_character_ else as.character(ref$estimand),
               scenario_id = if (is.null(ref$scenario_id)) NA_character_ else as.character(ref$scenario_id),
               stringsAsFactors = FALSE)
  } else {
    ref <- unlist(ref)
    data.frame(measure = names(ref), value = unname(as.numeric(ref)), estimand = NA_character_,
               scenario_id = NA_character_, stringsAsFactors = FALSE)
  }
}
