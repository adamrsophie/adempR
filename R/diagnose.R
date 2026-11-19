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
    a <- ifelse(is.na(est$alpha), 1 - est$level, est$alpha)
    both <- est$test_status == "valid" & est$interval_status == "valid" &
      est$alternative == "two.sided" & !is.na(est$null) & !is.na(est$level) &
      !is.na(est$p_value) & abs(a - (1 - est$level)) < 1e-12
    both[is.na(both)] <- FALSE
    if (any(both)) {
      b <- est[both, , drop = FALSE]
      excl <- b$null < b$lower | b$null > b$upper
      b$disagree <- (b$p_value <= a[both]) != excl
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
  structure(list(
    failures = failure_table(res), messages = msgs, generation_failures = gen_fail,
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
  keys <- data.frame(stage = "dgp_check", scenario_key = study$scenarios$keys,
                     component = "data", stringsAsFactors = FALSE)
  reg <- rng_registry(keys)
  base <- rng_base_seed(seed, "dgp_check")
  cache <- new.env(parent = emptyenv())
  out <- list()
  for (i in seq_along(study$scenarios$values)) {
    sc <- study$scenarios$values[[i]]
    key <- study$scenarios$keys[i]
    k <- rng_lookup(reg, "dgp_check", key, "data")
    seeds <- rng_substreams(rng_stream_seeds(base, k, cache)[[1]], seq_len(n))
    vals <- lapply(seq_len(n), function(r) {
      g <- with_rng_state(seeds[[r]], study$generator$fun(sc, list(
        replication = r, scenario_key = key, stage = "dgp_check")))
      d <- if (inherits(g, "adempr_data")) g$data else g
      summaries(d)
    })
    m <- do.call(rbind, vals)
    tg <- targets(sc)
    for (nm in names(tg)) {
      if (!(nm %in% colnames(m))) stopf("Summary `%s` is not returned by `summaries`.", nm)
      x <- m[, nm]
      se <- stats::sd(x) / sqrt(length(x))
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
#' * `"runtime"`: distribution of elapsed time per attempt by method.
#'
#' @param run A run or results.
#' @param type Plot type.
#' @param scenario Scenario (key, label or index) for `"estimates"`,
#'   `"est_se"` and `"zip"`; defaults to the first.
#' @param estimand Estimand; defaults to the first.
#' @param measure Measure for `"precision"`.
#' @param reference For `"zip"`: rank intervals against the true value
#'   (`"truth"`, a coverage zip plot) or against the declared null value
#'   (`"null"`, a power zip plot, where intervals that exclude the null are
#'   the rejections). The dotted lines give the Monte Carlo confidence
#'   interval of the proportion of intervals that contain the reference.
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
                                   "runtime"),
                     scenario = 1L, estimand = NULL, measure = "bias",
                     reference = c("truth", "null"), ...) {
  type <- match.arg(type)
  res <- sim_results(run)
  study <- res$manifest$study
  estimand <- estimand %||% names(study$estimands)[1]
  est <- res$estimates
  cols <- grDevices::hcl.colors(max(3L, length(study$methods)), "Dark 3")
  if (type %in% c("estimates", "est_se", "zip")) {
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
                              if (reference == "truth") "coverage" else "1 - rejection",
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
      lv <- unique(est$level[est$estimand == estimand & !is.na(est$level)])
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
