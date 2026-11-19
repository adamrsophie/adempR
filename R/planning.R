# Precision planning and pilots --------------------------------------------

#' Plan the number of replications
#'
#' For a proportion (coverage, rejection, usable result) with anticipated
#' value `p`, a Monte Carlo standard error `mcse` needs
#' `R >= p (1 - p) / mcse^2` independent replications; with `p` unknown
#' the conservative `p = 0.5` is used. For example `p = 0.95` and
#' `mcse = 0.005` give 1900 replications. A confidence interval half width
#' `half_width` is a different request: `R >= z^2 p (1 - p) / half_width^2`.
#' For a mean (bias) with anticipated standard deviation `sd` of the
#' contributions, `R >= sd^2 / mcse^2`.
#'
#' @param type `"proportion"` or `"mean"`.
#' @param p Anticipated proportion.
#' @param sd Anticipated standard deviation of the per replication
#'   contributions (for bias, the empirical SE).
#' @param mcse Target Monte Carlo standard error.
#' @param half_width Target half width of a normal approximation interval.
#' @param level Confidence level for `half_width`.
#' @return The required number of replications (an integer), with the
#'   inputs as attributes.
#' @examples
#' sim_plan_replications("proportion", p = 0.95, mcse = 0.005)  # 1900
#' sim_plan_replications("proportion", mcse = 0.005)            # worst case
#' sim_plan_replications("proportion", p = 0.95, half_width = 0.01)
#' sim_plan_replications("mean", sd = 0.2, mcse = 0.005)
#' @export
sim_plan_replications <- function(type = c("proportion", "mean"), p = NULL, sd = NULL,
                                  mcse = NULL, half_width = NULL, level = 0.95) {
  type <- match.arg(type)
  if (is.null(mcse) == is.null(half_width)) stopf("Give exactly one of `mcse` and `half_width`.")
  target <- mcse %||% half_width
  if (!is_scalar_number(target) || target <= 0) stopf("The precision target must be positive.")
  v <- if (type == "proportion") {
    if (is.null(p)) 0.25 else {
      if (!is_scalar_number(p) || p < 0 || p > 1) stopf("`p` must be in [0, 1].")
      p * (1 - p)
    }
  } else {
    if (!is_scalar_number(sd) || sd <= 0) stopf("`sd` must be a positive number for a mean.")
    sd^2
  }
  mult <- if (is.null(half_width)) 1 else stats::qnorm(1 - (1 - level) / 2)^2
  R <- ceiling(mult * v / target^2 - 1e-9)
  if (R > .Machine$integer.max) {
    stopf("The target needs %s replications, more than can be counted; relax the target.",
          format(R, big.mark = ","))
  }
  R <- as.integer(R)
  structure(max(R, 1L), type = type, p = p, sd = sd, mcse = mcse, half_width = half_width)
}

# Conservative p(1 - p) from a pilot: its maximum over the exact interval.
conservative_pq <- function(k, n, level = 0.95) {
  if (n == 0) return(0.25)
  a <- 1 - level
  lo <- if (k == 0) 0 else stats::qbeta(a / 2, k, n - k + 1)
  hi <- if (k == n) 1 else stats::qbeta(1 - a / 2, k + 1, n - k)
  if (lo <= 0.5 && hi >= 0.5) return(0.25)
  pm <- if (hi < 0.5) hi else lo
  pm * (1 - pm)
}

#' Run a pilot study
#'
#' Runs a small number of replications in a separate run directory labeled
#' `"pilot"`, then estimates the cost per task with its uncertainty, the
#' projected total cost and storage of the planned production run, failure
#' rates with exact upper bounds, and the replications needed to meet the
#' precision targets. For proportions the required replications use the
#' largest `p (1 - p)` within the pilot's exact interval, so a pilot that
#' happens to observe 0 or 100 percent cannot claim that almost no
#' replications are needed.
#'
#' Pilot results are not reused in the production run. A pilot directory
#' draws from its own family of random number streams, so a production run
#' with the same seed is independent of the pilot; reusing pilot data would
#' make the production sample depend on what the pilot showed. Every
#' measure named in `targets` is summarized; an unknown measure is an
#' error, and a target for which no requirement can be computed gives a
#' warning.
#'
#' @param study A [sim_study()].
#' @param replications Pilot replications per scenario.
#' @param seed Pilot seed.
#' @param planned Planned production replications per scenario (defaults
#'   to the largest requirement from the precision targets, if any).
#' @param targets Named vector of target MCSEs by measure, for example
#'   `c(coverage = 0.005, bias = 0.01)`; defaults to the protocol's
#'   `precision_targets`.
#' @param path Directory for the pilot run.
#' @param ... Passed to [sim_run()] (for example `backend`).
#' @return An object of class `adempr_pilot` with elements `run`,
#'   `summary`, `cost`, `precision` and `projection`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) {
#'     tt <- t.test(data)
#'     sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
#'                  lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95)
#'   }))
#' pilot <- sim_pilot(study, replications = 30, seed = 99,
#'                    targets = c(coverage = 0.01, bias = 0.01))
#' pilot
#' @export
sim_pilot <- function(study, replications = 50L, seed, planned = NULL, targets = NULL,
                      path = tempfile("adempR-pilot-"), ...) {
  stopifnot(inherits(study, "adempr_study"))
  run <- sim_run(study, replications = replications, seed = seed, path = path,
                 label = "pilot", ...)
  targets <- targets %||% study$protocol$precision_targets
  # Summarize with every measure that has a target; an unknown measure name
  # is an error from sim_measures().
  summ <- sim_summarise(run, measures = do.call(sim_measures, as.list(unique(c(
    "bias", "empse", "coverage", "rejection", "usable", names(targets))))))
  res <- sim_results(run)
  att <- res$attempts
  reps <- res$replications
  # cost per task: generation plus all methods (0 for data sets no method saw)
  per_task <- reps[c("scenario_key", "replication", "generation_elapsed")]
  per_task$methods_seconds <- 0
  if (!is.null(att) && nrow(att)) {
    ms <- stats::aggregate(att$elapsed, by = list(scenario_key = att$scenario_key,
                                                  replication = att$replication), FUN = sum)
    i <- match(paste(per_task$scenario_key, per_task$replication),
               paste(ms$scenario_key, ms$replication))
    per_task$methods_seconds[!is.na(i)] <- ms$x[i[!is.na(i)]]
  }
  per_task$task_seconds <- per_task$methods_seconds + per_task$generation_elapsed
  cost <- do.call(rbind, lapply(split(per_task, per_task$scenario_key), function(d) {
    data.frame(scenario_key = d$scenario_key[1],
               scenario_id = study$scenarios$labels[match(d$scenario_key[1], study$scenarios$keys)],
               mean_task_seconds = mean(d$task_seconds),
               se_task_seconds = stats::sd(d$task_seconds) / sqrt(nrow(d)),
               mean_generation_seconds = mean(d$generation_elapsed),
               mean_method_seconds = mean(d$methods_seconds), stringsAsFactors = FALSE)
  }))
  sizes <- file.info(list.files(file.path(path, "batches"), pattern = "\\.rds$",
                                full.names = TRUE))$size
  bytes_per_task <- sum(sizes) / max(1, nrow(per_task))
  s <- as.data.frame(summ)
  prec <- NULL
  if (length(targets)) {
    rows <- list()
    for (ms in names(targets)) {
      sub <- s[s$measure == ms, , drop = FALSE]
      for (i in seq_len(nrow(sub))) {
        need <- if (is.na(sub$value[i]) || sub$n_eligible[i] == 0) {
          NA_real_
        } else if (ms %in% c("coverage", "rejection", "usable", "noncoverage_lower",
                              "noncoverage_upper", "be_coverage")) {
          k <- round(sub$value[i] * sub$n_eligible[i])
          ceiling(conservative_pq(k, sub$n_eligible[i]) / targets[[ms]]^2)
        } else if (!is.na(sub$mcse[i]) && sub$n_eligible[i] > 1) {
          ceiling((sub$mcse[i]^2 * sub$n_eligible[i]) / targets[[ms]]^2)
        } else NA_real_
        rows[[length(rows) + 1L]] <- data.frame(
          scenario_id = sub$scenario_id[i], method = sub$method[i], estimand = sub$estimand[i],
          measure = ms, pilot_value = sub$value[i], pilot_n = sub$n_eligible[i],
          target_mcse = targets[[ms]], required_replications = need, stringsAsFactors = FALSE)
      }
    }
    prec <- bind_rows(rows)
    unplanned <- setdiff(names(targets),
                         prec$measure[!is.na(prec$required_replications)])
    if (length(unplanned)) {
      warnf("No replication requirement could be computed for target(s) %s: the pilot gave no eligible results or no Monte Carlo standard error for them.",
            paste(unplanned, collapse = ", "))
    }
  }
  if (is.null(planned)) {
    planned <- if (!is.null(prec) && any(!is.na(prec$required_replications))) {
      max(prec$required_replications, na.rm = TRUE)
    } else NA_integer_
  }
  n_tasks <- planned * length(study$scenarios)
  projection <- data.frame(
    planned_replications = planned, scenarios = length(study$scenarios),
    tasks = n_tasks,
    projected_seconds = if (is.na(planned)) NA_real_ else sum(cost$mean_task_seconds) * planned,
    projected_seconds_se = if (is.na(planned)) NA_real_ else
      sqrt(sum(cost$se_task_seconds^2)) * planned,
    projected_megabytes = if (is.na(planned)) NA_real_ else bytes_per_task * n_tasks / 1e6,
    stringsAsFactors = FALSE)
  structure(list(run = run, summary = summ, cost = cost, precision = prec,
                 projection = projection, failures = attr(summ, "failures")),
            class = "adempr_pilot")
}

#' @export
print.adempr_pilot <- function(x, ...) {
  cat(sprintf("<adempR pilot: %d replications per scenario>\n", x$run$manifest$replications))
  cat("\nCost per task (seconds):\n")
  print(x$cost[c("scenario_id", "mean_task_seconds", "se_task_seconds")], row.names = FALSE)
  if (!is.null(x$precision)) {
    cat("\nReplications needed for the precision targets:\n")
    print(x$precision[c("scenario_id", "method", "measure", "pilot_value", "target_mcse",
                        "required_replications")], row.names = FALSE)
  }
  cat("\nProjection for the production run:\n")
  print(x$projection, row.names = FALSE)
  f <- x$failures
  if (any(f$usable < f$generated)) {
    cat("\nMethods with unusable results in the pilot:\n")
    print(f[f$usable < f$generated, c("scenario_id", "method", "generated", "usable", "statuses")],
          row.names = FALSE)
  }
  cat("\nPilot results are not part of the production run.\n")
  invisible(x)
}

#' Allocate replications adaptively to reach a precision target
#'
#' Runs `initial` replications per scenario, then repeatedly estimates the
#' Monte Carlo standard error of `measure` in every scenario and gives each
#' scenario that misses `target` the replications it needs (from MCSE
#' proportional to `1 / sqrt(R)`, inflated by `inflate`), up to
#' `max_replications`, until every scenario meets the target or reaches
#' the cap. Scenarios that are already precise get no more replications, so
#' the budget goes where the Monte Carlo error is. The earlier replications
#' are kept (Morris's proposal for choosing the number of imputations
#' applied to scenarios). The allocation uses the estimated variability
#' only; it does not look at the value of the measure, so it does not
#' steer the estimate. Every round is recorded in the manifest.
#'
#' @param study A study.
#' @param measure Measure whose MCSE is targeted.
#' @param target Target MCSE.
#' @param seed,path Passed to [sim_run()].
#' @param initial Replications per scenario in the first round.
#' @param max_replications Cap per scenario.
#' @param max_rounds Maximum number of allocation rounds.
#' @param inflate Factor applied to the projected requirement.
#' @param ... Passed to [sim_run()].
#' @return The run, with `run$allocation`: one row per round and scenario
#'   giving the replications, the worst MCSE over methods and estimands, and
#'   whether the target is met.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(sd = c(1, 4)),
#'   generate = function(scenario, context) rnorm(20, 0, scenario$sd),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data))),
#'   measures = "bias")
#' run <- sim_run_adaptive(study, "bias", target = 0.05, seed = 1,
#'                         path = tempfile("adempR-run-"), initial = 20)
#' run$allocation
#' @export
sim_run_adaptive <- function(study, measure, target, seed, path, initial = 100L,
                             max_replications = 1e4, max_rounds = 5L, inflate = 1.1, ...) {
  stopifnot(inherits(study, "adempr_study"))
  if (!is_scalar_number(target) || target <= 0) stopf("`target` must be a positive number.")
  if (!is_count(initial, 2)) stopf("`initial` must be a whole number of at least 2.")
  if (!is_count(max_replications, initial)) stopf("`max_replications` must be at least `initial`.")
  run <- sim_run(study, replications = initial, seed = seed, path = path, ...)
  log <- list()
  for (round in seq_len(max_rounds + 1L)) {
    s <- as.data.frame(sim_summarise(run, measures = measure))
    s <- s[s$measure == measure, , drop = FALSE]
    worst <- tapply(s$mcse, s$scenario_key, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    have <- tapply(s$n_planned, s$scenario_key, max)
    keys <- names(worst)
    met <- !is.na(worst) & worst <= target
    capped <- have >= max_replications
    log[[round]] <- data.frame(round = round, scenario_id = study$scenarios$labels[match(keys, study$scenarios$keys)],
                               replications = as.integer(have), mcse = as.numeric(worst),
                               target_met = met, stringsAsFactors = FALSE)
    todo <- !met & !capped & !is.na(worst)
    if (!any(todo) || round > max_rounds) break
    need <- ceiling(have[todo] * (worst[todo] / target)^2 * inflate)
    need <- pmin(pmax(need, have[todo] + 1L), max_replications)
    run <- sim_run(study, replications = initial, seed = seed, path = path,
                   scenario_replications = stats::setNames(as.integer(need), keys[todo]), ...)
  }
  alloc <- bind_rows(log)
  rownames(alloc) <- NULL
  update_manifest(path, function(man) {
    man$adaptive <- list(measure = measure, target = target, allocation = alloc)
    man
  })
  run <- sim_open(path)
  run$allocation <- alloc
  run
}
