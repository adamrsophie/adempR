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
  R <- as.integer(ceiling(mult * v / target^2 - 1e-9))
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
#' Pilot results are not reused in the production run. Use a different
#' seed (or directory) for production; reusing pilot data would make the
#' production sample depend on what the pilot showed.
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
  summ <- sim_summarise(run, measures = sim_measures("bias", "empse", "coverage",
                                                     "rejection", "usable"))
  res <- sim_results(run)
  targets <- targets %||% study$protocol$precision_targets
  att <- res$attempts
  reps <- res$replications
  # cost per task: generation plus all methods
  per_task <- stats::aggregate(att$elapsed, by = list(scenario_key = att$scenario_key,
                                                      replication = att$replication), FUN = sum)
  names(per_task)[3] <- "methods_seconds"
  per_task <- merge(per_task, reps[c("scenario_key", "replication", "generation_elapsed")],
                    by = c("scenario_key", "replication"))
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
