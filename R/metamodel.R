# Descriptive meta-model ---------------------------------------------------

#' Describe how much each factor moves a performance measure
#'
#' Fits a linear model of the cell values of one measure on the scenario
#' factors (treated as categorical, main effects and optionally two-way
#' interactions), weighted by the inverse squared Monte Carlo standard
#' error, and reports each term's share of the explained sum of squares.
#' The result is descriptive: it summarizes the cells that were run, it
#' does not predict performance in scenarios that were not run, and a
#' misspecified model attributes its own lack of fit to the factors (Morris
#' 2024). The cell estimates are returned beside the model and should be
#' reported with it.
#'
#' @param run A run or results.
#' @param measure Measure.
#' @param method Method; defaults to the first.
#' @param estimand Estimand; defaults to the first.
#' @param factors Factors to include; defaults to every factor that varies.
#' @param interactions Include two-way interactions.
#' @return A list of class `adempr_factor_influence` with `anova` (term,
#'   degrees of freedom, sum of squares, share of the explained sum of
#'   squares), `cells` (the cell estimates used) and `r_squared`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L), sd = c(1, 3)),
#'   generate = function(scenario, context) rnorm(scenario$n, 0, scenario$sd),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(mean = function(data) sim_estimate("mean", mean(data))),
#'   measures = "empse")
#' run <- sim_run(study, replications = 50, seed = 1, path = tempfile("adempR-run-"))
#' sim_factor_influence(run, "empse")
#' @export
sim_factor_influence <- function(run, measure, method = NULL, estimand = NULL, factors = NULL,
                                 interactions = FALSE) {
  res <- sim_results(run)
  study <- res$manifest$study
  method <- method %||% names(study$methods)[1]
  estimand <- estimand %||% names(study$estimands)[1]
  s <- as.data.frame(sim_summarise(res, measures = measure))
  s <- s[s$method == method & s$estimand == estimand & s$measure == measure &
           is.finite(s$value), , drop = FALSE]
  sc <- as.data.frame(study$scenarios)
  facs <- setdiff(names(sc), c("scenario_key", "scenario_id"))
  facs <- facs[vapply(facs, function(f) !is.list(sc[[f]]) && length(unique(sc[[f]])) > 1L, logical(1))]
  factors <- factors %||% facs
  if (!length(factors) || !all(factors %in% facs)) {
    stopf("`factors` must name scenario factors that vary: %s.", paste(facs, collapse = ", "))
  }
  d <- merge(s[c("scenario_key", "value", "mcse")], sc[c("scenario_key", factors)], by = "scenario_key")
  for (f in factors) d[[f]] <- factor(d[[f]])
  rhs <- paste(factors, collapse = " + ")
  if (interactions && length(factors) > 1L) rhs <- paste0("(", rhs, ")^2")
  w <- if (all(is.finite(d$mcse) & d$mcse > 0)) 1 / d$mcse^2 else NULL
  form <- stats::as.formula(paste("value ~", rhs))
  fit <- stats::lm(form, data = d, weights = w)
  a <- stats::anova(fit)
  terms <- rownames(a)
  ss <- a[["Sum Sq"]]
  expl <- terms != "Residuals"
  tab <- data.frame(term = trimws(terms), df = a$Df, sum_sq = ss,
                    share = ifelse(expl, ss / sum(ss[expl]), NA_real_), stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  r2 <- if (sum(ss) > 0) sum(ss[expl]) / sum(ss) else NA_real_
  structure(list(anova = tab, cells = d, r_squared = r2, measure = measure, method = method,
                 estimand = estimand, weighted = !is.null(w)), class = "adempr_factor_influence")
}

#' @export
print.adempr_factor_influence <- function(x, ...) {
  cat(sprintf("Descriptive meta-model of %s for %s (%s), %d cells%s; R^2 = %.3f.\n",
              x$measure, x$method, x$estimand, nrow(x$cells),
              if (x$weighted) ", weighted by 1/MCSE^2" else "", x$r_squared))
  cat("It summarizes the cells run; it does not predict other scenarios. Report the cells too.\n")
  print(x$anova, row.names = FALSE, digits = 3)
  invisible(x)
}
