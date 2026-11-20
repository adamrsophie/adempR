# Descriptive meta-model ---------------------------------------------------

#' Describe how much each factor moves a performance measure
#'
#' Fits a linear model of the cell values of one measure on the scenario
#' factors (treated as categorical, main effects and optionally two-way
#' interactions), weighted by the inverse squared Monte Carlo standard
#' error, and reports each term's sum of squares given every term that
#' does not contain it (type II) and its share of the total of these. The
#' attribution does not depend on the order of the factors, also in
#' unbalanced, constrained or weighted designs, where the shares need not
#' add up to the explained part. When a cell has no positive finite Monte
#' Carlo standard error, all cells are weighted equally and the output says
#' so; cells without a finite value are excluded and counted, and
#' coefficients that the cells cannot estimate (aliasing) are counted.
#' Variation is judged on the retained cells: a level without any retained
#' cell is reported, a factor left with one level is dropped, and the
#' result is marked `reduced`, since it then describes the retained cells
#' and not the full design; with fewer than two retained cells or no
#' varying factor there is no model (`unavailable`).
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
#'   degrees of freedom, type II sum of squares, share), `cells` (the cell
#'   estimates used), `r_squared`, `weighting`, `excluded` (cells without a
#'   finite value), `aliased` (inestimable coefficients), `status`
#'   (`"ok"`, `"reduced"` when exclusions removed cells, levels or factors,
#'   or `"unavailable"` when no model can be fitted), `dropped` (factors
#'   that no longer vary), `removed_levels` and `note`.
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
  n_cells <- length(study$scenarios$keys)
  sc <- as.data.frame(study$scenarios)
  # factor names come from the design, not from the columns the data frame
  # adds (the scenario weight)
  facs <- unique(unlist(lapply(study$scenarios$values, names)))
  facs <- facs[vapply(facs, function(f) !is.list(sc[[f]]) && length(unique(sc[[f]])) > 1L, logical(1))]
  factors <- factors %||% facs
  if (!length(factors) || !all(factors %in% facs)) {
    stopf("`factors` must name scenario factors that vary: %s.", paste(facs, collapse = ", "))
  }
  d <- merge(s[c("scenario_key", "value", "mcse")], sc[c("scenario_key", factors)], by = "scenario_key")
  # The cells left after exclusions (failed cells, no finite value) can
  # lack whole factor levels. Variation is judged again on them: levels
  # that are gone are reported, factors left with one level are dropped,
  # and the analysis is labeled as reduced, because it no longer answers
  # the question about the full design.
  removed <- bind_rows(lapply(factors, function(f) {
    gone <- setdiff(unique(as.character(sc[[f]])), unique(as.character(d[[f]])))
    if (length(gone)) data.frame(factor = f, level = gone, stringsAsFactors = FALSE)
  }))
  if (is.null(removed)) removed <- data.frame(factor = character(), level = character())
  keep_f <- factors[vapply(factors, function(f) length(unique(d[[f]])) > 1L, logical(1))]
  dropped <- setdiff(factors, keep_f)
  status <- if (length(dropped) || nrow(removed) || nrow(d) < n_cells) "reduced" else "ok"
  note <- c(if (nrow(d) < n_cells) sprintf("%d of %d cells excluded (no finite value)", n_cells - nrow(d), n_cells),
            if (nrow(removed)) paste("levels without any retained cell:",
                                     paste(paste0(removed$factor, "=", removed$level), collapse = ", ")),
            if (length(dropped)) paste("factors that no longer vary, dropped:", paste(dropped, collapse = ", ")))
  if (nrow(d) < 2L || !length(keep_f)) {
    return(structure(list(anova = NULL, cells = d, r_squared = NA_real_, measure = measure,
                          method = method, estimand = estimand, weighted = FALSE,
                          weighting = NA_character_, excluded = n_cells - nrow(d), aliased = NA_integer_,
                          status = "unavailable", dropped = dropped, removed_levels = removed,
                          note = paste(c(note, "too few retained cells or no varying factor; no model"),
                                       collapse = "; ")),
                     class = "adempr_factor_influence"))
  }
  factors <- keep_f
  # syntactic stand-ins, so any factor name (with spaces or symbols) works
  # in a formula
  vars <- paste0(".f", seq_along(factors))
  for (k in seq_along(factors)) d[[vars[k]]] <- factor(d[[factors[k]]])
  labels <- if (interactions && length(vars) > 1L) {
    c(vars, apply(utils::combn(vars, 2L), 2L, paste, collapse = ":"))
  } else vars
  ok_w <- is.finite(d$mcse) & d$mcse > 0
  w <- if (nrow(d) && all(ok_w)) 1 / d$mcse^2 else rep(1, nrow(d))
  weighting <- if (nrow(d) && all(ok_w)) "inverse squared MCSE" else
    sprintf("equal (%d cell(s) without a positive finite MCSE)", sum(!ok_w))
  fit_of <- function(lab) {
    f <- stats::reformulate(if (length(lab)) lab else "1", response = "value")
    m <- stats::lm(f, data = d, weights = w)
    c(rss = stats::deviance(m), rank = m$rank)
  }
  # Each term's sum of squares given every term that does not contain it
  # (type II): it does not depend on the order of the factors, also when
  # the design is unbalanced or weighted. Its degrees of freedom are the
  # parameters it adds, which an incomplete design can reduce.
  parts <- strsplit(labels, ":", fixed = TRUE)
  tt <- vapply(seq_along(labels), function(t) {
    contains <- vapply(parts, function(p) all(parts[[t]] %in% p), logical(1))
    without <- labels[!contains]
    a <- fit_of(without)
    b <- fit_of(c(without, labels[t]))
    c(a[["rss"]] - b[["rss"]], b[["rank"]] - a[["rank"]])
  }, numeric(2))
  ss <- tt[1, ]
  df_t <- tt[2, ]
  fit <- stats::lm(stats::reformulate(labels, response = "value"), data = d, weights = w)
  aliased <- sum(is.na(stats::coef(fit)))
  rename <- function(x) {
    for (k in seq_along(vars)) x <- gsub(paste0("\\", vars[k], "\\b"), factors[k], x)
    x
  }
  tss <- fit_of(character())[["rss"]]
  rss_full <- stats::deviance(fit)
  tab <- data.frame(term = c(rename(labels), "Residuals"), df = c(df_t, stats::df.residual(fit)),
                    sum_sq = c(ss, rss_full), share = c(ss / sum(ss), NA_real_),
                    stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  r2 <- if (tss > 0) 1 - rss_full / tss else NA_real_
  d <- d[setdiff(names(d), vars)]
  structure(list(anova = tab, cells = d, r_squared = r2, measure = measure, method = method,
                 estimand = estimand, weighted = nrow(d) > 0 && all(ok_w), weighting = weighting,
                 excluded = n_cells - nrow(d), aliased = aliased, status = status,
                 dropped = dropped, removed_levels = removed,
                 note = if (length(note)) paste(c(note, "the model describes the retained cells only, not the full design"),
                                                collapse = "; ") else NA_character_),
            class = "adempr_factor_influence")
}

#' @export
print.adempr_factor_influence <- function(x, ...) {
  if (identical(x$status, "unavailable")) {
    cat(sprintf("Descriptive meta-model of %s for %s (%s): unavailable. %s.\n",
                x$measure, x$method, x$estimand, x$note))
    return(invisible(x))
  }
  if (identical(x$status, "reduced")) cat(sprintf("Reduced analysis: %s.\n", x$note))
  cat(sprintf("Descriptive meta-model of %s for %s (%s), %d cells, weighting %s; R^2 = %.3f.\n",
              x$measure, x$method, x$estimand, nrow(x$cells),
              x$weighting %||% (if (x$weighted) "inverse squared MCSE" else "equal"), x$r_squared))
  if (isTRUE(x$excluded > 0)) cat(sprintf("%d cell(s) without a finite value excluded.\n", x$excluded))
  if (isTRUE(x$aliased > 0)) cat(sprintf("%d coefficient(s) aliased (not estimable from these cells).\n", x$aliased))
  cat("Sums of squares of each term given the terms that do not contain it (type II); shares\n",
      "do not depend on factor order but need not add to R^2 in unbalanced designs.\n", sep = "")
  cat("It summarizes the cells run; it does not predict other scenarios. Report the cells too.\n")
  print(x$anova, row.names = FALSE, digits = 3)
  invisible(x)
}
