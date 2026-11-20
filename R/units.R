# Independent units -------------------------------------------------------
#
# With an antithetic generator, replications 2k - 1 and 2k are dependent
# and form one independent unit (pair k); otherwise every replication is
# its own unit. Every summary path (sim_summarise(), streaming, grid
# averages, decompositions, control variates, weighted averages) takes its
# MCSE from the units through the helpers here, so the analysis route does
# not change the uncertainty.

# Contributions of the measures whose value is the mean of one
# contribution per eligible replication, with the replications they come
# from; NULL for other measures. The eligibility masks are those of
# measure_cell().
linear_contributions <- function(m, rows, study, estimand) {
  if (!is.character(m) || length(m) != 1L) return(NULL)
  truth <- rows$truth
  usable <- rows$valid
  has_truth <- !is.na(truth)
  int_ok <- rows$estimate_status %in% USABLE_FOR_INTERVAL & rows$interval_status == "valid"
  test_ok <- rows$estimate_status %in% USABLE_FOR_INTERVAL & rows$test_status == "valid"
  e <- rows$estimate - truth
  out <- function(x, sel) list(x = as.numeric(x), replication = rows$replication[sel])
  switch(m,
    bias = out(e[usable & has_truth], usable & has_truth),
    mae = out(abs(e[usable & has_truth]), usable & has_truth),
    mse = out(e[usable & has_truth]^2, usable & has_truth),
    mean_se = {
      sel <- usable & rows$se_status == "valid"
      out(rows$se[sel], sel)
    },
    width = out((rows$upper - rows$lower)[int_ok], int_ok),
    coverage = {
      sel <- int_ok & has_truth
      out(rows$lower[sel] <= truth[sel] & truth[sel] <= rows$upper[sel], sel)
    },
    noncoverage_lower = {
      sel <- int_ok & has_truth
      out(truth[sel] < rows$lower[sel], sel)
    },
    noncoverage_upper = {
      sel <- int_ok & has_truth
      out(truth[sel] > rows$upper[sel], sel)
    },
    rejection = out(reject_indicator(rows[test_ok, , drop = FALSE], study$alpha), test_ok),
    usable = out(usable, rep(TRUE, nrow(rows))),
    usable_coverage = {
      if (any(int_ok & !has_truth)) return(NULL)
      out(int_ok & has_truth & rows$lower <= truth & truth <= rows$upper, rep(TRUE, nrow(rows)))
    },
    NULL)
}

# Measures whose contributions are events (0 or 1).
is_proportion_measure <- function(m) {
  is.character(m) && length(m) == 1L && m %in% c(PROPORTION_MEASURES, "usable_coverage")
}

# Sums and counts of contributions per unit.
unit_sums <- function(x, units) {
  s <- rowsum(as.numeric(x), units, reorder = FALSE)[, 1]
  nk <- rowsum(rep(1, length(units)), units, reorder = FALSE)[, 1]
  list(s = unname(s), n = unname(nk))
}

# Delete-one-unit estimates of a mean of contributions, in closed form
# from the unit sums and counts. The jackknife of these equals
# sd(unit means) / sqrt(units) when every unit has the same count.
unit_loo_means <- function(s, nk) (sum(s) - s) / (sum(nk) - nk)

# Replace the MCSE and interval of kernel result `r` by those of a
# delete-one-unit jackknife with leave-one-out values `th`. For an event
# proportion whose unit means are all equal (the plug-in MCSE is then 0,
# which is not certainty), the interval bounds the mean of independent
# unit means in [0, 1]: it lies between the share of units at 1 and the
# share above 0, so exact binomial limits for those two shares, each at
# half the error rate, hold jointly with probability at least `level`. A
# raw binomial interval over the dependent replications would not.
unit_jackknife <- function(r, th, level, proportion = FALSE, unit_means = NULL,
                           what = "antithetic pairs") {
  note <- function(txt) paste(c(if (!is.na(r$note)) r$note, txt), collapse = "; ")
  G <- length(th)
  r$n_units <- G
  if (G < 2L) {
    r$mcse <- r$lower <- r$upper <- NA_real_
    r$note <- note(sprintf("%s: fewer than two, no MCSE", what))
    return(r)
  }
  if (anyNA(th) || !all(is.finite(th))) {
    r$mcse <- r$lower <- r$upper <- NA_real_
    r$note <- note(sprintf("%s: the measure is undefined without some unit, no MCSE", what))
    return(r)
  }
  r$mcse <- sqrt((G - 1) / G * sum((th - mean(th))^2))
  z <- stats::qnorm(1 - (1 - level) / 2)
  r$lower <- r$value - z * r$mcse
  r$upper <- r$value + z * r$mcse
  txt <- sprintf("MCSE from %d %s (delete-one-unit jackknife)", G, what)
  if (proportion) {
    r$lower <- max(0, r$lower)
    r$upper <- min(1, r$upper)
    if (!is.null(unit_means) && length(unit_means) == G && all(unit_means == unit_means[1])) {
      r$lower <- cp_interval(sum(unit_means >= 1), G, level)[1, "lower"]
      r$upper <- cp_interval(sum(unit_means > 0), G, level)[1, "upper"]
      r$status <- "boundary"
      txt <- sprintf("all %d %s have the same mean, so the plug-in MCSE is 0 but the uncertainty is not; the interval bounds the mean of independent unit means in [0, 1] by exact limits for the shares of units at 1 and above 0",
                     G, what)
    }
  }
  r$note <- note(txt)
  r
}

# MCSE of a linear measure from its contributions and their units.
linear_unit_result <- function(r, x, units, level, proportion = FALSE, what = "antithetic pairs") {
  us <- unit_sums(x, units)
  unit_jackknife(r, unit_loo_means(us$s, us$n), level, proportion, us$s / us$n, what)
}
