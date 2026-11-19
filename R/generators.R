# Data-generating building blocks -----------------------------------------
#
# A small validated collection. Each function checks its inputs and never
# repairs them: an impossible probability, a non positive definite
# correlation matrix or an inconsistent dimension is an error. Use these
# inside generator functions; they draw from the random number stream that
# the execution layer assigns.

check_prob <- function(p, what, open = FALSE) {
  if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1) || (open && any(p <= 0 | p >= 1))) {
    stopf("`%s` must contain probabilities in %s.", what, if (open) "(0, 1)" else "[0, 1]")
  }
  invisible(p)
}

#' Validate a correlation or covariance matrix
#'
#' Checks that a matrix is square, symmetric, has a positive diagonal (unit
#' diagonal for a correlation matrix) and is positive semidefinite
#' (eigenvalues of the matrix scaled to unit diagonal above `-tol` times
#' the largest). The matrix is never
#' modified: a failing matrix is an error, because silently repairing it
#' would change the data-generating mechanism.
#'
#' @param m A matrix.
#' @param type `"correlation"` or `"covariance"`.
#' @param tol Relative eigenvalue tolerance.
#' @return `m`, invisibly, if valid.
#' @examples
#' sim_check_matrix(matrix(c(1, 0.5, 0.5, 1), 2))
#' try(sim_check_matrix(matrix(c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3)))
#' @export
sim_check_matrix <- function(m, type = c("correlation", "covariance"), tol = 1e-10) {
  type <- match.arg(type)
  if (!is.matrix(m) || !is.numeric(m) || nrow(m) != ncol(m) || anyNA(m) || any(!is.finite(m))) {
    stopf("The %s matrix must be a finite numeric square matrix.", type)
  }
  if (max(abs(m - t(m))) > 1e-12 * max(1, max(abs(m)))) stopf("The %s matrix is not symmetric.", type)
  d <- diag(m)
  if (type == "correlation") {
    if (any(abs(d - 1) > 1e-12)) stopf("A correlation matrix must have a unit diagonal.")
    if (any(abs(m) > 1 + 1e-12)) stopf("Correlations must lie in [-1, 1].")
  } else if (any(d <= 0)) {
    stopf("A covariance matrix must have a positive diagonal.")
  }
  # Checked on the correlation scale, so that a large variance elsewhere in
  # the matrix cannot hide an indefinite block within the tolerance.
  s <- 1 / sqrt(d)
  ev <- eigen(m * outer(s, s), symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tol * max(abs(ev))) {
    stopf("The %s matrix is not positive semidefinite (smallest eigenvalue %s). It is not repaired automatically.",
          type, format(min(ev), digits = 4))
  }
  invisible(m)
}

#' Multivariate normal draws
#'
#' @param n Number of draws.
#' @param mean Mean vector.
#' @param sigma Covariance matrix (validated with [sim_check_matrix()]).
#' @return An `n` by `length(mean)` matrix.
#' @examples
#' set.seed(1)
#' x <- sim_rmvnorm(1000, c(0, 0), matrix(c(1, 0.6, 0.6, 1), 2))
#' cor(x)
#' @export
sim_rmvnorm <- function(n, mean, sigma) {
  if (!is_count(n)) stopf("`n` must be a whole number.")
  sim_check_matrix(sigma, "covariance")
  if (length(mean) != nrow(sigma)) stopf("`mean` and `sigma` have incompatible dimensions.")
  ev <- eigen(sigma, symmetric = TRUE)
  root <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), length(mean)) %*% t(ev$vectors)
  z <- matrix(stats::rnorm(n * length(mean)), n)
  sweep(z %*% root, 2, mean, `+`)
}

#' Calibrate an intercept to a target marginal mean
#'
#' Finds the intercept `a` such that the marginal mean of
#' `inverse_link(a + lp)` equals `target`, where `lp` is a large sample of
#' the linear predictor drawn by `lp_sampler`. This calibrates, for
#' example, a logistic model to a required marginal prevalence, or a
#' missingness model to a required proportion missing. The solver is
#' [stats::uniroot()] on the sample average; the achieved value is reported
#' with its Monte Carlo standard error so the calibration can be checked.
#' A conditional coefficient is not a marginal effect, and a prevalence on
#' the latent scale is not the observed prevalence; this function targets
#' the observed scale marginal mean.
#'
#' @param target Target marginal mean in (0, 1).
#' @param lp_sampler Function of `n` returning `n` draws of the linear
#'   predictor without intercept.
#' @param n Size of the calibration sample.
#' @param seed Seed for the calibration sample (the caller's random number
#'   state is left unchanged).
#' @param inverse_link Inverse link function.
#' @param interval Search interval for the intercept.
#' @return A list with `intercept`, `achieved`, `achieved_se` and `n`.
#' @examples
#' cal <- sim_calibrate_intercept(0.2, function(n) 0.8 * rnorm(n), n = 1e5, seed = 1)
#' cal
#' @export
sim_calibrate_intercept <- function(target, lp_sampler, n = 1e5, seed = 1L,
                                    inverse_link = stats::plogis, interval = c(-30, 30)) {
  check_prob(target, "target", open = TRUE)
  check_function(lp_sampler, "lp_sampler")
  lp <- with_rng_state(rng_base_seed(seed), lp_sampler(n))
  if (!is.numeric(lp) || length(lp) != n || any(!is.finite(lp))) {
    stopf("`lp_sampler` must return `n` finite numbers.")
  }
  f <- function(a) mean(inverse_link(a + lp)) - target
  if (f(interval[1]) * f(interval[2]) > 0) {
    stopf("The target %s cannot be reached for intercepts in [%s, %s].", target,
          interval[1], interval[2])
  }
  a <- stats::uniroot(f, interval, tol = 1e-10)$root
  p <- inverse_link(a + lp)
  list(intercept = a, achieved = mean(p), achieved_se = stats::sd(p) / sqrt(n),
       target = target, n = n)
}

#' Survival times from a Weibull proportional hazards model
#'
#' Draws event times with hazard `lambda * shape * t^(shape - 1) *
#' exp(lp)` by inversion, then applies independent exponential censoring
#' with rate `censoring_rate` (0 for none) and administrative censoring at
#' `horizon`.
#'
#' @param lp Linear predictor, one value per subject.
#' @param lambda,shape Weibull scale and shape (positive).
#' @param censoring_rate Rate of independent exponential censoring.
#' @param horizon Administrative censoring time.
#' @return A data frame with `time`, `status` (1 event, 0 censored) and
#'   `event_time` (the uncensored time, protected information that methods
#'   should not see; pass it through [sim_data()] metadata if needed).
#' @examples
#' set.seed(1)
#' d <- sim_rweibull_ph(rnorm(200), lambda = 0.1, shape = 1.5, horizon = 5)
#' mean(d$status)
#' @export
sim_rweibull_ph <- function(lp, lambda, shape, censoring_rate = 0, horizon = Inf) {
  if (!is.numeric(lp) || anyNA(lp)) stopf("`lp` must be numeric without missing values.")
  if (!is_scalar_number(lambda) || lambda <= 0 || !is_scalar_number(shape) || shape <= 0) {
    stopf("`lambda` and `shape` must be positive numbers.")
  }
  if (!is_scalar_number(censoring_rate) || censoring_rate < 0) {
    stopf("`censoring_rate` must be nonnegative.")
  }
  if (!is_scalar_number(horizon) || horizon <= 0) stopf("`horizon` must be positive.")
  n <- length(lp)
  t_event <- (-log(stats::runif(n)) / (lambda * exp(lp)))^(1 / shape)
  t_cens <- if (censoring_rate > 0) stats::rexp(n, censoring_rate) else rep(Inf, n)
  time <- pmin(t_event, t_cens, horizon)
  data.frame(time = time, status = as.integer(t_event <= pmin(t_cens, horizon)),
             event_time = t_event)
}

#' Impose missing values
#'
#' Sets values of `var` to `NA` completely at random (`"MCAR"`, with
#' probability `prop`) or at random given other fully observed variables
#' (`"MAR"`, logistic in `weights` applied to `predictors`, with the
#' intercept calibrated per data set so the expected proportion missing in
#' that data set equals `prop`). Missing not at random mechanisms can be
#' expressed by including `var` itself in `predictors`; the label is then
#' `"MNAR"`.
#'
#' @param data A data frame.
#' @param var Variable to make incomplete.
#' @param prop Expected proportion missing.
#' @param mechanism `"MCAR"`, `"MAR"` or `"MNAR"`.
#' @param predictors Variables driving missingness (MAR, MNAR).
#' @param weights Coefficients for `predictors` (on standardized
#'   predictors).
#' @return `data` with missing values and an attribute `missing` holding the
#'   missingness indicator.
#' @examples
#' set.seed(2)
#' d <- data.frame(x = rnorm(500), y = rnorm(500))
#' d2 <- sim_add_missing(d, "y", prop = 0.3, mechanism = "MAR", predictors = "x", weights = 1)
#' mean(is.na(d2$y))
#' @export
sim_add_missing <- function(data, var, prop, mechanism = c("MCAR", "MAR", "MNAR"),
                            predictors = character(), weights = numeric()) {
  mechanism <- match.arg(mechanism)
  if (!is.data.frame(data) || !(var %in% names(data))) stopf("`var` must be a column of `data`.")
  check_prob(prop, "prop")
  n <- nrow(data)
  if (mechanism == "MCAR") {
    miss <- stats::runif(n) < prop
  } else {
    if (!length(predictors) || length(predictors) != length(weights)) {
      stopf("MAR and MNAR need `predictors` and matching `weights`.")
    }
    if (mechanism == "MAR" && var %in% predictors) {
      stopf("A mechanism that depends on `%s` itself is MNAR, not MAR.", var)
    }
    if (mechanism == "MNAR" && !(var %in% predictors)) {
      stopf("An MNAR mechanism here must include `%s` among the predictors.", var)
    }
    X <- scale(as.matrix(data[predictors]))
    if (anyNA(X)) stopf("Missingness predictors must be fully observed and nonconstant.")
    lp <- drop(X %*% weights)
    if (prop <= 0 || prop >= 1) {
      miss <- rep(prop >= 1, n)
    } else {
      a <- stats::uniroot(function(a) mean(stats::plogis(a + lp)) - prop, c(-50, 50),
                          tol = 1e-10)$root
      miss <- stats::runif(n) < stats::plogis(a + lp)
    }
  }
  data[[var]][miss] <- NA
  attr(data, "missing") <- miss
  data
}

#' Clustered normal outcomes with a random intercept
#'
#' Generates `n_clusters` clusters of size `cluster_size` (a single size or
#' one per cluster) with outcome `y = mu + u_j + e_ij`, where the random
#' intercept variance and residual variance are set from the intraclass
#' correlation `icc` and the total variance.
#'
#' @param n_clusters Number of clusters.
#' @param cluster_size Cluster sizes.
#' @param mu Mean.
#' @param icc Intraclass correlation in [0, 1).
#' @param total_var Total variance.
#' @return A data frame with `cluster`, `y` and `u` (the latent random
#'   effect, which methods should not see).
#' @examples
#' set.seed(3)
#' d <- sim_rclustered(20, 10, mu = 1, icc = 0.2)
#' head(d)
#' @export
sim_rclustered <- function(n_clusters, cluster_size, mu = 0, icc = 0, total_var = 1) {
  if (!is_count(n_clusters, 1)) stopf("`n_clusters` must be a positive whole number.")
  if (!(length(cluster_size) %in% c(1L, n_clusters)) || any(cluster_size < 1) ||
      any(cluster_size != round(cluster_size))) {
    stopf("`cluster_size` must be positive whole numbers, one or one per cluster.")
  }
  if (!is_scalar_number(icc) || icc < 0 || icc >= 1) stopf("`icc` must be in [0, 1).")
  if (!is_scalar_number(total_var) || total_var <= 0) stopf("`total_var` must be positive.")
  sizes <- rep_len(cluster_size, n_clusters)
  u <- stats::rnorm(n_clusters, 0, sqrt(icc * total_var))
  cl <- rep(seq_len(n_clusters), sizes)
  data.frame(cluster = cl, y = mu + u[cl] + stats::rnorm(sum(sizes), 0, sqrt((1 - icc) * total_var)),
             u = u[cl])
}

#' Contaminate values with outliers
#'
#' Replaces a random proportion `prop` of `x` by draws from a normal
#' distribution with standard deviation `scale` times that of the clean
#' distribution (scale contamination), to study robustness.
#'
#' @param x Numeric vector.
#' @param prop Contamination proportion.
#' @param scale Scale inflation of contaminated values.
#' @param sd Standard deviation of the clean values.
#' @return `x` with contaminated entries; attribute `contaminated` marks
#'   them.
#' @examples
#' set.seed(4)
#' x <- sim_contaminate(rnorm(100), prop = 0.1, scale = 10)
#' sum(attr(x, "contaminated"))
#' @export
sim_contaminate <- function(x, prop, scale = 10, sd = 1) {
  check_prob(prop, "prop")
  if (!is_scalar_number(scale) || scale <= 0) stopf("`scale` must be positive.")
  hit <- stats::runif(length(x)) < prop
  x[hit] <- stats::rnorm(sum(hit), 0, scale * sd)
  attr(x, "contaminated") <- hit
  x
}

#' Calibrate an exponential censoring rate to a target censored proportion
#'
#' With event times `T` drawn by `event_sampler`, independent exponential
#' censoring with rate `r` and administrative censoring at `horizon`, the
#' expected censored proportion is `1 - E[exp(-r T) 1(T <= horizon)]`.
#' The rate is found with [stats::uniroot()] on a large calibration sample
#' drawn from its own seed, and the achieved proportion is reported with
#' its Monte Carlo standard error. Pass the rate to [sim_rweibull_ph()] as
#' `censoring_rate`.
#'
#' @param target Target censored proportion in (0, 1).
#' @param event_sampler Function of `n` returning `n` event times.
#' @param horizon Administrative censoring time.
#' @param n Size of the calibration sample.
#' @param seed Seed for the calibration sample (the caller's random number
#'   state is left unchanged).
#' @param upper Upper end of the search interval for the rate.
#' @return A list with `rate`, `achieved`, `achieved_se`, `minimum` (the
#'   censored proportion without random censoring) and `n`.
#' @examples
#' cal <- sim_calibrate_censoring(0.3, function(n) {
#'   sim_rweibull_ph(rnorm(n), lambda = 0.1, shape = 1.5)$event_time
#' }, horizon = 10, seed = 1)
#' cal$rate
#' @export
sim_calibrate_censoring <- function(target, event_sampler, horizon = Inf, n = 1e5, seed = 1L,
                                    upper = 1e3) {
  check_prob(target, "target", open = TRUE)
  check_function(event_sampler, "event_sampler")
  if (!is_scalar_number(horizon) || horizon <= 0) stopf("`horizon` must be positive.")
  t <- with_rng_state(rng_base_seed(seed), event_sampler(n))
  if (!is.numeric(t) || length(t) != n || anyNA(t) || any(t < 0)) {
    stopf("`event_sampler` must return `n` nonnegative event times.")
  }
  inside <- t <= horizon
  cens <- function(r) 1 - mean(exp(-r * t) * inside)
  minimum <- cens(0)
  if (target < minimum) {
    stopf("Administrative censoring alone censors %.3f, above the target %s.", minimum, target)
  }
  if (cens(upper) < target) stopf("The target %s cannot be reached with rates up to %s.", target, upper)
  r <- if (target == minimum) 0 else stats::uniroot(function(r) cens(r) - target, c(0, upper),
                                                    tol = 1e-12)$root
  ind <- 1 - exp(-r * t) * inside
  list(rate = r, achieved = mean(ind), achieved_se = stats::sd(ind) / sqrt(n),
       minimum = minimum, target = target, n = n)
}

#' Draws from a Gaussian copula with given margins
#'
#' Draws correlated normal variables with correlation matrix `corr`,
#' transforms them to uniforms and applies each margin's quantile
#' function. The correlation applies to the normal scores, not to the
#' variables themselves (Spearman correlations equal
#' `6 / pi * asin(corr / 2)` for continuous margins).
#'
#' @param n Number of draws.
#' @param margins Named list of quantile functions, one per variable, for
#'   example `list(age = function(p) qnorm(p, 50, 10), y = function(p)
#'   qpois(p, 3))`.
#' @param corr Correlation matrix of the normal scores (validated with
#'   [sim_check_matrix()]).
#' @return A data frame with one column per margin.
#' @examples
#' set.seed(1)
#' d <- sim_rcopula(1000, list(x = function(p) qexp(p), y = function(p) qbinom(p, 1, 0.3)),
#'                  matrix(c(1, 0.5, 0.5, 1), 2))
#' cor(d, method = "spearman")
#' @export
sim_rcopula <- function(n, margins, corr) {
  if (!is_count(n, 1)) stopf("`n` must be a positive whole number.")
  if (!is.list(margins) || !length(margins) || is.null(names(margins)) ||
      !all(vapply(margins, is.function, logical(1)))) {
    stopf("`margins` must be a named list of quantile functions.")
  }
  k <- length(margins)
  corr <- as.matrix(corr)
  if (!identical(dim(corr), c(k, k))) stopf("`corr` must be a %d by %d matrix.", k, k)
  sim_check_matrix(corr, "correlation")
  z <- sim_rmvnorm(n, rep(0, k), corr)
  u <- stats::pnorm(z)
  out <- lapply(seq_len(k), function(j) {
    v <- margins[[j]](u[, j])
    if (length(v) != n) stopf("Margin `%s` must return one value per probability.", names(margins)[j])
    v
  })
  names(out) <- names(margins)
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Plasmode data: resampled real covariates with simulated outcomes
#'
#' A plasmode simulation keeps the joint distribution of covariates from a
#' real data set by resampling its rows, and generates the outcome (or
#' treatment) from a specified model so the true effect is known. The rows
#' drawn are stored in the attribute `rows`, so a data set can be traced to
#' the source.
#'
#' @param data The source data frame (covariates).
#' @param n Number of rows to draw.
#' @param outcome Function of the resampled data returning a vector (added
#'   as column `y`) or a data frame or named list of new columns; `NULL`
#'   returns the resampled covariates only.
#' @param replace Sample rows with replacement (bootstrap); `FALSE` draws a
#'   subsample, which needs `n <= nrow(data)`.
#' @return A data frame with attribute `rows`.
#' @examples
#' set.seed(1)
#' pm <- sim_plasmode(mtcars, n = 50, outcome = function(d) {
#'   rbinom(nrow(d), 1, plogis(-1 + 0.5 * d$am + 0.02 * (d$hp - 150)))
#' })
#' head(pm)
#' @export
sim_plasmode <- function(data, n = nrow(data), outcome = NULL, replace = TRUE) {
  if (!is.data.frame(data) || !nrow(data)) stopf("`data` must be a nonempty data frame.")
  if (!is_count(n, 1)) stopf("`n` must be a positive whole number.")
  if (!replace && n > nrow(data)) stopf("Without replacement `n` cannot exceed %d rows.", nrow(data))
  rows <- sample.int(nrow(data), n, replace = replace)
  d <- data[rows, , drop = FALSE]
  rownames(d) <- NULL
  if (!is.null(outcome)) {
    check_function(outcome, "outcome")
    y <- outcome(d)
    if (is.list(y)) {
      if (is.null(names(y)) || any(!nzchar(names(y)))) stopf("`outcome` must return named columns.")
      for (nm in names(y)) d[[nm]] <- y[[nm]]
    } else {
      if (length(y) != n) stopf("`outcome` must return one value per row.")
      d$y <- y
    }
  }
  attr(d, "rows") <- rows
  d
}
