# Mergeable moment accumulators ---------------------------------------------
#
# Streaming means, variances and covariances that can be merged exactly
# across batches (Chan, Golub and LeVeque 1979 pairwise update), so
# summaries can be built without holding all contributions in memory.

#' Mergeable summary accumulators
#'
#' `sim_accumulator()` creates an accumulator for one or two numeric
#' streams; `sim_accumulate()` adds values; `sim_merge()` combines two
#' accumulators exactly; `sim_accumulated()` returns the count, means,
#' variances, covariance and the Monte Carlo standard error of the mean.
#' The updates use the numerically stable pairwise formulas of Chan, Golub
#' and LeVeque, so batch order does not matter beyond rounding.
#'
#' @param x,y Numeric vectors (for `sim_accumulate()`, `y` is optional and
#'   paired with `x`).
#' @param acc,a,b Accumulators.
#' @return An accumulator (class `adempr_accumulator`) or, for
#'   `sim_accumulated()`, a list of summaries.
#' @examples
#' x <- rnorm(1000)
#' a <- sim_accumulate(sim_accumulator(), x[1:400])
#' b <- sim_accumulate(sim_accumulator(), x[401:1000])
#' m <- sim_accumulated(sim_merge(a, b))
#' all.equal(c(m$mean_x, m$var_x), c(mean(x), var(x)))
#' @export
sim_accumulator <- function() {
  structure(list(n = 0, mx = 0, my = 0, m2x = 0, m2y = 0, cxy = 0, md = 0, m2d = 0,
                 paired = NA),
            class = "adempr_accumulator")
}

#' @rdname sim_accumulator
#' @export
sim_accumulate <- function(acc, x, y = NULL) {
  stopifnot(inherits(acc, "adempr_accumulator"))
  if (!is.numeric(x) || !all(is.finite(x))) stopf("`x` must be finite numbers.")
  if (!is.null(y) && (!is.numeric(y) || length(y) != length(x) || !all(is.finite(y)))) {
    stopf("`y` must be finite numbers paired with `x`.")
  }
  if (!length(x)) return(acc)
  b <- sim_accumulator()
  # a double count does not overflow at 2^31
  b$n <- as.double(length(x))
  b$mx <- mean(x)
  b$m2x <- sum((x - b$mx)^2)
  b$paired <- !is.null(y)
  if (!is.null(y)) {
    b$my <- mean(y)
    b$m2y <- sum((y - b$my)^2)
    b$cxy <- sum((x - b$mx) * (y - b$my))
    # the difference is accumulated directly: subtracting large variances
    # and covariances would lose a small difference
    d <- x - y
    b$md <- mean(d)
    b$m2d <- sum((d - b$md)^2)
  }
  sim_merge(acc, b)
}

#' @rdname sim_accumulator
#' @export
sim_merge <- function(a, b) {
  stopifnot(inherits(a, "adempr_accumulator"), inherits(b, "adempr_accumulator"))
  if (a$n == 0) return(b)
  if (b$n == 0) return(a)
  if (!identical(a$paired, b$paired)) stopf("Cannot merge paired and unpaired accumulators.")
  n <- a$n + b$n
  dx <- b$mx - a$mx
  dy <- b$my - a$my
  out <- a
  out$n <- n
  out$mx <- a$mx + dx * b$n / n
  out$my <- a$my + dy * b$n / n
  out$m2x <- a$m2x + b$m2x + dx^2 * a$n * b$n / n
  out$m2y <- a$m2y + b$m2y + dy^2 * a$n * b$n / n
  out$cxy <- a$cxy + b$cxy + dx * dy * a$n * b$n / n
  dd <- b$md - a$md
  out$md <- a$md + dd * b$n / n
  out$m2d <- a$m2d + b$m2d + dd^2 * a$n * b$n / n
  out
}

#' @rdname sim_accumulator
#' @export
sim_accumulated <- function(acc) {
  stopifnot(inherits(acc, "adempr_accumulator"))
  n <- acc$n
  vx <- if (n > 1) acc$m2x / (n - 1) else NA_real_
  out <- list(n = n, mean_x = if (n) acc$mx else NA_real_, var_x = vx,
              mcse_mean_x = if (n > 1) sqrt(vx / n) else NA_real_)
  if (isTRUE(acc$paired)) {
    out$mean_y <- acc$my
    out$var_y <- if (n > 1) acc$m2y / (n - 1) else NA_real_
    out$cov_xy <- if (n > 1) acc$cxy / (n - 1) else NA_real_
    out$mean_diff <- acc$md
    out$mcse_mean_diff <- if (n > 1) sqrt(acc$m2d / (n - 1) / n) else NA_real_
  }
  out
}
