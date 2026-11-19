# Random number streams --------------------------------------------------
#
# Design (documented in vignette("reproducibility")):
# * One base seed per run is converted to a L'Ecuyer-CMRG state with the
#   normal kind "Inversion" and the sample kind "Rejection".
# * A run owns an append-only registry that maps a scientific identity
#   (stage, scenario key, component) to a stream index k. Stream k is the
#   base state advanced by parallel::nextRNGStream() k times.
# * Replication r of a stream uses substream r, the stream state advanced by
#   parallel::nextRNGSubStream() r times. Substream 0 is never used.
# * A fresh registry allocates indices in sorted key order, so the display
#   order of scenarios and methods is irrelevant. Extending a registry only
#   appends; existing indices are never renumbered.

RNG_KIND <- c("L'Ecuyer-CMRG", "Inversion", "Rejection")

rng_save <- function() {
  list(
    kind = RNGkind(),
    seed = if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
  )
}

rng_restore <- function(state) {
  suppressWarnings(RNGkind(state$kind[1], state$kind[2], state$kind[3]))
  if (is.null(state$seed)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state$seed, envir = globalenv())
  }
  invisible(NULL)
}

rng_set_state <- function(seed) {
  suppressWarnings(RNGkind(RNG_KIND[1], RNG_KIND[2], RNG_KIND[3]))
  assign(".Random.seed", seed, envir = globalenv())
}

# Evaluate `expr` with the global RNG set to `seed`; the caller's RNG kind
# and state (including a missing .Random.seed) are restored afterwards.
with_rng_state <- function(seed, expr) {
  state <- rng_save()
  on.exit(rng_restore(state), add = TRUE)
  rng_set_state(seed)
  force(expr)
}

# Evaluate a function that must be deterministic (exact and numerical
# truths, generator checks). The caller's generator is restored, and a call
# that draws random numbers stops with an error instead of silently
# advancing the caller's stream or giving a different value at every call.
deterministic_call <- function(what, fun, ...,
                               hint = " Use a Monte Carlo truth (sim_truth_mc()) for a value that needs simulation.") {
  state <- rng_save()
  on.exit(rng_restore(state), add = TRUE)
  set.seed(1L, kind = RNG_KIND[1], normal.kind = RNG_KIND[2], sample.kind = RNG_KIND[3])
  before <- get(".Random.seed", envir = globalenv())
  out <- fun(...)
  after <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv())
  }
  if (!identical(before, after)) {
    stopf("%s draws random numbers, but it must be deterministic.%s", what, hint)
  }
  out
}

rng_base_seed <- function(seed, domain = NULL) {
  if (!is_count(seed)) {
    stopf("`seed` must be a single nonnegative whole number.")
  }
  state <- rng_save()
  on.exit(rng_restore(state), add = TRUE)
  set.seed(seed, kind = RNG_KIND[1], normal.kind = RNG_KIND[2],
           sample.kind = RNG_KIND[3])
  base <- get(".Random.seed", envir = globalenv())
  if (!is.null(domain)) base[2:7] <- domain_state(seed, domain)
  base
}

# Work outside the production run (pilots, DGP checks) uses a separate
# family of streams: the same seed then gives independent random numbers,
# because stream indices restart at 1 in every registry. The L'Ecuyer-CMRG
# state of the family is taken from an md5 hash of the seed and the domain
# name (192 bits, reduced to the valid ranges of the two components).
domain_state <- function(seed, domain) {
  h1 <- md5_string(paste0("adempR rng domain\r", domain, "\r", seed))
  h <- paste0(h1, md5_string(paste0(h1, "\r2")))
  words <- vapply(0:5, function(i) {
    hx <- substr(h, 8 * i + 1, 8 * i + 8)
    strtoi(substr(hx, 1, 4), 16L) * 65536 + strtoi(substr(hx, 5, 8), 16L)
  }, numeric(1))
  m <- c(rep(4294967087, 3), rep(4294944443, 3))
  x <- words %% (m - 1) + 1
  as.integer(ifelse(x >= 2^31, x - 2^32, x))
}

#' Random number stream registry
#'
#' `rng_registry()` creates or extends the append-only registry that maps
#' scientific identities to L'Ecuyer-CMRG stream indices. It is used
#' internally by [sim_run()] and stored in the run manifest; it is exported
#' so that the allocation can be inspected and tested.
#'
#' @param keys A data frame with character columns `stage`, `scenario_key`
#'   and `component`.
#' @param registry An existing registry to extend, or `NULL` for a fresh one.
#' @return A data frame with columns `index`, `stage`, `scenario_key` and
#'   `component`, ordered by `index`. Rows of `registry` are kept unchanged
#'   and new identities are appended in sorted order.
#' @examples
#' keys <- data.frame(stage = "dgp", scenario_key = c("b", "a"),
#'                    component = "data")
#' reg <- rng_registry(keys)
#' reg
#' # adding a key never renumbers existing ones
#' rng_registry(rbind(keys, data.frame(stage = "dgp", scenario_key = "0",
#'                                     component = "data")), reg)
#' @export
rng_registry <- function(keys, registry = NULL) {
  need <- c("stage", "scenario_key", "component")
  if (!is.data.frame(keys) || !all(need %in% names(keys))) {
    stopf("`keys` must be a data frame with columns %s.", paste(need, collapse = ", "))
  }
  keys <- unique(keys[need])
  keys[] <- lapply(keys, as.character)
  if (anyNA(keys) || any(!nzchar(as.matrix(keys)))) {
    stopf("RNG stream identities must be complete (no missing or empty fields).")
  }
  if (is.null(registry)) {
    registry <- data.frame(index = integer(), stage = character(),
                           scenario_key = character(), component = character(),
                           stringsAsFactors = FALSE)
  }
  check_registry(registry)
  id_new <- paste(keys$stage, keys$scenario_key, keys$component, sep = "\r")
  id_old <- paste(registry$stage, registry$scenario_key, registry$component,
                  sep = "\r")
  add <- keys[!(id_new %in% id_old), , drop = FALSE]
  if (nrow(add)) {
    add <- add[order(add$stage, add$scenario_key, add$component, method = "radix"), ,
               drop = FALSE]
    start <- if (nrow(registry)) max(registry$index) else 0L
    add$index <- as.integer(start + seq_len(nrow(add)))
    registry <- rbind(registry, add[c("index", need)])
  }
  rownames(registry) <- NULL
  registry
}

# A registry maps each complete identity to its own positive integer index.
check_registry <- function(registry) {
  need <- c("index", "stage", "scenario_key", "component")
  if (!is.data.frame(registry) || !all(need %in% names(registry))) {
    stopf("An RNG registry must be a data frame with columns %s.", paste(need, collapse = ", "))
  }
  idx <- registry$index
  if (!is.numeric(idx) || anyNA(idx) || any(idx < 1 | idx != round(idx))) {
    stopf("RNG registry indices must be positive whole numbers.")
  }
  if (anyDuplicated(idx)) stopf("The RNG registry assigns one index to several identities.")
  ids <- registry[c("stage", "scenario_key", "component")]
  if (anyNA(ids) || any(!nzchar(as.matrix(ids)))) {
    stopf("The RNG registry has an incomplete identity.")
  }
  if (anyDuplicated(paste(ids$stage, ids$scenario_key, ids$component, sep = "\r"))) {
    stopf("The RNG registry contains duplicate identities.")
  }
  invisible(registry)
}

rng_lookup <- function(registry, stage, scenario_key, component) {
  check_registry(registry)
  hit <- which(registry$stage == stage & registry$scenario_key == scenario_key &
                 registry$component == component)
  if (length(hit) != 1L) {
    stopf("No RNG stream is registered for stage `%s`, scenario `%s`, component `%s`.",
          stage, scenario_key, component)
  }
  registry$index[hit]
}

# Stream states for indices 1..max(k), computed by iteration and cached in
# the supplied environment.
rng_stream_seeds <- function(base, k, cache = NULL) {
  if (is.null(cache)) cache <- new.env(parent = emptyenv())
  if (is.null(cache$streams)) cache$streams <- list()
  have <- length(cache$streams)
  need <- max(k)
  if (need > have) {
    cur <- if (have == 0L) base else cache$streams[[have]]
    for (i in seq.int(have + 1L, need)) {
      cur <- parallel::nextRNGStream(cur)
      cache$streams[[i]] <- cur
    }
  }
  cache$streams[k]
}

# Substream states for replications `reps` of a stream.
rng_substreams <- function(stream_seed, reps) {
  reps <- as.integer(reps)
  if (!length(reps)) return(list())
  if (any(reps < 1L)) stopf("Replication numbers must be positive.")
  out <- vector("list", length(reps))
  ord <- order(reps)
  cur <- stream_seed
  pos <- 0L
  for (j in ord) {
    while (pos < reps[j]) {
      cur <- parallel::nextRNGSubStream(cur)
      pos <- pos + 1L
    }
    out[[j]] <- cur
  }
  out
}

#' Look up the random number state of a task
#'
#' Returns the L'Ecuyer-CMRG state assigned to replication `replication` of
#' a registered stream. Useful for audits and for replaying a single task
#' outside the package.
#'
#' @param seed The run's base seed.
#' @param registry A registry from [rng_registry()] (for example
#'   `run$manifest$rng_registry`).
#' @param stage,scenario_key,component The stream identity.
#' @param replication Replication number (substream), a single positive
#'   whole number.
#' @param domain `NULL` for the production run, or the RNG domain recorded
#'   in the manifest (`manifest$rng_domain`, for example `"pilot"`).
#' @return An integer vector usable as `.Random.seed`.
#' @examples
#' reg <- rng_registry(data.frame(stage = "dgp", scenario_key = "s1",
#'                                component = "data"))
#' rng_task_seed(1, reg, "dgp", "s1", "data", replication = 3)
#' @export
rng_task_seed <- function(seed, registry, stage, scenario_key, component,
                          replication, domain = NULL) {
  if (!is_count(replication) || replication < 1) {
    stopf("`replication` must be a single positive whole number.")
  }
  k <- rng_lookup(registry, stage, scenario_key, component)
  base <- rng_base_seed(seed, domain)
  rng_substreams(rng_stream_seeds(base, k)[[1]], replication)[[1]]
}

# Streams for parallel or resampling work inside one method fit. They are
# derived from the method's own task state through a hash into a separate
# region of the L'Ecuyer-CMRG cycle, so they never coincide with a
# registered stream, and they are the same at every replay of the task.
inner_seed_fun <- function(task_seed) {
  force(task_seed)
  function(k) {
    if (!is_count(k, 1)) stopf("`k` must be a positive whole number.")
    base <- task_seed
    base[2:7] <- domain_state(paste(task_seed, collapse = ","), "inner")
    out <- vector("list", k)
    cur <- base
    for (i in seq_len(k)) {
      cur <- parallel::nextRNGStream(cur)
      out[[i]] <- cur
    }
    out
  }
}
