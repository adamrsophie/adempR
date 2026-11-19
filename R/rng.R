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

rng_base_seed <- function(seed) {
  if (!is_count(seed)) {
    stopf("`seed` must be a single nonnegative whole number.")
  }
  state <- rng_save()
  on.exit(rng_restore(state), add = TRUE)
  set.seed(seed, kind = RNG_KIND[1], normal.kind = RNG_KIND[2],
           sample.kind = RNG_KIND[3])
  get(".Random.seed", envir = globalenv())
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
  if (is.null(registry)) {
    registry <- data.frame(index = integer(), stage = character(),
                           scenario_key = character(), component = character(),
                           stringsAsFactors = FALSE)
  }
  id_new <- paste(keys$stage, keys$scenario_key, keys$component, sep = "\r")
  id_old <- paste(registry$stage, registry$scenario_key, registry$component,
                  sep = "\r")
  if (anyDuplicated(id_old)) stopf("The RNG registry contains duplicate identities.")
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

rng_lookup <- function(registry, stage, scenario_key, component) {
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
#' @param replication Replication number (substream).
#' @return An integer vector usable as `.Random.seed`.
#' @examples
#' reg <- rng_registry(data.frame(stage = "dgp", scenario_key = "s1",
#'                                component = "data"))
#' rng_task_seed(1, reg, "dgp", "s1", "data", replication = 3)
#' @export
rng_task_seed <- function(seed, registry, stage, scenario_key, component,
                          replication) {
  k <- rng_lookup(registry, stage, scenario_key, component)
  base <- rng_base_seed(seed)
  rng_substreams(rng_stream_seeds(base, k)[[1]], replication)[[1]]
}
