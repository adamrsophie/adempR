# Scenario design ----------------------------------------------------------

# Scenario keys index the random number streams, so their version stays
# fixed when the encoding of other objects (functions, formulas,
# environments) is revised: the same scenario keeps the same streams.
SCENARIO_KEY_VERSION <- "2"

scenario_key_of <- function(values) {
  values <- values[order(names(values), method = "radix")]
  md5_string(c(SCENARIO_KEY_VERSION, canonical(values)))
}

scenario_label_of <- function(values) {
  parts <- vapply(names(values), function(n) {
    v <- values[[n]]
    txt <- if (is.atomic(v) && length(v) == 1L) {
      if (is.double(v)) format(v, digits = 6) else as.character(v)
    } else {
      paste0("<", class(v)[1], ":", short_hash(md5_string(canonical(v))), ">")
    }
    paste0(n, "=", txt)
  }, character(1))
  paste(parts, collapse = ", ")
}

new_scenarios <- function(values, design, excluded = NULL, constraints = list(),
                          distribution = NULL, weights = NULL) {
  if (!length(values)) stopf("The scenario design contains no scenarios.")
  # every constructor (factorial, table, sampled, combined) passes here, so
  # factor names are checked once for all: identity columns of the public
  # tables cannot be overwritten by a factor
  for (v in values) {
    nm <- names(v)
    if (is.null(nm) || any(!nzchar(nm))) stopf("All scenario factors must be named.")
    if (anyDuplicated(nm)) stopf("Duplicated scenario factor name(s): %s.", paste(unique(nm[duplicated(nm)]), collapse = ", "))
    bad <- intersect(nm, RESERVED_FACTORS)
    if (length(bad)) stopf("`%s` is a reserved name and cannot be a scenario factor.", bad[1])
    if (!is.null(weights) && "weight" %in% nm) {
      stopf("`weight` is reserved for the design weights and cannot be a scenario factor of a weighted design.")
    }
  }
  keys <- vapply(values, scenario_key_of, character(1))
  dup <- duplicated(keys)
  if (any(dup)) {
    stopf("Duplicate scenario configurations: %s.",
          paste(unique(vapply(values[dup], scenario_label_of, character(1))),
                collapse = "; "))
  }
  labels <- vapply(values, scenario_label_of, character(1))
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != length(values) ||
        any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
      stopf("`weights` must be nonnegative finite numbers, one per scenario, with a positive sum.")
    }
    weights <- weights / sum(weights)
  }
  structure(list(
    values = values,
    keys = unname(keys),
    labels = unname(labels),
    design = design,
    excluded = excluded %||% data.frame(scenario = character(), reason = character(),
                                        stringsAsFactors = FALSE),
    constraints = constraints,
    distribution = distribution,
    weights = weights
  ), class = "adempr_scenarios")
}

#' Constraints on scenario combinations
#'
#' A constraint is a predicate evaluated on each candidate scenario (a named
#' list of factor values). `type = "exclude"` removes scenarios where the
#' predicate is `TRUE` and records the reason; such scenarios are a
#' deliberate design choice and are listed in reports. `type = "invalid"`
#' declares configurations that are logically impossible; any candidate
#' scenario matching it is an error, because the grid itself is wrong.
#'
#' @param when A function of one argument (the scenario list) returning a
#'   single `TRUE` or `FALSE`.
#' @param reason A short explanation, kept in the protocol and the report.
#' @param type `"exclude"` or `"invalid"`.
#' @return An object of class `adempr_constraint`.
#' @examples
#' sim_constraint(function(s) s$n < 50 && s$p > 10,
#'                reason = "more predictors than events per variable allow")
#' @export
sim_constraint <- function(when, reason, type = c("exclude", "invalid")) {
  check_function(when, "when")
  if (!is_scalar_character(reason)) stopf("`reason` must be a nonempty string.")
  structure(list(when = when, reason = reason, type = match.arg(type)),
            class = "adempr_constraint")
}

apply_constraints <- function(values, constraints) {
  excluded <- list()
  keep <- rep(TRUE, length(values))
  for (con in constraints) {
    if (!inherits(con, "adempr_constraint")) {
      stopf("Constraints must be created with sim_constraint().")
    }
    for (i in seq_along(values)) {
      hit <- con$when(values[[i]])
      if (!is.logical(hit) || length(hit) != 1L || is.na(hit)) {
        stopf("Constraint \"%s\" returned %s for scenario %s; it must return TRUE or FALSE.",
              con$reason, paste(deparse(hit), collapse = ""), scenario_label_of(values[[i]]))
      }
      if (hit) {
        if (con$type == "invalid") {
          stopf("Scenario %s is an invalid configuration: %s.",
                scenario_label_of(values[[i]]), con$reason)
        }
        if (keep[i]) {
          excluded[[length(excluded) + 1L]] <- data.frame(
            scenario = scenario_label_of(values[[i]]), reason = con$reason,
            stringsAsFactors = FALSE)
        }
        keep[i] <- FALSE
      }
    }
  }
  list(values = values[keep],
       excluded = bind_rows(excluded) %||%
         data.frame(scenario = character(), reason = character(),
                    stringsAsFactors = FALSE))
}

RESERVED_FACTORS <- c("scenario_key", "scenario_id", "replication")

check_factor_values <- function(name, v) {
  if (!nzchar(name)) stopf("All scenario factors must be named.")
  if (name %in% RESERVED_FACTORS) {
    stopf("`%s` is a reserved name and cannot be a scenario factor.", name)
  }
  if (is.function(v) || is.environment(v)) {
    stopf("Scenario factor `%s` must be a vector or a list, not a %s.", name, typeof(v))
  }
  if (!length(v)) stopf("Scenario factor `%s` has no levels.", name)
  if (is.atomic(v) && anyNA(v)) {
    stopf("Scenario factor `%s` contains missing values.", name)
  }
  if (is.factor(v)) v <- as.character(v)
  v
}

#' Define scenarios
#'
#' `sim_scenarios()` builds a full factorial design from named factors,
#' optionally removing combinations with [sim_constraint()].
#' `sim_scenarios_table()` uses an explicit table of scenarios.
#' `sim_scenarios_sample()` draws scenarios reproducibly from a declared
#' distribution. `sim_scenarios_bind()` combines designs, which expresses
#' nested factors (a different grid for each level of a factor).
#'
#' Every scenario gets a key computed from its content (the md5 hash of its
#' canonical factor values), so scientific identity does not depend on row
#' order. The key is used for random number streams, storage and
#' resumption.
#'
#' @param ... For `sim_scenarios()`, named factors: each is an atomic
#'   vector of levels or a list of levels (for matrix or other non scalar
#'   settings). For `sim_scenarios_bind()`, scenario designs to combine.
#' @param .constraints A list of [sim_constraint()] objects.
#' @param .weights Optional nonnegative weights, one per retained scenario,
#'   for prespecified global summaries.
#' @param .max_scenarios Refuse to build grids larger than this.
#' @return An object of class `adempr_scenarios`.
#' @examples
#' sc <- sim_scenarios(n = c(30L, 100L), mu = c(0, 0.5), sigma = 1)
#' sc
#' as.data.frame(sc)
#'
#' sim_scenarios(n = c(20L, 200L), p = c(2L, 10L),
#'   .constraints = list(sim_constraint(function(s) s$n < s$p * 5,
#'                                      reason = "too few observations per predictor")))
#' @export
sim_scenarios <- function(..., .constraints = list(), .weights = NULL,
                          .max_scenarios = 1e5) {
  factors <- list(...)
  if (!length(factors)) stopf("Supply at least one named scenario factor.")
  nm <- names(factors) %||% rep("", length(factors))
  if (anyDuplicated(nm)) stopf("Scenario factor names must be unique.")
  factors <- Map(check_factor_values, nm, factors)
  sizes <- vapply(factors, length, integer(1))
  total <- prod(sizes)
  if (total > .max_scenarios) {
    stopf("The factorial grid has %s scenarios, more than `.max_scenarios` = %s. Use constraints, a smaller design or raise the limit deliberately.",
          format(total, big.mark = ","), format(.max_scenarios, big.mark = ","))
  }
  idx <- expand.grid(lapply(sizes, seq_len), KEEP.OUT.ATTRS = FALSE)
  values <- lapply(seq_len(nrow(idx)), function(i) {
    out <- lapply(seq_along(factors), function(j) {
      v <- factors[[j]]
      if (is.list(v)) v[[idx[i, j]]] else v[idx[i, j]]
    })
    names(out) <- names(factors)
    out
  })
  res <- apply_constraints(values, .constraints)
  new_scenarios(res$values, design = "factorial", excluded = res$excluded,
                constraints = .constraints, weights = .weights)
}

#' @rdname sim_scenarios
#' @param table A data frame with one row per scenario. List columns are
#'   allowed.
#' @export
sim_scenarios_table <- function(table, .constraints = list(), .weights = NULL) {
  if (!is.data.frame(table) || !nrow(table)) {
    stopf("`table` must be a data frame with at least one row.")
  }
  for (n in names(table)) check_factor_values(n, table[[n]])
  values <- lapply(seq_len(nrow(table)), function(i) {
    out <- lapply(names(table), function(n) {
      v <- table[[n]]
      if (is.factor(v)) v <- as.character(v)
      if (is.list(v)) v[[i]] else v[i]
    })
    names(out) <- names(table)
    out
  })
  res <- apply_constraints(values, .constraints)
  new_scenarios(res$values, design = "table", excluded = res$excluded,
                constraints = .constraints, weights = .weights)
}

#' @rdname sim_scenarios
#' @param sampler A function of no arguments that returns one scenario as a
#'   named list, using the random number generator.
#' @param n Number of scenarios to draw.
#' @param seed Seed for drawing the scenarios. The draw uses its own
#'   L'Ecuyer-CMRG stream and leaves the caller's random number state
#'   unchanged.
#' @param distribution A text description of the scenario distribution,
#'   stored in the protocol. Performance averaged over sampled scenarios
#'   estimates performance averaged over this distribution, which differs
#'   from performance conditional on any single scenario.
#' @export
sim_scenarios_sample <- function(sampler, n, seed, distribution,
                                 .constraints = list()) {
  check_function(sampler, "sampler")
  if (!is_count(n, 1)) stopf("`n` must be a positive whole number.")
  if (!is_scalar_character(distribution)) {
    stopf("`distribution` must describe the scenario distribution.")
  }
  base <- rng_base_seed(seed)
  values <- with_rng_state(base, lapply(seq_len(n), function(i) {
    s <- sampler()
    if (!is.list(s) || is.null(names(s)) || any(!nzchar(names(s)))) {
      stopf("`sampler` must return a named list.")
    }
    s
  }))
  res <- apply_constraints(values, .constraints)
  new_scenarios(res$values, design = "sampled", excluded = res$excluded,
                constraints = .constraints,
                distribution = list(description = distribution, seed = seed, n = n,
                                    sampler_fingerprint = fingerprint(sampler)))
}

#' @rdname sim_scenarios
#' @export
sim_scenarios_bind <- function(...) {
  designs <- list(...)
  if (!all(vapply(designs, inherits, logical(1), "adempr_scenarios"))) {
    stopf("All arguments must be scenario designs.")
  }
  values <- do.call(c, lapply(designs, `[[`, "values"))
  excluded <- bind_rows(lapply(designs, `[[`, "excluded"))
  new_scenarios(values, design = "combined", excluded = excluded,
                constraints = do.call(c, lapply(designs, `[[`, "constraints")))
}

#' @export
length.adempr_scenarios <- function(x) length(x$values)

#' @export
as.data.frame.adempr_scenarios <- function(x, ...) {
  nms <- unique(unlist(lapply(x$values, names)))
  cols <- lapply(nms, function(n) {
    vals <- lapply(x$values, function(v) if (n %in% names(v)) v[[n]] else NA)
    scalar <- all(vapply(vals, function(v) is.atomic(v) && length(v) == 1L, logical(1)))
    if (scalar) unlist(vals) else I(vals)
  })
  names(cols) <- nms
  df <- data.frame(scenario_key = x$keys, scenario_id = x$labels,
                   stringsAsFactors = FALSE)
  for (n in nms) df[[n]] <- cols[[n]]
  if (!is.null(x$weights)) df$weight <- x$weights
  df
}

#' @export
print.adempr_scenarios <- function(x, ...) {
  cat(sprintf("<adempR scenarios: %d %s scenario%s>\n", length(x), x$design,
              if (length(x) == 1L) "" else "s"))
  n_show <- min(length(x), 10L)
  for (i in seq_len(n_show)) cat("  ", x$labels[i], "\n", sep = "")
  if (length(x) > n_show) cat(sprintf("  ... and %d more\n", length(x) - n_show))
  if (nrow(x$excluded)) {
    cat(sprintf("Excluded by design: %d\n", nrow(x$excluded)))
  }
  if (!is.null(x$distribution)) {
    cat("Sampled from: ", x$distribution$description, "\n", sep = "")
  }
  invisible(x)
}

scenario_index <- function(scenarios, key) {
  i <- match(key, scenarios$keys)
  if (anyNA(i)) stopf("Unknown scenario key: %s.", paste(key[is.na(i)], collapse = ", "))
  i
}

# Fingerprint of values that index random number streams (common random
# number groups, shared truth factors); versioned like scenario keys.
stream_fingerprint <- function(...) {
  md5_string(c(SCENARIO_KEY_VERSION, canonical(list(...))))
}
