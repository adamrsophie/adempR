# Blinded evaluation ------------------------------------------------------

#' Export coded data sets for blinded analysis
#'
#' Regenerates the data sets of a run from their stored streams and writes
#' each to `dir/data/<code>.rds`, where the code is random and reveals
#' neither the scenario nor the replication. Only the data that methods
#' see are written (no truth, no metadata). The key linking codes to
#' scenarios, replications and data checksums is written to
#' `dir/key.rds`; keep it away from the analysts. Analysts return one row
#' per data set and estimand, which [sim_external_method()] turns into a
#' method of the study.
#'
#' @param run A run.
#' @param dir Output directory (must not exist or be empty).
#' @param scenarios Scenario labels or keys to export; defaults to all.
#' @param replications Replications to export; defaults to all planned.
#' @param seed Seed of the random codes; the caller's random state is left
#'   untouched.
#' @return Invisibly, the key: a data frame with `code`, `scenario_key`,
#'   `replication` and `data_md5`.
#' @export
sim_export_blinded <- function(run, dir, scenarios = NULL, replications = NULL, seed = 1L) {
  stopifnot(inherits(run, "adempr_run"))
  if (!is_scalar_character(dir)) stopf("`dir` must be a path.")
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE))) {
    stopf("%s is not empty.", dir)
  }
  man <- run$manifest
  study <- man$study
  keys <- if (is.null(scenarios)) study$scenarios$keys else
    vapply(scenarios, function(s) resolve_scenario(study, s), character(1))
  reps <- replications %||% seq_len(man$replications)
  if (!all(vapply(reps, is_count, logical(1), min = 1))) stopf("`replications` must be positive whole numbers.")
  tasks <- expand.grid(replication = as.integer(reps), scenario_key = keys, stringsAsFactors = FALSE)
  codes <- with_rng_state(rng_base_seed(seed, "blind"), {
    hex <- c(0:9, letters[1:6])
    unique_codes <- character()
    while (length(unique_codes) < nrow(tasks)) {
      unique_codes <- unique(c(unique_codes, vapply(seq_len(nrow(tasks) - length(unique_codes)),
        function(i) paste(sample(hex, 12, replace = TRUE), collapse = ""), character(1))))
    }
    unique_codes
  })
  dir.create(file.path(dir, "data"), recursive = TRUE, showWarnings = FALSE)
  ctx <- list(seed = man$seed, domain = man$rng_domain, registry = man$rng_registry,
              cache = new.env(parent = emptyenv()), dgp_stage = "dgp", method_stage = "method")
  md5 <- rep(NA_character_, nrow(tasks))
  with_packages(study$generator$packages, {
    for (i in seq_len(nrow(tasks))) {
      k <- tasks$scenario_key[i]
      r <- tasks$replication[i]
      seeds <- task_seeds(ctx, k, r, study, character())[[1]]
      out <- execute_task_seeded(study, list(scenario_key = k, replication = r, seeds = seeds,
                                             methods = character()), keep_data = TRUE)
      if (!identical(out$replication$generation_status, "ok")) next
      md5[i] <- out$replication$data_md5
      saveRDS(out$data, file.path(dir, "data", paste0(codes[i], ".rds")), version = 3)
    }
  })
  key <- data.frame(code = codes, scenario_key = tasks$scenario_key,
                    replication = tasks$replication, data_md5 = md5, stringsAsFactors = FALSE)
  key <- key[!is.na(key$data_md5), , drop = FALSE]
  rownames(key) <- NULL
  saveRDS(key, file.path(dir, "key.rds"), version = 3)
  writeLines(c("Each file in data/ is one data set; its name is a random code.",
               "Return a table with columns code, estimand, estimate and, where available,",
               "se, lower, upper, level, p_value, null and note: one row per data set and estimand.",
               "key.rds links codes to scenarios and must not be given to the analysts."),
             file.path(dir, "README.txt"))
  invisible(key)
}

#' Turn externally computed results into a method
#'
#' Wraps results computed outside R or by blinded analysts on data sets
#' exported with [sim_export_blinded()]. The returned method identifies
#' each data set by its checksum (the data passed to a method are exactly
#' the exported data), looks up the rows for its code and returns them as
#' estimates. Add it to the study (for example with [sim_amend()] or a new
#' study revision) and extend the run: the results are then evaluated like
#' any other method. A data set without a returned row is an error of that
#' method on that data set and is classified by the failure policy.
#'
#' @param results Data frame with columns `code`, `estimand`, `estimate`
#'   and optionally `se`, `lower`, `upper`, `level`, `p_value`, `null`,
#'   `note`.
#' @param key The key written by [sim_export_blinded()] (data frame or path
#'   to `key.rds`).
#' @param ... Passed to [sim_method()] (for example `description`).
#' @return A [sim_method()].
#' @export
sim_external_method <- function(results, key, ...) {
  if (is.character(key)) key <- readRDS(key)
  if (!is.data.frame(results) || !all(c("code", "estimand", "estimate") %in% names(results))) {
    stopf("`results` must be a data frame with columns code, estimand and estimate.")
  }
  if (!is.data.frame(key) || !all(c("code", "data_md5") %in% names(key))) {
    stopf("`key` must come from sim_export_blinded().")
  }
  unknown <- setdiff(unique(results$code), key$code)
  if (length(unknown)) stopf("Unknown code(s) in `results`: %s.", paste(utils::head(unknown, 5), collapse = ", "))
  cols <- intersect(c("se", "lower", "upper", "level", "p_value", "null", "note"), names(results))
  table <- results[c("code", "estimand", "estimate", cols)]
  lookup <- key[c("code", "data_md5")]
  fun <- function(data) {
    h <- adempR_data_md5(data)
    code <- lookup$code[match(h, lookup$data_md5)]
    if (is.na(code)) stop("This data set was not exported for external analysis.")
    rows <- table[table$code == code, , drop = FALSE]
    if (!nrow(rows)) stop("No external result was returned for this data set.")
    lapply(seq_len(nrow(rows)), function(i) {
      args <- as.list(rows[i, c("estimand", "estimate", cols), drop = FALSE])
      do.call(adempR::sim_estimate, args)
    })
  }
  sim_method(fun, estimands = unique(results$estimand), ...)
}

#' Checksum of a data set as seen by methods
#'
#' The checksum stored in the replication table (`data_md5`) and used by
#' [sim_external_method()] to recognize exported data sets.
#'
#' @param data A data set.
#' @return An MD5 hex string.
#' @keywords internal
#' @export
adempR_data_md5 <- function(data) md5_bytes(serialize(data, NULL, version = 3))

# Replace method names by neutral labels in every table of a results object
# and drop what identifies a method (code, description, version, artifacts).
blind_results <- function(res, seed = 1L) {
  study <- res$manifest$study
  mn <- names(study$methods)
  perm <- with_rng_state(rng_base_seed(seed, "blind-report"), sample(seq_along(mn)))
  labels <- paste0("method_", blind_letters(length(mn)))
  map <- stats::setNames(labels[order(perm)], mn)
  ren <- function(df) {
    if (!is.null(df) && "method" %in% names(df)) df$method <- unname(map[df$method])
    df
  }
  res$estimates <- ren(res$estimates)
  res$attempts <- ren(res$attempts)
  res$diagnostics <- ren(res$diagnostics)
  res$artifacts <- NULL
  ms <- lapply(study$methods, function(m) {
    m$description <- NULL
    m$version <- NULL
    m$settings_source <- NULL
    m
  })
  names(ms) <- map[names(ms)]
  study$methods <- ms[sort(names(ms))]
  man <- res$manifest
  fps <- man$fingerprints$methods
  names(fps) <- map[names(fps)]
  man$fingerprints$methods <- fps
  man$study <- study
  res$manifest <- man
  res$blinding <- data.frame(label = unname(map), method = names(map), stringsAsFactors = FALSE)
  res
}

blind_letters <- function(k) {
  out <- character(k)
  for (i in seq_len(k)) {
    n <- i
    s <- ""
    while (n > 0) {
      s <- paste0(LETTERS[(n - 1) %% 26 + 1], s)
      n <- (n - 1) %/% 26
    }
    out[i] <- s
  }
  out
}
