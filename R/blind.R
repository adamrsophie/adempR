# Blinded evaluation ------------------------------------------------------

#' Export coded data sets for blinded analysis
#'
#' Regenerates the data sets of a run from their stored streams and writes
#' each to `dir/data/<code>.rds`, where the code is random and reveals
#' neither the scenario nor the replication. Only the data that methods
#' see are written (no truth, no metadata). The key linking codes to
#' scenarios, replications and data checksums is written to `key`, outside
#' `dir`, so the directory can be handed to the analysts as it is.
#' Analysts return one row per data set and estimand, which
#' [sim_external_method()] turns into a method of the study.
#'
#' By default every committed data set that was generated successfully is
#' exported: the replications of each scenario (including extra
#' replications of adaptively extended scenarios) and, for a sequential
#' run, those up to the analyzed look. A regenerated data set whose
#' checksum differs from the one recorded when the run made it (changed
#' generator inputs) is not exported. Data sets left out are listed in the
#' attribute `omitted` of the key, with the reason.
#'
#' @param run A run.
#' @param dir Output directory (must not exist or be empty).
#' @param key File for the key; defaults to `<dir>-key.rds` next to `dir`.
#' @param scenarios Scenario labels or keys to export; defaults to all.
#' @param replications Replications to export; defaults to all committed
#'   ones.
#' @param seed Seed of the random codes. The default draws a fresh seed,
#'   so the codes cannot be rebuilt from the package; it is stored in the
#'   key. The caller's random state is left untouched.
#' @return Invisibly, the key: a data frame with `code`, `scenario_key`,
#'   `replication` and `data_md5`, with attributes `seed` and `omitted`.
#' @export
sim_export_blinded <- function(run, dir, key = NULL, scenarios = NULL, replications = NULL,
                               seed = NULL) {
  stopifnot(inherits(run, "adempr_run"))
  if (!is_scalar_character(dir)) stopf("`dir` must be a path.")
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE))) {
    stopf("%s is not empty.", dir)
  }
  key_path <- key %||% file.path(dirname(normalizePath(dir, mustWork = FALSE)),
                                 paste0(basename(dir), "-key.rds"))
  if (!is_scalar_character(key_path)) stopf("`key` must be a path.")
  # normalizePath() leaves missing paths as given, with mixed separators on
  # Windows, so compare on one separator (and one case on Windows)
  canon <- function(x) {
    x <- gsub("[\\/]+", "/", normalizePath(x, winslash = "/", mustWork = FALSE))
    if (.Platform$OS.type == "windows") tolower(x) else x
  }
  if (startsWith(canon(key_path), paste0(sub("/$", "", canon(dir)), "/"))) {
    stopf("`key` must be outside `dir`, which is given to the analysts.")
  }
  seed <- seed %||% fresh_seed()
  man <- run$manifest
  study <- man$study
  keys <- if (is.null(scenarios)) study$scenarios$keys else
    vapply(scenarios, function(s) resolve_scenario(study, s), character(1))
  if (!is.null(replications) &&
      !all(vapply(replications, is_count, logical(1), min = 1))) {
    stopf("`replications` must be positive whole numbers.")
  }
  # the committed generations (per scenario counts, sequential prefix)
  reps_tab <- sim_results(run)$replications
  reps_tab <- reps_tab[reps_tab$scenario_key %in% keys, , drop = FALSE]
  omitted <- list()
  if (!is.null(replications)) {
    planned <- expand.grid(replication = as.integer(replications), scenario_key = keys,
                           stringsAsFactors = FALSE)
    have <- paste(reps_tab$scenario_key, reps_tab$replication)
    miss <- planned[!paste(planned$scenario_key, planned$replication) %in% have, , drop = FALSE]
    if (nrow(miss)) omitted[[length(omitted) + 1L]] <- data.frame(miss[c("scenario_key", "replication")],
      reason = "not generated in the run", stringsAsFactors = FALSE)
    reps_tab <- reps_tab[reps_tab$replication %in% replications, , drop = FALSE]
  }
  failed <- reps_tab$generation_status != "ok"
  if (any(failed)) omitted[[length(omitted) + 1L]] <- data.frame(
    scenario_key = reps_tab$scenario_key[failed], replication = reps_tab$replication[failed],
    reason = "generation failed", stringsAsFactors = FALSE)
  tasks <- reps_tab[!failed, c("scenario_key", "replication", "data_md5"), drop = FALSE]
  tasks <- tasks[order(match(tasks$scenario_key, keys), tasks$replication), , drop = FALSE]
  codes <- with_rng_state(rng_base_seed(seed, "blind"), {
    hex <- c(0:9, letters[1:6])
    unique_codes <- character()
    while (length(unique_codes) < nrow(tasks)) {
      unique_codes <- unique(c(unique_codes, vapply(seq_len(nrow(tasks) - length(unique_codes)),
        function(i) paste0(sample(letters[1:6], 1L),
                           paste(sample(hex, 11, replace = TRUE), collapse = "")),
        character(1))))
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
      if (!identical(out$replication$generation_status, "ok")) {
        omitted[[length(omitted) + 1L]] <- data.frame(scenario_key = k, replication = r,
          reason = "regeneration failed", stringsAsFactors = FALSE)
        next
      }
      if (!is.na(tasks$data_md5[i]) && !identical(out$replication$data_md5, tasks$data_md5[i])) {
        omitted[[length(omitted) + 1L]] <- data.frame(scenario_key = k, replication = r,
          reason = "regenerated data differ from the recorded checksum", stringsAsFactors = FALSE)
        next
      }
      md5[i] <- out$replication$data_md5
      saveRDS(out$data, file.path(dir, "data", paste0(codes[i], ".rds")), version = 3)
    }
  })
  out <- data.frame(code = codes, scenario_key = tasks$scenario_key,
                    replication = tasks$replication, data_md5 = md5, stringsAsFactors = FALSE)
  out <- out[!is.na(out$data_md5), , drop = FALSE]
  rownames(out) <- NULL
  om <- bind_rows(omitted)
  attr(out, "seed") <- seed
  attr(out, "omitted") <- if (is.null(om)) data.frame(scenario_key = character(), replication = integer(),
                                                      reason = character()) else om
  dir.create(dirname(key_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, key_path, version = 3)
  writeLines(c("Each file in data/ is one data set; its name is a random code.",
               "Return a table with columns code, estimand, estimate and, where available,",
               "se, lower, upper, level, p_value, null and note: one row per data set and estimand."),
             file.path(dir, "README.txt"))
  invisible(out)
}

#' Turn externally computed results into a method
#'
#' Wraps results computed outside R or by blinded analysts on data sets
#' exported with [sim_export_blinded()]. The returned method finds the
#' code of each data set from its scenario and replication (so identical
#' data sets in different replications keep their own results), checks
#' that the data equal the exported data by their checksum, and returns
#' the rows for that code as estimates. The tables are validated when the
#' method is built: codes and data sets must be unique in the key, and
#' there must be at most one row per code and estimand. Add it to the study (for example with [sim_amend()] or a new
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
  if (!is.data.frame(key) ||
      !all(c("code", "scenario_key", "replication", "data_md5") %in% names(key))) {
    stopf("`key` must come from sim_export_blinded().")
  }
  if (anyDuplicated(key$code)) stopf("The key has duplicated codes.")
  task <- paste(key$scenario_key, key$replication, sep = "\r")
  if (anyDuplicated(task)) stopf("The key lists a scenario and replication more than once.")
  unknown <- setdiff(unique(results$code), key$code)
  if (length(unknown)) stopf("Unknown code(s) in `results`: %s.", paste(utils::head(unknown, 5), collapse = ", "))
  dup <- duplicated(results[c("code", "estimand")])
  if (any(dup)) {
    stopf("`results` has more than one row for code and estimand: %s.",
          paste(utils::head(paste(results$code[dup], results$estimand[dup]), 5), collapse = ", "))
  }
  if (!is.character(results$estimand) && !is.factor(results$estimand)) stopf("`estimand` must be text.")
  cols <- intersect(c("se", "lower", "upper", "level", "p_value", "null", "note"), names(results))
  for (cn in c("estimate", setdiff(cols, "note"))) {
    if (!is.numeric(results[[cn]]) && !all(is.na(results[[cn]]))) stopf("Column `%s` of `results` must be numeric.", cn)
  }
  if ("note" %in% cols && !is.character(results$note) && !all(is.na(results$note))) {
    stopf("Column `note` of `results` must be text.")
  }
  results$estimand <- as.character(results$estimand)
  table <- results[c("code", "estimand", "estimate", cols)]
  lookup <- data.frame(task = task, code = key$code, data_md5 = key$data_md5, stringsAsFactors = FALSE)
  fun <- function(data, context) {
    i <- match(paste(context$scenario_key, context$replication, sep = "\r"), lookup$task)
    if (is.na(i)) stop("This data set was not exported for external analysis.")
    if (!identical(adempR_data_md5(data), lookup$data_md5[i])) {
      stop("This data set differs from the one exported for external analysis.")
    }
    code <- lookup$code[i]
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
blind_results <- function(res, seed) {
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
  res$attempts_all <- ren(res$attempts_all)
  res$diagnostics <- ren(res$diagnostics)
  # free text (messages, warnings) can name methods
  for (tb in c("attempts", "attempts_all")) {
    if (!is.null(res[[tb]])) res[[tb]]$message <- res[[tb]]$warnings <- NA_character_
  }
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
  # fingerprints identify methods to anyone who has the study code
  man$fingerprints$methods <- NULL
  man$revisions <- NULL
  study$amendments <- NULL
  if (!is.null(man$sequential$looks)) man$sequential$looks <- ren(man$sequential$looks)
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

# Rename the method column of a table by a blinding key.
relabel_methods <- function(df, map) {
  if (is.data.frame(df) && "method" %in% names(df)) df$method <- unname(map[df$method])
  df
}

# A seed nobody chose: from the clock and the process, without touching
# the caller's random number generator.
fresh_seed <- function() {
  strtoi(substr(md5_string(c(format(Sys.time(), "%OS6"), Sys.getpid(), new_uuid())), 1, 7), 16L)
}
