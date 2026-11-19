# Execution of one task: one data set, all requested methods ---------------

method_component <- function(m, name) {
  if (is.null(m$coupling)) paste0("method:", name) else paste0("coupling:", m$coupling)
}

method_stream_components <- function(study) {
  unique(vapply(names(study$methods), function(n) {
    method_component(study$methods[[n]], n)
  }, character(1)))
}

# Evaluate `expr` capturing warnings and errors and timing.
run_captured <- function(expr_fun, timeout = NULL) {
  warns <- character()
  t0 <- proc.time()
  value <- NULL
  err <- NULL
  if (!is.null(timeout)) {
    setTimeLimit(elapsed = timeout, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  }
  value <- tryCatch(
    withCallingHandlers(expr_fun(), warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      err <<- e
      NULL
    }
  )
  if (!is.null(timeout)) setTimeLimit(elapsed = Inf)
  tt <- proc.time() - t0
  timed_out <- !is.null(err) &&
    grepl("elapsed time limit", conditionMessage(err), fixed = TRUE)
  list(value = value, error = err, timed_out = timed_out, warnings = warns,
       elapsed = unname(tt[["elapsed"]]),
       cpu = unname(tt[["user.self"]] + tt[["sys.self"]]))
}

collapse_warnings <- function(w) {
  if (!length(w)) return(NA_character_)
  w <- unique(w)
  paste(utils::head(w, 5L), collapse = " | ")
}

# `task` is a list with scenario_key, replication, seeds (list with `data`
# and one entry per method stream component) and methods (names to run).
execute_task_seeded <- function(study, task, dgp_stage = "dgp") {
  key <- task$scenario_key
  rep <- as.integer(task$replication)
  idx <- scenario_index(study$scenarios, key)
  scenario <- study$scenarios$values[[idx]]
  gctx <- list(replication = rep, scenario_key = key, stage = dgp_stage)
  gen <- with_rng_state(task$seeds$data, run_captured(function() {
    out <- study$generator$fun(scenario, gctx)
    check_generated(out)
    out
  }))
  repl <- data.frame(scenario_key = key, replication = rep,
                     generation_status = if (is.null(gen$error)) "ok" else "error",
                     generation_message = if (is.null(gen$error)) NA_character_ else
                       conditionMessage(gen$error),
                     generation_warnings = collapse_warnings(gen$warnings),
                     generation_elapsed = gen$elapsed, stringsAsFactors = FALSE)
  artifacts <- list()
  attempts <- list()
  estimates <- list()
  if (!is.null(gen$error)) {
    return(list(replication = repl, attempts = NULL, estimates = NULL,
                dataset_truth = NULL, artifacts = artifacts))
  }
  generated <- gen$value
  if (inherits(generated, "adempr_data")) {
    data <- generated$data
    dtruth <- generated$truth
    meta <- generated$meta
    check_generated(data)
  } else {
    data <- generated
    dtruth <- list()
    meta <- list()
  }
  dataset_truth <- if (length(dtruth)) {
    data.frame(scenario_key = key, replication = rep, estimand = names(dtruth),
               truth = unlist(dtruth, use.names = FALSE), stringsAsFactors = FALSE)
  } else NULL
  if (isTRUE(study$retention$data)) artifacts$data <- data
  if (isTRUE(study$retention$fits)) artifacts$fits <- list()

  for (mn in task$methods) {
    m <- study$methods[[mn]]
    comp <- method_component(m, mn)
    seed <- task$seeds[[comp]]
    pub_scenario <- if (m$oracle) scenario else scenario[intersect(m$scenario_fields, names(scenario))]
    mctx <- list(replication = rep, method = mn, scenario_key = key)
    if (m$oracle) {
      mctx$truth <- dtruth
      mctx$meta <- meta
    }
    retries <- study$failure$retries
    for (attempt in seq_len(retries + 1L)) {
      mdata <- isolate_data(data)
      fit <- with_rng_state(seed, run_captured(function() {
        call_flexible(m$fun, list(data = mdata, scenario = pub_scenario, context = mctx))
      }, m$timeout))
      status <- "ok"
      message <- NA_character_
      stage <- "fit"
      rows <- NULL
      warns <- fit$warnings
      ext_elapsed <- 0
      if (!is.null(fit$error)) {
        status <- if (fit$timed_out) "timeout" else "error"
        message <- conditionMessage(fit$error)
      } else {
        stage <- "extract"
        ext <- run_captured(function() {
          out <- if (is.null(m$extract)) fit$value else
            call_flexible(m$extract, list(fit = fit$value, data = mdata,
                                          scenario = pub_scenario, context = mctx))
          outputs_to_rows(out, m$estimands)
        })
        ext_elapsed <- ext$elapsed
        warns <- c(warns, ext$warnings)
        if (!is.null(ext$error)) {
          status <- "extraction_error"
          message <- conditionMessage(ext$error)
        } else {
          rows <- ext$value
        }
      }
      attempts[[length(attempts) + 1L]] <- data.frame(
        scenario_key = key, replication = rep, method = mn,
        attempt = as.integer(attempt), status = status, failed_stage =
          if (status == "ok") NA_character_ else stage,
        message = message, n_warnings = length(warns),
        warnings = collapse_warnings(warns),
        elapsed = fit$elapsed + ext_elapsed, cpu = fit$cpu,
        timeout = m$timeout %||% NA_real_, final = FALSE, stringsAsFactors = FALSE)
      if (status == "ok" || attempt > retries || status == "timeout") break
    }
    attempts[[length(attempts)]]$final <- TRUE
    if (is.null(rows)) rows <- do.call(rbind, lapply(m$estimands, missing_row))
    rows$attempt_status <- status
    rows$had_warning <- length(warns) > 0L
    rows <- cbind(data.frame(scenario_key = key, replication = rep, method = mn,
                             stringsAsFactors = FALSE), rows)
    estimates[[length(estimates) + 1L]] <- rows
    if (isTRUE(study$retention$fits) && status == "ok") artifacts$fits[[mn]] <- fit$value
  }
  list(replication = repl, attempts = bind_rows(attempts),
       estimates = bind_rows(estimates), dataset_truth = dataset_truth,
       artifacts = artifacts)
}

# Convenience used by validation and replay: computes seeds directly.
execute_task <- function(study, ctx, key, rep, methods) {
  seeds <- task_seeds(ctx, key, rep, study, methods)[[1]]
  execute_task_seeded(study, list(scenario_key = key, replication = rep,
                                  seeds = seeds, methods = methods),
                      dgp_stage = ctx$dgp_stage)
}

# Seeds for a set of replications of one scenario. Returns a list (one per
# replication) of lists with `data` and one entry per method component.
task_seeds <- function(ctx, key, reps, study, methods) {
  base <- rng_base_seed(ctx$seed)
  comps <- unique(vapply(methods, function(mn) method_component(study$methods[[mn]], mn),
                         character(1)))
  k_data <- rng_lookup(ctx$registry, ctx$dgp_stage, key, "data")
  data_seeds <- rng_substreams(rng_stream_seeds(base, k_data, ctx$cache)[[1]], reps)
  comp_seeds <- lapply(comps, function(cp) {
    k <- rng_lookup(ctx$registry, ctx$method_stage, key, cp)
    rng_substreams(rng_stream_seeds(base, k, ctx$cache)[[1]], reps)
  })
  names(comp_seeds) <- comps
  lapply(seq_along(reps), function(j) {
    s <- list(data = data_seeds[[j]])
    for (cp in comps) s[[cp]] <- comp_seeds[[cp]][[j]]
    s
  })
}
