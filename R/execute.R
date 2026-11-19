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
  # The message of R's time limit error is translated, so compare with both
  # the English text and its translation in the current language.
  limit_msgs <- unique(c("reached elapsed time limit",
                         gettext("reached elapsed time limit", domain = "R")))
  # Only the limit's own error is a timeout. An ordinary error raised after
  # the deadline (compiled code that ran past the cooperative limit) stays
  # an error; the attempt records the overrun separately.
  timed_out <- !is.null(err) && !is.null(timeout) &&
    any(vapply(limit_msgs, grepl, logical(1), x = conditionMessage(err), fixed = TRUE))
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
execute_task_seeded <- function(study, task, dgp_stage = "dgp", keep_data = FALSE) {
  key <- task$scenario_key
  rep <- as.integer(task$replication)
  idx <- scenario_index(study$scenarios, key)
  scenario <- study$scenarios$values[[idx]]
  gctx <- list(replication = rep, scenario_key = key, stage = dgp_stage)
  if (isTRUE(study$generator$antithetic)) gctx$antithetic <- rep %% 2L == 0L
  gpk <- study$generator$packages
  gen <- with_rng_state(task$seeds$data, run_captured(function() with_packages(gpk, {
    out <- study$generator$fun(scenario, gctx)
    check_generated(out)
    if (inherits(out, "adempr_data")) check_generated(out$data)
    out
  })))
  repl <- data.frame(scenario_key = key, replication = rep,
                     generation_status = if (is.null(gen$error)) "ok" else "error",
                     generation_message = if (is.null(gen$error)) NA_character_ else
                       conditionMessage(gen$error),
                     generation_warnings = collapse_warnings(gen$warnings),
                     generation_elapsed = gen$elapsed, data_md5 = NA_character_,
                     stringsAsFactors = FALSE)
  artifacts <- list()
  attempts <- list()
  estimates <- list()
  diagnostics <- list()
  if (!is.null(gen$error)) {
    return(list(replication = repl, attempts = NULL, estimates = NULL,
                dataset_truth = NULL, artifacts = artifacts, diagnostics = NULL))
  }
  generated <- gen$value
  if (inherits(generated, "adempr_data")) {
    data <- generated$data
    dtruth <- generated$truth
    meta <- generated$meta
  } else {
    data <- generated
    dtruth <- list()
    meta <- list()
  }
  # checksum of the analysis data: identifies the data set across batches
  # and reveals a generator that returns the same data every time
  repl$data_md5 <- md5_bytes(serialize(data, NULL, version = 3))
  # A replication that already has committed fits is regenerated only to
  # run further methods on it (a retry, an added or changed method). Its
  # data must be the committed data set, or the methods would be compared
  # on different data under one replication number.
  if (length(task$expect_md5) && !is.na(task$expect_md5) &&
      !identical(task$expect_md5, repl$data_md5)) {
    stopf("Replication %d of scenario %s generated different data than the committed data set (checksum %s, now %s). The generator depends on something that changed since (a file, a global variable, a package version). Restore it, or use a new run directory.",
          rep, short_hash(key), short_hash(task$expect_md5), short_hash(repl$data_md5))
  }
  records <- NULL
  if (!is.null(study$generator$record)) {
    rc <- with_rng_state(task$seeds$data, run_captured(function() {
      with_packages(gpk, call_flexible(study$generator$record,
                                       list(data = isolate_data(data), meta = meta)))
    }))
    records <- callback_rows(rc, ".record_error",
                             data.frame(scenario_key = key, replication = rep,
                                        stringsAsFactors = FALSE))
  }
  dataset_truth <- if (length(dtruth)) {
    data.frame(scenario_key = key, replication = rep, estimand = names(dtruth),
               truth = unlist(dtruth, use.names = FALSE), stringsAsFactors = FALSE)
  } else NULL
  if (isTRUE(study$retention$data)) artifacts$data <- data
  if (isTRUE(study$retention$fits)) artifacts$fits <- list()

  old_cuda <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA)
  restore_cuda <- function() {
    if (is.na(old_cuda)) Sys.unsetenv("CUDA_VISIBLE_DEVICES") else Sys.setenv(CUDA_VISIBLE_DEVICES = old_cuda)
  }
  on.exit(restore_cuda(), add = TRUE)
  for (mn in task$methods) {
    m <- study$methods[[mn]]
    comp <- method_component(m, mn)
    seed <- task$seeds[[comp]]
    pub_scenario <- if (m$oracle) scenario else scenario[intersect(m$scenario_fields, names(scenario))]
    mctx <- list(replication = rep, method = mn, scenario_key = key,
                 inner_seeds = inner_seed_fun(seed))
    if (m$oracle) {
      mctx$truth <- dtruth
      mctx$meta <- meta
    }
    # a method that declared gpu = TRUE sees only its assigned device while
    # it runs; the variable is restored afterwards, and on any exit
    gpu_dev <- if (isTRUE(m$gpu) && !is.null(task$gpu)) as.character(task$gpu) else NA_character_
    if (!is.na(gpu_dev)) {
      mctx$gpu <- gpu_dev
      Sys.setenv(CUDA_VISIBLE_DEVICES = gpu_dev)
    }
    retries <- study$failure$retries
    for (attempt in seq_len(retries + 1L)) {
      mdata <- isolate_data(data)
      # Fit and extraction share the method's stream: random numbers drawn in
      # `extract` continue the stream of the fit and never touch the caller's
      # generator.
      both <- with_rng_state(seed, {
        fit <- run_captured(function() {
          with_packages(m$packages,
            call_flexible(m$fun, list(data = mdata, scenario = pub_scenario, context = mctx)))
        }, m$timeout)
        # the limit covers the whole attempt: extraction gets what the fit
        # left of it
        ext_limit <- if (is.null(m$timeout)) NULL else max(m$timeout - fit$elapsed, 0.01)
        ext <- if (is.null(fit$error)) run_captured(function() {
          out <- if (is.null(m$extract)) fit$value else with_packages(m$packages,
            call_flexible(m$extract, list(fit = fit$value, data = mdata,
                                          scenario = pub_scenario, context = mctx)))
          outputs_to_rows(out, m$estimands, study$estimands)
        }, ext_limit) else NULL
        list(fit = fit, ext = ext)
      })
      fit <- both$fit
      ext <- both$ext
      status <- "ok"
      message <- NA_character_
      stage <- "fit"
      rows <- NULL
      warns <- fit$warnings
      ext_elapsed <- 0
      if (!is.null(fit$error)) {
        message <- conditionMessage(fit$error)
        status <- if (fit$timed_out) "timeout" else if (is_memory_error(message)) "memory" else "error"
      } else {
        stage <- "extract"
        ext_elapsed <- ext$elapsed
        warns <- c(warns, ext$warnings)
        if (!is.null(ext$error)) {
          status <- if (ext$timed_out) "timeout" else "extraction_error"
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
        timeout = m$timeout %||% NA_real_,
        # the limit is cooperative: compiled code can run past it unnoticed
        overrun = !is.null(m$timeout) && status != "timeout" && fit$elapsed + ext_elapsed > m$timeout,
        final = FALSE, fallback = FALSE, gpu = gpu_dev, stringsAsFactors = FALSE)
      if (status == "ok" || attempt > retries || status == "timeout") break
    }
    attempts[[length(attempts)]]$final <- TRUE
    used_fallback <- FALSE
    if (status != "ok" && !is.null(m$fallback)) {
      fb <- with_rng_state(seed, run_captured(function() {
        out <- with_packages(m$packages,
          call_flexible(m$fallback, list(data = isolate_data(data), scenario = pub_scenario,
                                         context = mctx)))
        outputs_to_rows(out, m$estimands, study$estimands)
      }, m$timeout))
      # the fallback is recorded as an attempt of its own: its outcome,
      # warnings and cost
      fb_status <- if (is.null(fb$error)) "ok" else if (fb$timed_out) "timeout" else
        if (is_memory_error(conditionMessage(fb$error))) "memory" else "error"
      attempts[[length(attempts) + 1L]] <- data.frame(
        scenario_key = key, replication = rep, method = mn,
        attempt = as.integer(attempt) + 1L, status = fb_status,
        failed_stage = if (fb_status == "ok") NA_character_ else "fallback",
        message = if (is.null(fb$error)) NA_character_ else conditionMessage(fb$error),
        n_warnings = length(fb$warnings), warnings = collapse_warnings(fb$warnings),
        elapsed = fb$elapsed, cpu = fb$cpu, timeout = m$timeout %||% NA_real_,
        overrun = !is.null(m$timeout) && fb_status != "timeout" && fb$elapsed > m$timeout,
        final = FALSE, fallback = TRUE, gpu = gpu_dev, stringsAsFactors = FALSE)
      if (is.null(fb$error)) {
        rows <- fb$value
        used_fallback <- TRUE
        # the warnings policy applies to the output that is used
        warns <- fb$warnings
      }
    }
    if (status == "ok" && !is.null(m$diagnostics)) {
      dg <- with_rng_state(seed, run_captured(function() {
        with_packages(m$packages, call_flexible(m$diagnostics, list(fit = fit$value, data = mdata)))
      }))
      if (!is.null(dg$error) || length(unlist(dg$value))) {
        diagnostics[[length(diagnostics) + 1L]] <- callback_rows(dg, ".diagnostics_error",
          data.frame(scenario_key = key, replication = rep, method = mn,
                     stringsAsFactors = FALSE))
      }
    }
    restore_cuda()
    if (is.null(rows)) rows <- do.call(rbind, lapply(m$estimands, missing_row))
    rows$fallback <- used_fallback
    rows$attempt_status <- status
    rows$had_warning <- length(warns) > 0L
    rows <- cbind(data.frame(scenario_key = key, replication = rep, method = mn,
                             stringsAsFactors = FALSE), rows)
    estimates[[length(estimates) + 1L]] <- rows
    if (isTRUE(study$retention$fits) && status == "ok") artifacts$fits[[mn]] <- fit$value
  }
  list(replication = repl, attempts = bind_rows(attempts),
       estimates = bind_rows(estimates), dataset_truth = dataset_truth,
       artifacts = artifacts, diagnostics = bind_rows(diagnostics), records = records,
       data = if (keep_data) data, meta = if (keep_data) {
         if (inherits(generated, "adempr_data")) generated[c("truth", "meta")] else list()
       })
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
  base <- rng_base_seed(ctx$seed, ctx$domain)
  comps <- unique(vapply(methods, function(mn) method_component(study$methods[[mn]], mn),
                         character(1)))
  k_data <- rng_lookup(ctx$registry, ctx$dgp_stage, data_stream_key(study, key), "data")
  # antithetic pairs share the data substream of their first replication
  data_reps <- if (isTRUE(study$generator$antithetic)) 2L * ((reps + 1L) %/% 2L) - 1L else reps
  data_seeds <- rng_substreams(rng_stream_seeds(base, k_data, ctx$cache)[[1]], data_reps)
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

# Rows of a `record` or `diagnostics` callback: one per named number. An
# error, unnamed or duplicated names, or a value that is not a number gives
# a single error row (name `err_name`) with the reason, never silent `NA`s.
callback_rows <- function(cb, err_name, id) {
  v <- if (is.null(cb$error)) unlist(cb$value) else NULL
  reason <- if (!is.null(cb$error)) conditionMessage(cb$error) else
    if (!length(v)) "returned nothing" else
    if (is.null(names(v)) || any(!nzchar(names(v))) || anyDuplicated(names(v))) "values need unique names" else
    if (!(is.numeric(v) || is.logical(v))) "values must be numbers" else NULL
  if (!is.null(reason)) {
    return(cbind(id, data.frame(name = err_name, value = NA_real_, message = reason,
                                stringsAsFactors = FALSE)))
  }
  cbind(id[rep(1L, length(v)), , drop = FALSE],
        data.frame(name = names(v), value = as.numeric(v), message = NA_character_,
                   stringsAsFactors = FALSE), row.names = NULL)
}

# Memory exhaustion is a failure class of its own: it says more about the
# computing environment than about the method.
is_memory_error <- function(msg) {
  grepl("cannot allocate|memory exhausted|vector memory limit|std::bad_alloc", msg,
        ignore.case = TRUE)
}
