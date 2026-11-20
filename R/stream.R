# Out-of-core summaries ----------------------------------------------------

#' Summarize a run one batch at a time
#'
#' [sim_summarise()] reads every committed estimate into memory. For runs
#' too large for that, `sim_summarise_stream()` reads one committed batch
#' at a time and applies the same selection (stored fingerprints, the latest
#' retry round, deduplication of repeated deliveries with the same conflict
#' check, data set truths of the selected generation) and the same
#' eligibility rules, and adds each eligible contribution to mergeable
#' accumulators ([sim_accumulator()]) that carry the moments each kernel
#' needs. The values, MCSEs, intervals and statuses come from the same
#' kernels as in [sim_summarise()], so they agree up to rounding, and every
#' planned cell is returned, with status `"no_eligible"` when nothing
#' could be used. Memory grows with the number of tasks (short keys to
#' resolve superseded and duplicated rows, and the data set truths) and the
#' number of cells, not with the size of the stored results.
#'
#' @param run A run.
#' @param measures Any of `"bias"`, `"mse"`, `"empse"`, `"mean_se"` and
#'   `"coverage"`.
#' @return A data frame with `scenario_key`, `scenario_id`, `method`,
#'   `estimand`, `measure`, `value`, `mcse`, `lower`, `upper`, `n`,
#'   `status` and `note`.
#' @export
sim_summarise_stream <- function(run, measures = c("bias", "mse", "empse", "coverage")) {
  stopifnot(inherits(run, "adempr_run"))
  measures <- match.arg(measures, c("bias", "mse", "empse", "mean_se", "coverage"), several.ok = TRUE)
  man <- run$manifest
  study <- man$study
  level <- study$level
  p <- run_paths(run$path)
  keys <- study$scenarios$keys
  R <- planned_replications(man, keys)
  fps <- unlist(man$fingerprints$methods)
  side <- list.files(p$batches, pattern = "\\.md5$", full.names = TRUE)
  files <- sub("\\.md5$", ".rds", side)
  inv <- file.path(p$batches, "INVENTORY")
  listed <- if (file.exists(inv)) unique(readLines(inv, warn = FALSE)) else character()
  gone <- setdiff(listed[nzchar(listed)], sub("\\.md5$", "", basename(side)))
  if (length(gone)) {
    warning(sprintf("committed batch(es) %s are missing from the run directory; resume the run to rerun their tasks",
                    paste(gone, collapse = ", ")), call. = FALSE)
  }
  read_batch <- function(f) {
    side <- sub("\\.rds$", ".md5", f)
    ok <- file.exists(f) && identical(file_md5(f), readLines(side, warn = FALSE)[1])
    if (!ok) {
      warning(sprintf("corrupt or missing file %s ignored", basename(f)), call. = FALSE)
      return(NULL)
    }
    o <- readRDS(f)
    if (!identical(o$schema_version, SCHEMA_VERSION)) {
      stopf("A committed batch uses result schema %s; expected %s.", o$schema_version, SCHEMA_VERSION)
    }
    o
  }
  keep <- function(df) {
    if (is.null(df) || !nrow(df)) return(NULL)
    if (is.null(df$round)) df$round <- 0L
    df <- df[df$scenario_key %in% keys & df$method %in% names(fps), , drop = FALSE]
    df <- df[df$replication <= R[df$scenario_key], , drop = FALSE]
    df[df$method_fp == fps[df$method], , drop = FALSE]
  }
  task_id <- function(df) paste(df$scenario_key, df$replication, df$method, sep = "\r")
  sci <- c(RAW_ESTIMATE_COLUMNS, "attempt_status", "had_warning", "fallback")
  # pass 1: the latest retry round of every task, the signature of every
  # delivered row (repeated deliveries must be identical, as in
  # sim_results()), and the generations and data set truths
  last <- new.env(hash = TRUE, parent = emptyenv())
  sig <- new.env(hash = TRUE, parent = emptyenv())
  reps <- list()
  dtruth <- list()
  for (f in files) {
    o <- read_batch(f)
    if (is.null(o)) next
    e <- keep(o$estimates)
    if (!is.null(e)) {
      ids <- task_id(e)
      top <- tapply(e$round, ids, max)
      for (i in names(top)) if (is.null(last[[i]]) || top[[i]] > last[[i]]) assign(i, top[[i]], envir = last)
      rk <- paste(ids, e$estimand, e$round, sep = "\r")
      sg <- do.call(paste, c(lapply(e[intersect(sci, names(e))], col_text), sep = "\r"))
      for (j in seq_along(rk)) {
        prev <- sig[[rk[j]]]
        if (is.null(prev)) assign(rk[j], sg[j], envir = sig)
        else if (!identical(prev, sg[j])) {
          stopf("Conflicting duplicate estimate results: the same task (scenario_key=%s, replication=%s, method=%s, estimand=%s) was committed twice with different scientific content. The run directory is inconsistent; inspect the batches before continuing.",
                e$scenario_key[j], e$replication[j], e$method[j], e$estimand[j])
        }
      }
    }
    r <- o$replications
    if (!is.null(r) && nrow(r)) {
      if (is.null(r$round)) r$round <- 0L
      reps[[length(reps) + 1L]] <- r[intersect(c("scenario_key", "replication", "generator_fp", "round",
                                                 "generation_status", "generation_message", "data_md5"),
                                               names(r))]
    }
    d <- o$dataset_truth
    if (!is.null(d) && nrow(d)) {
      if (is.null(d$round)) d$round <- 0L
      dtruth[[length(dtruth) + 1L]] <- d
    }
  }
  dtruth <- select_dataset_truth(bind_rows(dtruth), select_generations(bind_rows(reps)))
  truth <- sim_truth(run)
  truth <- truth[truth$truth_status == "ok", c("scenario_key", "estimand", "truth", "truth_se",
                                               "truth_kind"), drop = FALSE]
  ds_est <- names(study$estimands)[vapply(study$estimands, function(e) e$truth$kind == "dataset", logical(1))]
  # pass 2: accumulate eligible contributions, once per task and estimand
  seen <- new.env(hash = TRUE, parent = emptyenv())
  acc <- new.env(hash = TRUE, parent = emptyenv())
  # With antithetic pairs, the pair (the independent unit) can span
  # batches, so the sums of every pair are kept (one entry per pair and
  # cell, like the task keys kept above) and the MCSE comes from them at
  # the end, as in sim_summarise().
  anti <- isTRUE(study$generator$antithetic)
  add <- function(cell, what, x, reps = NULL) {
    k <- paste(cell, what, sep = "\r")
    a <- acc[[k]] %||% list(acc = sim_accumulator(), pos = 0, neg = 0, other = 0,
                            tmin = Inf, tmax = -Inf, levels = numeric(),
                            units = if (anti) new.env(hash = TRUE, parent = emptyenv()), shift = NA_real_)
    fin <- is.finite(x)
    a$acc <- sim_accumulate(a$acc, x[fin])
    a$pos <- a$pos + sum(x == Inf, na.rm = TRUE)
    a$neg <- a$neg + sum(x == -Inf, na.rm = TRUE)
    a$other <- a$other + sum(is.na(x))
    if (anti && any(fin)) {
      # sums of x - shift and its square per pair; the shift (the first
      # value seen) keeps the squares of large values well conditioned
      if (is.na(a$shift)) a$shift <- x[fin][1]
      xs <- x[fin] - a$shift
      u <- as.character(rep_units(study, reps[fin]))
      s1 <- rowsum(xs, u, reorder = FALSE)
      s2 <- rowsum(xs^2, u, reorder = FALSE)
      nk <- rowsum(rep(1, length(u)), u, reorder = FALSE)
      for (j in seq_len(nrow(s1))) {
        id <- rownames(s1)[j]
        old <- a$units[[id]] %||% c(0, 0, 0)
        assign(id, old + c(s1[j, 1], s2[j, 1], nk[j, 1]), envir = a$units)
      }
    }
    assign(k, a, envir = acc)
  }
  # pair based MCSE of a streamed measure
  unit_stream <- function(r, a, m) {
    if (!anti || is.null(a) || !r$status %in% c("ok", "boundary") || !is.finite(r$value) ||
        a$pos + a$neg + a$other > 0) return(r)
    u <- do.call(rbind, mget(ls(a$units), envir = a$units))
    if (is.null(u)) return(r)
    if (m == "empse") {
      # leave-one-pair-out SD from the remaining sums
      n <- sum(u[, 3]) - u[, 3]
      s1 <- sum(u[, 1]) - u[, 1]
      s2 <- sum(u[, 2]) - u[, 2]
      th <- sqrt(pmax(0, (s2 - s1^2 / n) / (n - 1)))
      th[n < 2] <- NA_real_
      return(unit_jackknife(r, th, level))
    }
    means <- u[, 1] / u[, 3] + a$shift
    unit_jackknife(r, unit_loo_means(u[, 1], u[, 3]) + a$shift, level,
                   proportion = m == "coverage", unit_means = means)
  }
  note_cell <- function(cell, truth, lv) {
    k <- paste(cell, ".meta", sep = "\r")
    a <- acc[[k]] %||% list(tmin = Inf, tmax = -Inf, levels = numeric())
    if (length(truth)) {
      a$tmin <- min(a$tmin, truth)
      a$tmax <- max(a$tmax, truth)
    }
    a$levels <- sort(unique(c(a$levels, lv[!is.na(lv)])))
    assign(k, a, envir = acc)
  }
  for (f in files) {
    o <- read_batch(f)
    if (is.null(o)) next
    e <- keep(o$estimates)
    if (is.null(e)) next
    ids <- task_id(e)
    e <- e[e$round == unlist(mget(ids, envir = last)), , drop = FALSE]
    rk <- paste(task_id(e), e$estimand, sep = "\r")
    fresh <- !vapply(rk, exists, logical(1), envir = seen, inherits = FALSE) & !duplicated(rk)
    e <- e[fresh, , drop = FALSE]
    for (k in rk[fresh]) assign(k, TRUE, envir = seen)
    if (!nrow(e)) next
    e <- merge(e, truth, by = c("scenario_key", "estimand"), all.x = TRUE, sort = FALSE)
    if (!is.null(dtruth) && nrow(dtruth)) {
      idx <- match(paste(e$scenario_key, e$replication, e$estimand),
                   paste(dtruth$scenario_key, dtruth$replication, dtruth$estimand))
      use <- e$estimand %in% ds_est & !is.na(idx)
      e$truth[use] <- dtruth$truth[idx[use]]
      e$truth_se[use] <- 0
      e$truth_kind[use] <- "dataset"
    }
    miss <- is.na(e$truth_kind)
    e$truth_kind[miss] <- vapply(e$estimand[miss], function(id) study$estimands[[id]]$truth$kind, character(1))
    e <- classify_rows(e, study$estimands, study$failure)
    cell <- paste(e$scenario_key, e$method, e$estimand, sep = "\r")
    # the eligibility masks of measure_cell()
    usable <- e$valid
    has_truth <- !is.na(e$truth)
    se_ok <- e$se_status == "valid"
    int_ok <- e$estimate_status %in% USABLE_FOR_INTERVAL & e$interval_status == "valid"
    err <- e$estimate - e$truth
    for (cl in unique(cell)) {
      i <- cell == cl
      note_cell(cl, e$truth[i & usable & has_truth], e$level[i & int_ok & has_truth])
      if ("bias" %in% measures) add(cl, "bias", err[i & usable & has_truth], e$replication[i & usable & has_truth])
      if ("mse" %in% measures) add(cl, "mse", err[i & usable & has_truth]^2, e$replication[i & usable & has_truth])
      if ("empse" %in% measures) {
        add(cl, "empse", e$estimate[i & usable], e$replication[i & usable])
        add(cl, "empse_err", err[i & usable & has_truth], e$replication[i & usable & has_truth])
      }
      if ("mean_se" %in% measures) add(cl, "mean_se", e$se[i & usable & se_ok], e$replication[i & usable & se_ok])
      if ("coverage" %in% measures) {
        j <- i & int_ok & has_truth
        add(cl, "coverage", as.numeric(e$lower[j] <= e$truth[j] & e$truth[j] <= e$upper[j]), e$replication[j])
      }
    }
  }
  # every planned cell, also those with no eligible contribution
  out <- list()
  for (k in keys) for (mn in names(study$methods)) for (en in study$methods[[mn]]$estimands) {
    cl <- paste(k, mn, en, sep = "\r")
    meta <- acc[[paste(cl, ".meta", sep = "\r")]]
    varying <- identical(study$estimands[[en]]$truth$kind, "dataset")
    for (m in measures) {
      what <- if (m == "empse" && varying) "empse_err" else m
      a <- acc[[paste(cl, what, sep = "\r")]]
      st <- if (is.null(a)) list(n = 0, mean_x = NA_real_, var_x = NA_real_, m4_x = NA_real_) else sim_accumulated(a$acc)
      n <- st$n + if (is.null(a)) 0 else a$pos + a$neg + a$other
      r <- if (!measure_applies(m, study$estimands[[en]]$type)) {
        kernel_result(status = "not_applicable",
                      note = sprintf("measure does not apply to a %s estimand", study$estimands[[en]]$type))
      } else switch(m,
        empse = {
          r <- k_empse_moments(n, st$var_x, st$m4_x, if (is.null(a)) 0 else a$pos + a$neg + a$other, level)
          if (varying) r$note <- paste(c("SD of the errors (data set truth)",
                                         if (!is.na(r$note)) r$note), collapse = "; ")
          r
        },
        coverage = {
          r <- k_prop_counts(if (st$n) round(st$mean_x * st$n) else 0, n, level)
          lv <- if (is.null(meta)) numeric() else meta$levels
          r$note <- paste(c(if (length(lv) == 1L) paste("nominal level", format(lv)),
                            if (length(lv) > 1L) paste("intervals at several nominal levels:",
                                                       paste(format(lv), collapse = ", ")),
                            if (!is.na(r$note)) r$note), collapse = "; ")
          if (!nzchar(r$note)) r$note <- NA_character_
          r
        },
        k_mean_moments(n, st$mean_x, st$var_x, if (is.null(a)) 0 else a$pos,
                       if (is.null(a)) 0 else a$neg, level))
      r <- unit_stream(r, a, m)
      out[[length(out) + 1L]] <- data.frame(scenario_key = k,
        scenario_id = study$scenarios$labels[match(k, keys)], method = mn, estimand = en,
        measure = m, value = r$value, mcse = r$mcse, lower = r$lower, upper = r$upper,
        n = as.integer(r$n), status = r$status, note = r$note, stringsAsFactors = FALSE)
    }
  }
  out <- bind_rows(out)
  rownames(out) <- NULL
  out
}
