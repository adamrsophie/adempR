# Out-of-core summaries ----------------------------------------------------

#' Summarize a run one batch at a time
#'
#' [sim_summarise()] reads every committed estimate into memory. For runs
#' too large for that, `sim_summarise_stream()` reads one committed batch
#' at a time, applies the same selection (stored fingerprints, the latest
#' retry round, deduplication of repeated deliveries) and validity rules,
#' and adds each usable row to mergeable accumulators ([sim_accumulator()]).
#' Memory then grows with the number of tasks (one short key per task, to
#' resolve superseded rows) and the number of cells, not with the size of
#' the stored results. Values agree with [sim_summarise()] up to rounding.
#'
#' @param run A run.
#' @param measures Any of `"bias"`, `"mse"`, `"empse"`, `"mean_se"` and
#'   `"coverage"`.
#' @return A data frame with `scenario_key`, `scenario_id`, `method`,
#'   `estimand`, `measure`, `value`, `mcse` and `n`.
#' @export
sim_summarise_stream <- function(run, measures = c("bias", "mse", "empse", "coverage")) {
  stopifnot(inherits(run, "adempr_run"))
  measures <- match.arg(measures, c("bias", "mse", "empse", "mean_se", "coverage"), several.ok = TRUE)
  man <- run$manifest
  study <- man$study
  p <- run_paths(run$path)
  keys <- study$scenarios$keys
  R <- planned_replications(man, keys)
  fps <- unlist(man$fingerprints$methods)
  files <- sub("\\.md5$", ".rds", list.files(p$batches, pattern = "\\.md5$", full.names = TRUE))
  read_batch <- function(f) {
    side <- sub("\\.rds$", ".md5", f)
    ok <- file.exists(f) && identical(file_md5(f), readLines(side, warn = FALSE)[1])
    if (!ok) {
      warning(sprintf("corrupt or missing file %s ignored", basename(f)), call. = FALSE)
      return(NULL)
    }
    readRDS(f)
  }
  keep <- function(df) {
    if (is.null(df) || !nrow(df)) return(NULL)
    if (is.null(df$round)) df$round <- 0L
    df <- df[df$scenario_key %in% keys & df$method %in% names(fps), , drop = FALSE]
    df <- df[df$replication <= R[df$scenario_key], , drop = FALSE]
    df[df$method_fp == fps[df$method], , drop = FALSE]
  }
  task_id <- function(df) paste(df$scenario_key, df$replication, df$method, sep = "\r")
  # pass 1: the latest retry round of every task and the per data set truths
  last <- new.env(hash = TRUE, parent = emptyenv())
  dtruth <- list()
  for (f in files) {
    o <- read_batch(f)
    e <- keep(o$estimates)
    if (!is.null(e)) {
      ids <- task_id(e)
      for (i in unique(ids)) {
        r <- max(e$round[ids == i])
        if (is.null(last[[i]]) || r > last[[i]]) assign(i, r, envir = last)
      }
    }
    if (!is.null(o$dataset_truth) && nrow(o$dataset_truth)) dtruth[[length(dtruth) + 1L]] <- o$dataset_truth
  }
  dtruth <- bind_rows(dtruth)
  truth <- sim_truth(run)
  truth <- truth[truth$truth_status == "ok", c("scenario_key", "estimand", "truth"), drop = FALSE]
  ds_est <- names(study$estimands)[vapply(study$estimands, function(e) e$truth$kind == "dataset", logical(1))]
  # pass 2: accumulate usable rows, once per task and estimand
  seen <- new.env(hash = TRUE, parent = emptyenv())
  acc <- new.env(hash = TRUE, parent = emptyenv())
  add <- function(cell, what, x) {
    k <- paste(cell, what, sep = "\r")
    a <- acc[[k]] %||% sim_accumulator()
    assign(k, sim_accumulate(a, x), envir = acc)
  }
  for (f in files) {
    o <- read_batch(f)
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
    }
    e$truth_kind <- vapply(e$estimand, function(id) study$estimands[[id]]$truth$kind, character(1))
    e <- classify_rows(e, study$estimands, study$failure)
    cell <- paste(e$scenario_key, e$method, e$estimand, sep = "\r")
    ok_pt <- e$valid & is.finite(e$estimate)
    ok_tr <- ok_pt & !is.na(e$truth)
    ok_ci <- e$interval_status == "valid" & e$estimate_status %in% USABLE_FOR_INTERVAL & !is.na(e$truth)
    ok_se <- ok_pt & e$se_status == "valid"
    for (cl in unique(cell)) {
      i <- cell == cl
      if (any(c("bias", "mse") %in% measures) && any(i & ok_tr)) {
        err <- e$estimate[i & ok_tr] - e$truth[i & ok_tr]
        add(cl, "bias", err)
        add(cl, "mse", err^2)
      }
      if ("empse" %in% measures && any(i & ok_pt)) add(cl, "empse", e$estimate[i & ok_pt])
      if ("mean_se" %in% measures && any(i & ok_se)) add(cl, "mean_se", e$se[i & ok_se])
      if ("coverage" %in% measures && any(i & ok_ci)) {
        j <- i & ok_ci
        add(cl, "coverage", as.numeric(e$lower[j] <= e$truth[j] & e$truth[j] <= e$upper[j]))
      }
    }
  }
  out <- list()
  for (k in ls(acc)) {
    parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    m <- parts[4]
    if (!m %in% measures) next
    s <- sim_accumulated(acc[[k]])
    n <- s$n
    val <- if (m == "empse") sqrt(s$var_x) else s$mean_x
    mc <- switch(m, empse = if (n > 1) val / sqrt(2 * (n - 1)) else NA_real_,
                 coverage = sqrt(val * (1 - val) / n), s$mcse_mean_x)
    out[[length(out) + 1L]] <- data.frame(scenario_key = parts[1],
      scenario_id = study$scenarios$labels[match(parts[1], keys)], method = parts[2],
      estimand = parts[3], measure = m, value = val, mcse = mc, n = as.integer(n),
      stringsAsFactors = FALSE)
  }
  out <- bind_rows(out)
  if (is.null(out)) return(out)
  out <- out[order(match(out$scenario_key, keys), out$method, out$estimand,
                   match(out$measure, measures)), , drop = FALSE]
  rownames(out) <- NULL
  out
}
