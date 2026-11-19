# Run directory storage ----------------------------------------------------
#
# Layout of a run directory:
#   manifest.rds        study, seed, RNG registry, fingerprints, revisions
#   manifest.json       human readable summary of the manifest
#   LOCK/owner.json     coordinator lock (one active coordinator per run)
#   batches/<id>.rds    committed results of a batch of tasks
#   batches/<id>.md5    checksum sidecar; its presence is the commit point
#   truth/              computed truth values
#   tmp/                files being written; never read as results
#   quarantine/         corrupt or uncommitted files moved aside
#
# Commit protocol: the batch is serialized to tmp/, its md5 is written to a
# temporary sidecar, the data file is renamed into batches/, and finally the
# sidecar is renamed into batches/. Renames within one directory tree on a
# local file system are atomic on POSIX systems and on NTFS. A data file
# without a sidecar is uncommitted; a sidecar whose checksum does not match
# marks a corrupt file. Network and object stores are not covered by this
# contract.

SCHEMA_VERSION <- "2.0.0"

fault_point <- function(name) {
  if (identical(getOption("adempR.fault"), name)) {
    if (.Platform$OS.type == "unix") {
      tools::pskill(Sys.getpid(), tools::SIGKILL)
      Sys.sleep(5)
    }
    stop("adempR fault injection: ", name, call. = FALSE)
  }
  invisible(NULL)
}

run_paths <- function(path) {
  list(root = path,
       manifest = file.path(path, "manifest.rds"),
       manifest_json = file.path(path, "manifest.json"),
       lock = file.path(path, "LOCK"),
       batches = file.path(path, "batches"),
       truth = file.path(path, "truth"),
       tmp = file.path(path, "tmp"),
       quarantine = file.path(path, "quarantine"))
}

ensure_dirs <- function(p) {
  for (d in c(p$root, p$batches, p$truth, p$tmp, p$quarantine)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

file_md5 <- function(f) unname(tools::md5sum(f))

# Atomically write an R object with a checksum sidecar.
atomic_save <- function(obj, dir, name, tmpdir) {
  tmp_data <- file.path(tmpdir, paste0(name, ".rds"))
  tmp_side <- file.path(tmpdir, paste0(name, ".md5"))
  saveRDS(obj, tmp_data)
  fault_point("after_data_write")
  writeLines(file_md5(tmp_data), tmp_side)
  final_data <- file.path(dir, paste0(name, ".rds"))
  final_side <- file.path(dir, paste0(name, ".md5"))
  if (!file.rename(tmp_data, final_data)) stopf("Could not commit %s.", final_data)
  fault_point("before_sidecar")
  if (!file.rename(tmp_side, final_side)) stopf("Could not commit %s.", final_side)
  fault_point("after_commit")
  invisible(final_data)
}

# Read committed objects. Corrupt and uncommitted files are skipped and
# reported; with `repair = TRUE` (only under the coordinator lock) they are
# moved to the quarantine directory so their tasks are rerun.
read_committed <- function(dir, quarantine, repair = FALSE) {
  side <- list.files(dir, pattern = "\\.md5$", full.names = TRUE)
  data_files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  committed <- sub("\\.md5$", ".rds", side)
  problems <- character()
  orphan <- setdiff(data_files, committed)
  for (f in orphan) {
    if (repair) file.rename(f, file.path(quarantine, basename(f)))
    problems <- c(problems, sprintf("uncommitted file %s %s", basename(f),
                                    if (repair) "moved to quarantine" else "ignored"))
  }
  out <- list()
  for (i in seq_along(side)) {
    f <- committed[i]
    expected <- tryCatch(readLines(side[i], warn = FALSE)[1], error = function(e) NA)
    ok <- file.exists(f) && !is.na(expected) && identical(file_md5(f), expected)
    obj <- if (ok) tryCatch(readRDS(f), error = function(e) NULL) else NULL
    if (is.null(obj)) {
      if (repair) {
        if (file.exists(f)) file.rename(f, file.path(quarantine, basename(f)))
        file.rename(side[i], file.path(quarantine, basename(side[i])))
      }
      problems <- c(problems, sprintf("corrupt or missing file %s %s", basename(f),
                                      if (repair) "moved to quarantine; its tasks will be rerun" else
                                        "ignored; resume the run to repair it"))
      next
    }
    out[[basename(f)]] <- obj
  }
  if (length(problems)) warning(paste(problems, collapse = "\n"), call. = FALSE)
  list(objects = out, problems = problems)
}

# Locking ----------------------------------------------------------------

pid_alive <- function(pid) {
  if (.Platform$OS.type != "unix") return(NA)
  isTRUE(tools::pskill(pid, 0L))
}

lock_owner <- function(p) {
  f <- file.path(p$lock, "owner.json")
  if (!file.exists(f)) return(NULL)
  tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
}

acquire_lock <- function(p) {
  host <- Sys.info()[["nodename"]]
  # A random token identifies this coordinator: process ids repeat across
  # containers that share a file system and a host name.
  # tempfile() names do not touch the R random number stream of the caller.
  token <- paste0(Sys.getpid(), "-", basename(tempfile("")), "-", format(Sys.time(), "%OS6"))
  if (!dir.create(p$lock, showWarnings = FALSE)) {
    owner <- lock_owner(p)
    stale <- !is.null(owner) && identical(owner$host, host) &&
      identical(pid_alive(owner$pid), FALSE)
    if (!stale || !take_over_stale_lock(p, owner)) {
      stopf("The run directory is locked by %s. Only one coordinator may run at a time. If that process is no longer running, call sim_unlock(\"%s\").",
            if (is.null(owner)) "an unknown process" else
              sprintf("process %s on %s since %s", owner$pid, owner$host, owner$time),
            p$root)
    }
    message(sprintf("Recovered a stale lock left by process %s (no longer running).",
                    owner$pid))
  }
  jsonlite::write_json(list(pid = Sys.getpid(), host = host, time = now_utc(), token = token),
                       file.path(p$lock, "owner.json"), auto_unbox = TRUE)
  token
}

# Two coordinators can both see the same dead owner. Each moves the lock
# aside under a unique name (an atomic rename) and checks that what it moved
# is the dead owner's lock; one that moved a fresh lock puts it back. Only
# then is the lock created again, and dir.create() admits one winner.
take_over_stale_lock <- function(p, owner) {
  aside <- paste0(p$lock, ".stale-", Sys.getpid(), "-", basename(tempfile("")))
  if (!file.rename(p$lock, aside)) return(FALSE)
  moved <- tryCatch(jsonlite::fromJSON(file.path(aside, "owner.json")), error = function(e) NULL)
  if (!identical(moved[c("pid", "host", "time")], owner[c("pid", "host", "time")])) {
    # ponytail: if a third process created a lock in between, the moved
    # live lock cannot be restored; its owner then fails at release, not here.
    file.rename(aside, p$lock)
    return(FALSE)
  }
  unlink(aside, recursive = TRUE)
  dir.create(p$lock, showWarnings = FALSE)
}

release_lock <- function(p, token) {
  owner <- lock_owner(p)
  if (!is.null(owner) && identical(owner$token, token)) {
    unlink(p$lock, recursive = TRUE)
  }
  invisible(TRUE)
}

#' Remove a stale run lock
#'
#' A run directory admits one coordinator at a time. On Unix, a lock left
#' by a process that is no longer running on the same host is recovered
#' automatically. A lock held on another host, or any lock on Windows
#' (where adempR cannot check whether a process is alive), must be removed
#' with `sim_unlock()` after making sure no coordinator is running.
#' `sim_unlock()` refuses to remove a lock whose owner is a running process
#' on this host unless `force = TRUE`.
#'
#' @param path Run directory.
#' @param force Remove the lock even when its owner is running on this host.
#' @return `TRUE`, invisibly.
#' @examples
#' path <- tempfile("adempR-run-")
#' dir.create(file.path(path, "LOCK"), recursive = TRUE)
#' sim_unlock(path)
#' @export
sim_unlock <- function(path, force = FALSE) {
  lk <- file.path(path, "LOCK")
  if (!dir.exists(lk)) return(invisible(TRUE))
  owner <- lock_owner(list(lock = lk))
  if (!force && !is.null(owner) && identical(owner$host, Sys.info()[["nodename"]]) &&
      isTRUE(pid_alive(owner$pid))) {
    stopf("The lock is held by process %s on this host, which is still running. Stop it first, or use `force = TRUE`.",
          owner$pid)
  }
  unlink(lk, recursive = TRUE)
  invisible(TRUE)
}

# Environment and provenance ---------------------------------------------

environment_info <- function(packages = character()) {
  pk <- unique(c("adempR", packages))
  vers <- vapply(pk, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  si <- Sys.info()
  list(
    r_version = R.version.string,
    platform = R.version$platform,
    os = paste(si[["sysname"]], si[["release"]]),
    machine = si[["machine"]],
    rng_kind = RNG_KIND,
    packages = as.list(vers),
    # File names only: a library path can contain the user name.
    blas = tryCatch(basename(extSoftVersion()[["BLAS"]]), error = function(e) NA_character_),
    lapack = tryCatch(basename(La_library()), error = function(e) NA_character_),
    threads = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                   "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"),
                                 unset = NA)),
    git = git_info(),
    time = now_utc()
  )
}

# Version control state of the working directory, read with the git command
# line tool when it is installed. NULL outside a repository or without git.
git_info <- function(dir = getOption("adempR.git_dir", getwd())) {
  git <- Sys.which("git")
  if (!nzchar(git)) return(NULL)
  run <- function(...) {
    out <- tryCatch(suppressWarnings(system2(git, c("-C", shQuote(dir), ...),
                                             stdout = TRUE, stderr = FALSE)),
                    error = function(e) character())
    if (!is.null(attr(out, "status")) && attr(out, "status") != 0) character() else out
  }
  sha <- run("rev-parse", "HEAD")
  if (length(sha) != 1L) return(NULL)
  status <- run("status", "--porcelain", "--untracked-files=no")
  remote <- run("remote", "get-url", "origin")
  list(commit = sha,
       branch = run("rev-parse", "--abbrev-ref", "HEAD")[1] %||% NA_character_,
       dirty = length(status) > 0L,
       remote = if (length(remote)) strip_credentials(remote[1]) else NA_character_)
}

# Never store a token or password embedded in a remote URL.
strip_credentials <- function(url) sub("^([a-z+]+://)[^@/]*@", "\\1", url)

# Source references of stored functions hold the absolute path of the file
# they were sourced from and the working directory. They are removed before
# a study is written to a manifest or bundle; fingerprints do not use them.
strip_sources <- function(x) {
  if (is.function(x)) {
    return(if (is.primitive(x)) x else utils::removeSource(x))
  }
  if (is.list(x) && !is.data.frame(x) && length(x)) x[] <- lapply(x, strip_sources)
  x
}

# Variables carried by the enclosing environments of the functions in a
# study, other than the global environment and package namespaces. A
# serialized closure includes them, so an export lists them.
closure_contents <- function(x, out = character()) {
  if (is.function(x) && !is.primitive(x)) {
    e <- environment(x)
    if (!is.null(e) && !identical(e, globalenv()) && !isNamespace(e) &&
        !identical(e, baseenv()) && !identical(e, emptyenv())) {
      out <- c(out, ls(e, all.names = TRUE))
    }
    return(unique(out))
  }
  if (is.list(x) && !is.data.frame(x)) for (el in x) out <- closure_contents(el, out)
  unique(out)
}

# Environment fields compared across sessions (not the time).
environment_signature <- function(env) {
  env[c("r_version", "platform", "os", "packages", "blas", "lapack")]
}

write_manifest <- function(p, manifest) {
  tmp <- file.path(p$tmp, "manifest.rds")
  saveRDS(manifest, tmp)
  file.rename(tmp, p$manifest)
  js <- list(
    schema_version = manifest$schema_version, study = manifest$study$name,
    label = manifest$label, seed = manifest$seed,
    created = manifest$created, updated = now_utc(),
    replications = manifest$replications,
    fingerprints = manifest$fingerprints,
    rng_registry = manifest$rng_registry,
    revisions = manifest$revisions, environment = manifest$environment
  )
  jsonlite::write_json(js, p$manifest_json, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null", na = "null")
}

read_manifest <- function(path) {
  p <- run_paths(path)
  if (!file.exists(p$manifest)) stopf("No adempR run found at %s.", path)
  m <- readRDS(p$manifest)
  if (!identical(m$schema_version, SCHEMA_VERSION)) {
    stopf("The run at %s uses result schema %s; this version of adempR reads schema %s.",
          path, m$schema_version, SCHEMA_VERSION)
  }
  m
}
