# Run directory storage ----------------------------------------------------
#
# Layout of a run directory:
#   manifest.rds        study, seed, RNG registry, fingerprints, revisions
#   manifest.rds.md5    checksum of the manifest
#   manifest.rds.bak    previous manifest, read when the checksum fails
#   manifest.json       human readable summary of the manifest
#   LOCK/owner.json     coordinator lock (one active coordinator per run)
#   batches/<id>.rds    committed results of a batch of tasks
#   batches/<id>.md5    checksum sidecar; its presence is the commit point
#   batches/INVENTORY   names of all committed batches (detects lost batches)
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
  # unique temporary names: two writers never rename each other's files
  stem <- basename(tempfile(paste0(name, "-")))
  tmp_data <- file.path(tmpdir, paste0(stem, ".rds"))
  tmp_side <- file.path(tmpdir, paste0(stem, ".md5"))
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

# Publish an object under `name` unless one is already committed there, and
# return the committed object. Writers claim the name with an exclusive
# directory, so two concurrent writers never replace each other's value and
# every reader reports the first value committed. A claim left by a dead
# process, or held past `wait` seconds, is removed.
# ponytail: the claim is held only for one atomic_save, so the timeout is
# a crash guard, not a lease; add a lease if writes ever become slow.
publish_first <- function(obj, dir, name, tmpdir, wait = 30) {
  f <- file.path(dir, paste0(name, ".rds"))
  side <- file.path(dir, paste0(name, ".md5"))
  committed <- function() file.exists(f) && file.exists(side) &&
    identical(file_md5(f), readLines(side, warn = FALSE)[1])
  claim <- file.path(dir, paste0(name, ".claim"))
  deadline <- Sys.time() + wait
  repeat {
    if (committed()) return(readRDS(f))
    if (dir.create(claim, showWarnings = FALSE)) {
      on.exit(unlink(claim, recursive = TRUE), add = TRUE)
      writeLines(as.character(Sys.getpid()), file.path(claim, "pid"))
      if (committed()) return(readRDS(f))
      atomic_save(obj, dir, name, tmpdir)
      return(obj)
    }
    pid <- suppressWarnings(as.integer(tryCatch(
      readLines(file.path(claim, "pid"), warn = FALSE)[1], error = function(e) NA)))
    if ((!is.na(pid) && identical(pid_alive(pid), FALSE)) || Sys.time() > deadline) {
      unlink(claim, recursive = TRUE)
      deadline <- Sys.time() + wait
    }
    Sys.sleep(0.05)
  }
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
  inv <- file.path(dir, "INVENTORY")
  listed <- if (file.exists(inv)) unique(readLines(inv, warn = FALSE)) else character()
  listed <- listed[nzchar(listed)]
  present <- sub("\\.rds$", "", names(out))
  gone <- setdiff(listed, sub("\\.md5$", "", basename(side)))
  for (g in gone) {
    problems <- c(problems, sprintf("committed batch %s is missing from the run directory; %s",
                                    g, if (repair) "its tasks will be rerun" else
                                      "resume the run to rerun its tasks"))
  }
  if (repair) {
    # a batch committed by a process that stopped before appending it to the
    # inventory is added here, so its later loss is detected too
    keep <- union(intersect(listed, present), present)
    if (!identical(keep, listed)) {
      tmp <- paste0(inv, ".tmp")
      writeLines(keep, tmp)
      file.rename(tmp, inv)
    }
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

# Two coordinators can both see the same dead owner. A takeover is done
# only while holding a second lock, LOCK.takeover, created with an atomic
# dir.create(): one contender at a time checks, under that lock, that the
# run lock still belongs to the dead owner, removes it and creates it again.
# A coordinator that is not taking over cannot replace an existing run lock
# (its dir.create() fails), so no live lock is ever moved or removed.
take_over_stale_lock <- function(p, owner) {
  mutex <- paste0(p$lock, ".takeover")
  host <- Sys.info()[["nodename"]]
  if (!dir.create(mutex, showWarnings = FALSE)) {
    # a takeover lock left by a contender that died is itself stale
    mo <- tryCatch(jsonlite::fromJSON(file.path(mutex, "owner.json")), error = function(e) NULL)
    if (is.null(mo) || !identical(mo$host, host) || !identical(pid_alive(mo$pid), FALSE)) {
      return(FALSE)
    }
    unlink(mutex, recursive = TRUE)
    if (!dir.create(mutex, showWarnings = FALSE)) return(FALSE)
  }
  on.exit(unlink(mutex, recursive = TRUE), add = TRUE)
  jsonlite::write_json(list(pid = Sys.getpid(), host = host), file.path(mutex, "owner.json"),
                       auto_unbox = TRUE)
  current <- lock_owner(p)
  if (is.null(current) || !identical(current[c("pid", "host", "time")], owner[c("pid", "host", "time")]) ||
      !identical(pid_alive(current$pid), FALSE)) {
    return(FALSE)
  }
  unlink(p$lock, recursive = TRUE)
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

# Functions are exported with only the bindings they need. A function
# whose enclosing environment is not a package namespace is copied into a
# new environment that holds exactly the variables its body and formals
# name, found where R would find them (enclosing environments and the
# global environment), with functions among them pruned the same way. The
# exported study therefore runs in a clean process and carries no other
# objects of the session. A function that looks variables up by name (get,
# eval, environment and similar) cannot be pruned faithfully; it keeps its
# whole enclosing environment and is reported.
DYNAMIC_LOOKUP <- c("get", "get0", "mget", "exists", "environment", "parent.frame",
                    "sys.function", "sys.call", "eval", "evalq", "ls", "as.list")

export_functions <- function(x, memo = NULL) {
  if (is.null(memo)) {
    memo <- new.env(parent = emptyenv())
    memo$pairs <- list(); memo$carried <- character(); memo$whole <- character()
    memo$unshipped <- character()
    out <- export_functions(x, memo)
    attr(out, "carried") <- unique(memo$carried)
    attr(out, "whole") <- unique(memo$whole)
    attr(out, "unshipped") <- unique(memo$unshipped)
    return(out)
  }
  if (is.function(x)) return(export_function(x, memo))
  if (is.list(x) && !is.data.frame(x) && length(x)) {
    a <- attributes(x)
    x[] <- lapply(x, export_functions, memo = memo)
    attributes(x) <- a
  }
  x
}

# Names a piece of code can look up: every symbol, except the field names
# after `$` and `@`, which are not variables.
code_names <- function(x) {
  if (is.symbol(x)) return(if (nzchar(as.character(x))) as.character(x) else character())
  if (is.call(x)) {
    fn <- x[[1L]]
    if (is.symbol(fn) && as.character(fn) %in% c("$", "@") && length(x) == 3L) {
      return(c(as.character(fn), code_names(x[[2L]])))
    }
    return(unlist(lapply(as.list(x), code_names)))
  }
  if (is.pairlist(x) || is.list(x)) return(unlist(lapply(as.list(x), code_names)))
  character()
}

export_function <- function(f, memo) {
  if (is.primitive(f)) return(f)
  for (pr in memo$pairs) if (identical(pr[[1]], f)) return(pr[[2]])
  e <- environment(f)
  if (is.null(e) || isNamespace(e) || identical(e, baseenv()) || identical(e, emptyenv())) return(f)
  names_used <- unique(c(code_names(body(f)), unlist(lapply(formals(f), code_names))))
  if (any(DYNAMIC_LOOKUP %in% names_used) && !identical(e, globalenv())) {
    memo$whole <- c(memo$whole, ls(e, all.names = TRUE))
    memo$carried <- c(memo$carried, ls(e, all.names = TRUE))
    # names that fall through to the global environment do not travel with
    # a whole environment
    for (v in names_used) {
      h <- e
      while (!identical(h, globalenv()) && !isNamespace(h) && !identical(h, emptyenv()) &&
             !identical(h, baseenv()) && !exists(v, envir = h, inherits = FALSE)) h <- parent.env(h)
      if (identical(h, globalenv()) && exists(v, envir = h, inherits = FALSE)) {
        memo$unshipped <- c(memo$unshipped, v)
      }
    }
    return(f)
  }
  # the new environment's parent is the first namespace on the chain (a
  # function defined in a package's environment), else the global
  # environment, so package functions resolve as before
  top <- e
  while (!identical(top, globalenv()) && !isNamespace(top) &&
         !identical(top, emptyenv()) && !identical(top, baseenv())) top <- parent.env(top)
  ne <- new.env(parent = if (isNamespace(top)) top else globalenv())
  g <- utils::removeSource(f)
  environment(g) <- ne
  memo$pairs[[length(memo$pairs) + 1L]] <- list(f, g)
  for (v in names_used) {
    h <- e
    repeat {
      if (isNamespace(h) || identical(h, emptyenv()) || identical(h, baseenv())) { h <- NULL; break }
      if (exists(v, envir = h, inherits = FALSE)) break
      if (identical(h, globalenv())) { h <- NULL; break }
      h <- parent.env(h)
    }
    if (is.null(h)) next
    val <- get(v, envir = h, inherits = FALSE)
    if (is.function(val)) val <- export_function(val, memo)
    assign(v, val, envir = ne)
    memo$carried <- c(memo$carried, v)
  }
  g
}

# Environment fields compared across sessions (not the time).
environment_signature <- function(env) {
  env[c("r_version", "platform", "os", "packages", "blas", "lapack")]
}

# The manifest is written to a temporary file and renamed into place, with
# an md5 sidecar; the previous manifest is kept as manifest.rds.bak. This
# survives process death at any point, not power loss (nothing is flushed).
write_manifest <- function(p, manifest) {
  if (!dir.exists(p$tmp)) dir.create(p$tmp, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile("manifest-", tmpdir = p$tmp, fileext = ".rds")
  saveRDS(manifest, tmp)
  md5 <- file_md5(tmp)
  tmp_md5 <- sub("\\.rds$", ".md5", tmp)
  writeLines(md5, tmp_md5)
  side <- paste0(p$manifest, ".md5")
  # The backup is the previous manifest only when that one verifies: a
  # damaged primary never replaces a good backup. The backup carries its
  # own checksum.
  if (manifest_verifies(p$manifest)) {
    ok <- file.copy(p$manifest, paste0(p$manifest, ".bak"), overwrite = TRUE) &&
      file.copy(side, paste0(p$manifest, ".bak.md5"), overwrite = TRUE)
    if (!isTRUE(ok)) stopf("Could not back up the manifest of %s.", p$root)
  }
  if (!file.rename(tmp, p$manifest)) {
    unlink(c(tmp, tmp_md5))
    stopf("Could not commit the manifest of %s.", p$root)
  }
  if (!file.rename(tmp_md5, side)) stopf("Could not commit the manifest checksum of %s.", p$root)
  js <- list(
    schema_version = manifest$schema_version, study = manifest$study$name,
    label = manifest$label, seed = manifest$seed,
    created = manifest$created, updated = now_utc(),
    replications = manifest$replications,
    fingerprints = manifest$fingerprints,
    rng_registry = manifest$rng_registry,
    revisions = manifest$revisions, environment = manifest$environment,
    sequential = manifest$sequential
  )
  tmpj <- tempfile("manifest-", tmpdir = p$tmp, fileext = ".json")
  jsonlite::write_json(js, tmpj, auto_unbox = TRUE, pretty = TRUE,
                       digits = NA, null = "null", na = "null")
  if (!file.rename(tmpj, p$manifest_json)) {
    unlink(tmpj)
    stopf("Could not write %s.", p$manifest_json)
  }
  invisible(TRUE)
}

# A manifest file that exists and matches its checksum sidecar (a file
# without a sidecar is accepted, as written by older versions).
manifest_verifies <- function(f) {
  if (!file.exists(f)) return(FALSE)
  side <- paste0(f, ".md5")
  !file.exists(side) || identical(readLines(side, warn = FALSE)[1], file_md5(f))
}

# Change a stored manifest under the run lock.
update_manifest <- function(path, f) {
  p <- run_paths(path)
  token <- acquire_lock(p)
  on.exit(release_lock(p, token), add = TRUE)
  write_manifest(p, f(read_manifest(path)))
  invisible(TRUE)
}

read_manifest <- function(path) {
  p <- run_paths(path)
  if (!file.exists(p$manifest)) stopf("No adempR run found at %s.", path)
  m <- if (manifest_verifies(p$manifest)) tryCatch(readRDS(p$manifest), error = function(e) NULL)
  if (is.null(m)) {
    bak <- paste0(p$manifest, ".bak")
    m <- if (manifest_verifies(bak)) tryCatch(readRDS(bak), error = function(e) NULL)
    if (is.null(m)) stopf("The manifest of %s is damaged and has no readable backup.", path)
    warnf("The manifest of %s failed its checksum; using the previous manifest (manifest.rds.bak). Results committed after it are still read from their batches.", path)
  }
  if (!identical(m$schema_version, SCHEMA_VERSION)) {
    stopf("The run at %s uses result schema %s; this version of adempR reads schema %s.",
          path, m$schema_version, SCHEMA_VERSION)
  }
  m
}
