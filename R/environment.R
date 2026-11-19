# Software environment: package sources and lockfiles ------------------------

# Packages whose code a study runs: the declared ones and those named with
# `::` or `:::` in its functions.
env_packages <- function(study) {
  unique(c(study_packages(study), tryCatch(referenced_packages(study), error = function(e) character())))
}

# Base and recommended packages ship with R; their version is the R version.
BASE_PACKAGES <- c("base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods",
                   "parallel", "splines", "stats", "stats4", "tcltk", "tools", "utils")

package_cache <- new.env(parent = emptyenv())

# Name, version and source of every installed package in the dependency
# closure (Depends, Imports, LinkingTo) of `pk`. Only DESCRIPTION fields are
# read, so nothing is loaded and no network is used. The result is cached
# per session for a given set of packages and library paths.
package_records <- function(pk) {
  pk <- sort(unique(setdiff(c("adempR", pk), BASE_PACKAGES)))
  key <- paste(c(pk, .libPaths()), collapse = "\r")
  if (!is.null(package_cache[[key]])) return(package_cache[[key]])
  fields <- c("Package", "Version", "Depends", "Imports", "LinkingTo", "Repository", "RemoteType",
              "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteRef", "RemoteSha", "RemoteUrl")
  seen <- list()
  todo <- pk
  while (length(todo)) {
    p <- todo[1]
    todo <- todo[-1]
    if (!is.null(seen[[p]]) || p %in% BASE_PACKAGES) next
    d <- suppressWarnings(tryCatch(utils::packageDescription(p, fields = fields),
                                   error = function(e) NULL))
    if (is.null(d) || !is.list(d)) {
      seen[[p]] <- list(Package = p, Version = NA_character_, Source = "not installed")
      next
    }
    d <- d[!vapply(d, function(v) length(v) == 0L || is.na(v), logical(1))]
    deps <- unlist(lapply(d[intersect(c("Depends", "Imports", "LinkingTo"), names(d))], function(v) {
      v <- trimws(sub("\\(.*$", "", strsplit(gsub("\n", " ", v), ",")[[1]]))
      v[nzchar(v) & v != "R"]
    }))
    src <- if (!is.null(d$RemoteType) && !identical(d$RemoteType, "standard")) {
      switch(tolower(d$RemoteType), github = "GitHub", gitlab = "GitLab", bitbucket = "Bitbucket",
             git2r = , git = "git", url = "URL", local = "Local", d$RemoteType)
    } else if (!is.null(d$Repository)) "Repository" else "unknown"
    rec <- c(list(Package = p, Version = d$Version, Source = src),
             d[intersect(c("Repository", "RemoteType", "RemoteHost", "RemoteUsername",
                           "RemoteRepo", "RemoteRef", "RemoteSha", "RemoteUrl"), names(d))])
    # a local path can contain the user name; keep its last component only
    if (!is.null(rec$RemoteUrl) && identical(src, "Local")) rec$RemoteUrl <- basename(rec$RemoteUrl)
    seen[[p]] <- rec
    todo <- c(todo, setdiff(deps, names(seen)))
  }
  out <- seen[order(names(seen))]
  assign(key, out, envir = package_cache)
  out
}

# Repository URLs in effect, with the CRAN placeholder resolved.
repository_urls <- function() {
  r <- getOption("repos")
  if (is.null(r) || !length(r)) return(list(CRAN = "https://cloud.r-project.org"))
  r[r == "@CRAN@"] <- "https://cloud.r-project.org"
  if (is.null(names(r)) || any(!nzchar(names(r)))) names(r) <- paste0("repo", seq_along(r))
  as.list(r)
}

# Whether an renv project is active, and the checksum of its lockfile (not
# its path, which can contain the user name).
renv_status <- function() {
  proj <- Sys.getenv("RENV_PROJECT", unset = "")
  lock <- if (nzchar(proj)) file.path(proj, "renv.lock") else ""
  list(active = nzchar(proj),
       lockfile_md5 = if (nzchar(lock) && file.exists(lock)) unname(tools::md5sum(lock)) else NA_character_)
}

#' Lockfile and check of the software environment of a run
#'
#' Every session that works on a run records its R version, platform,
#' numerical libraries, thread settings, collation locale and time zone,
#' whether an renv project was active (with the checksum of its lockfile),
#' and the version and source (repository, GitHub commit, local build) of
#' every package in the dependency closure of adempR and of the packages
#' the study uses. `sim_lockfile()` writes that record as a lockfile in the
#' format of the renv package, so the package versions can be restored
#' with `renv::restore(lockfile = )` before regenerating results.
#' `sim_check_environment()` compares the packages installed now with a
#' run or a lockfile and reports every difference. [sim_export()] puts the
#' lockfile in the bundle, and its `regenerate.R` runs this check first.
#'
#' A lockfile pins R packages only. The R version, system libraries (BLAS,
#' LAPACK, compilers) and thread settings also affect numerical results;
#' they are recorded in the manifest, and a container image (for example
#' from the Rocker project) pins them.
#'
#' @param run A run, results, a run directory, or (for
#'   `sim_check_environment()`) the path of a lockfile.
#' @param file Path of the lockfile to write.
#' @param session `"created"` (the session that created the run) or
#'   `"latest"` (the last session that worked on it).
#' @return `sim_lockfile()`: the path, invisibly. `sim_check_environment()`:
#'   a data frame with one row per recorded package (`package`, `recorded`,
#'   `installed`, `source`, `match`), with attribute `r_version` giving the
#'   recorded and current R versions; a warning lists differences.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 2, seed = 1, path = tempfile("adempR-run-"))
#' lock <- sim_lockfile(run, tempfile(fileext = ".lock"))
#' head(sim_check_environment(lock))
#' @export
sim_lockfile <- function(run, file, session = c("created", "latest")) {
  session <- match.arg(session)
  env <- run_environment(run, session)
  write_lockfile(env, file)
  invisible(file)
}

#' @rdname sim_lockfile
#' @export
sim_check_environment <- function(run) {
  if (is.character(run) && length(run) == 1L && file.exists(run) && !dir.exists(run)) {
    lock <- jsonlite::read_json(run)
    recs <- lock$Packages
    r_rec <- lock$R$Version
  } else {
    env <- run_environment(run, "created")
    recs <- env$package_sources
    if (is.null(recs)) recs <- lapply(names(env$packages), function(p) list(Package = p, Version = env$packages[[p]]))
    r_rec <- env$r_number %||% env$r_version
  }
  pk <- vapply(recs, function(r) r$Package, character(1))
  recorded <- vapply(recs, function(r) as.character(r$Version %||% NA_character_), character(1))
  installed <- vapply(pk, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))
  out <- data.frame(package = pk, recorded = recorded, installed = unname(installed),
                    source = vapply(recs, function(r) as.character(r$Source %||% NA_character_), character(1)),
                    stringsAsFactors = FALSE)
  # "0.2-20" and "0.2.20" are the same version
  same_version <- function(a, b) {
    ok <- !is.na(a) & !is.na(b)
    ok[ok] <- vapply(which(ok), function(i) {
      isTRUE(tryCatch(package_version(a[i]) == package_version(b[i]), error = function(e) a[i] == b[i]))
    }, logical(1))
    ok
  }
  out$match <- same_version(out$installed, out$recorded)
  rownames(out) <- NULL
  r_now <- paste(R.version$major, R.version$minor, sep = ".")
  attr(out, "r_version") <- c(recorded = as.character(r_rec), current = r_now)
  bad <- out[!out$match, , drop = FALSE]
  r_differs <- !grepl(r_now, as.character(r_rec), fixed = TRUE)
  if (nrow(bad) || r_differs) {
    warning(paste0("The software environment differs from the recorded one",
                   if (r_differs) sprintf(" (R %s recorded, %s now)", r_rec, r_now),
                   if (nrow(bad)) paste0(": ", paste(sprintf("%s %s recorded, %s installed", bad$package,
                                                            bad$recorded, ifelse(is.na(bad$installed), "none", bad$installed)),
                                                    collapse = "; ")),
                   ". Results may differ numerically; renv::restore() with the lockfile installs the recorded versions."),
            call. = FALSE)
  }
  out
}

run_environment <- function(run, session) {
  if (is.character(run)) run <- sim_open(run)
  man <- if (inherits(run, "adempr_results") || inherits(run, "adempr_run")) run$manifest else
    stopf("`run` must be a run, results or a run directory.")
  if (session == "latest" && length(man$sessions)) man$sessions[[length(man$sessions)]] else man$environment
}

write_lockfile <- function(env, file) {
  recs <- env$package_sources
  if (is.null(recs)) {
    recs <- stats::setNames(lapply(names(env$packages), function(p) {
      list(Package = p, Version = env$packages[[p]], Source = "unknown")
    }), names(env$packages))
  }
  recs <- recs[!vapply(recs, function(r) is.na(r$Version %||% NA), logical(1))]
  repos <- env$repositories %||% list(CRAN = "https://cloud.r-project.org")
  lock <- list(
    R = list(Version = env$r_number %||% sub("^R version ([0-9.]+).*$", "\\1", env$r_version),
             Repositories = unname(Map(function(n, u) list(Name = n, URL = u), names(repos), unlist(repos)))),
    Packages = recs)
  jsonlite::write_json(lock, file, auto_unbox = TRUE, pretty = TRUE)
  invisible(file)
}
