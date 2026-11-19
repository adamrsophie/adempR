# Reports and reproducibility bundles -------------------------------------

html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

html_table <- function(df, digits = 4) {
  if (is.null(df) || !NROW(df)) return("<p><em>None.</em></p>")
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  cells <- lapply(df, function(col) {
    if (is.double(col)) fmt_num(col, digits) else html_escape(col)
  })
  head <- paste0("<tr>", paste0("<th>", html_escape(names(df)), "</th>", collapse = ""), "</tr>")
  body <- vapply(seq_len(nrow(df)), function(i) {
    paste0("<tr>", paste0("<td>", vapply(cells, `[`, character(1), i), "</td>", collapse = ""),
           "</tr>")
  }, character(1))
  paste0("<div class=\"tw\"><table>", head, paste(body, collapse = ""), "</table></div>")
}

html_plot <- function(fun, width = 900, height = 500) {
  if (!isTRUE(capabilities("png"))) return("<p><em>PNG device unavailable.</em></p>")
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  ok <- tryCatch({
    grDevices::png(f, width = width, height = height, res = 110)
    fun()
    TRUE
  }, error = function(e) conditionMessage(e), finally = grDevices::dev.off())
  if (!isTRUE(ok)) return(paste0("<p><em>Plot unavailable: ", html_escape(ok), "</em></p>"))
  b64 <- jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
  paste0("<img alt=\"plot\" src=\"data:image/png;base64,", b64, "\">")
}

fun_source <- function(f) {
  if (is.null(f)) return("")
  paste(deparse(f, control = c("keepNA", "keepInteger")), collapse = "\n")
}

# Verified bibliographic metadata (checked against Crossref on 2026-11-19).
REFERENCES <- data.frame(
  key = c("morris2019", "white2024", "williams2024", "siepe2024", "pawel2024",
          "pawel2025", "naimi2025", "howard2021"),
  text = c(
    "Morris TP, White IR, Crowther MJ (2019). Using simulation studies to evaluate statistical methods. Statistics in Medicine 38(11):2074-2102.",
    "White IR, Pham TM, Quartagno M, Morris TP (2024). How to check a simulation study. International Journal of Epidemiology 53(1):dyad134.",
    "Williams C, Yang Y, Lagisz M, Morrison K, Ricolfi L, Warton DI, Nakagawa S (2024). Transparent reporting items for simulation studies evaluating statistical methods: foundations for reproducibility and reliability. Methods in Ecology and Evolution 15(11):1926-1939.",
    "Siepe BS, Bartos F, Morris TP, Boulesteix AL, Heck DW, Pawel S (2024). Simulation studies for methodological research in psychology: a standardized template for planning, preregistration, and reporting. Psychological Methods.",
    "Pawel S, Kook L, Reeve K (2024). Pitfalls and potentials in simulation studies: questionable research practices in comparative simulation studies allow for spurious claims of superiority of any method. Biometrical Journal 66(1):2200091.",
    "Pawel S, Bartos F, Siepe BS, Lohmann A (2026). Handling missingness, failures, and non-convergence in simulation studies: a review of current practices and recommendations. The American Statistician 80(1):31-48.",
    "Naimi AI, Benkeser D, Rudolph JE (2025). Computing true parameter values in simulation studies using Monte Carlo integration. Epidemiology 36(5):690-693.",
    "Howard SR, Ramdas A, McAuliffe J, Sekhon J (2021). Time-uniform, nonparametric, nonasymptotic confidence sequences. The Annals of Statistics 49(2):1055-1080."
  ),
  doi = c("10.1002/sim.8086", "10.1093/ije/dyad134", "10.1111/2041-210X.14415",
          "10.1037/met0000695", "10.1002/bimj.202200091", "10.1080/00031305.2025.2540002",
          "10.1097/EDE.0000000000001873", "10.1214/20-AOS1991"),
  stringsAsFactors = FALSE
)

#' Write an HTML report of a run
#'
#' Writes a self-contained HTML file (plots are embedded) documenting the
#' protocol and its completeness, amendments and run revisions, the actual
#' scenario design with exclusions, estimands and truth computation,
#' methods with their code and fingerprints, failure handling, performance
#' measure definitions and denominators, results with Monte Carlo
#' uncertainty, failures, diagnostics, the computing environment and the
#' random number scheme. No rendering software is needed. Text is a fixed
#' template filled from the results; the report draws no scientific
#' conclusions and never declares a best method.
#'
#' @param run A run or results.
#' @param file Output HTML file.
#' @param overwrite Overwrite an existing file.
#' @param plots Include plots.
#' @return The file path, invisibly.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = c(10L, 40L)),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) {
#'     sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)))
#'   }))
#' run <- sim_run(study, replications = 20, seed = 1, path = tempfile("adempR-run-"))
#' f <- sim_report(run, tempfile(fileext = ".html"), plots = FALSE)
#' file.exists(f)
#' @export
sim_report <- function(run, file, overwrite = FALSE, plots = TRUE) {
  if (!is_scalar_character(file)) stopf("`file` must be a path.")
  if (file.exists(file) && !overwrite) {
    stopf("%s exists; use `overwrite = TRUE` to replace it.", file)
  }
  res <- sim_results(run)
  man <- res$manifest
  study <- man$study
  summ <- sim_summarise(res)
  diag <- sim_diagnose(res)
  s <- as.data.frame(summ)
  h <- character()
  add <- function(...) h <<- c(h, paste0(...))
  add("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">",
      "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
      "<title>", html_escape(study$name), " simulation report</title><style>",
      "body{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:16px;",
      "color:#1d2733;background:#fff;line-height:1.45}",
      "table{border-collapse:collapse;font-size:13px}td,th{border:1px solid #d5dbe1;",
      "padding:3px 6px;text-align:left}th{background:#eef2f6}.tw{overflow-x:auto}",
      "h1,h2{color:#17324d}img{max-width:100%}code,pre{background:#f4f6f8}",
      "pre{padding:8px;overflow-x:auto}.note{color:#5a6773}</style></head><body>")
  add("<h1>Simulation study report: ", html_escape(study$name), "</h1>")
  add("<p class=\"note\">Generated ", now_utc(), " by adempR ",
      as.character(utils::packageVersion("adempR")), ". Run label: ", html_escape(man$label),
      ". Study fingerprint: <code>", html_escape(man$study_fingerprint), "</code>.</p>")
  frozen <- study$status$frozen
  add("<p>", if (is.null(frozen)) "The protocol was not frozen before this run." else
    paste0("The protocol was frozen on ", html_escape(frozen$time), " with hash <code>", html_escape(frozen$hash), "</code>",
           if (!is.null(frozen$git_commit)) paste0(" at git commit <code>", html_escape(frozen$git_commit), "</code>",
             if (isTRUE(frozen$git_dirty)) " (with uncommitted changes)" else "") else "",
           ". A frozen hash proves identity of the specification, not that it was fixed before results were seen; a commit or registration pushed to a public archive (for example a Zenodo or Software Heritage snapshot, or an OSF registration) fixes it in time."),
    " ", if (is.null(study$protocol$registration)) "No external registration is recorded." else
      paste0("External registration: ", html_escape(paste(unlist(study$protocol$registration), collapse = ", ")), "."),
    "</p>")
  add("<h2>1. Aims and protocol</h2><p>", html_escape(study$aims %||% "No aims stated."), "</p>")
  p <- study$protocol
  prot <- data.frame(field = c("purpose", "scope", "intended users", "DGM rationale",
                               "estimand rationale", "method rationale", "measure rationale",
                               "replication rationale", "failure handling", "retention", "reporting"),
                     value = c(p$purpose, p$scope %||% "", p$users %||% "", p$dgm_rationale %||% "",
                               p$estimand_rationale %||% "", p$method_rationale %||% "",
                               p$measure_rationale %||% "", p$replication_rationale %||% "",
                               p$failure_handling %||% "", p$retention %||% "", p$reporting %||% ""),
                     stringsAsFactors = FALSE)
  add(html_table(prot))
  add("<h3>Documentation completeness</h3>", html_table(sim_completeness(study)))
  add("<h2>2. Data-generating mechanisms</h2><p>",
      html_escape(study$generator$description %||% "No description."), " Design: ",
      html_escape(study$scenarios$design), ", ", length(study$scenarios), " scenario(s).</p>")
  sc <- res$scenarios
  sc[] <- lapply(sc, function(col) if (is.list(col)) vapply(col, function(v) {
    paste(deparse(v), collapse = " ")
  }, character(1)) else col)
  add(html_table(sc[setdiff(names(sc), "scenario_key")]))
  add("<h3>Scenarios excluded by design</h3>", html_table(study$scenarios$excluded))
  if (!is.null(study$scenarios$distribution)) {
    add("<p>Scenarios were sampled from: ", html_escape(study$scenarios$distribution$description),
        ". Averages over scenarios estimate performance averaged over this distribution.</p>")
  }
  add("<details><summary>Generator code</summary><pre>",
      html_escape(fun_source(study$generator$fun)), "</pre></details>")
  add("<h2>3. Estimands and truth</h2>")
  ed <- do.call(rbind, lapply(study$estimands, function(e) {
    data.frame(id = e$id, quantity = e$quantity, population = e$population %||% "",
               conditioning = e$conditioning %||% "", effect = e$effect, scale = e$scale,
               units = e$units %||% "", truth = e$truth$kind, stringsAsFactors = FALSE)
  }))
  add(html_table(ed))
  tr <- res$truth
  tr$scenario <- study$scenarios$labels[match(tr$scenario_key, study$scenarios$keys)]
  add(html_table(tr[c("scenario", "estimand", "truth_kind", "truth", "truth_se", "truth_error",
                      "truth_n", "truth_method")]))
  add("<h2>4. Methods</h2>")
  md <- do.call(rbind, lapply(names(study$methods), function(mn) {
    m <- study$methods[[mn]]
    data.frame(method = mn, estimands = paste(m$estimands, collapse = ", "),
               visible_scenario_fields = paste(m$scenario_fields, collapse = ", "),
               oracle = m$oracle, version = m$version %||% "",
               timeout = m$timeout %||% NA_real_, coupling = m$coupling %||% "",
               fingerprint = short_hash(res$manifest$fingerprints$methods[[mn]] %||% ""),
               description = m$description %||% "", stringsAsFactors = FALSE)
  }))
  add(html_table(md))
  for (mn in names(study$methods)) {
    add("<details><summary>Code of ", html_escape(mn), "</summary><pre>",
        html_escape(fun_source(study$methods[[mn]]$fun)), "</pre></details>")
  }
  add("<h2>5. Performance measures and failure handling</h2>")
  mi <- MEASURE_INFO[MEASURE_INFO$measure %in% unique(s$measure), , drop = FALSE]
  add(html_table(mi[c("measure", "eligibility", "direction", "mcse_method")]))
  f <- lapply(unclass(study$failure), html_escape)
  add("<p>Failure policy: nonconverged estimates are ", f$nonconverged,
      "; outputs with warnings are ", if (f$warnings == "record") "kept (warnings recorded)" else "invalid",
      "; boundary estimates are ", if (f$boundary == "record") "kept (flag recorded)" else "invalid",
      "; infinite estimates are ", if (f$infinite == "include") "kept as extended real values" else "invalid",
      "; infrastructure retries: ", f$retries, ". Data sets are never redrawn after a method fails.",
      " Performance is reported conditional on each measure's eligibility rule, next to the probability of a usable result.</p>")
  add("<h2>6. Execution</h2><p>Replications planned per scenario: ", html_escape(man$replications),
      ". Base seed: ", html_escape(man$seed), ". Random numbers: L'Ecuyer-CMRG streams (normal kind Inversion, sample kind Rejection) from an append-only registry of ",
      nrow(man$rng_registry), " streams keyed by stage, scenario and component; replication r uses substream r.</p>")
  sq <- man$sequential
  if (!is.null(sq)) {
    lk <- sq$looks
    last <- lk[lk$look == max(lk$look), c("scenario_id", "method", "estimand", "measure", "value",
                                           "n_eligible", "lower", "upper", "half_width"),
               drop = FALSE]
    add("<h3>Sequential stopping</h3><p>", html_escape(sequential_note(sq)),
        " Looks: ", html_escape(paste(sq$design$looks, collapse = ", ")),
        "; target half width ", html_escape(sq$design$half_width), "; familywise alpha ",
        html_escape(sq$design$alpha), " over ", sq$design$cells, " cells.</p>",
        html_table(last))
  }
  env <- man$environment
  g <- env$git
  git_rows <- if (is.null(g)) {
    data.frame(item = "git", value = "not a git repository, or git not installed")
  } else {
    data.frame(item = c("git commit", "git branch", "uncommitted changes", "git remote"),
               value = c(g$commit, g$branch %||% "", if (isTRUE(g$dirty)) "yes" else "no",
                         g$remote %||% ""))
  }
  add(html_table(rbind(
    data.frame(item = c("R", "platform", "OS", "BLAS", "LAPACK", names(env$packages)),
               value = c(env$r_version, env$platform, env$os, env$blas %||% "",
                         env$lapack %||% "", unlist(env$packages)),
               stringsAsFactors = FALSE),
    git_rows)))
  sess <- man$sessions %||% list()
  if (length(sess)) {
    sig0 <- environment_signature(env)
    differ <- vapply(sess, function(s) !identical(environment_signature(s), sig0), logical(1))
    commits <- unique(unlist(lapply(sess, function(s) s$git$commit)))
    add("<p>Sessions that worked on this run: ", length(sess), ". Sessions with a different software environment from the one that created the directory: ",
        sum(differ), ". Git commits used: ",
        if (length(commits)) html_escape(paste(substr(commits, 1, 12), collapse = ", ")) else "none recorded",
        ".</p>")
  }
  add("<h2>7. Failures and validity</h2>", html_table(attr(summ, "failures")[
    c("scenario_id", "method", "estimand", "generated", "generation_failed", "attempted",
      "usable", "statuses", "se_valid", "interval_valid", "test_valid", "with_warnings")]))
  if (!is.null(diag$messages)) add("<h3>Most frequent errors and warnings</h3>", html_table(diag$messages))
  add("<h2>8. Results</h2><p>Each value is followed by its Monte Carlo standard error and interval; n is the number of eligible replications out of the generated ones. Status other than ok explains why a value is unavailable or needs care.</p>")
  add(html_table(s[c("scenario_id", "method", "estimand", "label", "value", "mcse", "lower",
                     "upper", "n_eligible", "n_generated", "mcse_total", "truth_sens_lower",
                     "truth_sens_upper", "status", "note")]))
  rt <- attr(summ, "runtime")
  if (!is.null(rt)) {
    rt$scenario_id <- study$scenarios$labels[match(rt$scenario_key, study$scenarios$keys)]
    add("<h3>Computational cost (all attempts, including failures)</h3>",
        html_table(rt[c("scenario_id", "method", "attempts", "failed", "timeouts", "total_seconds",
                        "mean_seconds", "median_seconds", "p90_seconds", "max_seconds",
                        "mean_seconds_successful")]))
  }
  if (plots) {
    add("<h2>9. Diagnostics</h2>")
    for (i in seq_along(study$scenarios$keys)) {
      add("<h3>", html_escape(study$scenarios$labels[i]), "</h3>")
      add(html_plot(function() sim_plot(res, "estimates", scenario = i)))
      if (any(res$estimates$interval_status == "valid" &
              res$estimates$scenario_key == study$scenarios$keys[i] & !is.na(res$estimates$truth))) {
        add(html_plot(function() sim_plot(res, "zip", scenario = i)))
      }
    }
    add(html_plot(function() sim_plot(res, "failures")))
    for (ms in intersect(c("bias", "coverage", "rmse"), s$measure)) {
      add(html_plot(function() sim_plot(res, "precision", measure = ms),
                    height = 120 + 18 * sum(s$measure == ms)))
    }
  }
  add("<h3>Outlying values (flagged, not removed)</h3>", html_table(diag$outliers))
  add("<h2>10. Amendments, revisions and deviations</h2>")
  am <- if (length(study$amendments)) do.call(rbind, lapply(study$amendments, function(a) {
    data.frame(time = a$time, reason = a$reason, components = paste(a$components, collapse = ", "),
               after_freeze = a$after_freeze, stringsAsFactors = FALSE)
  })) else NULL
  add("<h3>Protocol amendments</h3>", html_table(am))
  rv <- if (length(man$revisions)) do.call(rbind, lapply(man$revisions, function(r) {
    data.frame(time = r$time, changes = paste(r$changes, collapse = "; "), stringsAsFactors = FALSE)
  })) else NULL
  add("<h3>Run revisions</h3>", html_table(rv))
  add("<h2>11. Limitations</h2><ul>",
      "<li>Performance measures are estimated from a finite number of replications; Monte Carlo standard errors quantify this uncertainty only under independent replications.</li>",
      "<li>Measures are conditional on each method producing the required output; when failures are informative, conditional performance and pairwise comparisons describe the successful subset, not the full data-generating distribution.</li>",
      if (any(s$truth_kind == "mc")) "<li>Some truths are Monte Carlo approximations; outer MCSEs are conditional on them, and sensitivity ranges show their effect.</li>" else "",
      if (any(s$truth_kind == "unknown")) "<li>Some truths are unknown, so truth based measures are unavailable for them.</li>" else "",
      "<li>Results describe the scenarios studied; they do not extend to other data-generating mechanisms.</li>",
      "<li>Exact replay assumes the same software environment and numerical libraries.</li></ul>")
  add("<h2>References</h2><ol>", paste0("<li>", html_escape(REFERENCES$text),
      " <a href=\"https://doi.org/", REFERENCES$doi, "\">doi:", REFERENCES$doi, "</a></li>",
      collapse = ""), "</ol></body></html>")
  writeLines(h, file, useBytes = TRUE)
  invisible(file)
}

major_version <- function(v) sub("\\..*$", "", v)

csv_safe <- function(df) {
  if (is.null(df)) return(data.frame())
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(col) {
    if (is.list(col)) {
      col <- vapply(col, function(v) paste(deparse(v), collapse = " "), character(1))
    }
    if (is.character(col) || is.factor(col)) col <- csv_neutralize(as.character(col))
    col
  })
  df
}

# A spreadsheet treats a cell starting with =, +, -, @, tab or carriage
# return as a formula. Such text cells (labels, error messages) get a
# leading single quote; numbers written as text are left alone.
csv_neutralize <- function(x) {
  risky <- !is.na(x) & grepl("^[=+@\t\r-]", x) &
    is.na(suppressWarnings(as.numeric(x)))
  x[risky] <- paste0("'", x[risky])
  x
}

#' Export and import a reproducibility bundle
#'
#' `sim_export()` writes a self-contained directory with the study
#' specification (including code), all result tables as CSV files and as an
#' R object, truth values, the random number registry, the environment,
#' the manifest with revisions, an HTML report, a README explaining
#' regeneration and reanalysis, a reanalysis script, and md5 checksums of
#' every file. Paths inside the bundle are relative, numerical library
#' names are recorded without their directories, and source references
#' (which hold the paths of sourced files) are removed from stored
#' functions. The git remote URL, when recorded, names the repository and
#' its owner. A function defined inside another function carries the
#' objects of its enclosing environment; `sim_export()` lists them in a
#' message so private objects are not shared unknowingly. Text cells that
#' a spreadsheet would read as formulas are written with a leading single
#' quote. Generated data are copied only with `include_data = TRUE`,
#' because they may be large or private.
#'
#' `overwrite = TRUE` replaces an earlier bundle entirely and refuses any
#' other non-empty directory. `sim_import()` verifies the checksums and
#' refuses a bundle whose schema major version differs from this version
#' of adempR.
#'
#' A bundle contains R code (the study's functions). The checksums detect
#' corruption, not tampering: import bundles only from trusted sources, as
#' you would source an R script.
#'
#' `sim_import()` verifies the checksums and returns results that work with
#' [sim_summarise()], [sim_compare()], [sim_diagnose()], [sim_plot()] and
#' [sim_report()].
#'
#' @param run A run or results.
#' @param dir Bundle directory.
#' @param overwrite Allow writing into an existing nonempty directory.
#' @param report Include the HTML report.
#' @param include_data Copy retained data and fit artifacts.
#' @return `sim_export()` returns `dir` invisibly; `sim_import()` returns
#'   an object of class `adempr_results`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 10L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' run <- sim_run(study, replications = 10, seed = 1, path = tempfile("adempR-run-"))
#' b <- sim_export(run, tempfile("adempR-bundle-"), report = FALSE)
#' res <- sim_import(b)
#' sim_summarise(res)
#' @export
sim_export <- function(run, dir, overwrite = FALSE, report = TRUE, include_data = FALSE) {
  if (!is_scalar_character(dir)) stopf("`dir` must be a path.")
  if (dir.exists(dir) && length(list.files(dir, all.files = TRUE, no.. = TRUE))) {
    if (!overwrite) stopf("%s exists and is not empty; use `overwrite = TRUE`.", dir)
    # Only an earlier bundle is replaced, and it is replaced entirely, so no
    # stale file from the earlier export survives into the new checksums.
    if (!file.exists(file.path(dir, "CHECKSUMS.md5"))) {
      stopf("%s is not empty and is not an adempR bundle; refusing to overwrite it.", dir)
    }
    unlink(dir, recursive = TRUE)
  }
  res <- sim_results(run)
  man <- res$manifest
  study <- man$study
  dir.create(file.path(dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "code"), showWarnings = FALSE)
  summ <- sim_summarise(res)
  tabs <- list(scenarios = res$scenarios, replications = res$replications,
               attempts = res$attempts, estimates = res$estimates, truth = res$truth,
               summary = as.data.frame(summ), failures = attr(summ, "failures"),
               runtime = attr(summ, "runtime"), rng_registry = man$rng_registry,
               excluded_scenarios = study$scenarios$excluded)
  for (nm in names(tabs)) {
    utils::write.csv(csv_safe(tabs[[nm]]), file.path(dir, "tables", paste0(nm, ".csv")),
                     row.names = FALSE, na = "NA")
  }
  writeLines(fun_source(study$generator$fun), file.path(dir, "code", "generator.R"))
  for (mn in names(study$methods)) {
    writeLines(c(fun_source(study$methods[[mn]]$fun),
                 if (!is.null(study$methods[[mn]]$extract))
                   c("", "# extract", fun_source(study$methods[[mn]]$extract))),
               file.path(dir, "code", paste0("method-", mn, ".R")))
  }
  for (en in names(study$estimands)) {
    tf <- study$estimands[[en]]$truth$fun
    if (!is.null(tf)) writeLines(fun_source(tf), file.path(dir, "code", paste0("truth-", en, ".R")))
  }
  carried <- closure_contents(study)
  if (length(carried)) {
    message(sprintf("Functions in this study are closures whose enclosing environments are exported with them, including: %s. Check that none of these objects is private.",
                    paste(utils::head(carried, 20), collapse = ", ")))
  }
  exported <- res
  exported$path <- NULL
  exported$manifest$study <- strip_sources(exported$manifest$study)
  if (!is.null(exported$replications)) exported$replications$artifact <- NA_character_
  if (include_data && !is.null(res$path)) {
    src <- file.path(res$path, "artifacts")
    if (dir.exists(src)) {
      dir.create(file.path(dir, "artifacts"), showWarnings = FALSE)
      file.copy(list.files(src, full.names = TRUE), file.path(dir, "artifacts"))
      exported$replications$artifact <- res$replications$artifact
    }
  }
  saveRDS(exported, file.path(dir, "results.rds"))
  jsonlite::write_json(list(
    schema_version = SCHEMA_VERSION, study = study$name, aims = study$aims,
    label = man$label, seed = man$seed, replications = man$replications,
    study_fingerprint = man$study_fingerprint, fingerprints = man$fingerprints,
    revisions = man$revisions, amendments = study$amendments,
    frozen = study$status$frozen, protocol = unclass(study$protocol),
    environment = man$environment, sessions = man$sessions, sequential = man$sequential,
    exported = now_utc()),
    file.path(dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null",
    na = "null", digits = NA)
  if (report) sim_report(exported, file.path(dir, "report.html"), overwrite = TRUE)
  writeLines(c(
    paste0("# Reproducibility bundle: ", study$name), "",
    paste0("Exported ", now_utc(), " with adempR ", utils::packageVersion("adempR"), "."), "",
    "## Contents", "",
    "- `results.rds`: results object (study specification with code, all tables, manifest).",
    "- `tables/`: the same tables as CSV files (estimates with validity statuses, attempts, replications, truth, summary, failures, runtime, RNG registry).",
    "- `code/`: source of the generator, methods and truth functions.",
    "- `manifest.json`: seed, fingerprints, protocol, revisions, amendments, environment.",
    "- `report.html`: self-contained report (if exported).",
    "- `CHECKSUMS.md5`: md5 checksums of every file.", "",
    "## Reanalysis", "",
    "```r", "library(adempR)", "res <- sim_import(\".\")   # verifies the checksums",
    "sim_summarise(res)", "```", "",
    "## Regeneration", "",
    "The study can be rerun from `res$manifest$study` with the recorded seed:", "",
    "```r",
    paste0("run <- sim_run(res$manifest$study, replications = ", man$replications,
           ", seed = ", man$seed, ", path = \"rerun\")"), "```", "",
    "Exact regeneration requires the same R version, package versions and numerical",
    "libraries (see `manifest.json`). A package lockfile does not recreate missing",
    "system libraries."), file.path(dir, "README.md"))
  writeLines(c("library(adempR)", "res <- sim_import(\".\")", "print(sim_summarise(res))"),
             file.path(dir, "reanalyse.R"))
  files <- list.files(dir, recursive = TRUE, all.files = TRUE)
  files <- setdiff(files, "CHECKSUMS.md5")
  sums <- tools::md5sum(file.path(dir, files))
  writeLines(paste(unname(sums), files, sep = "  "), file.path(dir, "CHECKSUMS.md5"))
  invisible(dir)
}

#' @rdname sim_export
#' @export
sim_import <- function(dir) {
  ck <- file.path(dir, "CHECKSUMS.md5")
  if (!file.exists(ck)) stopf("No CHECKSUMS.md5 in %s; not an adempR bundle.", dir)
  lines <- readLines(ck)
  if (!any(grepl("  results\\.rds$", lines))) {
    stopf("CHECKSUMS.md5 does not list results.rds; not a valid adempR bundle.")
  }
  sums <- sub("  .*$", "", lines)
  files <- sub("^[0-9a-f]+  ", "", lines)
  actual <- unname(tools::md5sum(file.path(dir, files)))
  bad <- files[is.na(actual) | actual != sums]
  if (length(bad)) stopf("Checksum mismatch or missing file(s) in the bundle: %s.",
                         paste(bad, collapse = ", "))
  mj <- file.path(dir, "manifest.json")
  version <- if (file.exists(mj)) jsonlite::read_json(mj)$schema_version else NULL
  if (!is_scalar_character(version)) stopf("The bundle has no schema version in manifest.json.")
  if (!identical(major_version(version), major_version(SCHEMA_VERSION))) {
    stopf("The bundle has schema version %s; this version of adempR reads schema %s.x.",
          version, major_version(SCHEMA_VERSION))
  }
  res <- readRDS(file.path(dir, "results.rds"))
  if (!inherits(res, "adempr_results")) stopf("results.rds does not hold adempR results.")
  res
}
