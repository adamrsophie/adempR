# Protocol export ---------------------------------------------------------

#' Export the study protocol as a preregistration document
#'
#' Writes the protocol in the section order of the ADEMP-PreReg template
#' for simulation studies (Siepe et al. 2024): general information, aims,
#' data-generating mechanism, estimands and other targets, methods,
#' performance measures, and other details (replications and Monte Carlo
#' error, failure handling, software, random numbers). Every entry is taken
#' from the study object, so the document matches what [sim_run()] will
#' execute; fields that were not specified are marked as such, which shows
#' what is still missing before registration (for example on OSF).
#'
#' @param study A study.
#' @param file Output Markdown file.
#' @param replications Planned replications per scenario, if decided.
#' @param overwrite Overwrite an existing file.
#' @param code Include generator, truth and method code.
#' @return The file path, invisibly.
#' @references Siepe BS, Bartoš F, Morris TP, Boulesteix AL, Heck DW,
#'   Pawel S (2024). Simulation studies for methodological research in
#'   psychology: a standardized template for planning, preregistration,
#'   and reporting. Psychological Methods.
#' @export
sim_protocol_export <- function(study, file, replications = NULL, overwrite = FALSE, code = TRUE) {
  stopifnot(inherits(study, "adempr_study"))
  if (!is_scalar_character(file)) stopf("`file` must be a path.")
  if (file.exists(file) && !overwrite) stopf("%s exists; use `overwrite = TRUE` to replace it.", file)
  p <- study$protocol
  na <- "*Not specified.*"
  txt <- function(x) {
    if (is.null(x) || !length(x)) return(na)
    v <- unlist(x)
    if (!is.null(names(v)) && any(nzchar(names(v)))) v <- paste0(names(v), " = ", v)
    paste(v, collapse = "; ")
  }
  md_table <- function(df) {
    if (is.null(df) || !nrow(df)) return(na)
    df[] <- lapply(df, function(col) {
      v <- if (is.list(col)) vapply(col, function(e) paste(deparse(e), collapse = " "), character(1)) else
        as.character(col)
      v[is.na(v)] <- ""
      gsub("|", "\\|", gsub("\n", " ", v, fixed = TRUE), fixed = TRUE)
    })
    c(paste0("| ", paste(names(df), collapse = " | "), " |"),
      paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
      apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  code_block <- function(f) c("```r", fun_source(f), "```")
  frozen <- study$status$frozen
  out <- c(
    paste0("# Simulation study protocol: ", study$name), "",
    "Section order follows the ADEMP-PreReg template (Siepe et al. 2024).", "",
    "## 0. General information", "",
    paste0("- Study name: ", study$name),
    paste0("- Purpose: ", txt(p$purpose)),
    paste0("- Scope: ", txt(p$scope)),
    paste0("- Intended users: ", txt(p$users)),
    paste0("- Phase of methodological research: ", txt(p$extra$phase)),
    paste0("- Registration: ", txt(p$registration)),
    paste0("- Protocol status: ", if (is.null(frozen)) "not frozen" else
      paste0("frozen on ", frozen$time, ", hash `", frozen$hash, "`")),
    paste0("- Study fingerprint: `", sim_fingerprint(study), "`"),
    paste0("- Amendments: ", length(study$amendments)), "",
    "## 1. Aims", "", study$aims %||% na, "",
    "## 2. Data-generating mechanism", "",
    study$generator$description %||% na, "",
    paste0("Rationale: ", txt(p$dgm_rationale)), "",
    paste0("Design: ", study$scenarios$design, ", ", length(study$scenarios), " scenario(s)."), "",
    md_table(as.data.frame(study$scenarios)[setdiff(names(as.data.frame(study$scenarios)), "scenario_key")]), "")
  if (!is.null(study$scenarios$excluded) && NROW(study$scenarios$excluded)) {
    out <- c(out, "Scenarios excluded by design:", "", md_table(study$scenarios$excluded), "")
  }
  if (code) out <- c(out, "Generator:", "", code_block(study$generator$fun), "")
  ed <- do.call(rbind, lapply(study$estimands, function(e) {
    data.frame(id = e$id, quantity = e$quantity, population = e$population %||% "",
               effect = e$effect, scale = e$scale, truth = e$truth$kind, stringsAsFactors = FALSE)
  }))
  out <- c(out, "## 3. Estimands and other targets", "",
           paste0("Rationale: ", txt(p$estimand_rationale)), "", md_table(ed), "")
  if (code) for (e in study$estimands) {
    if (is.function(e$truth$fun)) out <- c(out, paste0("Truth of `", e$id, "` (", e$truth$kind, "):"), "",
                                           code_block(e$truth$fun), "")
  }
  md <- do.call(rbind, lapply(names(study$methods), function(mn) {
    m <- study$methods[[mn]]
    data.frame(method = mn, role = m$role %||% "method", estimands = paste(m$estimands, collapse = ", "),
               packages = paste(m$packages, collapse = ", "),
               settings_source = m$settings_source %||% "",
               description = m$description %||% "", stringsAsFactors = FALSE)
  }))
  out <- c(out, "## 4. Methods", "", paste0("Rationale: ", txt(p$method_rationale)), "", md_table(md), "")
  if (code) for (mn in names(study$methods)) {
    out <- c(out, paste0("Method `", mn, "`:"), "", code_block(study$methods[[mn]]$fun), "")
  }
  mi <- MEASURE_INFO[MEASURE_INFO$measure %in% names(study$measures), c("measure", "eligibility", "mcse_method")]
  out <- c(out, "## 5. Performance measures", "", paste0("Rationale: ", txt(p$measure_rationale)), "",
           md_table(mi), "",
           paste0("Nominal level of intervals: ", study$level, "; significance level: ", study$alpha, "."), "",
           "## 6. Other", "",
           "### Number of repetitions and Monte Carlo error", "",
           paste0("- Planned replications per scenario: ", if (is.null(replications)) na else replications),
           paste0("- Rationale: ", txt(p$replication_rationale)),
           paste0("- Precision targets (Monte Carlo SE): ", txt(p$precision_targets)), "",
           "### Handling of failures and missing results", "",
           paste0("- Plan: ", txt(p$failure_handling)),
           paste0("- Retries: ", study$failure$retries, "; warnings: ", txt(study$failure$warnings),
                  "; fallback: ", txt(study$failure$fallback)),
           "- Data sets are never redrawn after a method fails.", "",
           "### Software, random numbers and retention", "",
           paste0("- R ", R.version$major, ".", R.version$minor, ", adempR ",
                  as.character(utils::packageVersion("adempR"))),
           "- Random numbers: L'Ecuyer-CMRG streams, one per scenario, replication and stage, registered in the run manifest.",
           paste0("- Retention: ", txt(p$retention)),
           paste0("- Reporting: ", txt(p$reporting)), "")
  extra <- p$extra[setdiff(names(p$extra), "phase")]
  if (length(extra)) {
    out <- c(out, "### Additional fields", "",
             paste0("- ", names(extra), ": ", vapply(extra, txt, character(1))), "")
  }
  writeLines(out, file)
  invisible(file)
}
