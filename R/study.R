# Study specification, protocol, validation -------------------------------

PURPOSES <- c("unspecified", "exploratory", "confirmatory", "method_comparison",
              "design_power")

#' Study protocol metadata
#'
#' Scientific documentation that accompanies an executable study: aims,
#' scope, justification of the data-generating mechanisms, estimands,
#' methods and performance measures, planned replications and precision,
#' failure handling, data retention and reporting. Fields left empty are
#' reported by [sim_completeness()]; they never block an exploratory run.
#'
#' A locally frozen protocol is not a preregistration. Record an actual
#' registration (registry, identifier, date) in `registration`.
#'
#' @param purpose One of `"unspecified"`, `"exploratory"`,
#'   `"confirmatory"`, `"method_comparison"`, `"design_power"`.
#' @param scope,users,dgm_rationale,estimand_rationale,method_rationale,measure_rationale,replication_rationale,failure_handling,retention,reporting
#'   Free text.
#' @param precision_targets Optional named numeric vector of target Monte
#'   Carlo standard errors by measure (for example `c(coverage = 0.005)`).
#' @param registration Optional list with `registry`, `id`, `date` and
#'   `url` of an external registration.
#' @param ... Further named fields, kept verbatim and shown in the report:
#'   for example `phase` (the phase of methodological research, I to IV,
#'   Heinze et al. 2023), `neutral` (whether the authors have a stake in
#'   any method), `inventor_methods`, `primary_measure`, or OCTAVE fields
#'   of a clinical trial simulation.
#' @return An object of class `adempr_protocol`.
#' @examples
#' sim_protocol(purpose = "method_comparison",
#'              dgm_rationale = "sample sizes typical of pilot trials",
#'              precision_targets = c(coverage = 0.005))
#' @export
sim_protocol <- function(purpose = "unspecified", scope = NULL, users = NULL,
                         dgm_rationale = NULL, estimand_rationale = NULL,
                         method_rationale = NULL, measure_rationale = NULL,
                         replication_rationale = NULL, precision_targets = NULL,
                         failure_handling = NULL, retention = NULL,
                         reporting = NULL, registration = NULL, ...) {
  purpose <- match.arg(purpose, PURPOSES)
  if (!is.null(precision_targets) &&
      (!is.numeric(precision_targets) || is.null(names(precision_targets)) ||
       any(!(precision_targets > 0)))) {
    stopf("`precision_targets` must be a named vector of positive numbers.")
  }
  if (!is.null(registration) && !is.list(registration)) {
    stopf("`registration` must be a list with registry, id, date and url.")
  }
  structure(list(purpose = purpose, scope = scope, users = users,
                 dgm_rationale = dgm_rationale,
                 estimand_rationale = estimand_rationale,
                 method_rationale = method_rationale,
                 measure_rationale = measure_rationale,
                 replication_rationale = replication_rationale,
                 precision_targets = precision_targets,
                 failure_handling = failure_handling, retention = retention,
                 reporting = reporting, registration = registration,
                 extra = list(...)),
            class = "adempr_protocol")
}

#' Data retention policy
#'
#' Estimates, intervals, tests, statuses and timings are always kept.
#' Optionally the generated data sets and the raw method outputs are kept as
#' artifacts, which makes replay possible without regeneration and allows
#' later diagnostics of fitted objects. [sim_capabilities()] reports which
#' analyses remain possible under the chosen policy.
#'
#' @param data Keep every generated data set.
#' @param fits Keep every raw method output (before extraction).
#' @return An object of class `adempr_retention`.
#' @examples
#' sim_retention(data = TRUE)
#' @export
sim_retention <- function(data = FALSE, fits = FALSE) {
  if (!is.logical(data) || length(data) != 1L || is.na(data) ||
      !is.logical(fits) || length(fits) != 1L || is.na(fits)) {
    stopf("`data` and `fits` must be TRUE or FALSE.")
  }
  structure(list(data = data, fits = fits), class = "adempr_retention")
}

#' Define a simulation study
#'
#' Bundles the ADEMP components (aims, data-generating mechanisms,
#' estimands, methods and performance measures) with the protocol, the
#' failure policy and the retention policy into one validated, immutable
#' specification. Execution state lives in run directories, never in the
#' study object.
#'
#' @param name Study identifier.
#' @param aims Text stating the aims.
#' @param scenarios A design from [sim_scenarios()].
#' @param generate A function `generate(scenario, context)` or a
#'   [sim_generator()].
#' @param estimands A list of [sim_estimand()] objects.
#' @param methods A named list of [sim_method()] objects or plain functions
#'   (wrapped with default settings).
#' @param measures Performance measures from [sim_measures()] or a
#'   character vector of measure names.
#' @param protocol A [sim_protocol()].
#' @param failure A [sim_failure_policy()].
#' @param retention A [sim_retention()].
#' @param alpha Significance level for rejection rates when a method does
#'   not report its own.
#' @param level Confidence level for Monte Carlo intervals of performance
#'   measures.
#' @return An object of class `adempr_study`.
#' @examples
#' study <- sim_study(
#'   name = "normal_mean",
#'   aims = "Coverage of the t interval for a normal mean",
#'   scenarios = sim_scenarios(n = c(10L, 50L), mu = 0),
#'   generate = function(scenario, context) rnorm(scenario$n, scenario$mu),
#'   estimands = list(sim_estimand("mean",
#'     truth = sim_truth_exact(function(scenario) scenario$mu))),
#'   methods = list(t_test = function(data) {
#'     tt <- t.test(data)
#'     sim_estimate("mean", mean(data), se = sd(data) / sqrt(length(data)),
#'                  lower = tt$conf.int[1], upper = tt$conf.int[2], level = 0.95,
#'                  p_value = tt$p.value, null = 0)
#'   })
#' )
#' study
#' @export
sim_study <- function(name, aims = NULL, scenarios, generate, estimands, methods,
                      measures = sim_measures(), protocol = sim_protocol(),
                      failure = sim_failure_policy(), retention = sim_retention(),
                      alpha = 0.05, level = 0.95) {
  check_id(name, "name")
  if (!inherits(scenarios, "adempr_scenarios")) {
    stopf("`scenarios` must be created with sim_scenarios() or a related function.")
  }
  if (is.function(generate)) generate <- sim_generator(generate)
  if (!inherits(generate, "adempr_generator")) {
    stopf("`generate` must be a function or a sim_generator().")
  }
  if (inherits(estimands, "adempr_estimand")) estimands <- list(estimands)
  if (!is.list(estimands) || !length(estimands) ||
      !all(vapply(estimands, inherits, logical(1), "adempr_estimand"))) {
    stopf("`estimands` must be a list of sim_estimand() objects.")
  }
  ids <- vapply(estimands, `[[`, character(1), "id")
  if (anyDuplicated(ids)) stopf("Estimand identifiers must be unique.")
  names(estimands) <- ids
  if (!is.list(methods) || !length(methods)) stopf("`methods` must be a nonempty named list.")
  mnames <- names(methods)
  if (is.null(mnames) || any(!nzchar(mnames)) || anyDuplicated(mnames)) {
    stopf("`methods` must be a list with unique names.")
  }
  for (m in mnames) check_id(m, "method name")
  methods <- lapply(methods, function(m) {
    if (is.function(m)) m <- sim_method(m)
    if (!inherits(m, "adempr_method")) stopf("Each method must be a function or a sim_method().")
    if (is.null(m$estimands)) m$estimands <- ids
    m
  })
  if (is.character(measures)) measures <- do.call(sim_measures, as.list(measures))
  if (!inherits(measures, "adempr_measures")) stopf("`measures` must come from sim_measures().")
  if (!inherits(protocol, "adempr_protocol")) stopf("`protocol` must come from sim_protocol().")
  if (!inherits(failure, "adempr_failure_policy")) stopf("`failure` must come from sim_failure_policy().")
  if (!inherits(retention, "adempr_retention")) stopf("`retention` must come from sim_retention().")
  if (!is_scalar_number(alpha) || alpha <= 0 || alpha >= 1) stopf("`alpha` must be in (0, 1).")
  if (!is_scalar_number(level) || level <= 0 || level >= 1) stopf("`level` must be in (0, 1).")
  structure(list(name = name, aims = aims, scenarios = scenarios,
                 generator = generate, estimands = estimands, methods = methods,
                 measures = measures, protocol = protocol, failure = failure,
                 retention = retention, alpha = alpha, level = level,
                 status = list(label = "draft", frozen = NULL),
                 amendments = list()),
            class = "adempr_study")
}

#' @export
print.adempr_study <- function(x, ...) {
  cat(sprintf("<adempR study `%s` (%s)>\n", x$name, x$status$label))
  if (!is.null(x$aims)) cat("Aims: ", x$aims, "\n", sep = "")
  cat(sprintf("Scenarios: %d (%s design)\n", length(x$scenarios), x$scenarios$design))
  cat("Estimands: ", paste(names(x$estimands), collapse = ", "), "\n", sep = "")
  cat("Methods:   ", paste(names(x$methods), collapse = ", "), "\n", sep = "")
  cat("Measures:  ", paste(names(x$measures), collapse = ", "), "\n", sep = "")
  if (!is.null(x$status$frozen)) {
    cat(sprintf("Frozen %s, hash %s\n", x$status$frozen$time,
                short_hash(x$status$frozen$hash)))
  }
  if (length(x$amendments)) cat(sprintf("Amendments: %d\n", length(x$amendments)))
  invisible(x)
}

# Stage fingerprints. A whole-study fingerprint identifies the protocol; the
# stage fingerprints decide what must be recomputed.
study_fingerprints <- function(study) {
  list(
    generator = generator_fingerprint(study$generator),
    scenarios = stats::setNames(study$scenarios$keys, study$scenarios$keys),
    design = scenario_design_fingerprint(study$scenarios),
    methods = lapply(study$methods, method_fingerprint),
    truth = lapply(study$estimands, function(e) truth_fingerprint(e$truth)),
    estimands = lapply(study$estimands, function(e) {
      fingerprint(e[setdiff(names(e), "truth")])
    }),
    evaluation = fingerprint(lapply(study$measures, function(m) {
      if (is.character(m)) m else unclass(m)
    }), unclass(study$failure)[setdiff(names(study$failure), "retries")],
    study$alpha, study$level),
    # Settings that change how fits are executed but not what a method is:
    # recorded so that a change is visible, never used to select results.
    execution = execution_settings(study)
  )
}

# Everything about the scenario set that is not a scenario's identity:
# weights (keyed by scenario, so order does not matter), the sampling
# distribution, the design description, constraints and exclusions.
scenario_design_fingerprint <- function(sc) {
  w <- sc$weights
  if (!is.null(w)) {
    w <- stats::setNames(w, sc$keys)
    w <- w[order(names(w))]
  }
  fingerprint(w, sc$distribution, sc$design, sc$constraints, sc$excluded)
}

execution_settings <- function(study) {
  list(retries = study$failure$retries,
       timeout = lapply(study$methods, function(m) m$timeout))
}

#' @rdname sim_freeze
#' @export
sim_fingerprint <- function(study) {
  stopifnot(inherits(study, "adempr_study"))
  fp <- study_fingerprints(study)
  fingerprint(fp, unclass(study$protocol), study$aims, study$name,
              unclass(study$retention))
}

#' Validate a study before running it
#'
#' Performs deterministic checks: identifiers, method and estimand
#' declarations, compatibility of declared effect types and scales, method
#' scenario fields, generator checks for every scenario (impossible
#' parameter values, non positive definite matrices, probabilities outside
#' the unit interval), and exact or numerical truth for every scenario.
#' With `smoke = TRUE`, one data set is also generated for each scenario and
#' every method is run once on it, using a dedicated validation stream that
#' does not affect any later run.
#'
#' @param study A [sim_study()].
#' @param smoke Run the smoke test.
#' @param seed Seed for the smoke test.
#' @param error If `TRUE`, stop when problems are found.
#' @param large_n Optional named list of scenario fields to override for a
#'   large sample sanity run, for example `list(n = 1e5)`. One data set per
#'   scenario is generated with the overrides and every method is run once;
#'   an estimate more than 4 of its standard errors (or, without a standard
#'   error, more than 10% of the absolute truth plus 0.01) away from an
#'   exact or numerical truth is reported as a warning. A consistent method
#'   should be close to its estimand in a large sample; a large gap often
#'   means the estimand and the method target different quantities.
#' @return An object of class `adempr_validation`: a data frame of checks
#'   with columns `check`, `status` (`"ok"`, `"warning"`, `"error"`) and
#'   `message`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' sim_validate(study, smoke = TRUE, seed = 1)
#' @export
sim_validate <- function(study, smoke = FALSE, seed = 1L, error = TRUE, large_n = NULL) {
  stopifnot(inherits(study, "adempr_study"))
  out <- list()
  add <- function(check, status, message) {
    out[[length(out) + 1L]] <<- data.frame(check = check, status = status,
                                           message = message, stringsAsFactors = FALSE)
  }
  est_ids <- names(study$estimands)
  scen_fields <- unique(unlist(lapply(study$scenarios$values, names)))
  for (pkg in study_packages(study)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      add(paste0("package:", pkg), "error", "declared package is not installed")
    }
  }
  for (mn in names(study$methods)) {
    m <- study$methods[[mn]]
    unknown <- setdiff(m$estimands, est_ids)
    if (length(unknown)) {
      add(paste0("method:", mn), "error",
          sprintf("targets unknown estimand(s) %s", paste(unknown, collapse = ", ")))
    }
    bad_fields <- setdiff(m$scenario_fields, scen_fields)
    if (length(bad_fields)) {
      add(paste0("method:", mn), "error",
          sprintf("declares unknown scenario field(s) %s", paste(bad_fields, collapse = ", ")))
    }
    for (en in names(m$declares)) {
      if (!(en %in% est_ids)) {
        add(paste0("method:", mn), "error", sprintf("declares metadata for unknown estimand %s", en))
        next
      }
      d <- m$declares[[en]]
      e <- study$estimands[[en]]
      for (fld in intersect(names(d), c("effect", "scale"))) {
        if (!identical(d[[fld]], e[[fld]])) {
          add(paste0("method:", mn), "error",
              sprintf("declares %s `%s` for estimand `%s`, whose %s is `%s`; the method does not estimate this estimand",
                      fld, d[[fld]], en, fld, e[[fld]]))
        }
      }
    }
    if (m$oracle) add(paste0("method:", mn), "warning", "oracle method: sees privileged information")
  }
  unused <- setdiff(est_ids, unlist(lapply(study$methods, `[[`, "estimands")))
  if (length(unused)) {
    add("estimands", "warning", sprintf("no method targets %s", paste(unused, collapse = ", ")))
  }
  if (!is.null(study$generator$check)) {
    for (i in seq_along(study$scenarios$values)) {
      msg <- tryCatch({
        deterministic_call("The generator check", study$generator$check,
                           study$scenarios$values[[i]])
        NULL
      }, error = function(e) conditionMessage(e))
      if (!is.null(msg)) {
        add(paste0("scenario:", study$scenarios$labels[i]), "error",
            paste("invalid configuration:", msg))
      }
    }
  }
  for (en in est_ids) {
    tr <- study$estimands[[en]]$truth
    if (tr$kind %in% c("exact", "numeric")) {
      for (i in seq_along(study$scenarios$values)) {
        msg <- tryCatch({
          compute_truth_one(tr, study$scenarios$values[[i]])
          NULL
        }, error = function(e) conditionMessage(e))
        if (!is.null(msg)) {
          add(paste0("truth:", en), "error",
              sprintf("scenario %s: %s", study$scenarios$labels[i], msg))
        }
      }
    }
    if (tr$kind == "unknown") {
      add(paste0("truth:", en), "warning", "truth unknown: bias, MSE and coverage are unavailable")
    }
  }
  if (smoke) {
    smk <- tryCatch(with_packages(study_packages(study), smoke_test(study, seed)), error = function(e) {
      data.frame(check = "smoke", status = "error", message = conditionMessage(e),
                 stringsAsFactors = FALSE)
    })
    out[[length(out) + 1L]] <- smk
  }
  if (!is.null(large_n)) {
    if (!is.list(large_n) || is.null(names(large_n)) || !all(nzchar(names(large_n)))) {
      stopf("`large_n` must be a named list of scenario fields, for example list(n = 1e5).")
    }
    out[[length(out) + 1L]] <- with_packages(study_packages(study),
                                             large_n_check(study, seed, large_n))
  }
  res <- bind_rows(out) %||% data.frame(check = character(), status = character(),
                                         message = character(), stringsAsFactors = FALSE)
  if (!nrow(res)) {
    res <- data.frame(check = "all", status = "ok", message = "no problems found",
                      stringsAsFactors = FALSE)
  }
  class(res) <- c("adempr_validation", "data.frame")
  if (error && any(res$status == "error")) {
    print(res)
    stopf("Study validation failed with %d error(s).", sum(res$status == "error"))
  }
  res
}

large_n_check <- function(study, seed, overrides) {
  big <- study
  big$scenarios$values <- lapply(study$scenarios$values, function(v) {
    v[names(overrides)] <- overrides
    v
  })
  keys <- data.frame(stage = "large_n", scenario_key = study$scenarios$keys,
                     component = "data", stringsAsFactors = FALSE)
  mkeys <- expand.grid(stage = "large_n_method", scenario_key = study$scenarios$keys,
                       component = method_stream_components(study), stringsAsFactors = FALSE)
  ctx <- list(seed = seed, registry = rng_registry(rbind(keys, mkeys)),
              cache = new.env(parent = emptyenv()), dgp_stage = "large_n",
              method_stage = "large_n_method")
  rows <- list()
  for (i in seq_along(big$scenarios$values)) {
    sc <- big$scenarios$values[[i]]
    res <- execute_task(big, ctx, study$scenarios$keys[i], 1L, names(study$methods))
    est <- res$estimates
    lab <- paste0("large_n:", study$scenarios$labels[i])
    if (res$replication$generation_status != "ok" || is.null(est)) {
      rows[[length(rows) + 1L]] <- data.frame(check = lab, status = "warning",
        message = "no estimates from the large sample run", stringsAsFactors = FALSE)
      next
    }
    for (j in seq_len(nrow(est))) {
      tr <- study$estimands[[est$estimand[j]]]$truth
      if (!(tr$kind %in% c("exact", "numeric")) || !is.finite(est$estimate[j])) next
      truth <- compute_truth_one(tr, sc)$truth
      gap <- abs(est$estimate[j] - truth)
      se <- est$se[j]
      far <- if (is.finite(se) && se > 0) gap > 4 * se else gap > 0.1 * abs(truth) + 0.01
      rows[[length(rows) + 1L]] <- data.frame(check = paste0(lab, ":", est$method[j], ":",
                                                            est$estimand[j]),
        status = if (far) "warning" else "ok",
        message = sprintf("estimate %.4g, truth %.4g%s", est$estimate[j], truth,
                          if (far) ": far from the truth in a large sample (different estimand, or not consistent?)" else ""),
        stringsAsFactors = FALSE)
    }
  }
  bind_rows(rows)
}

smoke_test <- function(study, seed) {
  keys <- data.frame(stage = "validation", scenario_key = study$scenarios$keys,
                     component = "data", stringsAsFactors = FALSE)
  mkeys <- expand.grid(stage = "validation_method", scenario_key = study$scenarios$keys,
                       component = method_stream_components(study),
                       stringsAsFactors = FALSE)
  reg <- rng_registry(rbind(keys, mkeys))
  ctx <- list(seed = seed, registry = reg, cache = new.env(parent = emptyenv()),
              dgp_stage = "validation", method_stage = "validation_method")
  rows <- list()
  for (i in seq_along(study$scenarios$values)) {
    res <- execute_task(study, ctx, study$scenarios$keys[i], 1L, names(study$methods))
    att <- res$attempts
    gen <- res$replication
    msg <- if (gen$generation_status != "ok") {
      paste("generation failed:", gen$generation_message)
    } else {
      bad <- att[att$status != "ok", , drop = FALSE]
      if (nrow(bad)) {
        paste(sprintf("%s %s: %s", bad$method, bad$status, bad$message), collapse = "; ")
      } else {
        "all methods returned results"
      }
    }
    rows[[i]] <- data.frame(check = paste0("smoke:", study$scenarios$labels[i]),
                            status = if (gen$generation_status != "ok") "error" else
                              if (any(att$status != "ok")) "warning" else "ok",
                            message = msg, stringsAsFactors = FALSE)
  }
  bind_rows(rows)
}

#' @export
print.adempr_validation <- function(x, ...) {
  cat(sprintf("<adempR validation: %d error(s), %d warning(s)>\n",
              sum(x$status == "error"), sum(x$status == "warning")))
  df <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  print(df, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Completeness of the scientific documentation
#'
#' Lists which protocol items are documented and which are missing, grouped
#' by ADEMP component and following the reporting items of Morris et al.
#' (2019), Siepe et al. (2024) and Williams et al. (2024). Missing items are
#' information, not errors; exploratory runs are never blocked.
#'
#' @param study A [sim_study()].
#' @return A data frame with columns `section`, `item`, `status` and
#'   `detail`.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' sim_completeness(study)
#' @export
sim_completeness <- function(study) {
  stopifnot(inherits(study, "adempr_study"))
  p <- study$protocol
  has <- function(x) !is.null(x) && length(x) && !all(is.na(x)) &&
    !(is.character(x) && !any(nzchar(x)))
  item <- function(section, name, present, detail = "") {
    data.frame(section = section, item = name,
               status = if (present) "documented" else "missing",
               detail = detail, stringsAsFactors = FALSE)
  }
  est_doc <- vapply(study$estimands, function(e) has(e$population) && e$effect != "unspecified",
                    logical(1))
  meth_doc <- vapply(study$methods, function(m) has(m$description) || has(m$version), logical(1))
  rows <- list(
    item("Aims", "aims", has(study$aims)),
    item("Aims", "purpose", p$purpose != "unspecified", p$purpose),
    item("Aims", "scope and intended users", has(p$scope) || has(p$users)),
    item("Data-generating mechanisms", "description", has(study$generator$description)),
    item("Data-generating mechanisms", "rationale for factors and levels", has(p$dgm_rationale)),
    item("Data-generating mechanisms", "excluded scenarios explained",
         TRUE, sprintf("%d excluded", nrow(study$scenarios$excluded))),
    item("Estimands", "population and effect type for every estimand", all(est_doc),
         paste(names(study$estimands)[!est_doc], collapse = ", ")),
    item("Estimands", "truth provider for every estimand",
         all(vapply(study$estimands, function(e) e$truth$kind != "unknown", logical(1)))),
    item("Estimands", "rationale", has(p$estimand_rationale)),
    item("Methods", "description or version for every method", all(meth_doc),
         paste(names(study$methods)[!meth_doc], collapse = ", ")),
    item("Methods", "rationale for the methods compared", has(p$method_rationale)),
    item("Performance measures", "rationale", has(p$measure_rationale)),
    item("Performance measures", "planned replications justified", has(p$replication_rationale)),
    item("Performance measures", "precision targets", has(p$precision_targets)),
    item("Failures", "failure handling stated", has(p$failure_handling)),
    item("Reproducibility", "retention policy stated", has(p$retention)),
    item("Reporting", "reporting plan", has(p$reporting)),
    item("Registration", "external registration", has(p$registration),
         "a frozen protocol hash is not a registration"),
    item("Aims", "phase of methodological research (I to IV)", has(p$extra$phase),
         "sim_protocol(phase = ...)"),
    # reported by adempR for every run
    item("Performance measures", "Monte Carlo standard error for every measure", TRUE,
         "reported by adempR"),
    item("Failures", "failures and exclusions reported", TRUE, "flow table in the report"),
    item("Reproducibility", "code archived", TRUE,
         "code in the report and in sim_export() bundles")
  )
  bind_rows(rows)
}

#' Freeze and amend a protocol
#'
#' `sim_freeze()` marks a study as frozen before a production run and
#' records its fingerprint and time. A frozen hash proves the identity of
#' the specification; it does not prove that it was fixed before results
#' were seen. `sim_amend()` changes components of a study and records the
#' reason, the time and the fingerprints before and after, so deviations
#' from the frozen protocol are reported. [sim_run()] refuses a frozen
#' study whose fingerprint differs from the frozen (or last amended) hash,
#' so every change after freezing goes through `sim_amend()`.
#' `sim_fingerprint()` returns the whole-study fingerprint; it covers the
#' scenario weights, sampling distribution, constraints and exclusions as
#' well as the scenarios themselves.
#'
#' @param study A [sim_study()].
#' @param label `"pilot"` or `"production"`.
#' @param ... Components to replace, as named arguments of [sim_study()]
#'   (for example `methods = ...`).
#' @param reason Required reason for the amendment.
#' @return The modified study.
#' @examples
#' study <- sim_study("demo", scenarios = sim_scenarios(n = 20L),
#'   generate = function(scenario, context) rnorm(scenario$n),
#'   estimands = list(sim_estimand("mean", sim_truth_exact(function(scenario) 0))),
#'   methods = list(avg = function(data) sim_estimate("mean", mean(data))))
#' frozen <- sim_freeze(study)
#' amended <- sim_amend(frozen, alpha = 0.01, reason = "stricter level agreed")
#' amended$amendments
#' @export
sim_freeze <- function(study, label = c("production", "pilot")) {
  stopifnot(inherits(study, "adempr_study"))
  label <- match.arg(label)
  g <- git_info()
  study$status <- list(label = label,
                       frozen = list(hash = sim_fingerprint(study), time = now_utc(),
                                     git_commit = g$commit, git_dirty = g$dirty))
  study
}

#' @rdname sim_freeze
#' @export
sim_amend <- function(study, ..., reason) {
  stopifnot(inherits(study, "adempr_study"))
  if (missing(reason) || !is_scalar_character(reason)) {
    stopf("An amendment needs a `reason`.")
  }
  changes <- list(...)
  allowed <- c("aims", "scenarios", "generate", "estimands", "methods", "measures",
               "protocol", "failure", "retention", "alpha", "level")
  bad <- setdiff(names(changes), allowed)
  if (!length(changes) || is.null(names(changes)) || length(bad)) {
    stopf("Amend named components among: %s.", paste(allowed, collapse = ", "))
  }
  before <- sim_fingerprint(study)
  args <- list(name = study$name, aims = study$aims, scenarios = study$scenarios,
               generate = study$generator, estimands = unname(study$estimands),
               methods = study$methods, measures = study$measures,
               protocol = study$protocol, failure = study$failure,
               retention = study$retention, alpha = study$alpha, level = study$level)
  args[names(changes)] <- changes
  new <- do.call(sim_study, args)
  new$amendments <- c(study$amendments, list(list(
    time = now_utc(), reason = reason, components = names(changes),
    fingerprint_before = before, fingerprint_after = sim_fingerprint(new),
    after_freeze = !is.null(study$status$frozen))))
  new$status <- study$status
  if (!is.null(study$status$frozen)) {
    new$status$frozen$amended_hash <- sim_fingerprint(new)
  }
  new
}
