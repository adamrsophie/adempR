# Regression specifications for the reliability contract of commit
# 814be78113bf0c52dfd7217550d75c2b881779d2.
# These tests describe required behavior and intentionally fail that commit.
# Run against the installed package; no repository files are changed.

.contract_rel_study <- function(generate, methods, n = 5L,
                             truth = function(scenario) 0) {
  adempR::sim_study(
    "contract_reliability", scenarios = adempR::sim_scenarios(n = n),
    generate = generate,
    estimands = list(adempR::sim_estimand("mean", adempR::sim_truth_exact(truth))),
    methods = methods, measures = "bias"
  )
}

.contract_rel_capture <- function(expr) {
  tryCatch(force(expr), error = function(e) e)
}

.contract_rel_error_info <- function(x) {
  if (inherits(x, "error")) conditionMessage(x) else ""
}

testthat::test_that("F1 failed regeneration cannot invalidate committed datasets", {
  testthat::skip_if_not_installed("adempR")
  withr::local_options(list(adempR.git_dir = tempdir()))
  flag <- tempfile("contract-rel-source-")
  writeLines("available", flag)
  withr::defer(unlink(flag), envir = environment())
  path <- tempfile("contract-rel-regeneration-")
  withr::defer(unlink(path, recursive = TRUE), envir = environment())

  generate <- local({
    source_file <- flag
    function(scenario, context) {
      if (!file.exists(source_file)) stop("data source temporarily unavailable")
      stats::rnorm(scenario$n)
    }
  })
  study <- .contract_rel_study(generate,
    list(a = function(data) adempR::sim_estimate("mean", mean(data))))
  original <- adempR::sim_run(study, 2, 1, path, progress = FALSE)
  reference <- adempR::sim_results(original)$estimates$estimate
  extended_study <- study
  extended_study$methods$b <- adempR::sim_method(
    function(data) adempR::sim_estimate("mean", stats::median(data)),
    estimands = "mean")

  unlink(flag)
  # A controlled extension error is permitted; corrupting already committed
  # scientific results is not. The extension's return value is immaterial.
  invisible(.contract_rel_capture(adempR::sim_run(
    extended_study, 2, 1, path, progress = FALSE)))
  retained <- .contract_rel_capture(adempR::sim_results(original))
  testthat::expect_false(inherits(retained, "error"),
    info = paste("Original results must remain readable:", .contract_rel_error_info(retained)))
  if (!inherits(retained, "error")) {
    testthat::expect_equal(retained$estimates$estimate, reference)
    testthat::expect_true(all(retained$replications$generation_status == "ok"))
  }

  writeLines("available again", flag)
  resumed <- .contract_rel_capture({
    rr <- adempR::sim_run(extended_study, 2, 1, path, progress = FALSE,
                          retry_failed = TRUE)
    adempR::sim_results(rr)
  })
  testthat::expect_false(inherits(resumed, "error"),
    info = paste("Restoring the source must allow recovery:", .contract_rel_error_info(resumed)))
  if (!inherits(resumed, "error")) {
    testthat::expect_equal(nrow(resumed$estimates), 4L)
    testthat::expect_equal(resumed$estimates$estimate[resumed$estimates$method == "a"],
                          reference)
  }
})

testthat::test_that("F2 collected runs allocate retry IDs above all imported rounds", {
  testthat::skip_if_not_installed("adempR")
  withr::local_options(list(adempR.git_dir = tempdir()))
  flag <- tempfile("contract-rel-service-")
  writeLines("unavailable", flag)
  withr::defer(unlink(flag), envir = environment())
  parent <- tempfile("contract-rel-shards-")
  into <- tempfile("contract-rel-collected-")
  withr::defer(unlink(c(parent, into), recursive = TRUE), envir = environment())
  method <- local({
    service_flag <- flag
    function(data) {
      if (file.exists(service_flag)) stop("transient service unavailable")
      adempR::sim_estimate("mean", mean(data))
    }
  })
  study <- .contract_rel_study(
    function(scenario, context) stats::rnorm(scenario$n),
    list(a = method), n = c(5L, 6L))
  for (k in 1:2) {
    adempR::sim_run_shard(study, k, 2, 2, 1, parent, progress = FALSE)
    adempR::sim_run_shard(study, k, 2, 2, 1, parent,
                          retry_failed = TRUE, progress = FALSE)
  }
  collected <- adempR::sim_collect_shards(study, 2, 2, 1, parent, into)
  imported <- adempR::sim_results(collected)
  imported_max <- max(imported$attempts_all$round)
  testthat::expect_equal(imported_max, 1L)

  unlink(flag)
  retried <- .contract_rel_capture(adempR::sim_resume(
    into, retry_failed = TRUE, progress = FALSE))
  testthat::expect_false(inherits(retried, "error"),
    info = .contract_rel_error_info(retried))
  if (!inherits(retried, "error")) {
    testthat::expect_true(isTRUE(retried$manifest$retry_round > imported_max),
      info = "A retry after collection must not reuse an imported round identifier.")
    current <- .contract_rel_capture(adempR::sim_results(retried))
    testthat::expect_false(inherits(current, "error"),
      info = paste("Successful retried results must reduce without conflicts:",
                   .contract_rel_error_info(current)))
    if (!inherits(current, "error")) {
      testthat::expect_true(all(current$attempts$status == "ok"))
      testthat::expect_true(all(current$attempts$round > imported_max))
      testthat::expect_equal(nrow(current$estimates), 4L)
      testthat::expect_true(any(current$attempts_all$round == imported_max))
    }
  }
})

testthat::test_that("F3 truth cache recovery cannot change truth under its old identity", {
  testthat::skip_if_not_installed("adempR")
  withr::local_options(list(adempR.git_dir = tempdir()))
  # A unique temporary global name allows a directly referenced global
  # without overwriting any pre-existing object in the caller's session.
  global_name <- make.names(basename(tempfile("contract_rel_truth_")))
  assign(global_name, 0, envir = globalenv())
  withr::defer(rm(list = global_name, envir = globalenv()), envir = environment())
  truth <- eval(substitute(function(scenario) VALUE,
                            list(VALUE = as.name(global_name))), globalenv())
  study <- .contract_rel_study(
    function(scenario, context) rep(0, scenario$n),
    list(a = function(data) adempR::sim_estimate("mean", 0)), truth = truth)
  path <- tempfile("contract-rel-truth-")
  withr::defer(unlink(path, recursive = TRUE), envir = environment())
  run <- adempR::sim_run(study, 2, 1, path, progress = FALSE)
  original <- adempR::sim_truth(run)
  testthat::expect_equal(original$truth, 0)

  assign(global_name, 100, envir = globalenv())
  unlink(list.files(file.path(path, "truth"), full.names = TRUE), recursive = TRUE)
  recovered <- .contract_rel_capture(adempR::sim_truth(run))
  # Both fixes are acceptable: reject recomputation with changed dependencies,
  # or recover a frozen original definition and obtain the original truth.
  testthat::expect_true(inherits(recovered, "error") ||
    isTRUE(all.equal(recovered$truth, original$truth)),
    info = "Cache recovery must refuse changed truth or restore the original value; it must not publish 100 under the identity of 0.")
  if (!inherits(recovered, "error") && isTRUE(all.equal(recovered$truth, original$truth))) {
    testthat::expect_identical(recovered$truth_fp, original$truth_fp)
  }
})

testthat::test_that("F4 changed generator package resolution cannot pool distinct DGPs", {
  testthat::skip_if_not_installed("adempR")
  testthat::skip_if_not_installed("xfun")
  withr::local_options(list(adempR.git_dir = tempdir()))
  # Use an explicit clean parent for the generator so a test runner's
  # private bindings cannot mask the declared package's exported attr().
  generate <- eval(quote(function(scenario, context) {
    as.numeric(length(attr(structure(1, value = 123), "val")))
  }), globalenv())
  study <- .contract_rel_study(adempR::sim_generator(generate, packages = "base"),
    list(a = function(data) adempR::sim_estimate("mean", mean(data))))
  path <- tempfile("contract-rel-package-")
  withr::defer(unlink(path, recursive = TRUE), envir = environment())
  run <- adempR::sim_run(study, 2, 1, path, progress = FALSE)
  original_values <- adempR::sim_results(run)$estimates$estimate
  testthat::expect_equal(original_values, c(1, 1))
  changed <- study
  changed$generator$packages <- "xfun"
  altered <- .contract_rel_capture(adempR::sim_run(changed, 4, 1, path, progress = FALSE))

  # Refusing the extension is the simplest correct policy. A future design
  # may instead introduce a distinct DGP identity, but must never pool
  # changed-generator data under the old identity.
  invariant <- inherits(altered, "error")
  if (!invariant) {
    result <- .contract_rel_capture(adempR::sim_results(altered))
    if (!inherits(result, "error")) {
      invariant <- !identical(altered$manifest$fingerprints$generator,
                              run$manifest$fingerprints$generator) &&
        isTRUE(all(result$replications$generator_fp ==
                     altered$manifest$fingerprints$generator)) &&
        isTRUE(all.equal(result$estimates$estimate, rep(0, 4)))
    }
  }
  testthat::expect_true(invariant,
    info = "A change in package-resolved generator behavior must refuse extension or create a distinct DGP identity; 1,1,0,0 under one identity is invalid.")
})
