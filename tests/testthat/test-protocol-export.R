test_that("the protocol is exported in ADEMP-PreReg order with every specified field", {
  st <- normal_study(n = c(10L, 20L), mu = 0,
                     protocol = sim_protocol(purpose = "method_comparison", phase = "II",
                                             dgm_rationale = "normal data | simple",
                                             precision_targets = c(bias = 0.01)))
  f <- tempfile(fileext = ".md")
  sim_protocol_export(st, f, replications = 500)
  md <- readLines(f)
  heads <- grep("^## ", md, value = TRUE)
  expect_equal(heads, c("## 0. General information", "## 1. Aims", "## 2. Data-generating mechanism",
                        "## 3. Estimands and other targets", "## 4. Methods",
                        "## 5. Performance measures", "## 6. Other"))
  expect_true(any(grepl("Phase of methodological research: II", md)))
  expect_true(any(grepl("Planned replications per scenario: 500", md)))
  expect_true(any(grepl("bias = 0.01", md)))
  expect_true(any(grepl("normal data \\| simple", md)))
  expect_true(any(grepl("t.test", md)))
  expect_true(any(grepl("Scope: \\*Not specified", md)))
  expect_error(sim_protocol_export(st, f), "exists")
  sim_protocol_export(st, f, overwrite = TRUE, code = FALSE)
  expect_false(any(grepl("t.test", readLines(f))))
})

test_that("the protocol lists custom measures, labels the levels and gives the whole failure policy", {
  st <- normal_study(n = 10L, mu = 0, level = 0.99,
                     measures = sim_measures("bias", sim_measure_interval_score()),
                     failure = sim_failure_policy(boundary = "invalid", retries = 2L),
                     methods = list(m = sim_method(function(data) list(fit = mean(data$y)),
                       extract = function(fit) sim_estimate("mean", fit$fit),
                       fallback = function(data) sim_estimate("mean", median(data$y)))))
  f <- tempfile(fileext = ".md")
  sim_protocol_export(st, f)
  md <- readLines(f)
  expect_true(any(grepl("\\| interval_score \\|", md)))
  expect_true(any(grepl("Monte Carlo confidence intervals of the performance measures: 0.99", md)))
  expect_false(any(grepl("Nominal level of intervals", md)))
  expect_true(any(grepl("estimates at a boundary invalid", md)))
  expect_true(any(grepl("nonconverged fits invalid", md)))
  expect_true(any(grepl("retries 2", md)))
  expect_true(any(grepl("Method `m`, extract:", md)))
  expect_true(any(grepl("Method `m`, fallback:", md)))
  expect_true(any(grepl("Measure `interval_score`, contribution:", md)))
})
