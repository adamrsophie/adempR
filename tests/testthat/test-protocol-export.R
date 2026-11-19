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
