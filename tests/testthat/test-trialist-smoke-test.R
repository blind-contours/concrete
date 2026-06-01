test_that("trialist smoke-test script is packaged and parseable", {
  smoke_script <- system.file("examples", "trialist-smoke-test.R", package = "concrete")
  expect_true(file.exists(smoke_script))
  expect_error(parse(smoke_script), regexp = NA)
})

test_that("trialist smoke-test script runs the default analysis", {
  testthat::skip_on_cran()
  smoke_script <- system.file("examples", "trialist-smoke-test.R", package = "concrete")
  smoke_env <- new.env(parent = globalenv())

  sys.source(smoke_script, envir = smoke_env)

  expect_true(exists("smoke_summary", envir = smoke_env))
  smoke_summary <- get("smoke_summary", envir = smoke_env)
  expect_s3_class(smoke_summary, "data.table")
  expect_true("cox_only" %in% smoke_summary$analysis)
  expect_true(all(smoke_summary$status == "ok"))
})
