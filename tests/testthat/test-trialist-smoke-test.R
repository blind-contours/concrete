test_that("trialist smoke-test script is packaged and parseable", {
  smoke_script <- system.file("examples", "trialist-smoke-test.R", package = "concrete")
  expect_true(file.exists(smoke_script))
  expect_error(parse(smoke_script), regexp = NA)
})
