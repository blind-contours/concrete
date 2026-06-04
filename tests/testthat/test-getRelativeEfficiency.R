test_that("getRelativeEfficiency validates inputs", {
  expect_error(getRelativeEfficiency(data.frame(a = 1), data.frame(a = 1)),
               "ConcreteOut")
})

test_that("getRelativeEfficiency computes the variance ratio", {
  mk <- function(se) {
    dt <- data.table::data.table(
      Intervention = "[A=1] - [A=0]", Estimand = "Risk Diff", Estimator = "tmle",
      Event = 1, Time = 365, `Pt Est` = 0.05, se = se
    )
    class(dt) <- union("ConcreteOut", class(dt))
    dt
  }
  adj <- mk(0.04)     # adjusted is more precise
  unadj <- mk(0.05)
  re <- getRelativeEfficiency(adj, unadj)
  expect_equal(re$RelEfficiency, (0.05^2) / (0.04^2), tolerance = 1e-8)
  expect_equal(re$VarReductionPct, 100 * (1 - 1 / re$RelEfficiency), tolerance = 1e-8)
  expect_gt(re$RelEfficiency, 1)   # adjustment helped
})
