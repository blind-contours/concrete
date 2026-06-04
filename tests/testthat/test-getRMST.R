test_that("trapezoidWeights integrate linear functions exactly", {
  w <- concrete:::trapezoidWeights(c(0, 1, 2))
  expect_equal(w, c(0.5, 1, 0.5))
  # integral of f(t) = t over [0, 2] is 2
  expect_equal(sum(w * c(0, 1, 2)), 2)
  # uneven grid: integral of f(t) = 1 over [0, 5] is 5
  grid <- c(0, 1, 3, 5)
  expect_equal(sum(concrete:::trapezoidWeights(grid) * rep(1, 4)), 5)
})

test_that("addWaldInference computes Wald p-values and CIs", {
  dt <- data.table::data.table(
    Estimand = c("Risk Diff", "Rel Risk", "Abs Risk"),
    `Pt Est` = c(0.10, 1.50, 0.30),
    se = c(0.05, 0.25, 0.04)
  )
  out <- concrete:::addWaldInference(dt, Signif = 0.05)
  # Risk Diff: z = 0.10 / 0.05 = 2
  expect_equal(out$pValue[1], 2 * stats::pnorm(-2), tolerance = 1e-8)
  # Rel Risk: null is 1, z = (1.5 - 1) / 0.25 = 2
  expect_equal(out$pValue[2], 2 * stats::pnorm(-2), tolerance = 1e-8)
  # Absolute risk is one-sample: no p-value
  expect_true(is.na(out$pValue[3]))
  expect_equal(out$`CI Low`[1], 0.10 - stats::qnorm(0.975) * 0.05, tolerance = 1e-8)
})

test_that("addWaldInference non-inferiority respects the confidence interval", {
  dt <- data.table::data.table(
    Estimand = "Risk Diff", `Pt Est` = 0.02, se = 0.03
  )
  # upper = smaller is better; CI Hi = 0.02 + 1.96*0.03 = 0.0788
  ni_pass <- concrete:::addWaldInference(dt, NIMargin = 0.10, NIDirection = "upper")
  expect_true(ni_pass$NonInferior)           # CI Hi 0.079 < margin 0.10
  ni_fail <- concrete:::addWaldInference(dt, NIMargin = 0.05, NIDirection = "upper")
  expect_false(ni_fail$NonInferior)          # CI Hi 0.079 > margin 0.05
})

test_that("getRMST is self-consistent and reuses the risk IC", {
  skip_on_cran()
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex)]
  data <- data[stats::complete.cases(data)]
  data[, arm := as.integer(trt == 2)]
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data <- data[, .(id, time, event, arm, age, sex)]

  args <- formatArguments(
    DataTable = data, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = c(730, 1460, 2190),
    TargetEvent = c(1, 2), CVArg = list(V = 2), UpdateMethod = "adaptive",
    EICStopRule = "absolute", MaxUpdateIter = 15, Verbose = FALSE
  )
  est <- suppressMessages(doConcrete(args))
  rmst <- getRMST(est, Horizon = 2190, Intervention = c(1, 2))

  expect_s3_class(rmst, "ConcreteOut")
  # all event types targeted -> event-free RMST is returned
  expect_true(any(rmst$Estimand == "RMST"))

  # RMST (event-free) == Horizon - sum_j LYL_j for each arm
  lyl <- rmst[rmst$Estimand == "Life Years Lost", ]
  for (arm in c("A=0", "A=1")) {
    implied <- 2190 - sum(lyl$`Pt Est`[lyl$Intervention == arm])
    got <- rmst$`Pt Est`[rmst$Estimand == "RMST" & rmst$Intervention == arm]
    expect_equal(got, implied, tolerance = 1e-6)
  }
  # life-years lost are non-negative and bounded by the horizon
  expect_true(all(lyl$`Pt Est` >= 0 & lyl$`Pt Est` <= 2190))
  # standard errors are finite and positive for the targeted estimands
  expect_true(all(is.finite(rmst$se) & rmst$se > 0))
})
