test_that("targetRMST directly targets RMST and is self-consistent", {
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
    ID = "id", Intervention = makeITT(), TargetTime = 2000, TargetEvent = c(1, 2),
    CVArg = list(V = 2), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 10, Verbose = FALSE
  )
  est <- suppressMessages(doConcrete(args))
  direct <- suppressWarnings(
    targetRMST(est, Horizon = 2000, Intervention = c(1, 2), MaxUpdateIter = 80))

  expect_s3_class(direct, "ConcreteOut")
  expect_true(all(c("RMST", "Life Years Lost", "RMST Diff") %in% direct$Estimand))

  # event-free RMST == Horizon - sum_j life-years-lost, per arm
  for (arm in c("A=0", "A=1")) {
    lyl <- direct$`Pt Est`[direct$Estimand == "Life Years Lost" &
                             direct$Intervention == arm]
    rmst <- direct$`Pt Est`[direct$Estimand == "RMST" & direct$Intervention == arm]
    expect_equal(rmst, 2000 - sum(lyl), tolerance = 1e-6)
  }

  # standard errors are finite and positive
  expect_true(all(is.finite(direct$se) & direct$se > 0))

  # direct targeting converges under the default hybrid rule on this example
  expect_true(all(attr(direct, "RMSTConverged")))
})

test_that("direct and pointwise RMST point estimates agree closely", {
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
    ID = "id", Intervention = makeITT(), TargetTime = c(1000, 1500, 2000),
    TargetEvent = c(1, 2), CVArg = list(V = 2), UpdateMethod = "adaptive",
    EICStopRule = "absolute", MaxUpdateIter = 10, Verbose = FALSE
  )
  est <- suppressMessages(doConcrete(args))

  direct <- suppressWarnings(targetRMST(est, Horizon = 2000, Intervention = c(1, 2)))
  pointwise <- getRMST(est, Horizon = 2000, Intervention = c(1, 2))

  d_rmst <- direct$`Pt Est`[direct$Estimand == "RMST" & direct$Intervention == "A=0"]
  p_rmst <- pointwise$`Pt Est`[pointwise$Estimand == "RMST" & pointwise$Intervention == "A=0"]
  # within ~3% on the same horizon (they differ only by integration grid)
  expect_lt(abs(d_rmst - p_rmst) / p_rmst, 0.03)
})
