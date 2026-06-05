test_that("senseCensoring validates its inputs", {
  expect_error(senseCensoring(list(), deltas = c(0, 1)), "ConcreteArgs")
})

test_that("senseCensoring imputation moves the estimate and is reference-safe", {
  skip_on_cran()
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex, bili)]
  data <- data[stats::complete.cases(data)][1:200]
  data[, arm := as.integer(trt == 2)]
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data <- data[, .(id, time, event, arm, age, sex, bili)]

  args <- formatArguments(
    DataTable = data, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = 3000, TargetEvent = 1,
    CVArg = list(V = 2), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 12, Verbose = FALSE)
  cens_before <- sum(args$DataTable$event == 0 & args$DataTable$time < 3000)

  sc <- suppressMessages(senseCensoring(args, deltas = c(0, 0.5, 1),
                                        Estimand = "Risk", Intervention = c(1, 2)))
  # number imputed grows with delta
  imp <- unique(sc[, c("delta", "n_imputed")])
  expect_equal(imp$n_imputed, c(0, floor(0.5 * cens_before), cens_before))
  # imputing censored as events monotonically increases each arm's risk
  byarm <- split(sc, sc$Event)[["1"]]
  riskA0 <- byarm$`Pt Est`[seq(1, nrow(byarm), by = 2)]   # arm-aligned rows per delta
  expect_true(riskA0[3] > riskA0[1])                       # delta=1 risk > delta=0 risk
  # the caller's args are NOT mutated (data.table reference-safety)
  expect_equal(sum(args$DataTable$event == 0 & args$DataTable$time < 3000), cens_before)
  expect_s3_class(attr(sc, "tippingPoint"), "data.table")
})
