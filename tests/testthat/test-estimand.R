test_that("makeEstimand records the five ICH E9(R1) attributes", {
  e <- makeEstimand(treatment = "active vs control", strategy = "hypothetical",
                    summary = "risk difference", intercurrent = "discontinuation")
  expect_s3_class(e, "ConcreteEstimand")
  expect_identical(e$`intercurrent event strategy`, "hypothetical")
  expect_identical(e$`population-level summary`, "risk difference")
  expect_output(print(e), "ICH E9")
})

test_that("applyIntercurrentEvent recodes the event type per strategy", {
  d <- data.table::data.table(
    id = 1:6, time = c(5, 3, 8, 2, 6, 4),
    event = c(1L, 2L, 0L, 9L, 1L, 9L), arm = c(1L, 1L, 0L, 0L, 1L, 0L))

  h <- applyIntercurrentEvent(d, "time", "event", Intercurrent = 9,
                              strategy = "hypothetical", Verbose = FALSE)
  expect_equal(h$event, c(1L, 2L, 0L, 0L, 1L, 0L))   # 9 -> 0 (censoring)

  cc <- applyIntercurrentEvent(d, "time", "event", Intercurrent = 9,
                               strategy = "composite", TargetEvent = 1, Verbose = FALSE)
  expect_equal(cc$event, c(1L, 2L, 0L, 1L, 1L, 1L))  # 9 -> 1 (event of interest)

  tp <- applyIntercurrentEvent(d, "time", "event", Intercurrent = 9,
                               strategy = "treatment policy", Verbose = FALSE)
  expect_equal(tp$event, d$event)                    # unchanged

  expect_s3_class(attr(h, "Estimand"), "ConcreteEstimand")
  expect_warning(
    applyIntercurrentEvent(d, "time", "event", Intercurrent = 99,
                           strategy = "hypothetical", Verbose = FALSE),
    "changed nothing")
})

test_that("doConcrete runs on a hypothetical-strategy recoded dataset", {
  skip_on_cran()
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex)]
  data <- data[stats::complete.cases(data)][1:150]
  data[, arm := as.integer(trt == 2)]
  # 1 = death, 2 = competing, 9 = intercurrent event (treatment discontinuation)
  set.seed(1)
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data[sample(.N, 20), event := 9L]
  data <- data[, .(id, time, event, arm, age, sex)]

  hyp <- applyIntercurrentEvent(data, "time", "event", Intercurrent = 9,
                                strategy = "hypothetical", Verbose = FALSE)
  args <- formatArguments(
    DataTable = hyp, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = 2000, TargetEvent = 1,
    CVArg = list(V = 2), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 10, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  out <- getOutput(est, Estimand = "RD", Intervention = c(1, 2), Simultaneous = FALSE)
  expect_true(all(is.finite(out[out$Estimator == "tmle", se])))
})
