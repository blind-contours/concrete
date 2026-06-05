test_that("cross-fitting predicts each subject out-of-fold", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex, bili)]
  data <- data[stats::complete.cases(data)][1:160]
  data[, arm := as.integer(trt == 2)]
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data <- data[, .(id, time, event, arm, age, sex, bili)]

  Model <- list(arm = "SL.glm",
                "0" = list(Cox = survival::Surv(time, event == 0) ~ arm + age + sex + bili),
                "1" = list(Cox = survival::Surv(time, event == 1) ~ arm + age + sex + bili),
                "2" = list(Cox = survival::Surv(time, event == 2) ~ arm + age + sex + bili))
  args <- formatArguments(
    DataTable = data, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = c(1000, 2000), TargetEvent = c(1, 2),
    CVArg = list(V = 2), Model = Model, CrossFit = TRUE, RenameCovs = FALSE, Verbose = FALSE)

  Data <- args$DataTable
  Regime <- concrete:::getRegime(args$Intervention, Data)
  CVFolds <- args$CVFolds
  est <- suppressMessages(concrete:::getCVInitialEstimate(
    Data = Data, Model = args$Model, CVFolds = CVFolds, MinNuisance = args$MinNuisance,
    TargetEvent = args$TargetEvent, TargetTime = args$TargetTime, Regime = Regime,
    ReturnModels = FALSE))
  cf_prop <- est[["A=1"]][["PropScore"]]   # P(A=1|W) for the treat-all regime

  for (fold in seq_along(CVFolds)) {
    va <- CVFolds[[fold]][["validation_set"]]
    tr <- CVFolds[[fold]][["training_set"]]
    g <- stats::glm(arm ~ age + sex + bili, data = as.data.frame(Data)[tr, ],
                    family = stats::binomial)
    p <- stats::predict(g, newdata = as.data.frame(Data)[va, ], type = "response")
    # each held-out subject's nuisance comes from a model fit without them
    expect_equal(as.numeric(cf_prop[va]), as.numeric(p), tolerance = 1e-8)
  }
  # and it is genuinely out-of-fold: differs from a full-sample fit
  gfull <- stats::glm(arm ~ age + sex + bili, data = as.data.frame(Data), family = stats::binomial)
  expect_gt(max(abs(cf_prop - stats::predict(gfull, type = "response"))), 0.01)
})

test_that("doConcrete runs with CrossFit = TRUE and gives finite inference", {
  skip_on_cran()
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex)]
  data <- data[stats::complete.cases(data)][1:160]
  data[, arm := as.integer(trt == 2)]
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data <- data[, .(id, time, event, arm, age, sex)]
  args <- formatArguments(
    DataTable = data, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = 2000, TargetEvent = 1,
    CVArg = list(V = 3), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 15, CrossFit = TRUE, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  out <- getOutput(est, Estimand = c("Risk", "RD"), Intervention = c(1, 2), Simultaneous = FALSE)
  tmle <- out[out$Estimator == "tmle", ]
  expect_true(all(is.finite(tmle$`Pt Est`)))
  expect_true(all(is.finite(tmle$se) & tmle$se > 0))
})

test_that("print.ConcreteEst does not error on a cross-fitted object", {
  skip_on_cran()
  data <- data.table::as.data.table(survival::pbc)
  data <- data[!is.na(trt), .(id, time, status, trt, age, sex)]
  data <- data[stats::complete.cases(data)][1:140]
  data[, arm := as.integer(trt == 2)]
  data[, event := data.table::fifelse(status == 2L, 1L,
          data.table::fifelse(status == 1L, 2L, 0L))]
  data <- data[, .(id, time, event, arm, age, sex)]
  args <- formatArguments(
    DataTable = data, EventTime = "time", EventType = "event", Treatment = "arm",
    ID = "id", Intervention = makeITT(), TargetTime = 2000, TargetEvent = 1,
    CVArg = list(V = 3), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 8, CrossFit = TRUE, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  expect_error(capture.output(print(est)), NA)   # must not error
})
