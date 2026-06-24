test_that("CensoringTV threads through formatArguments/doConcrete without breaking the baseline path", {
  skip_on_cran()
  skip_if_not_installed("SuperLearner")
  data <- data.table::as.data.table(survival::pbc)[1:150, c("time","status","trt","id","age","sex")]
  set.seed(1); data[, trt := sample(0:1, .N, TRUE)]
  tv <- data.table::rbindlist(lapply(c(400, 1200), function(v)
    data.table::data.table(id = data$id, time = v, lab = stats::rnorm(nrow(data)))))
  fa <- function(ctv) suppressMessages(formatArguments(DataTable = data, EventTime = "time",
    EventType = "status", Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = 2500,
    TargetEvent = NULL, CVArg = list(V = 2), MaxUpdateIter = 2, Model = NULL, CensoringTV = ctv))
  risk <- function(e) {
    o <- suppressMessages(getOutput(e)); o <- o[o$Estimator == "tmle" & o$Estimand == "Abs Risk", ]
    o[["Pt Est"]][1]
  }
  r0 <- expect_no_error(suppressMessages(doConcrete(fa(NULL))))
  rtv <- expect_no_error(suppressMessages(doConcrete(fa(tv))))
  expect_true(is.finite(risk(r0)) && is.finite(risk(rtv)))
  # tv overrides the IPCW -> the targeted risk should change
  expect_false(isTRUE(all.equal(risk(r0), risk(rtv))))
  # input validation: CensoringTV missing required columns
  expect_error(suppressMessages(doConcrete(fa(tv[, .(id, lab)]))))  # no 'time' column
})

test_that(".tvLOCF carries last observation forward and computes change-from-baseline", {
  ids <- 1:3
  tv <- data.frame(pid = c(1,1,2), time = c(1, 3, 2), x = c(10, 14, 20))
  m <- concrete:::.tvLOCF(ids, tv, "pid", "time", starts = c(0, 2, 4))
  expect_equal(m[["x_val"]][1, ], c(10, 10, 14))   # subj1: baseline before t=1, LOCF after
  expect_equal(m[["x_chg"]][1, ], c(0, 0, 4))       # change from baseline (10)
  expect_equal(m[["x_val"]][3, ], rep(stats::median(c(10,14,20)), 3))  # subj3 absent -> median fallback (not 0)
  expect_false(anyNA(m[["x_val"]]))                  # design matrix must be complete
})

test_that(".tvLOCF imputes missing measurements (drop NA + LOCF + median fallback) and flags missingness", {
  ids <- 1:3
  tv <- data.frame(pid = c(1, 1, 2), time = c(1, 3, 2), x = c(10, NA, 20))  # subj1 NA at t=3; subj3 absent
  m <- concrete:::.tvLOCF(ids, tv, "pid", "time", starts = c(0, 2, 4))
  expect_false(anyNA(m[["x_val"]]))                  # NA measurement must not propagate
  expect_equal(m[["x_val"]][1, ], c(10, 10, 10))     # NA at t=3 dropped -> LOCF the t=1 value
  expect_true("x_miss" %in% names(m))                # missingness indicator added when it varies
  expect_equal(m[["x_miss"]][3, ], c(1L, 1L, 1L))    # absent subject flagged imputed throughout
})

test_that(".tvCensoringInc does not pull post-horizon censoring into the last interval", {
  times <- c(0, 1, 2)
  base <- data.frame(A = c(0, 1, 0, 1), W = c(-1, 0, 1, 2))
  tv <- list()

  after <- concrete:::.tvCensoringInc(
    times = times,
    obsT = c(3, 3, 3, 3),
    censInd = c(1L, 1L, 1L, 1L),
    baseCov = base,
    tvMats = tv,
    SL.library = "SL.mean",
    n.folds = 1L
  )
  expect_lt(max(after), 1e-8)

  before <- concrete:::.tvCensoringInc(
    times = times,
    obsT = c(1.5, 1.5, 1.5, 1.5),
    censInd = c(1L, 1L, 1L, 1L),
    baseCov = base,
    tvMats = tv,
    SL.library = "SL.mean",
    n.folds = 1L
  )
  expect_gt(max(before), 0.1)
  expect_gt(max(before), 1e6 * max(after))
})
