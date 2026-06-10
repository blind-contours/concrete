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
  expect_equal(m[["x_val"]][3, ], c(0, 0, 0))        # subj3 absent -> 0
})
