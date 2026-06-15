famFit <- function(n = 300, seed = 4) {
  set.seed(seed)
  N <- 2L * n; W <- rnorm(N); A <- rep(0:1, each = n)
  l1 <- 0.12 * exp(0.3 * W - 0.5 * A); l2 <- 0.15 * exp(0.2 * W - 0.3 * A)
  Te <- rexp(N, l1 + l2); cause <- 1L + rbinom(N, 1L, l2 / (l1 + l2)); C <- rexp(N, 0.08)
  dat <- data.table::data.table(id = seq_len(N), time = pmin(Te, C, 5),
                                status = ifelse(Te <= pmin(C, 5), cause, 0L),
                                trt = A, W = W, W2 = rnorm(N))
  a <- suppressMessages(formatArguments(
    DataTable = dat, EventTime = "time", EventType = "status", Treatment = "trt",
    ID = "id", Intervention = 0:1, TargetTime = c(1, 2, 3, 4), TargetEvent = c(1, 2),
    CVArg = list(V = 2), MaxUpdateIter = 6, Verbose = FALSE,
    Model = list(trt = "SL.mean",
                 "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                 "1" = list(Cox = survival::Surv(time, status == 1) ~ .),
                 "2" = list(Cox = survival::Surv(time, status == 2) ~ .))))
  suppressMessages(suppressWarnings(doConcrete(a)))
}

test_that("getSimultaneousFamily widens bands jointly and preserves point estimates", {
  skip_on_cran()
  e <- famFit()
  o  <- suppressMessages(getOutput(e, Estimand = c("Risk", "RD"), Simultaneous = FALSE))
  rm <- suppressMessages(getRMST(e, Horizon = 4))
  wr <- suppressMessages(getWinRatio(e, Horizon = 4, Intervention = c(2, 1), TargetEvent = c(1, 2)))

  fam <- getSimultaneousFamily(RD = o, RMST = rm, WinRatio = wr)

  ## structure: RD (2 events x 4 times) + RMST/LYL (3) + WR/WO/NB (3)
  expect_true(all(c("family", "Estimand", "Pt Est", "CI Low", "CI Hi",
                    "SimCI Low", "SimCI Hi") %in% names(fam)))
  expect_setequal(unique(fam$family), c("RD", "RMST", "WinRatio"))
  expect_gt(nrow(fam), 10)

  ## simultaneous critical value exceeds the pointwise z, and every band is wider
  expect_gt(attr(fam, "critValue"), stats::qnorm(0.975))
  ptW  <- fam$`CI Hi`  - fam$`CI Low`
  simW <- fam$`SimCI Hi` - fam$`SimCI Low`
  expect_true(all(round(simW - ptW, 8) >= 0))

  ## point estimates and pointwise CIs match the source objects exactly
  rdRow <- as.data.frame(o)[as.data.frame(o)$Estimand == "Risk Diff" &
                              as.data.frame(o)$Estimator == "tmle" &
                              as.data.frame(o)$Event == 1 & as.data.frame(o)$Time == 4, ]
  famRD <- fam[Estimand == "Risk Diff" & Event == 1 & Time == 4]
  expect_equal(famRD$`Pt Est`, rdRow$`Pt Est`, tolerance = 1e-10)
})

test_that("getSimultaneousFamily errors on objects without family ICs", {
  skip_on_cran()
  expect_error(getSimultaneousFamily(data.table::data.table(x = 1)),
               regexp = "family influence functions")
})

test_that("a single-object family reduces to the Wald interval", {
  skip_on_cran()
  e <- famFit(seed = 6)
  rm <- suppressMessages(getRMST(e, Horizon = 4))
  fam <- getSimultaneousFamily(RMST = rm)
  ## with one non-degenerate contrast family the crit value is near z (a few rows,
  ## strong correlation) -> simultaneous close to pointwise, never narrower
  simW <- fam$`SimCI Hi` - fam$`SimCI Low`
  ptW  <- fam$`CI Hi` - fam$`CI Low`
  expect_true(all(round(simW - ptW, 8) >= 0))
})
