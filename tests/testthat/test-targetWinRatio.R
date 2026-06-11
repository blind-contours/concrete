simWRData <- function(n = 300, seed = 9) {
  set.seed(seed)
  N <- 2L * n
  W <- rnorm(N); A <- rep(0:1, each = n)
  l1 <- 0.10 * exp(0.3 * W - 0.5 * A); l2 <- 0.15 * exp(0.2 * W - 0.3 * A)
  Te <- rexp(N, l1 + l2); cause <- 1L + rbinom(N, 1L, l2 / (l1 + l2))
  C <- rexp(N, 0.08)
  data.table::data.table(
    id = seq_len(N), time = pmin(Te, C, 5),
    status = ifelse(Te <= pmin(C, 5), cause, 0L),
    trt = A, W = W, W2 = rnorm(N))
}

fitWREst <- function(dat, TT = seq(0.5, 4, by = 0.5)) {
  a <- suppressMessages(formatArguments(
    DataTable = dat, EventTime = "time", EventType = "status", Treatment = "trt",
    ID = "id", Intervention = 0:1, TargetTime = TT, TargetEvent = c(1, 2),
    CVArg = list(V = 2), MaxUpdateIter = 8, Verbose = FALSE,
    Model = list(trt = "SL.mean",
                 "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                 "1" = list(Cox = survival::Surv(time, status == 1) ~ .),
                 "2" = list(Cox = survival::Surv(time, status == 2) ~ .))))
  suppressMessages(suppressWarnings(doConcrete(a)))
}

test_that("targetWinRatio agrees with the plug-in on a dense grid", {
  skip_on_cran()
  e <- fitWREst(simWRData())
  i1 <- match("A=1", names(e)); i0 <- match("A=0", names(e))
  pl <- suppressMessages(getWinRatio(e, Horizon = 4, Intervention = c(i1, i0),
                                     TargetEvent = c(1, 2)))
  dr <- suppressMessages(targetWinRatio(e, Horizon = 4, Intervention = c(i1, i0),
                                        TargetEvent = c(1, 2)))
  expect_true(isTRUE(attr(dr, "WRConverged")))
  expect_identical(attr(dr, "Targeting"), "direct")
  expect_setequal(dr$Estimand,
                  c("Win Ratio", "Win Odds", "Net Benefit", "P(win)", "P(loss)", "P(tie)"))
  g <- function(o, est) as.data.frame(o)[as.data.frame(o)$Estimand == est, "Pt Est"][1]
  ## dense grid: the two estimators should be close (grid refinement + targeting)
  expect_lt(abs(g(dr, "Win Ratio") - g(pl, "Win Ratio")), 0.15)
  expect_lt(abs(g(dr, "Net Benefit") - g(pl, "Net Benefit")), 0.05)
  ## probabilities are coherent
  expect_equal(g(dr, "P(win)") + g(dr, "P(loss)") + g(dr, "P(tie)"), 1, tolerance = 1e-8)
  expect_true(all(as.data.frame(dr)$se > 0))
})

test_that("targetWinRatio solves its estimating equations", {
  skip_on_cran()
  e <- fitWREst(simWRData(seed = 10), TT = c(2, 4))   # sparse grid: targeting must work
  i1 <- match("A=1", names(e)); i0 <- match("A=0", names(e))
  dr <- suppressMessages(suppressWarnings(
    targetWinRatio(e, Horizon = 4, Intervention = c(i1, i0), TargetEvent = c(1, 2))))
  expect_true(isTRUE(attr(dr, "WRConverged")))
  wr <- as.data.frame(dr)[as.data.frame(dr)$Estimand == "Win Ratio", ]
  expect_true(wr[["CI Low"]] > 0 && wr[["CI Low"]] < wr[["Pt Est"]],
              info = "log-scale CI is ordered and positive")
})

test_that("targetWinRatio validates inputs", {
  expect_error(targetWinRatio(list()), regexp = "ConcreteEst")
  e <- fitWREst(simWRData(n = 150, seed = 11), TT = c(2, 4))
  expect_error(targetWinRatio(e, Intervention = 1), regexp = "treatment and control")
})
