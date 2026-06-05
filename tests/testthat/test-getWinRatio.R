test_that("getWinRatio returns coherent win statistics", {
  skip_on_cran()
  set.seed(1); n <- 500
  W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  # treatment lowers the hazard -> better survival -> win ratio > 1
  lam <- 6e-4 * exp(-0.6 * A + 0.4 * W + 0.2 * W2)
  Tt <- stats::rexp(n, lam); C <- stats::rexp(n, 3e-4)
  To <- pmin(Tt, C, 3000); ev <- ifelse(To >= 3000, 0L, ifelse(Tt <= C, 1L, 0L))
  d <- data.table::data.table(id = 1:n, time = To, event = ev, arm = A, W = W, W2 = W2)
  args <- formatArguments(
    DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm", ID = "id",
    Intervention = makeITT(), TargetTime = c(300, 600, 900, 1200, 1500), TargetEvent = 1,
    CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 20, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  wr <- getWinRatio(est, Horizon = 1500, Intervention = c(1, 2))

  expect_s3_class(wr, "ConcreteOut")
  expect_true(all(c("Win Ratio", "Win Odds", "Net Benefit",
                    "P(win)", "P(loss)", "P(tie)") %in% wr$Estimand))
  # win/loss/tie probabilities sum to 1
  probs <- wr[wr$Estimand %in% c("P(win)", "P(loss)", "P(tie)"), `Pt Est`]
  expect_equal(sum(probs), 1, tolerance = 1e-8)
  # all probabilities in [0, 1]
  expect_true(all(probs >= -1e-8 & probs <= 1 + 1e-8))
  # standard errors finite and positive for the comparative statistics
  comp <- wr[wr$Estimand %in% c("Win Ratio", "Win Odds", "Net Benefit"), ]
  expect_true(all(is.finite(comp$se) & comp$se > 0))
  # treatment improves survival here -> win ratio above 1, net benefit above 0
  expect_gt(wr$`Pt Est`[wr$Estimand == "Win Ratio"], 1)
  expect_gt(wr$`Pt Est`[wr$Estimand == "Net Benefit"], 0)
  # win ratio = P(win)/P(loss) internal consistency
  pw <- wr$`Pt Est`[wr$Estimand == "P(win)"]; pl <- wr$`Pt Est`[wr$Estimand == "P(loss)"]
  expect_equal(wr$`Pt Est`[wr$Estimand == "Win Ratio"], pw / pl, tolerance = 1e-8)
})

test_that("getWinRatio gives a win ratio near 1 when the arms are exchangeable", {
  skip_on_cran()
  set.seed(2); n <- 600
  W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  lam <- 6e-4 * exp(0.4 * W + 0.2 * W2)         # no treatment effect
  Tt <- stats::rexp(n, lam); C <- stats::rexp(n, 3e-4)
  To <- pmin(Tt, C, 3000); ev <- ifelse(To >= 3000, 0L, ifelse(Tt <= C, 1L, 0L))
  d <- data.table::data.table(id = 1:n, time = To, event = ev, arm = A, W = W, W2 = W2)
  args <- formatArguments(
    DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm", ID = "id",
    Intervention = makeITT(), TargetTime = c(300, 600, 900, 1200, 1500), TargetEvent = 1,
    CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 20, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  wr <- getWinRatio(est, Horizon = 1500, Intervention = c(1, 2))
  # CI for the win ratio should contain 1 under the null
  r <- wr[wr$Estimand == "Win Ratio", ]
  expect_lt(r$`CI Low`, 1); expect_gt(r$`CI Hi`, 1)
})
