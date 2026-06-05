test_that("breslowBaseHazIncrements reduces to Nelson-Aalen when eta is constant", {
  # with no covariate effect (eta = 0), the Breslow baseline equals the
  # Nelson-Aalen estimator: increment at each event time = d / (# at risk)
  time <- c(2, 4, 4, 6, 8, 8, 10)
  eventj <- c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE)
  Hazards <- data.table::data.table(Time = c(0, sort(unique(time))))
  bh <- concrete:::breslowBaseHazIncrements(rep(0, length(time)), time, eventj, Hazards)
  # manual Nelson-Aalen at each event time of type j
  na_inc <- function(t) sum(time == t & eventj) / sum(time >= t)
  for (t in c(2, 4, 6, 8)) {
    expect_equal(bh$BaseHaz[bh$Time == t], na_inc(t), tolerance = 1e-10)
  }
})

test_that("Coxnet and Cox give the same covariate-adjusted estimate (baseline fixed)", {
  skip_on_cran()
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 500
  W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  lam <- 6e-4 * exp(-0.6 * A + 0.5 * W + 0.4 * W2)
  Tt <- stats::rexp(n, lam); C <- stats::rexp(n, 3e-4)
  To <- pmin(Tt, C, 2500); ev <- ifelse(To >= 2500, 0L, ifelse(Tt <= C, 1L, 0L))
  d <- data.table::data.table(id = 1:n, time = To, event = ev, arm = A, W = W, W2 = W2)
  mk <- function(M) formatArguments(
    DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm", ID = "id",
    Intervention = makeITT(), TargetTime = 1200, TargetEvent = 1, CVArg = list(V = 5),
    Model = M, UpdateMethod = "adaptive", EICStopRule = "absolute", MaxUpdateIter = 20,
    RenameCovs = FALSE, Verbose = FALSE)
  rd <- function(M) {
    o <- as.data.table(getOutput(suppressMessages(doConcrete(mk(M))),
                                 Estimand = "RD", Intervention = c(1, 2), Simultaneous = FALSE))
    o[o$Estimator == "tmle" & o$Estimand == "Risk Diff", `Pt Est`]
  }
  rd_cox <- rd(list(arm = "SL.glm",
                    "0" = list(Cox = survival::Surv(time, event == 0) ~ arm + W + W2),
                    "1" = list(Cox = survival::Surv(time, event == 1) ~ arm + W + W2)))
  rd_coxnet <- rd(list(arm = "SL.glm", "0" = list(coxnet = "coxnet"), "1" = list(coxnet = "coxnet")))
  # with the baseline fixed, a lightly-penalized coxnet reconstructs essentially
  # the same conditional hazards as the Cox model -> same marginal estimate
  expect_equal(rd_coxnet, rd_cox, tolerance = 0.02)
})
