test_that("ensembleHazWeights returns a valid simplex weighting and favors lower loss", {
  set.seed(1)
  n <- 200; M <- 3
  # candidate 2 has uniformly lower cumulative hazard at observed times and higher
  # hazard at events -> lower NLL -> should get the most weight
  cumMat <- matrix(abs(rnorm(n * M, 1, 0.2)), n, M)
  cumMat[, 2] <- cumMat[, 2] * 0.5
  hazMat <- matrix(abs(rnorm(n * M, 0.1, 0.02)), n, M)
  hazMat[, 2] <- hazMat[, 2] * 2
  eventMask <- rep(c(TRUE, FALSE), length.out = n)
  w <- concrete:::ensembleHazWeights(cumMat, hazMat, eventMask)
  expect_length(w, M)
  expect_true(all(w >= -1e-8))
  expect_equal(sum(w), 1, tolerance = 1e-6)
  expect_equal(which.max(w), 2L)
})

test_that("the Cox hazard learner uses its own baseline so the correct model is selected", {
  skip_on_cran()
  set.seed(3); n <- 400
  W1 <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  l1 <- 5e-4 * exp(0.5 * A + 0.5 * W1 + 0.4 * W2); l2 <- 3e-4 * exp(-0.2 * A); lc <- 4e-4
  T1 <- stats::rexp(n, l1); T2 <- stats::rexp(n, l2); C <- stats::rexp(n, lc)
  To <- pmin(T1, T2, C, 2500)
  ev <- ifelse(To >= 2500, 0L, ifelse(T1 <= pmin(T2, C), 1L, ifelse(T2 <= C, 2L, 0L)))
  d <- data.table::data.table(id = 1:n, time = To, event = ev, arm = A, W1 = W1, W2 = W2)
  M <- list(arm = "SL.glm",
    "0" = list(trtonly = survival::Surv(time, event == 0) ~ arm,
               full = survival::Surv(time, event == 0) ~ arm + W1 + W2),
    "1" = list(trtonly = survival::Surv(time, event == 1) ~ arm,
               full = survival::Surv(time, event == 1) ~ arm + W1 + W2),
    "2" = list(trtonly = survival::Surv(time, event == 2) ~ arm,
               full = survival::Surv(time, event == 2) ~ arm + W1 + W2))
  args <- formatArguments(
    DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm", ID = "id",
    Intervention = makeITT(), TargetTime = 1500, TargetEvent = c(1, 2), CVArg = list(V = 5),
    Model = M, HazEnsemble = TRUE, RenameCovs = FALSE, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  sl <- attr(est, "InitFits")[["1"]]
  # with covariate effects real and the baseline fixed, the full model should win
  expect_gt(sl$SLCoef[["full"]], sl$SLCoef[["trtonly"]])
  expect_equal(sum(sl$SLCoef), 1, tolerance = 1e-6)
})

test_that("doConcrete runs with HazEnsemble = TRUE", {
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
    ID = "id", Intervention = makeITT(), TargetTime = 2000, TargetEvent = 1,
    CVArg = list(V = 3), UpdateMethod = "adaptive", EICStopRule = "absolute",
    MaxUpdateIter = 12, HazEnsemble = TRUE, Verbose = FALSE)
  est <- suppressMessages(doConcrete(args))
  out <- getOutput(est, Estimand = c("Risk", "RD"), Intervention = c(1, 2), Simultaneous = FALSE)
  expect_true(all(is.finite(out[out$Estimator == "tmle", se])))
})
