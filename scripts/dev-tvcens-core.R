## Core-level bias->fix: informative censoring driven by a time-varying covariate
## biases the targeted absolute risk from doConcrete() under baseline-only IPCW;
## supplying CensoringTV (the time-varying covariate for the censoring model)
## recovers the uncensored truth. Frailty U drives events + the covariate L(t) +
## censoring; baseline W is independent of U (so baseline-only G cannot capture it).
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
tau <- 5; visits <- c(1, 2, 3.5)

simCore <- function(n, seed) {
  set.seed(seed); W <- rnorm(n); W2 <- rnorm(n); U <- rnorm(n); trt <- rbinom(n, 1, 0.5)
  Tt <- rexp(n, 0.12 * exp(0.2 * W - 0.4 * trt + 1.0 * U))     # event time (frailty U)
  Lvis <- sapply(visits, function(v) U + rnorm(n, 0, 0.4))     # noisy U at visits
  # informative censoring on a fine grid, hazard rising with current L (LOCF)
  M <- 40L; gg <- seq(0, tau, length.out = M + 1L); dt <- tau / M; st <- gg[-(M + 1L)]
  vi <- vapply(st, function(g) sum(visits <= g), integer(1))
  Lmat <- matrix(0, n, M); for (j in 1:M) if (vi[j] >= 1) Lmat[, j] <- Lvis[, vi[j]]
  Cint <- rep(M + 1L, n); done <- logical(n)
  for (j in 1:M) { lamC <- 0.14 * exp(1.0 * Lmat[, j]); fire <- !done & runif(n) < 1 - exp(-lamC * dt)
    Cint[fire] <- j; done <- done | fire }
  Clat <- ifelse(Cint <= M, gg[pmin(Cint, M) + 1L], Inf)   # interval end (> 0)
  obsT <- pmin(Tt, Clat, tau); status <- ifelse(Tt <= pmin(Clat, tau), 1L, 0L)
  dat <- data.table(id = seq_len(n), time = obsT, status = status, trt = trt, W = W, W2 = W2)
  tv <- rbindlist(lapply(seq_along(visits), function(k)
    data.table(id = seq_len(n), time = visits[k], L = Lvis[, k])))
  list(dat = dat, tv = tv, truthA1 = mean(Tt[trt == 1] <= 4), truthA0 = mean(Tt[trt == 0] <= 4))
}

riskFit <- function(dat, ctv) {
  args <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time", EventType = "status",
    Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = 4, TargetEvent = 1,
    CVArg = list(V = 3), MaxUpdateIter = 30, Model = NULL, CensoringTV = ctv))
  o <- suppressMessages(getOutput(suppressMessages(doConcrete(args))))
  o <- o[o$Estimator == "tmle" & o$Estimand == "Abs Risk", ]
  c(A1 = o[["Pt Est"]][o$Intervention == "A=1"][1], A0 = o[["Pt Est"]][o$Intervention == "A=0"][1])
}

B <- 12L; n <- 2500L
R <- t(vapply(seq_len(B), function(b) {
  d <- simCore(n, 700 + b)
  r0 <- riskFit(d$dat, NULL); rt <- riskFit(d$dat, d$tv)
  c(truthA1 = d$truthA1, noTV_A1 = r0["A1"], tv_A1 = rt["A1"],
    truthA0 = d$truthA0, noTV_A0 = r0["A0"], tv_A0 = rt["A0"])
}, numeric(6)))
cat(sprintf("\ncore absolute-risk bias under informative censoring (%d reps, n=%d):\n", B, n))
cat(sprintf("  arm A=1: TRUTH=%.4f  no-tv=%.4f (bias %+.4f)  tv=%.4f (bias %+.4f)\n",
  mean(R[, "truthA1"]), mean(R[, "noTV_A1.A1"]), mean(R[, "noTV_A1.A1"] - R[, "truthA1"]),
  mean(R[, "tv_A1.A1"]), mean(R[, "tv_A1.A1"] - R[, "truthA1"])))
cat(sprintf("  arm A=0: TRUTH=%.4f  no-tv=%.4f (bias %+.4f)  tv=%.4f (bias %+.4f)\n",
  mean(R[, "truthA0"]), mean(R[, "noTV_A0.A0"]), mean(R[, "noTV_A0.A0"] - R[, "truthA0"]),
  mean(R[, "tv_A0.A0"]), mean(R[, "tv_A0.A0"] - R[, "truthA0"])))
