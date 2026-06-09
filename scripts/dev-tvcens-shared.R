## Validate the shared .tvCensoringInc(): baseline-only censoring G is biased under
## informative (time-varying-covariate-driven) censoring; adding the time-varying
## covariate recovers truth. Same DGP as scripts/dev-tvcens.R, now via the package
## internal estimator that the core + clinicalWinRatio will both use.
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE) }))
tau <- 1.5; M <- 30L; dt <- tau / M; times <- seq(0, tau, length.out = M + 1L)
starts <- times[-(M + 1L)]; visits <- c(0.2, 0.45, 0.75, 1.05)
visitIdx <- vapply(starts, function(g) sum(visits <= g), integer(1))

simData <- function(n, seed) {
  set.seed(seed); W <- rnorm(n); U <- rnorm(n)
  Tt <- rexp(n, 0.6 * exp(0.2 * W + 0.9 * U))
  Lvis <- sapply(visits, function(v) U + rnorm(n, 0, 0.35))
  Lmat <- matrix(0, n, M); for (j in 1:M) if (visitIdx[j] >= 1) Lmat[, j] <- Lvis[, visitIdx[j]]
  Cint <- rep(M + 1L, n); done <- logical(n)
  for (j in 1:M) { lamC <- 0.30 * exp(0.8 * Lmat[, j]); fire <- !done & runif(n) < 1 - exp(-lamC * dt)
    Cint[fire] <- j; done <- done | fire }
  eventInt <- ifelse(Tt <= tau, findInterval(Tt, times), M + 1L)
  obsT <- pmin(ifelse(Tt <= tau, Tt, tau), ifelse(Cint <= M, starts[pmin(Cint, M)], tau))
  censInd <- as.integer(Cint < eventInt & Cint <= M)
  base <- Lmat[, 1]
  list(W = W, Tt = Tt, Cint = Cint, obsT = obsT, censInd = censInd,
       tvMats = list(L_val = Lmat, L_chg = Lmat - base))
}

Gtau <- function(inc) apply(pmax(1 - inc, 0.02), 2, prod)     # G_i(tau) = prod(1 - hazard)
SL <- c("SL.mean", "SL.glm")
est <- function(d) {
  baseCov <- data.frame(W = d$W)
  incB <- .tvCensoringInc(times, d$obsT, d$censInd, baseCov, list(),       SL, 1L)   # baseline only
  incT <- .tvCensoringInc(times, d$obsT, d$censInd, baseCov, d$tvMats,     SL, 1L)   # + time-varying
  surv <- (d$Tt > tau) & (d$Cint > M)
  c(truth = mean(d$Tt > tau),
    ipcw_W  = mean(surv / Gtau(incB)),
    ipcw_WL = mean(surv / Gtau(incT)))
}

B <- 40L; n <- 4000L
R <- t(vapply(seq_len(B), function(b) est(simData(n, 3000 + b)), numeric(3)))
cat(sprintf("shared .tvCensoringInc (n=%d, %d reps):\n", n, B))
cat(sprintf("  TRUTH P(T>tau)                  = %.4f\n", mean(R[, "truth"])))
cat(sprintf("  IPCW baseline-only G(t|W)       = %.4f   bias %+.4f   <- BIASED\n",
            mean(R[, "ipcw_W"]),  mean(R[, "ipcw_W"]  - R[, "truth"])))
cat(sprintf("  IPCW time-varying G(t|W,L(t))   = %.4f   bias %+.4f   <- fixed\n",
            mean(R[, "ipcw_WL"]), mean(R[, "ipcw_WL"] - R[, "truth"])))
