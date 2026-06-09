## Prototype: time-varying covariates in the censoring model fix informative-
## censoring bias. Censoring is driven by a post-baseline measurement L(t) that is
## correlated with the outcome (via a latent frailty U, independent of baseline W).
## Target: theta = P(T > tau).  Estimators (IPCW):
##   - truth (uncensored latent T)
##   - naive IPCW, marginal G (assumes censoring ~ time only)
##   - IPCW with baseline-only G(t | W)            -> BIASED (misses L)
##   - IPCW with time-varying G(t | W, L(t))       -> should recover truth
## Same G that multiplies the win-ratio EIF clever covariates, so the fix transfers.
suppressWarnings(suppressMessages(library(stats)))
tau <- 1.5; M <- 30L; dt <- tau / M; grid <- seq(0, tau, length.out = M + 1L)
visits <- c(0.375, 0.75, 1.125)                 # tau/4, tau/2, 3tau/4
visitIdx <- vapply(grid[-(M + 1L)], function(g) sum(visits <= g), integer(1))  # LOCF visit per interval

simData <- function(n, seed) {
  set.seed(seed)
  W <- rnorm(n); U <- rnorm(n)                  # W (baseline) independent of U (frailty)
  T <- rexp(n, 0.6 * exp(0.2 * W + 0.9 * U))    # event time: frailty U drives it
  Lvis <- sapply(visits, function(v) U + rnorm(n, 0, 0.4))   # n x 3 noisy frailty readings
  Lmat <- matrix(0, n, M)                       # LOCF current L per interval (0 before first visit)
  for (j in 1:M) if (visitIdx[j] >= 1) Lmat[, j] <- Lvis[, visitIdx[j]]
  Cint <- rep(M + 1L, n); done <- logical(n)    # informative censoring: hazard rises with L(t)
  for (j in 1:M) { lamC <- 0.5 * exp(0.9 * Lmat[, j]); fire <- !done & runif(n) < 1 - exp(-lamC * dt)
    Cint[fire] <- j; done <- done | fire }
  eventInt <- ifelse(T <= tau, findInterval(T, grid), M + 1L)
  list(n = n, W = W, T = T, Cint = Cint, eventInt = eventInt, Lmat = Lmat)
}

## per-interval censoring long data (at risk until event/censor/tau; event = competing)
longCens <- function(d) {
  leaveInt <- pmin(d$eventInt, d$Cint); lastj <- pmin(leaveInt, M)
  si <- rep(seq_len(d$n), lastj); ji <- unlist(lapply(lastj, seq_len))
  data.frame(censEv = as.integer(ji == d$Cint[si] & d$Cint[si] < d$eventInt[si] & d$Cint[si] <= M),
             W = d$W[si], L = d$Lmat[cbind(si, ji)], t = grid[ji])
}
## predict G_i(tau) = prod_j (1 - censoring hazard) for each subject from a fitted model
Gtau <- function(model, d, useL) {
  Sj <- rep(seq_len(d$n), each = M); Jj <- rep(seq_len(M), times = d$n)
  pf <- data.frame(W = d$W[Sj], L = d$Lmat[cbind(Sj, Jj)], t = grid[Jj])
  p <- pmin(predict(model, pf, type = "response"), 0.95)
  apply(matrix(1 - p, nrow = M), 2, prod)
}

estimate <- function(d) {
  long <- longCens(d)
  m0  <- glm(censEv ~ poly(t, 2), binomial, long)              # marginal (time only)
  mW  <- glm(censEv ~ W + poly(t, 2), binomial, long)          # baseline only
  mWL <- glm(censEv ~ W + L + poly(t, 2), binomial, long)      # + time-varying L
  G0 <- pmax(Gtau(m0, d), 0.02); GW <- pmax(Gtau(mW, d), 0.02); GWL <- pmax(Gtau(mWL, d), 0.02)
  surv <- (d$T > tau) & (d$Cint > M)                           # observed event-free & uncensored to tau
  c(truth = mean(d$T > tau), naive = mean(surv / G0),
    ipcw_W = mean(surv / GW), ipcw_WL = mean(surv / GWL),
    pcens = mean(d$Cint <= M))
}

B <- 60L; n <- 4000L
R <- t(vapply(seq_len(B), function(b) estimate(simData(n, 1000 + b)), numeric(5)))
cat(sprintf("informative censoring driven by time-varying L(t) (n=%d, %d reps, ~%.0f%% censored)\n",
            n, B, 100 * mean(R[, "pcens"])))
cat(sprintf("  TRUTH  P(T>tau)               = %.4f\n", mean(R[, "truth"])))
cat(sprintf("  naive IPCW (G ~ time only)    = %.4f   bias %+.4f\n", mean(R[, "naive"]),  mean(R[, "naive"]  - R[, "truth"])))
cat(sprintf("  IPCW  G(t | W)  baseline-only = %.4f   bias %+.4f   <- BIASED\n", mean(R[, "ipcw_W"]),  mean(R[, "ipcw_W"]  - R[, "truth"])))
cat(sprintf("  IPCW  G(t | W, L(t)) time-var = %.4f   bias %+.4f   <- fixed\n", mean(R[, "ipcw_WL"]), mean(R[, "ipcw_WL"] - R[, "truth"])))
