## Win-ratio-level validation of time-varying censoring covariates in
## clinicalWinRatio(). DGP: illness-death (death > HFH) with a latent frailty U
## driving outcomes, a time-varying covariate L(t) = noisy U at visits, and
## informative censoring driven by L(t). Compare to brute-force uncensored truth:
##   - clinicalWinRatio WITHOUT censoring.tv  -> biased (informative censoring)
##   - clinicalWinRatio WITH    censoring.tv  -> recovers truth
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
tau <- 3; M <- 30L; dt <- tau / M; gstart <- seq(0, tau - dt, length.out = M); visits <- c(0.4, 0.9, 1.6, 2.3)

simWR <- function(n, seed) {
  set.seed(seed); W <- rnorm(n); U <- rnorm(n); arm <- rbinom(n, 1, 0.5)
  fW <- exp(0.3 * W + 0.7 * U); fT <- exp(-0.45 * arm); fE <- exp(-0.30 * arm)
  a01 <- 0.20 * fW * fE; a02 <- 0.10 * fW * fT; a12 <- 0.32 * fW * fT
  T01 <- rexp(n, a01); T02 <- rexp(n, a02); u0 <- pmin(T01, T02); hfh <- T01 < T02
  sH <- ifelse(hfh, u0, Inf); Dlat <- ifelse(hfh, u0 + rexp(n, a12), u0)   # latent death
  Lvis <- sapply(visits, function(v) U + rnorm(n, 0, 0.4))
  vi <- vapply(gstart, function(g) sum(visits <= g), integer(1))
  Lmat <- matrix(0, n, M); for (j in 1:M) if (vi[j] >= 1) Lmat[, j] <- Lvis[, vi[j]]
  Cint <- rep(M + 1L, n); done <- logical(n)
  for (j in 1:M) { lamC <- 0.18 * exp(0.8 * Lmat[, j]); fire <- !done & runif(n) < 1 - exp(-lamC * dt)
    Cint[fire] <- j; done <- done | fire }
  Clat <- ifelse(Cint <= M, gstart[pmin(Cint, M)], Inf)
  list(W = W, arm = arm, sH = sH, Dlat = Dlat, Clat = Clat, Lvis = Lvis)
}

cap <- function(x) ifelse(x <= tau, x, Inf)
bruteWR <- function(d, np = 2e6) {                      # uncensored truth: death > HFH
  it <- which(d$arm == 1); ic <- which(d$arm == 0)
  i <- sample(it, np, TRUE); j <- sample(ic, np, TRUE)
  xD <- cap(d$Dlat[i]); yD <- cap(d$Dlat[j]); xH <- cap(d$sH[i]); yH <- cap(d$sH[j])
  win <- xD > yD; loss <- xD < yD; dec <- xD != yD
  u <- !dec; win <- win | (u & xH > yH); loss <- loss | (u & xH < yH)
  mean(win) / mean(loss)
}

makeInputs <- function(d) {
  n <- length(d$W); termT <- pmin(d$Dlat, d$Clat, tau)
  died <- as.integer(d$Dlat <= d$Clat & d$Dlat <= tau)
  hosp <- ifelse(d$sH <= termT, d$sH, NA)
  dat <- data.frame(idv = seq_len(n), arm = d$arm, t_hosp = hosp,
                    t_term = termT, died = ifelse(termT >= tau, 0L, died), W = d$W)
  tv <- do.call(rbind, lapply(seq_along(visits), function(k)
    data.frame(idv = seq_len(n), time = visits[k], L = d$Lvis[, k])))
  list(dat = dat, tv = tv)
}
runWR <- function(dat, tv = NULL, idc = NULL)
  clinicalWinRatio(dat, arm = "arm", illness.time = "t_hosp", terminal.time = "t_term",
    terminal.status = "died", covariates = "W", horizon = tau, n.grid = 30L, n.folds = 1L,
    id = idc, censoring.tv = tv)$`Pt Est`[1]

## smoke: one run with and without censoring.tv
in1 <- makeInputs(simWR(1500L, 1L))
cat(sprintf("smoke: WR (no tv) = %.3f ; WR (tv) = %.3f\n",
            runWR(in1$dat), runWR(in1$dat, in1$tv, "idv")))

## bias check over reps
B <- 24L
R <- t(vapply(seq_len(B), function(b) {
  d <- simWR(1500L, 100 + b); ins <- makeInputs(d)
  c(truth = bruteWR(d), noTV = runWR(ins$dat), TV = runWR(ins$dat, ins$tv, "idv"))
}, numeric(3)))
cat(sprintf("\nwin-ratio bias under informative censoring (%d reps, n=1500/arm... 750/arm):\n", B))
cat(sprintf("  brute-force TRUTH WR      = %.4f\n", mean(R[, "truth"])))
cat(sprintf("  clinicalWinRatio no tv    = %.4f   bias %+.4f\n", mean(R[, "noTV"]), mean(R[, "noTV"] - R[, "truth"])))
cat(sprintf("  clinicalWinRatio + tv     = %.4f   bias %+.4f\n", mean(R[, "TV"]),   mean(R[, "TV"]   - R[, "truth"])))
