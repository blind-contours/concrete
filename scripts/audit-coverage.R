#!/usr/bin/env Rscript
# Across-the-board coverage audit on a competing-risks DGP with closed-form truth.
# One fit per rep; checks Risk, RD, RR, RMST (getRMST) and RMST (targetRMST) for
# coverage / bias, in a null and an alternative scenario.
#   Rscript scripts/audit-coverage.R [B] [n]
suppressWarnings(suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION"))
    devtools::load_all(".", quiet = TRUE) else library(concrete)
  library(data.table)
}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 60L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 400L
tau <- 1200; tau_max <- 2500; times <- c(300, 600, 900, tau)

scen <- list(
  null = list(l1_0 = 5e-4, b1 = 0.0, g1 = 0.4, l2_0 = 3e-4, b2 =  0.0, g2 = 0.3, lc = 4e-4),
  alt  = list(l1_0 = 5e-4, b1 = 0.5, g1 = 0.4, l2_0 = 3e-4, b2 = -0.2, g2 = 0.3, lc = 4e-4))
hazf <- function(p, A, W1, W2) {
  l1 <- p$l1_0 * exp(p$b1 * A + p$g1 * W1); l2 <- p$l2_0 * exp(p$b2 * A + p$g2 * W2)
  list(l1 = l1, l2 = l2, lt = l1 + l2)
}
truth <- function(p, Wn = 4e5) {
  W1 <- stats::rnorm(Wn); W2 <- stats::rnorm(Wn)
  h1 <- hazf(p, 1, W1, W2); h0 <- hazf(p, 0, W1, W2)
  F1 <- function(h, t) mean((h$l1 / h$lt) * (1 - exp(-h$lt * t)))   # cause-1 risk, arm
  rmst <- function(h) mean((1 - exp(-h$lt * tau)) / h$lt)
  list(risk1 = function(t) F1(h1, t), risk0 = function(t) F1(h0, t),
       rmst1 = rmst(h1), rmst0 = rmst(h0))
}
set.seed(20260607); TR <- lapply(scen, truth)

simData <- function(p, seed) {
  set.seed(seed)
  W1 <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5); h <- hazf(p, A, W1, W2)
  T1 <- stats::rexp(n, h$l1); T2 <- stats::rexp(n, h$l2); C <- stats::rexp(n, p$lc)
  To <- pmin(T1, T2, C, tau_max)
  ev <- ifelse(To >= tau_max, 0L, ifelse(T1 <= pmin(T2, C), 1L, ifelse(T2 <= C, 2L, 0L)))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W1 = W1, W2 = W2)
}
cov <- function(lo, hi, tr) is.finite(lo) & lo <= tr & tr <= hi
row <- function(est, tr, scn, lab) data.table(scenario = scn, Estimand = lab,
  bias = est[["Pt Est"]] - tr, cover = cov(est[["CI Low"]], est[["CI Hi"]], tr))

one <- function(s, sc) {
  p <- scen[[sc]]; tr <- TR[[sc]]; d <- simData(p, s)
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = times, TargetEvent = c(1, 2),
      CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
      MaxUpdateIter = 25, Verbose = FALSE)
    est <- suppressMessages(doConcrete(a))
    out <- as.data.table(getOutput(est, Estimand = c("Risk", "RD", "RR"),
                                   Intervention = c(1, 2), Simultaneous = FALSE))
    rows <- list()
    for (tm in times) {
      rk <- out[Estimator == "tmle" & Estimand == "Abs Risk" & Event == 1 & Time == tm & grepl("=1", Intervention)]
      if (nrow(rk)) rows[[length(rows)+1]] <- row(rk, tr$risk1(tm), sc, "Risk(A=1)")
      rd <- out[Estimator == "tmle" & Estimand == "Risk Diff" & Event == 1 & Time == tm]
      if (nrow(rd)) rows[[length(rows)+1]] <- row(rd, tr$risk1(tm) - tr$risk0(tm), sc, "RD")
      rr <- out[Estimator == "tmle" & Estimand == "Rel Risk" & Event == 1 & Time == tm]
      if (nrow(rr)) rows[[length(rows)+1]] <- row(rr, tr$risk1(tm) / tr$risk0(tm), sc, "RR")
    }
    rm <- as.data.table(getRMST(est, Horizon = tau, Intervention = c(1, 2)))
    r1 <- rm[Estimand == "RMST" & Intervention == "A=1"]; r0 <- rm[Estimand == "RMST" & Intervention == "A=0"]
    rd <- rm[grepl("RMST Diff", Estimand)]
    if (nrow(r1)) rows[[length(rows)+1]] <- row(r1, tr$rmst1, sc, "RMST(A=1)")
    if (nrow(r0)) rows[[length(rows)+1]] <- row(r0, tr$rmst0, sc, "RMST(A=0)")
    if (nrow(rd)) rows[[length(rows)+1]] <- row(rd, tr$rmst1 - tr$rmst0, sc, "RMST Diff")
    tr1 <- tryCatch({
      tm2 <- as.data.table(suppressWarnings(targetRMST(est, Horizon = tau, Intervention = c(1, 2))))
      tm2[Estimand == "RMST" & Intervention == "A=1"]
    }, error = function(e) NULL)
    if (!is.null(tr1) && nrow(tr1)) rows[[length(rows)+1]] <- row(tr1, tr$rmst1, sc, "targetRMST(A=1)")
    data.table::rbindlist(rows)
  }, error = function(e) NULL)
}
seeds <- 20260607L + seq_len(B)
envCores <- suppressWarnings(as.integer(Sys.getenv("VAL_CORES", unset = NA)))
uc <- if (!is.na(envCores) && envCores >= 1L) envCores else
  if (.Platform$OS.type == "unix") min(4L, max(1L, parallel::detectCores() - 1L)) else 1L
res <- data.table::rbindlist(c(
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "null", mc.cores = uc)),
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "alt",  mc.cores = uc))))
summ <- res[, .(Reps = .N, Bias = round(mean(bias), 4), Coverage = round(mean(cover), 3)),
            by = .(scenario, Estimand)][order(scenario, Estimand)]
print(summ)
data.table::fwrite(summ, "scripts/audit-coverage-summary.csv")
cat("\nTarget: Coverage ~ 0.95 for every estimand; bias ~ 0.\n")
