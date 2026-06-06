#!/usr/bin/env Rscript
# End-to-end validation of the estimator and the new inference on a DGP with
# closed-form truth (exponential competing-risk hazards). Two scenarios:
#   null : no treatment effect on cause 1 (true RD = 0)  -> type-I error + coverage
#   alt  : a real treatment effect on cause 1            -> bias + coverage + power
# Checks risk-difference coverage, Wald p-value calibration, and RMST coverage.
#   Rscript scripts/make-validation-sim.R [B] [n]
suppressWarnings(suppressMessages({library(concrete); library(data.table)}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 100L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 400L
tau <- 1200; tau_max <- 2500

scen <- list(
  null = list(l1_0 = 5e-4, b1 = 0.0, g1 = 0.4, l2_0 = 3e-4, b2 = 0.0, g2 = 0.3, lc = 4e-4),
  alt  = list(l1_0 = 5e-4, b1 = 0.5, g1 = 0.4, l2_0 = 3e-4, b2 = -0.2, g2 = 0.3, lc = 4e-4))

hazf <- function(p, A, W1, W2) {
  l1 <- p$l1_0 * exp(p$b1 * A + p$g1 * W1); l2 <- p$l2_0 * exp(p$b2 * A + p$g2 * W2)
  list(l1 = l1, l2 = l2, lt = l1 + l2)
}
truth <- function(p, Wn = 4e5) {
  W1 <- stats::rnorm(Wn); W2 <- stats::rnorm(Wn)
  h1 <- hazf(p, 1, W1, W2); h0 <- hazf(p, 0, W1, W2)
  # true marginal cause-1 absolute risk at an arbitrary time t, by arm
  F1 <- function(h, t) mean((h$l1 / h$lt) * (1 - exp(-h$lt * t)))
  rdAt <- function(t) F1(h1, t) - F1(h0, t)              # RD evaluated AT time t
  rmst <- function(a) { h <- hazf(p, a, W1, W2); mean((1 - exp(-h$lt * tau)) / h$lt) }
  list(rdAt = rdAt, RD = rdAt(tau), RMST1 = rmst(1), RMST0 = rmst(0))
}
set.seed(20260605)
TR <- lapply(scen, truth)
cat(sprintf("Truth: null RD=%.4f | alt RD=%.4f\n", TR$null$RD, TR$alt$RD))

simData <- function(p, seed) {
  set.seed(seed)
  W1 <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5); h <- hazf(p, A, W1, W2)
  T1 <- stats::rexp(n, h$l1); T2 <- stats::rexp(n, h$l2); C <- stats::rexp(n, p$lc)
  To <- pmin(T1, T2, C, tau_max)
  ev <- ifelse(To >= tau_max, 0L, ifelse(T1 <= pmin(T2, C), 1L, ifelse(T2 <= C, 2L, 0L)))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W1 = W1, W2 = W2)
}
cover <- function(lo, hi, truth) is.finite(lo) & lo <= truth & truth <= hi

one <- function(s, sc) {
  p <- scen[[sc]]; tr <- TR[[sc]]; d <- simData(p, s)
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = c(300,600,900,tau), TargetEvent = c(1, 2),
      CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
      MaxUpdateIter = 20, Verbose = FALSE)
    est <- suppressMessages(doConcrete(a))
    o <- as.data.table(getOutput(est, Estimand = "RD", Intervention = c(1, 2), Simultaneous = FALSE))
    rd <- o[Estimator == "tmle" & Estimand == "Risk Diff" & Event == 1]
    rd_truth <- vapply(rd$Time, tr$rdAt, numeric(1))   # truth matched to each target time
    rm <- as.data.table(getRMST(est, Horizon = tau, Intervention = c(1, 2)))
    rmst1 <- rm[Estimand == "RMST" & Intervention == "A=1"]
    data.table(scenario = sc, Time = rd$Time,
               rd_cov = cover(rd$`CI Low`, rd$`CI Hi`, rd_truth),
               rd_bias = rd$`Pt Est` - rd_truth,
               rd_reject = rd$pValue < 0.05,
               rmst_cov = cover(rmst1$`CI Low`, rmst1$`CI Hi`, tr$RMST1))
  }, error = function(e) NULL)
}
seeds <- 20260605L + seq_len(B)
envCores <- suppressWarnings(as.integer(Sys.getenv("VAL_CORES", unset = NA)))
useCores <- if (!is.na(envCores) && envCores >= 1L) envCores else
  if (.Platform$OS.type == "unix") min(6L, max(1L, parallel::detectCores() - 1L)) else 1L
runs <- data.table::rbindlist(c(
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "null", mc.cores = useCores)),
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "alt", mc.cores = useCores))))
summ <- runs[, .(Reps = .N, RD_coverage = mean(rd_cov), RD_bias = mean(rd_bias),
                 `RD_reject@.05` = mean(rd_reject), RMST_coverage = mean(rmst_cov)), by = scenario]
byTime <- runs[, .(Reps = .N, RD_coverage = mean(rd_cov), RD_bias = mean(rd_bias),
                   `RD_reject@.05` = mean(rd_reject)), by = .(scenario, Time)][order(scenario, Time)]
cat("\n== Pooled over target times ==\n"); print(summ)
cat("\n== By target time (RD truth matched per time) ==\n"); print(byTime)
data.table::fwrite(summ, "scripts/validation-sim-summary.csv")
data.table::fwrite(byTime, "scripts/validation-sim-bytime.csv")
cat("\nInterpretation: null RD_reject ~ 0.05 (type-I error); coverages ~ 0.95;",
    "alt RD_reject = power (grows with target time as the effect accrues).\n")
