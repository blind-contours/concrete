#!/usr/bin/env Rscript
# LIGHT smoke test: confirms the hazard-learner fixes hold and the new features
# run end-to-end. Small n, few reps, few iterations -- not the full validation.
#   Rscript scripts/smoke-validation.R [B] [n]
suppressWarnings(suppressMessages({library(concrete); library(data.table)}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 15L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 300L
tau <- 1200; tau_max <- 2500; MAXIT <- 12L
scen <- list(
  null = list(b1 = 0.0),     # no treatment effect on the (single) event -> true RD = 0
  alt  = list(b1 = 0.6))     # strong effect (the regime that exposed the Cox/coxnet bug)
hz <- function(b1, A, W, W2) 6e-4 * exp(b1 * A + 0.5 * W + 0.4 * W2)
truthRD <- function(b1, Wn = 3e5) {
  W <- stats::rnorm(Wn); W2 <- stats::rnorm(Wn)
  f <- function(a) mean(1 - exp(-hz(b1, a, W, W2) * tau))
  f(1) - f(0)
}
set.seed(11); TR <- sapply(scen, function(s) truthRD(s$b1))
cat(sprintf("Truth RD: null=%.4f  alt=%.4f\n", TR["null"], TR["alt"]))
simd <- function(b1, s) {
  set.seed(s); W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, .5)
  Tt <- stats::rexp(n, hz(b1, A, W, W2)); C <- stats::rexp(n, 3e-4)
  To <- pmin(Tt, C, tau_max); ev <- ifelse(To >= tau_max, 0L, ifelse(Tt <= C, 1L, 0L))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W = W, W2 = W2)
}
cov <- function(lo, hi, tr) is.finite(lo) & lo <= tr & tr <= hi
one <- function(s, sc) {
  d <- simd(scen[[sc]]$b1, s); tr <- TR[[sc]]
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = tau, TargetEvent = 1, CVArg = list(V = 5),
      UpdateMethod = "adaptive", EICStopRule = "absolute", MaxUpdateIter = MAXIT, Verbose = FALSE)
    o <- as.data.table(getOutput(suppressMessages(doConcrete(a)), Estimand = "RD",
                                 Intervention = c(1, 2), Simultaneous = FALSE))
    r <- o[Estimator == "tmle" & Estimand == "Risk Diff"]
    data.table(scenario = sc, bias = r$`Pt Est` - tr,
               cover = cov(r$`CI Low`, r$`CI Hi`, tr), reject = r$pValue < 0.05)
  }, error = function(e) NULL)
}
seeds <- 100L + seq_len(B)
uc <- if (.Platform$OS.type == "unix") min(4L, max(1L, parallel::detectCores() - 2L)) else 1L
res <- data.table::rbindlist(c(
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "null", mc.cores = uc)),
  Filter(Negate(is.null), parallel::mclapply(seeds, one, sc = "alt", mc.cores = uc))))
cat("\n=== RD smoke (n=", n, ", ", B, " reps, ", MAXIT, " iters) ===\n", sep = "")
print(res[, .(Reps = .N, Bias = round(mean(bias), 4), Coverage = round(mean(cover), 3),
              `Reject@.05` = round(mean(reject), 3)), by = scenario])

cat("\n=== feature end-to-end smoke (one dataset, must all run) ===\n")
d <- simd(0.6, 1)
a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
  ID = "id", Intervention = makeITT(), TargetTime = c(400, 800, tau), TargetEvent = 1,
  CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute", MaxUpdateIter = MAXIT,
  Verbose = FALSE)
chk <- function(label, expr) {
  r <- tryCatch({ expr; "ok" }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("  %-28s %s\n", label, r))
}
est <- suppressMessages(doConcrete(a))
chk("getOutput Risk/RD/RR", getOutput(est, c("Risk", "RD", "RR"), c(1, 2), Simultaneous = TRUE))
chk("getRMST",              getRMST(est, Horizon = tau, Intervention = c(1, 2)))
chk("targetRMST",          suppressWarnings(targetRMST(est, Horizon = tau, Intervention = c(1, 2))))
chk("getWinRatio",          getWinRatio(est, Horizon = tau, Intervention = c(1, 2)))
chk("getTmleDiagnostics",   getTmleDiagnostics(est, "components"))
chk("doConcrete CrossFit",  suppressMessages(doConcrete(formatArguments(a, CrossFit = TRUE))))
chk("doConcrete HazEnsemble", suppressMessages(doConcrete(formatArguments(a, HazEnsemble = TRUE))))
chk("senseCensoring",       suppressMessages(senseCensoring(a, deltas = c(0, 0.5), Estimand = "RD")))
cat("\nSmoke test done.\n")
