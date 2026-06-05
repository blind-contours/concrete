#!/usr/bin/env Rscript
# Coverage validity check for cross-fitting (CrossFit = TRUE) vs the standard
# in-sample nuisance fit, on a DGP with closed-form truth (exponential
# competing-risk hazards). Confirms cross-fitted influence-function inference is
# valid (near-nominal coverage); the largest *advantage* of cross-fitting appears
# with flexible ML learners under misspecification, which is heavier to simulate.
#   Rscript scripts/make-crossfit-coverage.R [B] [n]
suppressWarnings(suppressMessages({library(concrete); library(data.table)}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 60L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 350L
tau <- 1200; tau_max <- 2500
p <- list(l1_0 = 5e-4, b1 = 0.3, g1 = 0.4, l2_0 = 3e-4, b2 = -0.2, g2 = 0.3, lc = 4e-4)
haz <- function(A, W1, W2) {
  l1 <- p$l1_0 * exp(p$b1 * A + p$g1 * W1); l2 <- p$l2_0 * exp(p$b2 * A + p$g2 * W2)
  list(l1 = l1, l2 = l2, lt = l1 + l2)
}
truthRisk <- function(a, Wn = 4e5) {
  W1 <- stats::rnorm(Wn); W2 <- stats::rnorm(Wn); h <- haz(a, W1, W2)
  mean((h$l1 / h$lt) * (1 - exp(-h$lt * tau)))   # F_1(tau | a)
}
set.seed(20260605); TR <- c(A0 = truthRisk(0), A1 = truthRisk(1))
cat(sprintf("True F1(tau): A0=%.4f A1=%.4f\n", TR["A0"], TR["A1"]))
sim <- function(seed) {
  set.seed(seed)
  W1 <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5); h <- haz(A, W1, W2)
  T1 <- stats::rexp(n, h$l1); T2 <- stats::rexp(n, h$l2); C <- stats::rexp(n, p$lc)
  To <- pmin(T1, T2, C, tau_max)
  ev <- ifelse(To >= tau_max, 0L, ifelse(T1 <= pmin(T2, C), 1L, ifelse(T2 <= C, 2L, 0L)))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W1 = W1, W2 = W2)
}
cover <- function(pt, se, tr, z = 1.96) is.finite(se) & (pt - z*se <= tr) & (tr <= pt + z*se)
one <- function(seed, cf) {
  d <- sim(seed)
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = tau, TargetEvent = c(1, 2),
      CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
      MaxUpdateIter = 20, CrossFit = cf, Verbose = FALSE)
    o <- getOutput(suppressMessages(doConcrete(a)), Estimand = "Risk", Intervention = c(1, 2),
                   Simultaneous = FALSE)
    o <- as.data.table(o)[Estimator == "tmle" & Event == 1]
    data.table(method = if (cf) "cross-fit" else "standard",
               cov0 = cover(o[Intervention == "A=0", `Pt Est`], o[Intervention == "A=0", se], TR["A0"]),
               cov1 = cover(o[Intervention == "A=1", `Pt Est`], o[Intervention == "A=1", se], TR["A1"]),
               b0 = o[Intervention == "A=0", `Pt Est`] - TR["A0"],
               b1 = o[Intervention == "A=1", `Pt Est`] - TR["A1"])
  }, error = function(e) NULL)
}
seeds <- 20260605L + seq_len(B)
# Cap cores: cross-fitting refits the nuisance library V times per replicate, so
# too many parallel workers exhaust memory and get silently killed by mclapply.
useCores <- if (.Platform$OS.type == "unix") min(3L, max(1L, parallel::detectCores() - 1L)) else 1L
runM <- function(cf) data.table::rbindlist(Filter(Negate(is.null),
  parallel::mclapply(seeds, one, cf = cf, mc.cores = useCores)))
res <- rbind(runM(FALSE), runM(TRUE))
summ <- res[, .(Reps = .N, CovA0 = mean(cov0), CovA1 = mean(cov1),
                BiasA0 = mean(b0), BiasA1 = mean(b1)), by = method]
print(summ)
data.table::fwrite(summ, "scripts/crossfit-coverage-summary.csv")
cat("Done.\n")
