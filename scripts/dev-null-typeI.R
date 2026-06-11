## Null (type-I) validation for the core absolute-risk TMLE: treatment has NO
## effect on either cause. Checks per-arm risk coverage against closed-form truth,
## RD coverage at RD = 0, and the two-sided type-I error of the RD Wald test.
## This closes the main evidence gap: all prior null sims were win-ratio only.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; B <- 160L; n <- 500L; tau <- 5; targT <- 4

## DGP: two competing events + censoring, all depending on W but NOT on A.
## Closed-form truth: cause-1 CIF at t for rate l1(W), l2(W):
##   F1(t|W) = l1/(l1+l2) * (1 - exp(-(l1+l2) t)); marginal by MC integration.
l1f <- function(W) 0.20 * exp(0.4 * W)
l2f <- function(W) 0.10 * exp(-0.3 * W)
set.seed(99); Wbig <- rnorm(2e6)
truthF1 <- mean(l1f(Wbig)/(l1f(Wbig)+l2f(Wbig)) * (1 - exp(-(l1f(Wbig)+l2f(Wbig)) * targT)))

simNull <- function(n, seed) {
  set.seed(seed); N <- 2L * n
  W <- rnorm(N); W2 <- rnorm(N); A <- rep(0:1, each = n)
  l1 <- l1f(W); l2 <- l2f(W)
  Te <- rexp(N, l1 + l2); cause <- 1L + rbinom(N, 1L, l2 / (l1 + l2))
  C <- rexp(N, 0.10)
  time <- pmin(Te, C, tau)
  status <- ifelse(Te <= pmin(C, tau), cause, 0L)
  data.table(id = seq_len(N), time = time, status = status, trt = A, W = W, W2 = W2)
}

Mdl <- list(trt = c("SL.mean", "SL.glm"),
            "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
            "1" = list(Cox = survival::Surv(time, status == 1) ~ .),
            "2" = list(Cox = survival::Surv(time, status == 2) ~ .))

oneRep <- function(b) {
  on.exit(gc(FALSE))
  dat <- simNull(n, 3000 + b)
  a <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = targT, TargetEvent = c(1, 2),
        CVArg = list(V = 2), MaxUpdateIter = 15, Model = Mdl, Verbose = FALSE))
  o <- suppressMessages(getOutput(suppressMessages(doConcrete(a)), Estimand = c("Risk", "RD"),
                                  Simultaneous = FALSE))
  o <- as.data.table(o)[Estimator == "tmle" & Event == 1]
  r1 <- o[Estimand == "Abs Risk" & Intervention == "A=1"]
  r0 <- o[Estimand == "Abs Risk" & Intervention == "A=0"]
  rd <- o[Estimand == "Risk Diff"]
  c(est1 = r1[["Pt Est"]], cov1 = as.integer(r1[["CI Low"]] <= truthF1 & truthF1 <= r1[["CI Hi"]]),
    est0 = r0[["Pt Est"]], cov0 = as.integer(r0[["CI Low"]] <= truthF1 & truthF1 <= r0[["CI Hi"]]),
    rd = rd[["Pt Est"]], rdse = rd[["se"]],
    covRD = as.integer(rd[["CI Low"]] <= 0 & 0 <= rd[["CI Hi"]]),
    rej = as.integer(rd[["pValue"]] < 0.05))
}

R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b)
  tryCatch(oneRep(b), error = function(e) rep(NA_real_, 8)), mc.cores = MC))
colnames(R) <- c("est1","cov1","est0","cov0","rd","rdse","covRD","rej")
R <- R[stats::complete.cases(R), , drop = FALSE]
saveRDS(R, "/tmp/null-typeI.rds")
cat(sprintf("\n===== NULL DGP: core TMLE type-I / coverage (%d reps, n=%d/arm) =====\n", nrow(R), n))
cat(sprintf("  truth F1(%.1f)      = %.4f\n", targT, truthF1))
cat(sprintf("  A=1 risk: mean est  = %.4f  bias %+.4f   coverage %.3f\n",
            mean(R[,"est1"]), mean(R[,"est1"]) - truthF1, mean(R[,"cov1"])))
cat(sprintf("  A=0 risk: mean est  = %.4f  bias %+.4f   coverage %.3f\n",
            mean(R[,"est0"]), mean(R[,"est0"]) - truthF1, mean(R[,"cov0"])))
cat(sprintf("  RD (true 0): mean   = %+.4f  empSD %.4f  mean IF-SE %.4f (ratio %.2f)\n",
            mean(R[,"rd"]), sd(R[,"rd"]), mean(R[,"rdse"]), sd(R[,"rd"]) / mean(R[,"rdse"])))
cat(sprintf("  RD coverage at 0    = %.3f\n", mean(R[,"covRD"])))
cat(sprintf("  TYPE-I ERROR (a=.05)= %.3f   [acceptance band for B=%d: %.3f-%.3f]\n",
            mean(R[,"rej"]), nrow(R),
            0.05 - 1.96*sqrt(0.05*0.95/nrow(R)), 0.05 + 1.96*sqrt(0.05*0.95/nrow(R))))
