## Control cell for the strata-correction validation: SAME DGP but with SIMPLE
## (iid Bernoulli) randomization and adjusted models, no Strata argument.
## If SD/SE ratio here is also ~0.86, the residual gap in dev-strata-validation.R
## is common-mode small-sample IF conservatism in this DGP -- not a flaw in the
## covariate-adaptive correction. If it is ~1.0, the correction is incomplete.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; B <- 150L; n <- 600L; tau <- 5; targT <- 4
bzv <- c(-0.8, -0.3, 0.3, 0.8); trtEff <- -0.4
F1 <- function(l) 1 - exp(-l * targT)
truthRD <- mean(F1(0.20 * exp(bzv + trtEff))) - mean(F1(0.20 * exp(bzv)))

oneRep <- function(b) {
  on.exit(gc(FALSE)); set.seed(7000 + b)            # same seeds as the main sim
  Z <- sample(1:4, n, TRUE); W2 <- rnorm(n)
  Ablock <- integer(n)                               # drawn to keep the RNG stream
  for (s in 1:4) { ix <- which(Z == s)
    for (bb in split(ix, ceiling(seq_along(ix) / 4))) Ablock[bb] <- sample(rep(0:1, length.out = length(bb))) }
  A <- rbinom(n, 1, 0.5)                             # SIMPLE randomization
  T1 <- rexp(n, 0.20 * exp(bzv[Z] + trtEff * A)); C <- rexp(n, 0.08)
  dat <- data.table(id = seq_len(n), time = pmin(T1, C, tau),
                    status = as.integer(T1 <= pmin(C, tau)), trt = A,
                    Z = factor(Z), W2 = W2)
  a <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = targT, TargetEvent = 1,
        CVArg = list(V = 2), MaxUpdateIter = 12, Verbose = FALSE,
        Model = list(trt = "SL.mean",
                     "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                     "1" = list(Cox = survival::Surv(time, status == 1) ~ .))))
  o <- suppressMessages(getOutput(suppressMessages(suppressWarnings(doConcrete(a))),
                                  Estimand = "RD", Simultaneous = FALSE))
  o <- as.data.table(o)[Estimator == "tmle" & Estimand == "Risk Diff"]
  c(rd = o[["Pt Est"]][1], se = o[["se"]][1],
    cov = as.integer(o[["CI Low"]][1] <= truthRD & truthRD <= o[["CI Hi"]][1]))
}
R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b)
  tryCatch(oneRep(b), error = function(e) rep(NA_real_, 3)), mc.cores = MC))
R <- R[stats::complete.cases(R), , drop = FALSE]
cat(sprintf("\n--- SIMPLE randomization control (B=%d, n=%d, adjusted models) ---\n", nrow(R), n))
cat(sprintf("  truth RD = %+.4f  mean est %+.4f (bias %+.4f)\n",
            truthRD, mean(R[, 1]), mean(R[, 1]) - truthRD))
cat(sprintf("  empirical SD %.5f | mean iid SE %.5f | SD/SE ratio %.3f | coverage %.3f\n",
            sd(R[, 1]), mean(R[, 2]), sd(R[, 1]) / mean(R[, 2]), mean(R[, 3])))
