## Validate the covariate-adaptive (stratified permuted-block) variance
## correction. DGP: 4 prognostic strata, permuted blocks of 4 within stratum.
## Cell (a): hazard models UNADJUSTED for the stratum -> iid SE should be
##           conservative (over-cover); the corrected SE should be calibrated
##           (ratio empSD/meanSE ~ 1, coverage ~ .95).
## Cell (b): hazard models ADJUSTED for the stratum -> correction ~ 0, both agree.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; B <- 150L; n <- 600L; tau <- 5; targT <- 4

bzv <- c(-0.8, -0.3, 0.3, 0.8); trtEff <- -0.4
## truth for RD at targT by mega-MC (closed form per stratum, equal stratum probs)
F1 <- function(l) 1 - exp(-l * targT)
truthRD <- mean(F1(0.20 * exp(bzv + trtEff))) - mean(F1(0.20 * exp(bzv)))

simStrat <- function(n, seed) {
  set.seed(seed)
  Z <- sample(1:4, n, TRUE); W2 <- rnorm(n)
  A <- integer(n)
  for (s in 1:4) {                                  # permuted blocks of 4 within stratum
    ix <- which(Z == s)
    for (b in split(ix, ceiling(seq_along(ix) / 4))) A[b] <- sample(rep(0:1, length.out = length(b)))
  }
  T1 <- rexp(n, 0.20 * exp(bzv[Z] + trtEff * A)); C <- rexp(n, 0.08)
  time <- pmin(T1, C, tau); status <- as.integer(T1 <= pmin(C, tau))
  data.table(id = seq_len(n), time = time, status = status, trt = A, Z = factor(Z), W2 = W2)
}

fitRD <- function(dat, adjusted, useStrata) {
  haz <- if (adjusted) survival::Surv(time, status == 1) ~ . else survival::Surv(time, status == 1) ~ trt
  cen <- if (adjusted) survival::Surv(time, status == 0) ~ . else survival::Surv(time, status == 0) ~ trt
  a <- suppressMessages(formatArguments(DataTable = copy(dat), EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = targT, TargetEvent = 1,
        CVArg = list(V = 2), MaxUpdateIter = 12, Verbose = FALSE,
        Model = list(trt = "SL.mean", "0" = list(Cox = cen), "1" = list(Cox = haz)),
        Strata = if (useStrata) "Z" else NULL))
  o <- suppressMessages(getOutput(suppressMessages(suppressWarnings(doConcrete(a))),
                                  Estimand = "RD", Simultaneous = FALSE))
  o <- as.data.table(o)[Estimator == "tmle" & Estimand == "Risk Diff"]
  c(rd = o[["Pt Est"]][1], se = o[["se"]][1],
    cov = as.integer(o[["CI Low"]][1] <= truthRD & truthRD <= o[["CI Hi"]][1]))
}

runCell <- function(adjusted) {
  R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) tryCatch({
    on.exit(gc(FALSE)); dat <- simStrat(n, 7000 + b)
    c(fitRD(dat, adjusted, FALSE), fitRD(dat, adjusted, TRUE))
  }, error = function(e) rep(NA_real_, 6)), mc.cores = MC))
  colnames(R) <- c("rd0", "se0", "cov0", "rd1", "se1", "cov1")
  R <- R[stats::complete.cases(R), , drop = FALSE]
  cat(sprintf("\n--- models %s for stratum (B=%d, n=%d) ---\n",
              if (adjusted) "ADJUSTED" else "UNADJUSTED", nrow(R), n))
  cat(sprintf("  truth RD            = %+.4f   mean est %+.4f (bias %+.4f)\n",
              truthRD, mean(R[, "rd0"]), mean(R[, "rd0"]) - truthRD))
  cat(sprintf("  empirical SD of RD  = %.5f\n", sd(R[, "rd0"])))
  cat(sprintf("  iid SE:     mean %.5f  SD/SE ratio %.3f  coverage %.3f\n",
              mean(R[, "se0"]), sd(R[, "rd0"]) / mean(R[, "se0"]), mean(R[, "cov0"])))
  cat(sprintf("  strata SE:  mean %.5f  SD/SE ratio %.3f  coverage %.3f\n",
              mean(R[, "se1"]), sd(R[, "rd1"]) / mean(R[, "se1"]), mean(R[, "cov1"])))
  invisible(R)
}

Ru <- runCell(adjusted = FALSE)
Ra <- runCell(adjusted = TRUE)
saveRDS(list(unadjusted = Ru, adjusted = Ra), "/tmp/strata-validation.rds")
cat("\nExpectation: UNADJUSTED cell -> iid ratio < 1 (conservative), strata ratio ~ 1;",
    "\n             ADJUSTED cell   -> both ratios ~ 1, SEs nearly equal.\n")
