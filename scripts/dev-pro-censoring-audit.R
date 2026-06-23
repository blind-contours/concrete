## AUDIT: does the PRO block correctly IPCW-correct for censoring before the
## horizon? The conditional PRO net benefit and the (censoring-free) reach are
## invariant to the censoring rate; a correct estimator must reproduce them at
## every rate. Hierarchy: death > KCCQ (final-visit PRO, landmark = horizon).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
tau <- 4; delta <- 5; n <- 2000L; B <- 15L
bD <- function(W, A) 0.12 * exp(0.3 * W - 0.6 * A)
muY <- function(W, A) 58 + 6 * W + 7 * A; sdY <- 14
piObs <- function(W) plogis(1.2 + 0.4 * W)

## censoring-free truth: reach = P(both alive at tau); KCCQ net benefit among reachers
truth <- local({
  set.seed(1); Ng <- 2e5L; M <- 3e6L
  drw <- function(a){ W <- rnorm(Ng); tD <- rexp(Ng, bD(W, a)); Y <- muY(W, a) + rnorm(Ng, 0, sdY)
    list(al = tD > tau, Y = Y) }
  X <- drw(1); Z <- drw(0); rk <- mean(X$al) * mean(Z$al)
  ax <- which(X$al); ay <- which(Z$al); i <- sample(ax, M, TRUE); j <- sample(ay, M, TRUE)
  w <- mean(X$Y[i] > Z$Y[j] + delta); l <- mean(X$Y[i] < Z$Y[j] - delta)
  c(reach = rk, NB = w - l)                                  # conditional KCCQ net benefit
}); gc(FALSE)

run1 <- function(seed, crate) {
  set.seed(seed); A <- rep(0:1, each = n); W <- rnorm(2 * n)
  tD <- rexp(2 * n, bD(W, A)); C <- if (crate <= 0) rep(Inf, 2 * n) else rexp(2 * n, crate)
  obst <- pmin(tD, C, tau); reach <- tD > tau & C >= tau
  kccq <- ifelse(reach & runif(2 * n) < piObs(W), muY(W, A) + rnorm(2 * n, 0, sdY), NA)
  dat <- data.frame(arm = A, t_term = obst, died = as.integer(tD <= pmin(C, tau)),
                    W = W, W2 = rnorm(2 * n), kccq = kccq)
  o <- suppressMessages(suppressWarnings(clinicalPSNB(dat, arm = "arm",
        illness.time = character(0), terminal.time = "t_term", terminal.status = "died",
        covariates = c("W", "W2"), charter = "reach", horizon = tau, n.grid = 24L, n.folds = 1L,
        pro = list(marker = "kccq", landmark = tau, margin = delta,
                   direction = "higher.better", type = "continuous", label = "KCCQ"))))
  c(reach = o[Estimand == "Reach[KCCQ]", `Pt Est`], NB = o[Estimand == "NetBenefit[KCCQ]", `Pt Est`])
}

cat(sprintf("\n===== PRO censoring audit (death>KCCQ, n=%d/arm, %d reps) =====\n", n, B))
cat(sprintf("  truth: reach %.4f  NetBenefit[KCCQ] %.4f (invariant to censoring)\n", truth["reach"], truth["NB"]))
for (crate in c(0, 0.05, 0.15)) {
  R <- do.call(rbind, lapply(seq_len(B), function(b) tryCatch(run1(10 * crate * 100 + b, crate),
        error = function(e) c(NA, NA))))
  R <- R[stats::complete.cases(R), , drop = FALSE]
  pcens <- if (crate <= 0) 0 else 1 - exp(-crate * tau)
  cat(sprintf("  cens rate %.2f (~%.0f%% pre-horizon): reach %.4f (bias %+.4f) | NB[KCCQ] %.4f (bias %+.4f)\n",
      crate, 100 * pcens, mean(R[, "reach"]), mean(R[, "reach"]) - truth["reach"],
      mean(R[, "NB"]), mean(R[, "NB"]) - truth["NB"]))
}
cat("PRO-CENS-DONE\n")
