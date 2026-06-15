## Coverage validation of the multistate clinicalRMTIF analytic influence-function
## inference: bias, SE calibration (empirical SD vs mean IF-SE), and 95% CI
## coverage against a brute-force pairwise ground truth. K=2 (death > hosp).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; B <- 80L; n <- 700L; tau <- 5; targT <- 4

bD <- function(W, A) 0.10 * exp(0.3 * W - 0.5 * A)
bH <- function(W, A) 0.20 * exp(0.2 * W - 0.3 * A)

## brute-force truth: integral of (X better than Y) - (worse), worst-event ranking
set.seed(1); Nb <- 4e6
dr <- function(a) { W <- rnorm(Nb); list(tD = rexp(Nb, bD(W, a)), tH = rexp(Nb, bH(W, a))) }
X <- dr(1); Y <- dr(0)
rk <- function(d, t) ifelse(d$tD <= t, 2L, ifelse(d$tH <= t, 1L, 0L))
ts <- seq(0, targT, length.out = 600); dt <- ts[2] - ts[1]; zt <- numeric(length(ts))
for (k in seq_along(ts)) { rx <- rk(X, ts[k]); ry <- rk(Y, ts[k]); zt[k] <- mean(rx < ry) - mean(rx > ry) }
truth <- sum(zt) * dt - 0.5 * dt * (zt[1] + zt[length(zt)])

simOne <- function(seed) {
  set.seed(seed); A <- rep(0:1, each = n); W <- rnorm(2*n)
  tD <- rexp(2*n, bD(W, A)); tH <- rexp(2*n, bH(W, A)); C <- rexp(2*n, 0.05)
  dat <- data.frame(arm = A, t_hosp = ifelse(tH <= pmin(tD, C, tau), tH, NA),
                    t_term = pmin(tD, C, tau), died = as.integer(tD <= pmin(C, tau)),
                    W = W, W2 = rnorm(2*n))
  o <- suppressMessages(suppressWarnings(clinicalRMTIF(dat, arm = "arm", illness.time = "t_hosp",
        terminal.time = "t_term", terminal.status = "died", covariates = c("W", "W2"),
        horizon = targT, n.grid = 40, n.folds = 1)))
  r <- o[Estimand == "RMT-IF"]
  c(est = r[["Pt Est"]], se = r[["se"]],
    cov = as.integer(r[["CI Low"]] <= truth & truth <= r[["CI Hi"]]))
}
R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b)
  tryCatch({ on.exit(gc(FALSE)); simOne(4000 + b) }, error = function(e) rep(NA, 3)), mc.cores = MC))
R <- R[stats::complete.cases(R), , drop = FALSE]
cat(sprintf("\n===== clinicalRMTIF analytic-IF coverage (%d reps, n=%d/arm) =====\n", nrow(R), n))
cat(sprintf("  truth RMT-IF        = %.4f\n", truth))
cat(sprintf("  mean estimate       = %.4f   bias %+.4f\n", mean(R[,"est"]), mean(R[,"est"]) - truth))
cat(sprintf("  empirical SD        = %.4f\n", sd(R[,"est"])))
cat(sprintf("  mean IF-SE          = %.4f   (SD/SE ratio %.3f)\n", mean(R[,"se"]), sd(R[,"est"])/mean(R[,"se"])))
cat(sprintf("  95%% CI coverage     = %.3f\n", mean(R[,"cov"])))
cat("RMTIF-COV-DONE\n")
