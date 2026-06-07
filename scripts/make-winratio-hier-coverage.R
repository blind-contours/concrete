#!/usr/bin/env Rscript
# Coverage validation of the HIERARCHICAL (prioritized competing-risk) win ratio.
# Truth is computed by brute force: simulate the true marginal first-event law
# (T, J) under each arm, then apply the exact prioritized pairwise rule the
# estimator targets (death/cause-1 > cause-2). Confirms the curve-based estimator
# and its influence-function CIs recover that pairwise truth.
#   Rscript scripts/make-winratio-hier-coverage.R [B] [n]
suppressWarnings(suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION"))
    devtools::load_all(".", quiet = TRUE) else library(concrete)
  library(data.table)
}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 60L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 600L
tau <- 1500; tau_max <- 3000
gridBy <- suppressWarnings(as.numeric(Sys.getenv("GRID_BY", unset = "150")))
grid <- seq(gridBy, tau, by = gridBy)             # denser grid -> smaller discretization bias
prio <- c(1, 2)                                   # priority order: cause 1 > cause 2

# cause-specific hazards (two competing causes); treatment lowers severe cause 1
lam1 <- function(A, W, W2) 6e-4 * exp(-0.6 * A + 0.4 * W)
lam2 <- function(A, W, W2) 5e-4 * exp( 0.1 * A + 0.3 * W2)
lc   <- 3e-4

draw_arm <- function(a, M) {                      # true marginal (T, J) under arm a
  W <- stats::rnorm(M); W2 <- stats::rnorm(M)
  T1 <- stats::rexp(M, lam1(a, W, W2)); T2 <- stats::rexp(M, lam2(a, W, W2))
  Tt <- pmin(T1, T2); J <- ifelse(T1 < T2, 1L, 2L)
  J[Tt > tau] <- 0L                               # event-free by horizon
  list(T = Tt, J = J)
}
# rank: position in priority order (1 = highest priority / most severe); event-free = best
rank_of <- function(J) ifelse(J == 0L, length(prio) + 1L, match(J, prio))
# treated beats control under the prioritized rule (same rule the estimator targets)
pair_outcome <- function(t1, j1, t0, j0) {
  r1 <- rank_of(j1); r0 <- rank_of(j0)
  win  <- (r1 > r0) | (r1 == r0 & r1 <= length(prio) & t1 > t0)
  loss <- (r0 > r1) | (r1 == r0 & r1 <= length(prio) & t0 > t1)
  list(win = mean(win), loss = mean(loss))
}
set.seed(20260606); M <- 4e6
tr <- draw_arm(1, M); ct <- draw_arm(0, M)        # independent treated & control draws
po <- pair_outcome(tr$T, tr$J, ct$T, ct$J)
Pwin <- po$win; Ploss <- po$loss; Ptie <- 1 - Pwin - Ploss
TR <- list(`Win Ratio` = Pwin / Ploss,
           `Win Odds`  = (Pwin + Ptie/2) / (Ploss + Ptie/2),
           `Net Benefit` = Pwin - Ploss)
cat(sprintf("TRUE (brute force): P(win)=%.4f P(loss)=%.4f P(tie)=%.4f | WR=%.3f WO=%.3f NB=%.4f\n",
            Pwin, Ploss, Ptie, TR$`Win Ratio`, TR$`Win Odds`, TR$`Net Benefit`))

sim <- function(s) {                              # one trial-like data set with censoring
  set.seed(s)
  W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  T1 <- stats::rexp(n, lam1(A, W, W2)); T2 <- stats::rexp(n, lam2(A, W, W2))
  Tt <- pmin(T1, T2); J <- ifelse(T1 < T2, 1L, 2L)
  C <- stats::rexp(n, lc); To <- pmin(Tt, C, tau_max)
  ev <- ifelse(To >= tau_max, 0L, ifelse(Tt <= C, J, 0L))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W = W, W2 = W2)
}
one <- function(s) {
  d <- sim(s)
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = grid, TargetEvent = c(1, 2),
      CVArg = list(V = 5), UpdateMethod = "adaptive", EICStopRule = "absolute",
      MaxUpdateIter = 25, Verbose = FALSE)
    wr <- as.data.table(getWinRatio(suppressMessages(doConcrete(a)),
                                    Horizon = tau, Intervention = c(1, 2), TargetEvent = c(1, 2)))
    data.table::rbindlist(lapply(names(TR), function(lab) {
      r <- wr[Estimand == lab]
      data.table(Estimand = lab, bias = r$`Pt Est` - TR[[lab]],
                 cover = is.finite(r$`CI Low`) & r$`CI Low` <= TR[[lab]] & TR[[lab]] <= r$`CI Hi`)
    }))
  }, error = function(e) NULL)
}
seeds <- 20260606L + seq_len(B)
envCores <- suppressWarnings(as.integer(Sys.getenv("VAL_CORES", unset = NA)))
uc <- if (!is.na(envCores) && envCores >= 1L) envCores else
  if (.Platform$OS.type == "unix") min(4L, max(1L, parallel::detectCores() - 1L)) else 1L
res <- data.table::rbindlist(Filter(Negate(is.null), parallel::mclapply(seeds, one, mc.cores = uc)))
summ <- res[, .(Reps = .N, Bias = mean(bias), Coverage = mean(cover)), by = Estimand]
print(summ)
data.table::fwrite(summ, "scripts/winratio-hier-coverage-summary.csv")
cat("Interpretation: coverage ~ 0.95; small bias from target-grid discretization.\n")
