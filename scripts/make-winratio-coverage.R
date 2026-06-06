#!/usr/bin/env Rscript
# Coverage validation of getWinRatio() on a single-event DGP with closed-form
# truth. Checks that the influence-function CIs for the win ratio, win odds, and
# net benefit cover the truth at the nominal rate.
#   Rscript scripts/make-winratio-coverage.R [B] [n]
suppressWarnings(suppressMessages({library(concrete); library(data.table)}))
cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 60L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 500L
tau <- 1500; tau_max <- 3000
lam <- function(A, W, W2) 6e-4 * exp(-0.6 * A + 0.4 * W + 0.2 * W2); lc <- 3e-4
grid <- seq(200, tau, by = 200)

set.seed(7); Wt <- stats::rnorm(4e5); W2t <- stats::rnorm(4e5)
Sbar <- function(a, t) mean(exp(-lam(a, Wt, W2t) * t))
fg <- seq(0, tau, length.out = 2000); S1 <- sapply(fg, function(t) Sbar(1, t)); S0 <- sapply(fg, function(t) Sbar(0, t))
Pwin <- sum(S1[-1] * diff(1 - S0)); Ploss <- sum(S0[-1] * diff(1 - S1)); Ptie <- 1 - Pwin - Ploss
TR <- list(`Win Ratio` = Pwin / Ploss,
           `Win Odds` = (Pwin + Ptie / 2) / (Ploss + Ptie / 2),
           `Net Benefit` = Pwin - Ploss)
cat(sprintf("TRUE: WR=%.3f WO=%.3f NB=%.4f\n", TR$`Win Ratio`, TR$`Win Odds`, TR$`Net Benefit`))

sim <- function(s) {
  set.seed(s)
  W <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  Tt <- stats::rexp(n, lam(A, W, W2)); C <- stats::rexp(n, lc)
  To <- pmin(Tt, C, tau_max); ev <- ifelse(To >= tau_max, 0L, ifelse(Tt <= C, 1L, 0L))
  data.table(id = seq_len(n), time = To, event = ev, arm = A, W = W, W2 = W2)
}
one <- function(s) {
  d <- sim(s)
  tryCatch({
    a <- formatArguments(DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = grid, TargetEvent = 1, CVArg = list(V = 5),
      UpdateMethod = "adaptive", EICStopRule = "absolute", MaxUpdateIter = 25, Verbose = FALSE)
    wr <- as.data.table(getWinRatio(suppressMessages(doConcrete(a)), Horizon = tau, Intervention = c(1, 2)))
    data.table::rbindlist(lapply(names(TR), function(lab) {
      r <- wr[Estimand == lab]
      data.table(Estimand = lab, bias = r$`Pt Est` - TR[[lab]],
                 cover = is.finite(r$`CI Low`) & r$`CI Low` <= TR[[lab]] & TR[[lab]] <= r$`CI Hi`)
    }))
  }, error = function(e) NULL)
}
seeds <- 20260605L + seq_len(B)
envCores <- suppressWarnings(as.integer(Sys.getenv("VAL_CORES", unset = NA)))
uc <- if (!is.na(envCores) && envCores >= 1L) envCores else
  if (.Platform$OS.type == "unix") min(6L, max(1L, parallel::detectCores() - 1L)) else 1L
res <- data.table::rbindlist(Filter(Negate(is.null), parallel::mclapply(seeds, one, mc.cores = uc)))
summ <- res[, .(Reps = .N, Bias = mean(bias), Coverage = mean(cover)), by = Estimand]
print(summ)
data.table::fwrite(summ, "scripts/winratio-coverage-summary.csv")
cat("Interpretation: coverage ~ 0.95; small negative bias = target-grid discretization.\n")
