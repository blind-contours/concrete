## General-K analytic EIF coverage (oracle nuisances) at K=3 and K=4, validating
## the generalized assembly's inference. Coverage of the brute-force population WR.
suppressWarnings(suppressMessages(library(data.table)))
source("scripts/genwr-eif.R")
setGrid(40L)
args <- as.integer(commandArgs(trailingOnly = TRUE))
B  <- if (length(args) >= 1) args[1] else 120L
nA <- if (length(args) >= 2) args[2] else 1500L
censRate <- 0.12

runK <- function(K) {
  configK(K); configDGP(K)
  set.seed(50 + K)
  WRpop <- mean(vapply(1:5, function(b) bruteWR(simArm(12000L,1L,0), simArm(12000L,0L,0), 3e6)["WR"], numeric(1)))
  oneRep <- function(b) tryCatch({
    set.seed(4000 + 100*K + b)
    trt <- armSetup(simArm(nA, 1L, censRate), censRate)
    ctl <- armSetup(simArm(nA, 0L, censRate), censRate)
    e <- winRatioEIF(trt, ctl)
    c(e["WR"], e["seWR"], cov = as.integer(e["lo"] <= WRpop & WRpop <= e["hi"]))
  }, error = function(er) c(NA, NA, NA))
  R <- do.call(rbind, parallel::mclapply(seq_len(B), oneRep, mc.cores = 4L))
  R <- R[!is.na(R[,1]), , drop = FALSE]
  cat(sprintf("K=%d (%d tiers): pop WR=%.4f | reps=%d mean WR=%.4f emp.SD=%.4f EIF-SE=%.4f (ratio %.2f) coverage=%.3f\n",
      K, K, WRpop, nrow(R), mean(R[,1]), sd(R[,1]), mean(R[,2]), mean(R[,2])/sd(R[,1]), mean(R[,3])))
}
cat(sprintf("General-K oracle EIF coverage (B=%d, n=%d/arm, ~28%% cens):\n", B, nA))
runK(3L)
runK(4L)
