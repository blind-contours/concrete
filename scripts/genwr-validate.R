## Validate the K-general engine: oracle occupancy engine vs brute-force pairwise
## win ratio, at K=3 (regression) and K=4 (new). Light/single-threaded.
suppressWarnings(suppressMessages(library(data.table)))
source("scripts/genwr-engine.R")
setGrid(60L)

checkK <- function(K, nrep = 6L, nBrute = 6000L, nPop = 30000L, npairs = 3e6) {
  configK(K); configDGP(K)
  set.seed(100 + K)
  ePop <- engineWR(data.table(W = rnorm(nPop), a = 1L), data.table(W = rnorm(nPop), a = 0L))
  wr <- vapply(seq_len(nrep), function(b) {
    set.seed(900 + 10 * K + b)
    bruteWR(simArm(nBrute, 1L, 0), simArm(nBrute, 0L, 0), npairs)["WR"] }, numeric(1))
  cat(sprintf("K=%d: ENGINE pop WR=%.4f | BRUTE mean WR=%.4f (sd %.4f, %d reps n=%d) | dWR=%+.4f\n",
              K, ePop["WR"], mean(wr), sd(wr), nrep, nBrute, ePop["WR"] - mean(wr)))
}
checkK(3L)
checkK(4L)
