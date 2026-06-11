## Validate targetWinRatio() against getWinRatio() with brute-force ground truth.
## Cell A (sparse): K=2 hierarchy, TargetTime = c(2, 4) only -> the plug-in win
##   integral is a 2-point Riemann sum; the direct method integrates the full
##   event grid. Expect plug-in biased / under-covered, direct ~nominal.
## Cell B (dense): K=1 single event, dense grid -> both should be ~unbiased;
##   checks whether direct targeting improves the ~0.92 win-odds/net-benefit
##   coverage seen for the plug-in at n=500.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; B <- 120L; n <- 500L; tau <- 5; targT <- 4

## ---- DGP (shared): K=2 competing events, treatment protective ----
simDat <- function(n, seed, K2 = TRUE) {
  set.seed(seed); N <- 2L * n
  W <- rnorm(N); A <- rep(0:1, each = n)
  if (K2) { l1 <- 0.10 * exp(0.3 * W - 0.5 * A); l2 <- 0.15 * exp(0.2 * W - 0.3 * A) }
  else    { l1 <- 0.20 * exp(0.3 * W - 0.5 * A); l2 <- 0 * W }
  lt <- l1 + l2
  Te <- rexp(N, lt); cause <- 1L + rbinom(N, 1L, ifelse(lt > 0, l2 / lt, 0))
  C <- rexp(N, 0.08)
  time <- pmin(Te, C, tau); status <- ifelse(Te <= pmin(C, tau), cause, 0L)
  data.table(id = seq_len(N), time = time, status = status, trt = A, W = W, W2 = rnorm(N))
}

## ---- brute-force truth: pairwise hierarchical comparison at horizon targT ----
## first events (T, J) per arm, no censoring; rules: event-free beats any event;
## different priorities -> the HIGHER-priority (smaller J) first event loses;
## same priority -> later event time wins.
bruteTruth <- function(K2, Npairs = 4e6, seed = 1) {
  set.seed(seed)
  draw <- function(a) {
    W <- rnorm(Npairs)
    if (K2) { l1 <- 0.10 * exp(0.3 * W - 0.5 * a); l2 <- 0.15 * exp(0.2 * W - 0.3 * a) }
    else    { l1 <- 0.20 * exp(0.3 * W - 0.5 * a); l2 <- 0 * W }
    lt <- l1 + l2
    Te <- rexp(Npairs, lt); J <- 1L + rbinom(Npairs, 1L, ifelse(lt > 0, l2 / lt, 0))
    list(T = Te, J = J)
  }
  trt <- draw(1); ctl <- draw(0)
  ti <- trt$T; tj <- ctl$T; ji <- trt$J; jj <- ctl$J
  freeI <- ti > targT; freeJ <- tj > targT
  win <- (freeI & !freeJ) |
    (!freeI & !freeJ & ji > jj) |                      # treated's event lower priority
    (!freeI & !freeJ & ji == jj & ti > tj)             # same priority, treated later
  loss <- (freeJ & !freeI) |
    (!freeI & !freeJ & jj > ji) |
    (!freeI & !freeJ & ji == jj & tj > ti)
  c(Pwin = mean(win), Ploss = mean(loss), WR = mean(win) / mean(loss),
    NB = mean(win) - mean(loss))
}

fitBoth <- function(dat, K2, sparse) {
  TT <- if (sparse) c(2, 4) else seq(0.5, 4, by = 0.5)
  TE <- if (K2) c(1, 2) else 1
  Mdl <- c(list(trt = "SL.mean", "0" = list(Cox = survival::Surv(time, status == 0) ~ .),
                "1" = list(Cox = survival::Surv(time, status == 1) ~ .)),
           if (K2) list("2" = list(Cox = survival::Surv(time, status == 2) ~ .)))
  a <- suppressMessages(formatArguments(DataTable = dat, EventTime = "time", EventType = "status",
        Treatment = "trt", ID = "id", Intervention = 0:1, TargetTime = TT, TargetEvent = TE,
        CVArg = list(V = 2), MaxUpdateIter = 12, Verbose = FALSE, Model = Mdl))
  e <- suppressMessages(suppressWarnings(doConcrete(a)))
  trtIdx <- match("A=1", names(e)); ctlIdx <- match("A=0", names(e))
  pl <- suppressMessages(suppressWarnings(
    getWinRatio(e, Horizon = targT, Intervention = c(trtIdx, ctlIdx), TargetEvent = TE)))
  dr <- suppressMessages(suppressWarnings(
    targetWinRatio(e, Horizon = targT, Intervention = c(trtIdx, ctlIdx), TargetEvent = TE)))
  g <- function(o, est) { r <- as.data.frame(o); r <- r[r$Estimand == est, ]
    c(r[["Pt Est"]][1], r[["CI Low"]][1], r[["CI Hi"]][1]) }
  c(plWR = g(pl, "Win Ratio"), plNB = g(pl, "Net Benefit"),
    drWR = g(dr, "Win Ratio"), drNB = g(dr, "Net Benefit"),
    drCv = as.numeric(attr(dr, "WRConverged")))
}

runCell <- function(K2, sparse, label) {
  tr <- bruteTruth(K2)
  R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) tryCatch({
    on.exit(gc(FALSE)); fitBoth(simDat(n, 5000 + b, K2), K2, sparse)
  }, error = function(e) rep(NA_real_, 13)), mc.cores = MC))
  R <- R[stats::complete.cases(R), , drop = FALSE]
  covg <- function(lo, hi, tru) mean(R[, lo] <= tru & tru <= R[, hi])
  cat(sprintf("\n--- %s (B=%d, n=%d/arm) truth: WR=%.3f NB=%.3f ---\n",
              label, nrow(R), n, tr["WR"], tr["NB"]))
  cat(sprintf("  WR  plug-in: mean %.3f  bias %+.3f  coverage %.3f\n",
              mean(R[, "plWR1"]), mean(R[, "plWR1"]) - tr["WR"], covg("plWR2", "plWR3", tr["WR"])))
  cat(sprintf("  WR  direct : mean %.3f  bias %+.3f  coverage %.3f\n",
              mean(R[, "drWR1"]), mean(R[, "drWR1"]) - tr["WR"], covg("drWR2", "drWR3", tr["WR"])))
  cat(sprintf("  NB  plug-in: mean %+.3f  bias %+.3f  coverage %.3f\n",
              mean(R[, "plNB1"]), mean(R[, "plNB1"]) - tr["NB"], covg("plNB2", "plNB3", tr["NB"])))
  cat(sprintf("  NB  direct : mean %+.3f  bias %+.3f  coverage %.3f\n",
              mean(R[, "drNB1"]), mean(R[, "drNB1"]) - tr["NB"], covg("drNB2", "drNB3", tr["NB"])))
  cat(sprintf("  direct convergence rate: %.2f\n", mean(R[, "drCv"])))
  invisible(R)
}

Ra <- runCell(K2 = TRUE, sparse = TRUE, "Cell A: K=2 hierarchy, SPARSE 2-time grid")
Rb <- runCell(K2 = FALSE, sparse = FALSE, "Cell B: K=1, dense grid")
saveRDS(list(sparseK2 = Ra, denseK1 = Rb), "/tmp/targetwr-validation.rds")
cat("\nDONE\n")
