## Validate the productionized clinicalWinRatio(): K=2 (back-compat) + K=4 smokes,
## and a K=3 package-level coverage check (full path incl. SL + output assembly).
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
source("scripts/genwr-engine.R")

## build user-facing columns from the DGP (death > E2 > ... > E_K), horizon tau
makeData <- function(K, n, cens, seed) {
  configK(K); configDGP(K); set.seed(seed)
  D <- rbind(simArm(n, 1L, cens), simArm(n, 0L, cens))
  termT <- pmin(D$tD, D$C); died <- as.integer(is.finite(D$tD) & D$tD <= D$C)
  dat <- data.frame(arm = D$a, t_term = pmin(termT, tau),
                    died = ifelse(termT > tau, 0L, died), W = D$W)
  ill <- character(0)
  for (e in seq_len(K - 1L)) { col <- paste0("ill", e)
    dat[[col]] <- ifelse(is.finite(D[[paste0("t", e)]]) & D[[paste0("t", e)]] <= termT, D[[paste0("t", e)]], NA)
    ill <- c(ill, col) }
  list(dat = dat, ill = ill)
}
runWR <- function(dat, ill, folds = 1L) clinicalWinRatio(dat, arm = "arm", illness.time = ill,
  terminal.time = "t_term", terminal.status = "died", covariates = "W",
  horizon = 3, n.grid = 40L, n.folds = folds)

cat("--- smokes ---\n")
for (K in c(2L, 4L)) {
  md <- makeData(K, 1500L, 0.12, 10 + K)
  r <- runWR(md$dat, md$ill)
  cat(sprintf("K=%d (illness cols: %s): WR=%.3f CI=[%.3f,%.3f]  P(win+loss+tie)=%.3f\n",
      K, paste(md$ill, collapse = ","), r$`Pt Est`[1], r$`CI Low`[1], r$`CI Hi`[1],
      sum(r$`Pt Est`[4:6])))
}

cat("\n--- K=3 package coverage (B reps, n.folds=1) ---\n")
configK(3L); configDGP(3L); set.seed(3)
WRpop <- mean(vapply(1:5, function(b) bruteWR(simArm(12000L,1L,0), simArm(12000L,0L,0), 3e6)["WR"], numeric(1)))
B <- 60L
R <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) tryCatch({
  md <- makeData(3L, 1200L, 0.12, 2000 + b); r <- runWR(md$dat, md$ill)
  wr <- r[r$Estimand == "Win Ratio", ]
  c(wr$`Pt Est`, wr$se, as.integer(wr$`CI Low` <= WRpop & WRpop <= wr$`CI Hi`))
}, error = function(e) c(NA, NA, NA)), mc.cores = 4L))
R <- R[!is.na(R[,1]), , drop = FALSE]
cat(sprintf("pop WR=%.4f | reps=%d mean WR=%.4f emp.SD=%.4f EIF-SE=%.4f (ratio %.2f) coverage=%.3f\n",
    WRpop, nrow(R), mean(R[,1]), sd(R[,1]), mean(R[,2]), mean(R[,2])/sd(R[,1]), mean(R[,3])))
