## TRISCEND II-style hierarchical win ratio: can concrete analyze the win ratio
## "as done in TRISCEND II"? The TRISCEND II pivotal trial (Hahn et al., NEJM 2025)
## primary endpoint is a 7-tier hierarchical composite analyzed by win ratio
## (WR 2.02, 95% CI 1.56-2.62, driven mainly by the quality-of-life tiers):
##
##   1 death (any cause)              [terminal]
##   2 RV-assist device / transplant  [time-to-event, rare]
##   3 tricuspid-valve reintervention [time-to-event, rare]
##   4 HF hospitalization            [time-to-event]
##   5 KCCQ-OS improvement >= 10 pts [continuous PRO, landmark]
##   6 NYHA improvement >= 1 class    [ordinal PRO, landmark]
##   7 6-min walk improvement >= 30 m [continuous PRO, landmark]
##
## We validate concrete's clinicalWinRatio() over this exact hierarchy (4 hard
## tiers via the multistate engine + 3 PRO tiers) against a brute-force pairwise
## win ratio truth, and show (a) the engine recovers the WR and its per-tier
## decomposition, (b) the PRO tiers dominate (as in TRISCEND II), and (c) a
## mortality-weighted charter (PSNB) re-weights the composite.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
tau <- 1.0; tauL <- 1.0                              # 1-year horizon; PROs at 1 year
NF <- c("rvad","tvi","hfh")                          # non-fatal tiers, priority order
EVB <- c(rvad=1L, tvi=2L, hfh=4L)
nfRate <- c(rvad=0.012, tvi=0.040, hfh=0.16)         # rarer = higher priority
bAe <- c(rvad=-0.25, tvi=-0.5, hfh=-0.25)            # device lowers non-fatal rates
d0 <- 0.07; morb <- 1.6; bWd <- 0.35; bAd <- -0.25   # death: rises w/ morbidity, device protects
## PRO change-from-baseline at landmark (device improves; the WR drivers)
muKCCQ <- function(W,A) 5 + 4*W + 16*A;  sdKCCQ <- 18
muNYHA <- function(W,A) 0.4 + 0.2*W + 0.8*A; sdNYHA <- 0.9   # improvement (classes), >=0
muWALK <- function(W,A) 6 + 8*W + 32*A;  sdWALK <- 60
piObs  <- function(W) plogis(1.4 + 0.3*W)            # landmark PRO observation (MAR)

popcount <- function(s) sum(bitwAnd(bitwShiftR(s, 0:7), 1L))
## Gillespie full-history sim over [0,tau]: death rate rises with accrued morbidity
simOne <- function(W, A) {
  s <- 0L; t <- 0; tev <- c(rvad=Inf, tvi=Inf, hfh=Inf); tD <- Inf
  repeat {
    rD <- d0 * morb^popcount(s) * exp(bWd*W + bAd*A)
    rates <- c(D = rD)
    for (e in NF) if (bitwAnd(s, EVB[[e]])==0L) rates[e] <- nfRate[[e]]*exp(0.3*W + bAe[[e]]*A)
    tot <- sum(rates); t <- t + rexp(1, tot); if (t > tau) break
    pick <- names(rates)[which(runif(1)*tot <= cumsum(rates))[1]]
    if (pick=="D") { tD <- t; break }
    tev[[pick]] <- t; s <- bitwOr(s, EVB[[pick]])
  }
  c(tD=tD, tev)
}
simPop <- function(n, A) {
  W <- rnorm(n)
  H <- t(vapply(seq_len(n), function(i) simOne(W[i], A), numeric(4)))
  aliveL <- H[,"tD"] > tauL
  data.table(W=W, A=A, tD=H[,"tD"], t_rvad=H[,"rvad"], t_tvi=H[,"tvi"], t_hfh=H[,"hfh"],
             aliveL=aliveL,
             kccq=ifelse(aliveL, muKCCQ(W,A)+rnorm(n,0,sdKCCQ), NA),
             nyha=ifelse(aliveL, pmax(0, round(muNYHA(W,A)+rnorm(n,0,sdNYHA))), NA),
             walk=ifelse(aliveL, muWALK(W,A)+rnorm(n,0,sdWALK), NA))
}

## ---- brute-force 7-tier pairwise win ratio (uncensored complete-data truth) ----
cap <- function(x) ifelse(x <= tau, x, Inf)
truth <- local({
  set.seed(1); Np <- 60000L; M <- 4e6L
  X <- simPop(Np, 1L); Z <- simPop(Np, 0L)
  it <- sample(Np, M, TRUE); ic <- sample(Np, M, TRUE)
  win <- loss <- dec <- logical(M)
  Wk <- Lk <- numeric(7); labs <- c("death", NF, "KCCQ","NYHA","6MWD")
  hardCols <- c("tD","t_rvad","t_tvi","t_hfh")
  step <- function(x, y, k, margin=0) {            # later/larger wins; decided = won|lost
    u <- !dec
    w <- u & (x > y + margin); l <- u & (x < y - margin)   # Inf>Inf=FALSE => both-event-free ties
    Wk[k] <<- mean(w); Lk[k] <<- mean(l)
    win <<- win | w; loss <<- loss | l; dec <<- dec | w | l
  }
  for (k in seq_along(hardCols)) {                   # hard tiers: later/no event wins
    cc <- hardCols[k]; step(cap(X[[cc]][it]), cap(Z[[cc]][ic]), k, 0) }
  step(X$kccq[it], Z$kccq[ic], 5, 10)                # PRO tiers (only undecided=both reachers)
  step(X$nyha[it], Z$nyha[ic], 6, 0)
  step(X$walk[it], Z$walk[ic], 7, 30)
  list(WR = sum(Wk)/sum(Lk), Pwin = sum(Wk), Ploss = sum(Lk),
       Wk = setNames(Wk, labs), Lk = setNames(Lk, labs), reachKCCQ = mean(X$aliveL)*mean(Z$aliveL))
}); gc(FALSE)
cat(sprintf("\n===== TRISCEND II-style truth (brute-force pairwise, 7 tiers) =====\n"))
cat(sprintf("  WR = %.3f  (Pwin %.3f, Ploss %.3f)\n", truth$WR, truth$Pwin, truth$Ploss))
cat("  per-tier win - loss contribution (W^k - L^k):\n")
for (k in 1:7) cat(sprintf("    %-6s  W=%.3f L=%.3f  net=%+.3f\n",
  names(truth$Wk)[k], truth$Wk[k], truth$Lk[k], truth$Wk[k]-truth$Lk[k]))

proSpecs <- list(
  list(marker="kccq", landmark=tauL, margin=10, direction="higher.better", type="continuous", n.grid=30L, label="KCCQ"),
  list(marker="nyha", landmark=tauL, margin=0,  direction="higher.better", type="ordinal",            label="NYHA"),
  list(marker="walk", landmark=tauL, margin=30, direction="higher.better", type="continuous", n.grid=30L, label="6MWD"))

makeTrial <- function(seed, n1, n0) {
  set.seed(seed)
  X <- simPop(n1, 1L); Z <- simPop(n0, 0L); D <- rbind(X, Z)
  C <- rexp(nrow(D), 0.05)                            # light censoring
  obst <- pmin(D$tD, C, tau)
  obsNF <- function(tt) ifelse(tt <= pmin(C, tau, D$tD), tt, NA)
  mk <- function(v) ifelse(D$aliveL & C >= tau & runif(nrow(D)) < piObs(D$W), v, NA)  # 1-yr visit
  data.frame(arm=D$A, t_rvad=obsNF(D$t_rvad), t_tvi=obsNF(D$t_tvi), t_hfh=obsNF(D$t_hfh),
             t_term=obst, died=as.integer(D$tD <= pmin(C, tau)),
             W1=D$W, W2=rnorm(nrow(D)),
             kccq=mk(D$kccq), nyha=mk(D$nyha), walk=mk(D$walk))
}

runEngine <- function(dat) suppressMessages(suppressWarnings(clinicalWinRatio(
  dat, arm="arm", illness.time=c("t_rvad","t_tvi","t_hfh"), terminal.time="t_term",
  terminal.status="died", covariates=c("W1","W2"), horizon=tau, n.grid=24L, n.folds=1L,
  pro=proSpecs)))

## ---- single fit: WR + decomposition, at TRISCEND-like 2:1 allocation ----
cat(sprintf("\n===== concrete clinicalWinRatio, single fit (n=400 device / 200 control) =====\n"))
o <- runEngine(makeTrial(11, 400L, 200L))
print(o[, c("Estimand","Pt Est","se","CI Low","CI Hi","pValue")])

## ---- coverage over reps (validate WR point + inference vs truth) ----
B <- 30L
R <- do.call(rbind, lapply(seq_len(B), function(b) tryCatch({
  o <- runEngine(makeTrial(200+b, 400L, 200L)); wr <- as.data.frame(o)[o$Estimand=="Win Ratio",]
  c(wr=wr[["Pt Est"]], se=wr[["se"]],
    cov=as.integer(wr[["CI Low"]] <= truth$WR & truth$WR <= wr[["CI Hi"]]))
}, error=function(e) rep(NA,3))))
R <- R[stats::complete.cases(R),,drop=FALSE]
cat(sprintf("\n===== WR coverage vs brute-force truth (%d reps, n=400/200) =====\n", nrow(R)))
cat(sprintf("  truth WR %.3f | est %.3f bias %+.3f | empSD %.3f meanSE %.3f | cover %.3f\n",
    truth$WR, mean(R[,1]), mean(R[,1])-truth$WR, sd(R[,1]), mean(R[,2]), mean(R[,3])))

## ---- charter (PSNB): re-weight the composite toward mortality ----
cat(sprintf("\n===== mortality-weighted charter (PSNB), single fit =====\n"))
ch <- c(0.40, 0.15, 0.15, 0.15, 0.05, 0.05, 0.05)    # death 0.40 vs reach-driven default
ps <- suppressMessages(suppressWarnings(clinicalPSNB(
  makeTrial(11, 400L, 200L), arm="arm", illness.time=c("t_rvad","t_tvi","t_hfh"),
  terminal.time="t_term", terminal.status="died", covariates=c("W1","W2"),
  charter=ch, horizon=tau, n.grid=24L, n.folds=1L, pro=proSpecs)))
print(ps[Estimand %in% c("PSNB","PSWR"), c("Estimand","Pt Est","se","CI Low","CI Hi")])
cat("TRISCEND-DONE\n")
