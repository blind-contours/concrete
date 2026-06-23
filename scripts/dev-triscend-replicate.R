## Replicate the TRISCEND II primary result on a cohort CALIBRATED to the
## published 1-year numbers (Hahn et al., NEJM 2025; Quality-of-Life paper JACC
## 2024). We do NOT have the patient-level trial data; this generates a synthetic
## trial matching the published marginal inputs and asks whether concrete's
## clinicalWinRatio over the exact 7-tier hierarchy recovers the reported
## win ratio 2.02 (95% CI 1.56-2.62) and its quality-of-life-driven structure.
##
## Calibration targets (published):
##   1-yr all-cause mortality   12.6% device / 15.2% control
##   KCCQ-OS change             +22 device / +4 control      (margin >= 10)
##   NYHA improvement           93% vs 34% NYHA I/II          (margin >= 1 class)
##   6-min walk change          +11 device / -20 control m    (margin >= 30 m)
##   HFH                        no significant difference
##   RV-assist/transplant, TV reintervention  rare; device fewer reinterventions
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
n1 <- 267L; n0 <- 133L; tau <- 1.0                  # TRISCEND II: ~400 randomized 2:1, 1-year
NF <- c("rvad","tvi","hfh"); EVB <- c(rvad=1L, tvi=2L, hfh=4L)
## hazards calibrated to the published 1-year marginal rates
lamD <- function(W,A) 0.158 * 1.2^0 * exp(0.30*W - 0.20*A)           # death: 15.2% / 12.6%
lamRV <- function(W,A) 0.015 * exp(0.2*W)                            # RVAD/transplant ~1.5%
lamTV <- function(W,A) 0.051 * exp(0.2*W - 0.93*A)                   # reintervention 5% / 2%
lamHF <- function(W,A) 0.151 * exp(0.3*W - 0.05*A)                   # HFH ~14%, no sig diff
rate <- function(W,A,e) switch(e, rvad=lamRV(W,A), tvi=lamTV(W,A), hfh=lamHF(W,A))
## final-visit PRO change-from-baseline, calibrated to published deltas
muKCCQ <- function(W,A) 4 + 18*A + 3*W;  sdK <- 20                   # +22 / +4
muNYHA <- function(W,A) 0.35 + 0.95*A + 0.2*W; sdN <- 0.9           # ord improvement (classes)
muWALK <- function(W,A) -20 + 31*A + 5*W; sdL <- 70                  # +11 / -20 m
piObs  <- function(W) plogis(1.6 + 0.3*W)                            # ~90% have 1-yr visit

popcount <- function(s) sum(bitwAnd(bitwShiftR(s,0:7),1L))
simOne <- function(W,A){ s<-0L; t<-0; tev<-c(rvad=Inf,tvi=Inf,hfh=Inf); tD<-Inf
  repeat{ rD<-lamD(W,A)*1.2^popcount(s); rs<-c(D=rD)
    for(e in NF) if(bitwAnd(s,EVB[[e]])==0L) rs[e]<-rate(W,A,e)
    tot<-sum(rs); t<-t+rexp(1,tot); if(t>tau) break
    pick<-names(rs)[which(runif(1)*tot<=cumsum(rs))[1]]; if(pick=="D"){tD<-t;break}
    tev[[pick]]<-t; s<-bitwOr(s,EVB[[pick]]) }
  c(tD=tD, tev) }
simTrial <- function(seed){
  set.seed(seed); A<-c(rep(1L,n1),rep(0L,n0)); N<-n1+n0; W<-rnorm(N)
  H<-t(vapply(seq_len(N), function(i) simOne(W[i],A[i]), numeric(4)))
  C<-rexp(N,0.05); obst<-pmin(H[,"tD"],C,tau); reach<-H[,"tD"]>tau & C>=tau
  oN<-function(x) ifelse(x<=pmin(C,tau,H[,"tD"]), x, NA)
  mk<-function(v) ifelse(reach & runif(N)<piObs(W), v, NA)
  data.frame(arm=A, t_rvad=oN(H[,"rvad"]), t_tvi=oN(H[,"tvi"]), t_hfh=oN(H[,"hfh"]),
    t_term=obst, died=as.integer(H[,"tD"]<=pmin(C,tau)), W1=W, W2=rnorm(N),
    kccq=mk(muKCCQ(W,A)+rnorm(N,0,sdK)),
    nyha=mk(pmax(0,round(muNYHA(W,A)+rnorm(N,0,sdN)))),
    walk=mk(muWALK(W,A)+rnorm(N,0,sdL))) }
pro <- list(
  list(marker="kccq", landmark=tau, margin=10, direction="higher.better", type="continuous", label="KCCQ"),
  list(marker="nyha", landmark=tau, margin=0,  direction="higher.better", type="ordinal",            label="NYHA"),
  list(marker="walk", landmark=tau, margin=30, direction="higher.better", type="continuous", label="6MWD"))
fit <- function(dat) suppressMessages(suppressWarnings(clinicalWinRatio(dat, arm="arm",
  illness.time=c("t_rvad","t_tvi","t_hfh"), terminal.time="t_term", terminal.status="died",
  covariates=c("W1","W2"), horizon=tau, n.grid=24L, n.folds=1L, pro=pro)))

cat("\n===== ONE TRISCEND II-calibrated synthetic trial (n=267 device / 133 control) =====\n")
o <- fit(simTrial(101))
print(o[, c("Estimand","Pt Est","se","CI Low","CI Hi","pValue")])
cat(sprintf("\n  Published TRISCEND II: Win Ratio 2.02 (95%% CI 1.56-2.62)\n"))

## distribution of the win ratio across calibrated synthetic trials
B <- 100L
WR <- vapply(seq_len(B), function(b) tryCatch({
  d<-fit(simTrial(200+b)); as.data.frame(d)[d$Estimand=="Win Ratio","Pt Est"] },
  error=function(e) NA_real_), numeric(1))
WR <- WR[is.finite(WR)]
cat(sprintf("\n===== Win ratio across %d calibrated synthetic trials =====\n", length(WR)))
cat(sprintf("  mean WR %.2f | median %.2f | 2.5-97.5%% range %.2f-%.2f\n",
    mean(WR), median(WR), quantile(WR,0.025), quantile(WR,0.975)))
cat(sprintf("  published point 2.02 inside synthetic 95%% range: %s\n",
    if (quantile(WR,0.025) <= 2.02 & 2.02 <= quantile(WR,0.975)) "YES" else "no"))
cat("TRISCEND-REPLICATE-DONE\n")
