## End-to-end validation of a bottom PRO tier wired into clinicalPSNB /
## clinicalWinRatio, against a brute-force pairwise truth. Hierarchy: death > KCCQ
## (continuous PRO at a landmark, higher = better, margin delta). The PRO tier
## reach = P(both alive at horizon); within reach the comparison is Y_T vs Y_C+/-d.
## Checks: (1) reach well-behaved & matches truth; (2) clinicalWinRatio WR vs truth;
## (3) clinicalPSNB PSNB/PSWR coverage vs charter-weighted truth.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 1L; B <- 40L; n <- 800L; tau <- 4; tauL <- 4; delta <- 5   # final-visit PRO (landmark = horizon)
alpha <- c(0.5, 0.5)                                   # charter: death 0.5, KCCQ 0.5
bD <- function(W,A) 0.12*exp(0.3*W - 0.6*A)
muY <- function(W,A) 58 + 6*W + 7*A; sdY <- 14
piObs <- function(W) plogis(0.9 + 0.4*W)              # landmark KCCQ observation (MAR)

## ---- brute-force truth: death tier + KCCQ tier among both-alive pairs ----
truth <- local({
  set.seed(1); Ng <- 120000L; M <- 3e6L
  drw <- function(a){ W<-rnorm(Ng); tD<-rexp(Ng,bD(W,a)); Y<-muY(W,a)+rnorm(Ng,0,sdY)
    list(tD=tD, Y=Y) }
  X<-drw(1); Z<-drw(0)
  i<-sample.int(Ng,M,TRUE); j<-sample.int(Ng,M,TRUE)
  tDx<-X$tD[i]; tDy<-Z$tD[j]; Yx<-X$Y[i]; Yy<-Z$Y[j]
  Xa<-tDx>tau; Ya<-tDy>tau
  c1<-ifelse(Xa&!Ya,1L, ifelse(!Xa&Ya,-1L, ifelse(!Xa&!Ya, sign(tDx-tDy), 0L)))
  w1<-mean(c1==1); l1<-mean(c1==-1); D1<-w1-l1
  R2<-(c1==0L); r2<-mean(R2)                           # reach KCCQ = both alive at tau
  w2<-mean(R2 & (Yx > Yy+delta))/r2; l2<-mean(R2 & (Yx < Yy-delta))/r2; D2<-w2-l2
  Wk1<-w1; Lk1<-l1; Wk2<-r2*w2; Lk2<-r2*l2            # tier components
  WR<-(Wk1+Wk2)/(Lk1+Lk2)
  psnb<-alpha[1]*D1+alpha[2]*D2
  pswr<-(alpha[1]*w1+alpha[2]*w2)/(alpha[1]*l1+alpha[2]*l2)
  c(WR=WR, PSNB=psnb, PSWR=pswr, D1=D1, D2=D2, r2=r2, w2=w2, l2=l2)
}); gc(FALSE)

proSpec <- list(marker="kccq", landmark=tauL, margin=delta, direction="higher.better",
                type="continuous", n.grid=40L, label="KCCQ")

simOne <- function(seed) {
  set.seed(seed); A<-rep(0:1,each=n); W<-rnorm(2*n)
  tD<-rexp(2*n,bD(W,A)); C<-rexp(2*n,0.04); obst<-pmin(tD,C,tau)
  aliveL <- tD>tauL & C>=tauL
  kccq <- ifelse(aliveL & runif(2*n)<piObs(W), muY(W,A)+rnorm(2*n,0,sdY), NA)
  dat<-data.frame(arm=A, t_term=obst, died=as.integer(tD<=pmin(C,tau)),
                  W=W, W2=rnorm(2*n), kccq=kccq)
  ar<-list(data=dat, arm="arm", illness.time=character(0), terminal.time="t_term",
           terminal.status="died", covariates=c("W","W2"), horizon=tau, n.grid=30L,
           n.folds=1L, pro=proSpec)
  wr<-suppressMessages(suppressWarnings(do.call(clinicalWinRatio, ar)))
  ps<-suppressMessages(suppressWarnings(do.call(clinicalPSNB, c(ar, list(charter=alpha)))))
  WRr<-as.data.frame(wr)[wr$Estimand=="Win Ratio",]
  pn<-ps[Estimand=="PSNB"]; pw<-ps[Estimand=="PSWR"]; rk2<-ps[Estimand=="Reach[KCCQ]",`Pt Est`]
  c(wr=WRr[["Pt Est"]],
    pn=pn[["Pt Est"]], pnse=pn[["se"]], pncov=as.integer(pn[["CI Low"]]<=truth["PSNB"] & truth["PSNB"]<=pn[["CI Hi"]]),
    pw=pw[["Pt Est"]], pwse=pw[["se"]], pwcov=as.integer(pw[["CI Low"]]<=truth["PSWR"] & truth["PSWR"]<=pw[["CI Hi"]]),
    rk2=rk2)
}
R<-do.call(rbind, parallel::mclapply(seq_len(B), function(b)
  tryCatch({on.exit(gc(FALSE)); simOne(300+b)}, error=function(e) rep(NA,8)), mc.cores=MC))
R<-R[stats::complete.cases(R),,drop=FALSE]
cat(sprintf("\n===== PRO bottom tier end-to-end (%d reps, n=%d/arm, death>KCCQ) =====\n", nrow(R), n))
cat(sprintf("  truth: WR %.3f PSNB %.4f PSWR %.3f (D_death %.3f, D_KCCQ %.3f, reach_KCCQ %.3f)\n",
    truth["WR"], truth["PSNB"], truth["PSWR"], truth["D1"], truth["D2"], truth["r2"]))
cat(sprintf("  reach_KCCQ:  est %.4f bias %+.4f (truth %.4f)\n", mean(R[,"rk2"]), mean(R[,"rk2"])-truth["r2"], truth["r2"]))
cat(sprintf("  WinRatio:    est %.3f bias %+.3f\n", mean(R[,"wr"]), mean(R[,"wr"])-truth["WR"]))
cat(sprintf("  PSNB: mean %.4f bias %+.4f | empSD %.4f meanSE %.4f (SD/SE %.2f) | cover %.3f\n",
    mean(R[,"pn"]), mean(R[,"pn"])-truth["PSNB"], sd(R[,"pn"]), mean(R[,"pnse"]), sd(R[,"pn"])/mean(R[,"pnse"]), mean(R[,"pncov"])))
cat(sprintf("  PSWR: mean %.4f bias %+.4f | empSD %.4f meanSE %.4f (SD/SE %.2f) | cover %.3f\n",
    mean(R[,"pw"]), mean(R[,"pw"])-truth["PSWR"], sd(R[,"pw"]), mean(R[,"pwse"]), sd(R[,"pw"])/mean(R[,"pwse"]), mean(R[,"pwcov"])))
cat("PRO-ENGINE-DONE\n")
