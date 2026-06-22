## Validate clinicalWinRatio's covariate-adjusted CROSSOVER IPCW: does censoring
## treatment-switchers at the switch time + a crossover hazard recover the
## HYPOTHETICAL NO-SWITCHING win ratio? DGP: death > KCCQ (landmark = horizon).
## OMT subjects switch to the device with an INFORMATIVE hazard (sicker switch
## sooner); after switching they get the device death + KCCQ benefit. The ITT win
## ratio is diluted; the crossover-IPCW estimand should recover the no-switching
## truth (EVOQUE vs OMT had no one switched).
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet=TRUE); library(data.table) }))
data.table::setDTthreads(1L)
tau <- 3; delta <- 5; n <- 2000L; B <- 25L
lamE <- function(W) 0.10*exp(0.3*W)            # EVOQUE death hazard
lamO <- function(W) 0.20*exp(0.3*W)            # OMT (no-switch) death hazard
gamX <- function(W) 0.18*exp(0.6*W)            # OMT->device switch hazard (informative: sicker switch)
muKE <- function(W) 20 + 5*W; muKO <- function(W) 5 + 5*W; sdK <- 15

## no-switching truth: EVOQUE vs OMT with NO switching
truth <- local({ set.seed(1); Ng<-3e5L; M<-3e6L
  W1<-rnorm(Ng); tE<-rexp(Ng,lamE(W1)); YE<-muKE(W1)+rnorm(Ng,0,sdK)
  W0<-rnorm(Ng); tO<-rexp(Ng,lamO(W0)); YO<-muKO(W0)+rnorm(Ng,0,sdK)
  aE<-tE>tau; aO<-tO>tau                                   # alive at horizon
  i<-sample(Ng,M,TRUE); j<-sample(Ng,M,TRUE)
  c1<-ifelse(aE[i]&!aO[j],1L,ifelse(!aE[i]&aO[j],-1L,ifelse(!aE[i]&!aO[j],sign(tE[i]-tO[j]),0L)))
  R<-c1==0L; w1<-mean(c1==1);l1<-mean(c1==-1)
  w2<-mean(R & YE[i]>YO[j]+delta)/mean(R); l2<-mean(R & YE[i]<YO[j]-delta)/mean(R); r2<-mean(R)
  Wk1<-w1;Lk1<-l1;Wk2<-r2*w2;Lk2<-r2*l2; c(WR=(Wk1+Wk2)/(Lk1+Lk2)) }); gc(FALSE)

simOne <- function(seed){ set.seed(seed); N<-2L*n; A<-rep(0:1,each=n); W<-rnorm(N)
  ## EVOQUE arm
  tD <- numeric(N); Y <- numeric(N); sw <- rep(Inf,N)
  e <- A==1; tD[e]<-rexp(sum(e),lamE(W[e])); Y[e]<-muKE(W[e])+rnorm(sum(e),0,sdK)
  ## OMT arm: no-switch death tO; switch time s; if s<tO & s<tau -> device after s
  o<-which(A==0); Wo<-W[o]; tO<-rexp(length(o),lamO(Wo)); s<-rexp(length(o),gamX(Wo))
  switched <- s < pmin(tO, tau)
  tDo <- ifelse(switched, s + rexp(length(o), lamE(Wo)), tO)    # post-switch device survival
  Yo  <- ifelse(switched & s<tau, muKE(Wo)+rnorm(length(o),0,sdK), muKO(Wo)+rnorm(length(o),0,sdK))
  tD[o]<-tDo; Y[o]<-Yo; sw[o]<-ifelse(switched, s, Inf)
  C<-rexp(N,0.04); obst<-pmin(tD,C,tau); died<-as.integer(tD<=pmin(C,tau))
  reach<-tD>tau & C>=tau
  kccq<-ifelse(reach, Y, NA)                                   # KCCQ at horizon among survivors
  dat<-data.frame(arm=A,t_term=obst,died=died,W=W,W2=rnorm(N),kccq=kccq,xover=ifelse(is.finite(sw),sw,NA))
  pro<-list(list(marker="kccq",landmark=tau,margin=delta,direction="higher.better",type="continuous",label="KCCQ"))
  g<-function(o) as.data.frame(o)[o$Estimand=="Win Ratio","Pt Est"]
  itt<-g(suppressMessages(suppressWarnings(clinicalWinRatio(dat,arm="arm",illness.time=character(0),terminal.time="t_term",terminal.status="died",covariates=c("W","W2"),horizon=tau,n.grid=24L,n.folds=1L,SL.library=c("SL.mean","SL.glm"),pro=pro))))
  nos<-g(suppressMessages(suppressWarnings(clinicalWinRatio(dat,arm="arm",illness.time=character(0),terminal.time="t_term",terminal.status="died",covariates=c("W","W2"),horizon=tau,n.grid=24L,n.folds=1L,SL.library=c("SL.mean","SL.glm"),crossover="xover",pro=pro))))
  c(itt=itt, nos=nos)
}
R<-do.call(rbind, lapply(seq_len(B), function(b) tryCatch(simOne(100+b), error=function(e) c(NA,NA))))
R<-R[stats::complete.cases(R),,drop=FALSE]
cat(sprintf("\n===== crossover-IPCW win-ratio validation (%d reps, n=%d/arm) =====\n", nrow(R), n))
cat(sprintf("  NO-SWITCHING truth WR             = %.3f\n", truth["WR"]))
cat(sprintf("  ITT (with switching, diluted)     = %.3f  (bias vs no-switch %+.3f)\n", mean(R[,"itt"]), mean(R[,"itt"])-truth["WR"]))
cat(sprintf("  crossover-IPCW (no-switching est) = %.3f  (bias vs no-switch %+.3f)\n", mean(R[,"nos"]), mean(R[,"nos"])-truth["WR"]))
cat("XOVER-VALID-DONE\n")
