## Full validation suite for the RMT-IF estimators. Runs sequentially (one R
## process, MC=2) so the footprint never exceeds ~2 cores. Prints each block as
## it completes. Blocks:
##   A  clinicalRMTIF K=3 coverage vs brute force      [DECISIVE: structural bug test]
##   B  clinicalRMTIF K=2 bias diagnostic (grid x n)   [diagnoses the small bias]
##   C  getRMTIF first-event K=2 coverage              [first-event IF beyond K=1]
##   D  clinicalRMTIF K=2 cross-fitted (n.folds=5)     [EIF holds out-of-fold]
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; tau <- 5; targT <- 4
fl <- function(...) { cat(sprintf(...)); flush(stdout()) }

## ============================ Block A: K=3 multistate ======================
## priority: death > stroke > hosp. Markov DGP (independent latent exp times).
bD <- function(W,A) 0.08*exp(0.3*W - 0.5*A)
bS <- function(W,A) 0.12*exp(0.2*W - 0.3*A)   # stroke (higher-priority non-fatal)
bH <- function(W,A) 0.20*exp(0.2*W - 0.2*A)   # hosp   (lower-priority non-fatal)
set.seed(1); Nb <- 3e6
drk3 <- function(a){ W<-rnorm(Nb); list(tD=rexp(Nb,bD(W,a)), tS=rexp(Nb,bS(W,a)), tH=rexp(Nb,bH(W,a))) }
X3<-drk3(1); Y3<-drk3(0)
rank3 <- function(d,t) ifelse(d$tD<=t, 3L, ifelse(d$tS<=t & d$tS<d$tD, 2L,
                       ifelse(d$tH<=t & d$tH<d$tD, 1L, 0L)))
ts<-seq(0,targT,length.out=400); dt<-ts[2]-ts[1]; zt<-numeric(length(ts))
for(k in seq_along(ts)){rx<-rank3(X3,ts[k]);ry<-rank3(Y3,ts[k]);zt[k]<-mean(rx<ry)-mean(rx>ry)}
truthA <- sum(zt)*dt - 0.5*dt*(zt[1]+zt[length(zt)])
simA <- function(seed,n){ set.seed(seed); A<-rep(0:1,each=n); W<-rnorm(2*n)
  tD<-rexp(2*n,bD(W,A)); tS<-rexp(2*n,bS(W,A)); tH<-rexp(2*n,bH(W,A)); C<-rexp(2*n,0.04)
  obs<-pmin(tD,C,tau)
  dat<-data.frame(arm=A, t_stroke=ifelse(tS<obs, tS, NA), t_hosp=ifelse(tH<obs, tH, NA),
                  t_term=obs, died=as.integer(tD<=pmin(C,tau)), W=W, W2=rnorm(2*n))
  o<-suppressMessages(suppressWarnings(clinicalRMTIF(dat,arm="arm",
       illness.time=c("t_stroke","t_hosp"), terminal.time="t_term", terminal.status="died",
       covariates=c("W","W2"), horizon=targT, n.grid=40, n.folds=1)))
  r<-o[Estimand=="RMT-IF"]; c(est=r[["Pt Est"]], se=r[["se"]],
    cov=as.integer(r[["CI Low"]]<=truthA & truthA<=r[["CI Hi"]])) }
BA<-120L; nA<-700L
RA<-do.call(rbind,parallel::mclapply(seq_len(BA),function(b)tryCatch({on.exit(gc(FALSE));simA(100+b,nA)},
   error=function(e)rep(NA,3)),mc.cores=MC)); RA<-RA[stats::complete.cases(RA),,drop=FALSE]
fl("\n=== BLOCK A: clinicalRMTIF K=3 (death>stroke>hosp), %d reps, n=%d/arm ===\n",nrow(RA),nA)
fl("  truth %.4f | mean est %.4f (bias %+.4f) | empSD %.4f | meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
   truthA, mean(RA[,"est"]), mean(RA[,"est"])-truthA, sd(RA[,"est"]), mean(RA[,"se"]),
   sd(RA[,"est"])/mean(RA[,"se"]), mean(RA[,"cov"]))
fl("BLOCK-A-DONE\n")

## ============================ Block B: K=2 bias diagnostic =================
bD2<-function(W,A)0.10*exp(0.3*W-0.5*A); bH2<-function(W,A)0.20*exp(0.2*W-0.3*A)
set.seed(2)
drk2<-function(a){W<-rnorm(Nb);list(tD=rexp(Nb,bD2(W,a)),tH=rexp(Nb,bH2(W,a)))}
X2<-drk2(1);Y2<-drk2(0)
rank2<-function(d,t)ifelse(d$tD<=t,2L,ifelse(d$tH<=t & d$tH<d$tD,1L,0L))
zt<-numeric(length(ts)); for(k in seq_along(ts)){rx<-rank2(X2,ts[k]);ry<-rank2(Y2,ts[k]);zt[k]<-mean(rx<ry)-mean(rx>ry)}
truthB<-sum(zt)*dt-0.5*dt*(zt[1]+zt[length(zt)])
simB<-function(seed,n,ng){set.seed(seed);A<-rep(0:1,each=n);W<-rnorm(2*n)
  tD<-rexp(2*n,bD2(W,A));tH<-rexp(2*n,bH2(W,A));C<-rexp(2*n,0.05);obs<-pmin(tD,C,tau)
  dat<-data.frame(arm=A,t_hosp=ifelse(tH<obs,tH,NA),t_term=obs,died=as.integer(tD<=pmin(C,tau)),W=W,W2=rnorm(2*n))
  o<-suppressMessages(suppressWarnings(clinicalRMTIF(dat,arm="arm",illness.time="t_hosp",
     terminal.time="t_term",terminal.status="died",covariates=c("W","W2"),horizon=targT,n.grid=ng,n.folds=1)))
  r<-o[Estimand=="RMT-IF"];c(est=r[["Pt Est"]],se=r[["se"]],cov=as.integer(r[["CI Low"]]<=truthB&truthB<=r[["CI Hi"]]))}
cellB<-function(n,ng,B,seed0){R<-do.call(rbind,parallel::mclapply(seq_len(B),function(b)
  tryCatch({on.exit(gc(FALSE));simB(seed0+b,n,ng)},error=function(e)rep(NA,3)),mc.cores=MC))
  R<-R[stats::complete.cases(R),,drop=FALSE]
  fl("  n=%d grid=%d: bias %+.4f | SD/SE %.3f | cover %.3f (%d reps)\n",n,ng,mean(R[,"est"])-truthB,
     sd(R[,"est"])/mean(R[,"se"]),mean(R[,"cov"]),nrow(R))}
fl("\n=== BLOCK B: K=2 bias diagnostic (truth %.4f) ===\n",truthB)
cellB(700,40,70,200); cellB(1500,80,70,400)
fl("BLOCK-B-DONE\n")

## ============================ Block C: getRMTIF first-event ================
## competing risks (1=death,2=hosp); first-event ranking. brute-force truth.
l1f<-function(W,A)0.14*exp(0.3*W-0.5*A); l2f<-function(W,A)0.16*exp(0.2*W-0.3*A)
set.seed(3)
W<-rnorm(Nb);lt<-l1f(W,1)+l2f(W,1);Tt<-rexp(Nb,lt);Jt<-1L+rbinom(Nb,1,l2f(W,1)/lt)
W0<-rnorm(Nb);lt0<-l1f(W0,0)+l2f(W0,0);Tc<-rexp(Nb,lt0);Jc<-1L+rbinom(Nb,1,l2f(W0,0)/lt0)
fav<-c(`1`=2,`2`=1)  # death worse(2), hosp(1); event-free=0
rkfe<-function(T,J,t)ifelse(T<=t,fav[as.character(J)],0L)
zt<-numeric(length(ts));for(k in seq_along(ts)){rx<-rkfe(Tt,Jt,ts[k]);ry<-rkfe(Tc,Jc,ts[k]);zt[k]<-mean(rx<ry)-mean(rx>ry)}
truthC<-sum(zt)*dt-0.5*dt*(zt[1]+zt[length(zt)])
simC<-function(seed,n){set.seed(seed);A<-rep(0:1,each=n);W<-rnorm(2*n)
  lt<-l1f(W,A)+l2f(W,A);Te<-rexp(2*n,lt);cause<-1L+rbinom(2*n,1,l2f(W,A)/lt);C<-rexp(2*n,0.06)
  dat<-data.table(id=1:(2*n),time=pmin(Te,C,tau),status=ifelse(Te<=pmin(C,tau),cause,0L),trt=A,W=W,W2=rnorm(2*n))
  a<-suppressMessages(formatArguments(DataTable=dat,EventTime="time",EventType="status",Treatment="trt",
     ID="id",Intervention=0:1,TargetTime=c(1,2,3,4),TargetEvent=c(1,2),CVArg=list(V=2),MaxUpdateIter=8,Verbose=FALSE,
     Model=list(trt="SL.mean","0"=list(Cox=survival::Surv(time,status==0)~.),
       "1"=list(Cox=survival::Surv(time,status==1)~.),"2"=list(Cox=survival::Surv(time,status==2)~.))))
  o<-suppressMessages(getRMTIF(suppressMessages(suppressWarnings(doConcrete(a))),Horizon=4,Intervention=c(2,1),TargetEvent=c(1,2)))
  r<-o[Estimand=="RMT-IF"];c(est=r[["Pt Est"]],se=r[["se"]],cov=as.integer(r[["CI Low"]]<=truthC&truthC<=r[["CI Hi"]]))}
BC<-120L;nC<-600L
RC<-do.call(rbind,parallel::mclapply(seq_len(BC),function(b)tryCatch({on.exit(gc(FALSE));simC(700+b,nC)},
   error=function(e)rep(NA,3)),mc.cores=MC));RC<-RC[stats::complete.cases(RC),,drop=FALSE]
fl("\n=== BLOCK C: getRMTIF first-event K=2, %d reps, n=%d/arm (truth %.4f) ===\n",nrow(RC),nC,truthC)
fl("  mean est %.4f (bias %+.4f) | empSD %.4f | meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
   mean(RC[,"est"]),mean(RC[,"est"])-truthC,sd(RC[,"est"]),mean(RC[,"se"]),sd(RC[,"est"])/mean(RC[,"se"]),mean(RC[,"cov"]))
fl("BLOCK-C-DONE\n")

## ============================ Block D: K=2 cross-fitted ====================
BD<-50L;nD<-700L
RD<-do.call(rbind,parallel::mclapply(seq_len(BD),function(b)tryCatch({on.exit(gc(FALSE))
  set.seed(900+b);A<-rep(0:1,each=nD);W<-rnorm(2*nD)
  tD<-rexp(2*nD,bD2(W,A));tH<-rexp(2*nD,bH2(W,A));C<-rexp(2*nD,0.05);obs<-pmin(tD,C,tau)
  dat<-data.frame(arm=A,t_hosp=ifelse(tH<obs,tH,NA),t_term=obs,died=as.integer(tD<=pmin(C,tau)),W=W,W2=rnorm(2*nD))
  o<-suppressMessages(suppressWarnings(clinicalRMTIF(dat,arm="arm",illness.time="t_hosp",
     terminal.time="t_term",terminal.status="died",covariates=c("W","W2"),horizon=targT,n.grid=40,n.folds=5)))
  r<-o[Estimand=="RMT-IF"];c(est=r[["Pt Est"]],se=r[["se"]],cov=as.integer(r[["CI Low"]]<=truthB&truthB<=r[["CI Hi"]]))},
  error=function(e)rep(NA,3)),mc.cores=MC));RD<-RD[stats::complete.cases(RD),,drop=FALSE]
fl("\n=== BLOCK D: clinicalRMTIF K=2 CROSS-FITTED (n.folds=5), %d reps, n=%d/arm ===\n",nrow(RD),nD)
fl("  bias %+.4f | SD/SE %.3f | cover %.3f\n",mean(RD[,"est"])-truthB,sd(RD[,"est"])/mean(RD[,"se"]),mean(RD[,"cov"]))
fl("SUITE-DONE\n")
