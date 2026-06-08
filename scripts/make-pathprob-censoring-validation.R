## Path-probability one-step Theta_hat WITH random censoring + IPCW (G != 1).
## Validates the 1/G clever-covariate factor by coverage of the MC truth.
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
cli <- commandArgs(trailingOnly=TRUE); B <- if(length(cli)>=1) as.integer(cli[1]) else 40L
n <- if(length(cli)>=2) as.integer(cli[2]) else 2000L
useTrueG <- length(cli)>=4 && cli[4]=="trueG"
tau <- 1500; gl <- if(length(cli)>=3) as.integer(cli[3]) else 81L
grid <- seq(0, tau, length.out=gl); M <- length(grid)-1L; Mp1<-M+1L
starts<-grid[-Mp1]
b01<-7e-4; b02<-4e-4; b12<-1.0e-3; g01<-0.3; g02<-0.2; g12<-0.4
cc <- 3e-4                                                        # censoring rate (~36% by tau)
Gmin <- 0.05

simID <- function(n,seed){ set.seed(seed); W<-rnorm(n)
  a01<-b01*exp(g01*W);a02<-b02*exp(g02*W);a12<-b12*exp(g12*W)
  T01<-rexp(n,a01);T02<-rexp(n,a02);u0<-pmin(T01,T02);hfh<-T01<T02
  dpost<-ifelse(hfh,u0+rexp(n,a12),Inf); D<-ifelse(hfh,dpost,u0)  # overall death
  C<-rexp(n,cc); cap<-pmin(C,tau)
  d01<-as.integer(hfh&u0<cap); t01<-ifelse(d01==1,u0,Inf)
  d02<-as.integer(!hfh&u0<cap); t02<-ifelse(d02==1,u0,Inf)
  d12<-as.integer(hfh&u0<cap&dpost<cap); t12<-ifelse(d12==1,dpost,Inf)
  exit0<-pmin(u0,cap); entry1<-ifelse(d01==1,u0,Inf); exit1<-ifelse(d01==1,pmin(dpost,cap),Inf)
  obsT<-pmin(D,C,tau); censE<-as.integer(C<D & C<tau)
  data.table(W=W,d01=d01,t01=t01,d02=d02,t02=t02,d12=d12,t12=t12,exit0=exit0,entry1=entry1,exit1=exit1,obsT=obsT,censE=censE) }

oneStep <- function(d){ n<-nrow(d); W<-d$W; Cov<-data.frame(W=W); SLlib<-c("SL.mean","SL.glm")
  f01<-fitTransitionSL(rep(0,n),d$exit0,d$d01,Cov,grid,SL.library=SLlib,V=5)
  f02<-fitTransitionSL(rep(0,n),d$exit0,d$d02,Cov,grid,SL.library=SLlib,V=5)
  sub<-d$d01==1
  f12<-fitTransitionSL(d$entry1[sub],d$exit1[sub],d$d12[sub],data.frame(W=W[sub]),grid,SL.library=SLlib,V=5)
  fc <-fitTransitionSL(rep(0,n),d$obsT,d$censE,Cov,grid,SL.library=SLlib,V=5)      # censoring hazard
  I01<-predictTransitionSL(f01,Cov);I02<-predictTransitionSL(f02,Cov);I12<-predictTransitionSL(f12,Cov)
  IC <-predictTransitionSL(fc,Cov)
  cur<-multistateCurves(I01,I02,I12); S0<-cur$S0; pimat<-cur$pi; S12tau<-cur$S12toTau
  Glag<-if(useTrueG) matrix(exp(-cc*grid[1:M]), M, n) else pmax(rbind(1,exp(-apply(IC,2,cumsum)))[1:M,,drop=FALSE], Gmin)
  S0end<-S0[2:Mp1,,drop=FALSE]
  Pi<-colSums(pimat); Pigt<-apply(pimat,2,function(c) rev(cumsum(rev(c)))-c)
  ## IPCW clever covariates (divide by G)
  phi01<-((S0end*S12tau - Pigt)/S0end)/Glag; phi02<-(-Pigt/S0end)/Glag; phi12<-(-S12tau)/Glag
  ## martingale increments (Y already censoring-adjusted via exit0/exit1)
  Y0<-outer(starts,d$exit0,"<")
  inIv<-function(tt){o<-matrix(0L,M,n);k<-findInterval(tt,grid);ok<-is.finite(tt)&k>=1&k<=M;o[cbind(k[ok],which(ok))]<-1L;o}
  dM01<-inIv(d$t01)-Y0*I01; dM02<-inIv(d$t02)-Y0*I02
  Y1<-outer(starts,d$entry1,">=")&outer(starts,d$exit1,"<"); dM12<-inIv(d$t12)-Y1*I12
  xi<-colSums(phi01*dM01)+colSums(phi02*dM02)+colSums(phi12*dM12)
  os<-mean(Pi)+mean(xi); eif<-(Pi-os)+xi; c(os=os, se=sqrt(var(eif)/n)) }

## NOTE: truth = P(in state 1 at tau) in the FULL (uncensored) law. Recompute cleanly:
TRtruth<-{set.seed(99);N<-2e6;W<-rnorm(N);a01<-b01*exp(g01*W);a02<-b02*exp(g02*W);a12<-b12*exp(g12*W)
  T01<-rexp(N,a01);T02<-rexp(N,a02);u0<-pmin(T01,T02);hfh<-T01<T02;dpost<-ifelse(hfh,u0+rexp(N,a12),Inf)
  mean(hfh & u0<tau & dpost>=tau)}
cat(sprintf("TRUTH Theta=%.4f ; B=%d n=%d cc=%.0e (~%.0f%% cens)\n", TRtruth, B, n, cc, 100*(1-exp(-cc*tau))))
one<-function(seed) tryCatch({r<-oneStep(simID(n,seed)); data.table(os=r["os"],se=r["se"],
  cover=as.integer(r["os"]-1.96*r["se"]<=TRtruth & TRtruth<=r["os"]+1.96*r["se"]))}, error=function(e) NULL)
res<-data.table::rbindlist(Filter(Negate(is.null), parallel::mclapply(1:B+8000L, one, mc.cores=min(5L,parallel::detectCores()-1L))))
cat(sprintf("RESULT reps=%d mean(os)=%.4f bias=%.4f | coverage=%.3f | mean(se)=%.4f emp.sd=%.4f\n",
  nrow(res),mean(res$os),mean(res$os)-TRtruth,mean(res$cover),mean(res$se),sd(res$os)))
