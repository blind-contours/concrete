## Clinical win-ratio EIF variance: per-subject influence function assembled from
## validated pieces (D_a, D_Theta, level-1 single-event win on overall survival
## S^D=S0+p1, bilinear B), two-arm 1/pi weighting, delta method. Validate by
## coverage of the brute-force pairwise WR truth. v1: G==1 (admin censor at tau).
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
cli <- commandArgs(trailingOnly=TRUE); B <- if(length(cli)>=1) as.integer(cli[1]) else 1L
n <- if(length(cli)>=2) as.integer(cli[2]) else 3000L
tau <- 1500; gl <- if(length(cli)>=3) as.integer(cli[3]) else 41L
grid <- seq(0, tau, length.out = gl); M <- length(grid)-1L; Mp1 <- M+1L
ctrl<-list(b01=7e-4,b02=4e-4,b12=1.0e-3,g01=0.3,g02=0.2,g12=0.4)
trt <-list(b01=5e-4,b02=2.5e-4,b12=6e-4,g01=0.3,g02=0.2,g12=0.4)

simArm <- function(n,p,seed){ set.seed(seed); W<-rnorm(n)
  a01<-p$b01*exp(p$g01*W);a02<-p$b02*exp(p$g02*W);a12<-p$b12*exp(p$g12*W)
  T01<-rexp(n,a01);T02<-rexp(n,a02);u0<-pmin(T01,T02);hfh<-T01<T02
  s<-ifelse(hfh,u0,Inf);dpost<-ifelse(hfh,u0+rexp(n,a12),Inf);D<-ifelse(hfh,dpost,u0)
  est<-data.table(W=W,d01=as.integer(hfh&u0<tau),t01=ifelse(hfh&u0<tau,u0,Inf),
    d02=as.integer(!hfh&u0<tau),t02=ifelse(!hfh&u0<tau,u0,Inf),
    d12=as.integer(hfh&dpost<tau),t12=ifelse(hfh&dpost<tau,dpost,Inf),
    exit0=pmin(u0,tau),entry1=ifelse(hfh&u0<tau,s,Inf),exit1=ifelse(hfh&u0<tau,pmin(dpost,tau),Inf))
  bf<-data.table(D=ifelse(D<=tau,D,Inf),s=ifelse(hfh&s<=tau,s,Inf),state=ifelse(D<=tau,2L,ifelse(hfh&s<=tau,1L,0L)))
  list(est=est,bf=bf) }
bfWR <- function(bT,bC){ DT<-bT$D;DC<-bC$D;sT<-bT$s;sC<-bC$s;stT<-bT$state;stC<-bC$state
  w1<-DC<=tau&DT>DC;l1<-DT<=tau&DC>DT;t1<-!(w1|l1);rk<-function(st)ifelse(st==0L,3L,2L)
  w2<-t1&(rk(stT)>rk(stC)|(rk(stT)==rk(stC)&stT==1L&sT>sC));l2<-t1&(rk(stC)>rk(stT)|(rk(stT)==rk(stC)&stT==1L&sC>sT))
  Pw<-mean(w1|w2);Pl<-mean(l1|l2);c(Pwin=Pw,Ploss=Pl,WR=Pw/Pl) }

## ---- per-arm: fit hazards, curves, per-subject dM, and all building-block IFs ----
armFit <- function(est){
  n<-nrow(est);W<-est$W;Cov<-data.frame(W=W);SLlib<-c("SL.mean","SL.glm")
  f01<-fitTransitionSL(rep(0,n),est$exit0,est$d01,Cov,grid,SL.library=SLlib,V=5)
  f02<-fitTransitionSL(rep(0,n),est$exit0,est$d02,Cov,grid,SL.library=SLlib,V=5)
  sub<-est$d01==1
  f12<-fitTransitionSL(est$entry1[sub],est$exit1[sub],est$d12[sub],data.frame(W=W[sub]),grid,SL.library=SLlib,V=5)
  I01<-predictTransitionSL(f01,Cov);I02<-predictTransitionSL(f02,Cov);I12<-predictTransitionSL(f12,Cov)
  cur<-multistateCurves(I01,I02,I12); S0<-cur$S0; p1<-cur$p1; SD<-cur$SD; pimat<-cur$pi
  L12<-rbind(0,apply(I12,2,cumsum))
  ## martingale increments (M x n)
  Y0<-outer(grid[-Mp1],est$exit0,"<")
  inIv<-function(tt){o<-matrix(0L,M,n);k<-findInterval(tt,grid);ok<-is.finite(tt)&k>=1&k<=M;o[cbind(k[ok],which(ok))]<-1L;o}
  dM01<-inIv(est$t01)-Y0*I01; dM02<-inIv(est$t02)-Y0*I02; dM0<-dM01+dM02
  Y1<-outer(grid[-Mp1],est$entry1,">=")&outer(grid[-Mp1],est$exit1,"<"); dM12<-inIv(est$t12)-Y1*I12
  S0end<-S0[2:Mp1,,drop=FALSE]                                   # S0 at interval ends (M x n)
  ## --- IF of a = S0(tau) (per subject scalar) ---
  abar<-mean(S0[Mp1,]); Da<-(S0[Mp1,]-abar) - S0[Mp1,]*colSums(dM0/S0end)
  ## --- IF of Theta = p1(tau): use phi for Pi(tau) ---
  S12tau<-cur$S12toTau                                           # S12(t_k,tau) (M x n)
  Pi<-colSums(pimat); Pigt<-apply(pimat,2,function(c) rev(cumsum(rev(c)))-c)
  phi01<-(S0end*S12tau - Pigt)/S0end; phi02<- -Pigt/S0end; phi12<- -S12tau
  xiTheta<-colSums(phi01*dM01)+colSums(phi02*dM02)+colSums(phi12*dM12)
  Tbar<-mean(Pi); DTheta<-(Pi-Tbar)+xiTheta
  ## --- IF curves of S^D(t_j), j=1..Mp1 (M+1 x n): D_S0 (O(M)) + D_p1 (O(M^2)) ---
  DS0<-matrix(0,Mp1,n); cumdM0<-rbind(0,apply(dM0/S0end,2,cumsum))  # cumulative to grid j
  for(j in 1:Mp1) DS0[j,]<-(S0[j,]-mean(S0[j,])) - S0[j,]*cumdM0[j,]
  Dp1<-matrix(0,Mp1,n)
  for(j in 2:Mp1){ H<-j-1                                         # horizon = interval index H (t_j=grid[j])
    S12toH<-exp(-(matrix(L12[j,],M,n,byrow=TRUE)-L12[2:Mp1,,drop=FALSE]))  # S12(t_k,t_j)
    piH<-cur$S0mid*I01*S12toH; if(H<M) piH[(H+1):M,]<-0                   # enter by t_j, survive to t_j
    PiH<-colSums(piH); PigtH<-apply(piH,2,function(c) rev(cumsum(rev(c)))-c)
    p01<-(S0end*S12toH-PigtH)/S0end; p02<- -PigtH/S0end; p12<- -S12toH
    if(H<M){p01[(H+1):M,]<-0;p02[(H+1):M,]<-0;p12[(H+1):M,]<-0}           # only intervals u<=t_j
    xi<-colSums(p01*dM01)+colSums(p02*dM02)+colSums(p12*dM12)
    Dp1[j,]<-(PiH-mean(PiH))+xi }
  DSD<- -(DS0+Dp1)                                                # IF of F^D = 1 - S^D  (note: D_FD = -D_SD)
  list(n=n, SDbar=rowMeans(SD), DFD=DSD, a=abar, Da=Da, Theta=Tbar, DTheta=DTheta,
       hbar=rowMeans(pimat), pimat=pimat, S0mid=cur$S0mid, I01=I01, dM01=dM01,dM02=dM02,dM12=dM12,
       S0end=S0end, L12=L12, S0=S0)
}

## single-event win EIF on overall survival (mirror existing getWinRatio coef logic)
## win1 = sum SDmid_T * dFC ; returns per-subject Dwin1_T (treated) and Dwin1_C (control)
levelOneIF <- function(A, B){  # A=treated arm fit, B=control arm fit
  SDt<-A$SDbar; SDc<-B$SDbar    # length Mp1 (overall survival)
  Ft<-1-SDt; Fc<-1-SDc; dFc<-diff(Fc); dFt<-diff(Ft)           # increments length M
  SmT<-sqrt(SDt[1:M]*SDt[2:Mp1])                               # midpoint S^D_T
  ## P(win1)=sum SmT*dFc ; gradient wrt F_C(t_k): coef on dFc ; wrt S_T: SmT
  ## treated contributes via D[F^D_T]; control via D[F^D_C]
  ## use the single-event win gradient: Dwin1 = sum_k SmT[k]*dD_FC[k] + sum_k ( -nextSmT? )...
  ## Build exactly as integral: win1 = sum_k SmT[k]*(Fc[k+1]-Fc[k]); linear in F_C and in S_T.
  ## d/dF_C(t_j): appears in dFc[j-1] (coef +SmT[j-1]) and dFc[j] (coef -SmT[j]).
  ## d/dS^D_T enters via SmT; approximate gradient wrt F_T at grid as the mirror.
  ## treated-side gradient (vary S^D_T): win1 = sum_k SmT[k]*dFc[k]; SmT[k]=sqrt(SDt[k]SDt[k+1]).
  ##   dSmT[k]/dSDt[k]=0.5 SmT[k]/SDt[k]; /dSDt[k+1]=0.5 SmT[k]/SDt[k+1]. and dFD_T=-dSD_T.
  ## control-side gradient (vary F_C): coef_j = SmT[j-1]-SmT[j] (with edges).
  cFc<-numeric(Mp1)
  for(j in 1:Mp1){ a<- if(j-1>=1&&j-1<=M) SmT[j-1] else 0; b<- if(j>=1&&j<=M) SmT[j] else 0; cFc[j]<-a-b }
  Dwin1 <- as.numeric(crossprod(cFc, B$DFD))                    # control contributes via D[F^D_C]
  ## treated contributes via D[S^D_T]= -D[F^D_T]; coef on SDt[j]:
  cSDt<-numeric(Mp1)
  for(k in 1:M){ cSDt[k]<-cSDt[k]+0.5*SmT[k]/SDt[k]*dFc[k]; cSDt[k+1]<-cSDt[k+1]+0.5*SmT[k]/SDt[k+1]*dFc[k] }
  Dwin1 <- Dwin1 + as.numeric(crossprod(cSDt, -A$DFD))         # D[S^D_T] = -D[F^D_T]
  list(point=sum(SmT*dFc), Dt=cSDt, Dc=cFc)                    # store coef vectors for assembly
}

## ---- assemble clinical WR + per-subject IF, two-arm 1/pi weighting ----
oneRep <- function(seed_t, seed_c){
  AT<-armFit(simArm(n,trt,seed_t)$est); AC<-armFit(simArm(n,ctrl,seed_c)$est)
  piT<-0.5; piC<-0.5                                            # equal allocation (pooled n/arm)
  ## level 1 win/loss point + coefficient vectors
  L1w<-levelOneIF(AT,AC); L1l<-levelOneIF(AC,AT)
  ## level 2: W2 = a_T Theta_C + B_T ; B_T = sum hC[k]*tailT[k]
  hT<-AT$hbar; hC<-AC$hbar; tailT<-rev(cumsum(rev(hT)))-hT; tailC<-rev(cumsum(rev(hC)))-hC
  B_T<-sum(hC*tailT); B_C<-sum(hT*tailC)
  W2<-AT$a*AC$Theta+B_T; L2<-AC$a*AT$Theta+B_C
  Pwin<-L1w$point+W2; Ploss<-L1l$point+L2; WR<-Pwin/Ploss
  ## per-subject IFs (each subject contributes to its own arm; pooled vector length 2n with 1/pi)
  ## D_Pwin: level1(treated via Dt on S^D_T, control via Dc on F^D_C) + W2 + B
  ## treated subjects:
  D_Pwin_T <- as.numeric(crossprod(L1w$Dt, -AT$DFD))            # treated level-1 contribution (D S^D_T)
  D_Pwin_T <- D_Pwin_T + AC$Theta*AT$Da                        # product rule: Theta_C * D_aT
  ## B_T treated side: sum_k hC[k]*D[Theta_T^{>t_k}] -- reuse DTheta-style but moved limit; approx via tail of pimat IF
  ## (approximate D[Theta_T^{>t_k}] by the centered plug-in of tail pi + its xi using moved-limit phi; use cumulative)
  ## For tractability use: D[B_T treated] ~ sum over treated of contribution through h_T tail.
  ##   h_T[m] plug-in = pimat_T[m,]; weight on treated subject = sum_{k<m} hC[k] = cumsum(hC)[m-1].
  wC_lt<-c(0,cumsum(hC)[1:(M-1)])                              # sum_{k<m} hC[k]
  D_Pwin_T <- D_Pwin_T + (colSums(wC_lt*AT$pimat) - sum(wC_lt*hT))   # treated B via weighted pi plug-in IF (centered)
  ## control subjects:
  D_Pwin_C <- as.numeric(crossprod(L1w$Dc, AC$DFD))            # control level-1 (D F^D_C)
  D_Pwin_C <- D_Pwin_C + AT$a*AC$DTheta                        # product rule: a_T * D_Theta_C
  wT_gt<-rev(cumsum(rev(hT)))-hT                               # sum_{m>k} hT[m] (= tailT) weight on control HFH at k
  D_Pwin_C <- D_Pwin_C + (colSums(wT_gt*AC$pimat) - sum(wT_gt*hC))
  ## D_Ploss mirror
  D_Ploss_C <- as.numeric(crossprod(L1l$Dt, -AC$DFD)) + AT$Theta*AC$Da +
               (colSums(c(0,cumsum(hT)[1:(M-1)])*AC$pimat) - sum(c(0,cumsum(hT)[1:(M-1)])*hC))
  D_Ploss_T <- as.numeric(crossprod(L1l$Dc, AT$DFD)) + AC$a*AT$DTheta +
               (colSums((rev(cumsum(rev(hC)))-hC)*AT$pimat) - sum((rev(cumsum(rev(hC)))-hC)*hT))
  ## delta method per arm, pooled 1/pi weights
  Dwr_T <- (1/piT)*( D_Pwin_T/Ploss - Pwin/Ploss^2 * D_Ploss_T )
  Dwr_C <- (1/piC)*( D_Pwin_C/Ploss - Pwin/Ploss^2 * D_Ploss_C )
  Ntot<-AT$n+AC$n
  se <- sqrt( (sum(Dwr_T^2)+sum(Dwr_C^2)) / Ntot^2 )
  c(WR=WR, se=se)
}

if(B==1){
  TR<-{bfT<-simArm(1e6,trt,11)$bf;bfC<-simArm(1e6,ctrl,12)$bf;bfWR(bfT,bfC)}
  cat(sprintf("TRUTH WR=%.3f\n",TR["WR"]))
  r<-oneRep(1,2); cat(sprintf("EST WR=%.3f se=%.4f  95%%CI [%.3f,%.3f] covers:%s\n",
    r["WR"],r["se"],r["WR"]-1.96*r["se"],r["WR"]+1.96*r["se"],(r["WR"]-1.96*r["se"]<=TR["WR"])&(TR["WR"]<=r["WR"]+1.96*r["se"])))
} else {
  TR<-{bfT<-simArm(1e6,trt,11)$bf;bfC<-simArm(1e6,ctrl,12)$bf;bfWR(bfT,bfC)["WR"]}
  cat(sprintf("TRUTH WR=%.3f ; B=%d n=%d/arm grid=%d\n",TR,B,n,gl))
  res<-data.table::rbindlist(lapply(1:B,function(b){r<-tryCatch(oneRep(2*b,2*b+1),error=function(e)NULL)
    if(is.null(r))NULL else data.table(WR=r["WR"],se=r["se"],cover=as.integer(r["WR"]-1.96*r["se"]<=TR&TR<=r["WR"]+1.96*r["se"]))}))
  cat(sprintf("RESULT reps=%d mean(WR)=%.3f bias=%.3f | coverage=%.3f | mean(se)=%.4f emp.sd=%.4f\n",
    nrow(res),mean(res$WR),mean(res$WR)-TR,mean(res$cover),mean(res$se),sd(res$WR)))
}
