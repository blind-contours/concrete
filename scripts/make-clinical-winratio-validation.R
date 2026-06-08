## Clinical (death-priority) win ratio in a Markov illness-death model:
## estimator (multistate curves) vs brute-force pairwise truth on full histories.
## v1: no random censoring (admin at tau).
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
tau <- 1500; grid <- seq(0, tau, length.out = 121); M <- length(grid) - 1L; Mp1 <- M + 1L
ctrl <- list(b01=7e-4, b02=4e-4, b12=1.0e-3, g01=0.3, g02=0.2, g12=0.4)
trt  <- list(b01=5e-4, b02=2.5e-4, b12=6e-4, g01=0.3, g02=0.2, g12=0.4)

simArm <- function(n, p, seed) {
  set.seed(seed); W <- rnorm(n)
  a01<-p$b01*exp(p$g01*W); a02<-p$b02*exp(p$g02*W); a12<-p$b12*exp(p$g12*W)
  T01<-rexp(n,a01); T02<-rexp(n,a02); u0<-pmin(T01,T02); hfh<-T01<T02
  s<-ifelse(hfh,u0,Inf); dpost<-ifelse(hfh,u0+rexp(n,a12),Inf); D<-ifelse(hfh,dpost,u0)
  est <- data.table(W=W,
    d01=as.integer(hfh & u0<tau), t01=ifelse(hfh & u0<tau,u0,Inf),
    d02=as.integer(!hfh & u0<tau), t02=ifelse(!hfh & u0<tau,u0,Inf),
    d12=as.integer(hfh & dpost<tau), t12=ifelse(hfh & dpost<tau,dpost,Inf),
    exit0=pmin(u0,tau), entry1=ifelse(hfh & u0<tau,s,Inf), exit1=ifelse(hfh & u0<tau,pmin(dpost,tau),Inf))
  bf <- data.table(D=ifelse(D<=tau,D,Inf), s=ifelse(hfh & s<=tau,s,Inf),
                   state=ifelse(D<=tau,2L,ifelse(hfh & s<=tau,1L,0L)))
  list(est=est, bf=bf)
}

bfWinRatio <- function(bT, bC) {                       # paired (element-wise) clinical rule
  DT<-bT$D; DC<-bC$D; sT<-bT$s; sC<-bC$s; stT<-bT$state; stC<-bC$state
  win1 <- DC<=tau & DT>DC; loss1 <- DT<=tau & DC>DT; tie1 <- !(win1|loss1)
  rk <- function(st) ifelse(st==0L,3L,2L)              # state0 best; both alive in tie1
  win2 <- tie1 & (rk(stT)>rk(stC) | (rk(stT)==rk(stC) & stT==1L & sT>sC))
  loss2<- tie1 & (rk(stC)>rk(stT) | (rk(stT)==rk(stC) & stT==1L & sC>sT))
  Pw<-mean(win1|win2); Pl<-mean(loss1|loss2); c(Pwin=Pw,Ploss=Pl,WR=Pw/Pl)
}

armMarginal <- function(est) {
  n<-nrow(est); W<-est$W; Cov<-data.frame(W=W); SLlib<-c("SL.mean","SL.glm")
  f01<-fitTransitionSL(rep(0,n),est$exit0,est$d01,Cov,grid,SL.library=SLlib,V=5)
  f02<-fitTransitionSL(rep(0,n),est$exit0,est$d02,Cov,grid,SL.library=SLlib,V=5)
  sub<-est$d01==1
  f12<-fitTransitionSL(est$entry1[sub],est$exit1[sub],est$d12[sub],data.frame(W=W[sub]),grid,SL.library=SLlib,V=5)
  cur<-multistateCurves(predictTransitionSL(f01,Cov),predictTransitionSL(f02,Cov),predictTransitionSL(f12,Cov))
  SDmid <- sqrt(rowMeans(cur$SD)[1:M]*rowMeans(cur$SD)[2:Mp1])   # midpoint overall survival
  list(SD=rowMeans(cur$SD), SDmid=SDmid, a=mean(cur$S0[Mp1,]), Theta=mean(cur$p1[Mp1,]), h=rowMeans(cur$pi))
}
assembleWR <- function(mT, mC) {
  dC<-mC$SD[1:M]-mC$SD[2:Mp1]; dT<-mT$SD[1:M]-mT$SD[2:Mp1]
  winL1 <- sum(mT$SDmid*dC); lossL1 <- sum(mC$SDmid*dT)          # midpoint single-event win
  hT<-mT$h; hC<-mC$h; tailT<-rev(cumsum(rev(hT)))-hT; tailC<-rev(cumsum(rev(hC)))-hC
  B_T<-sum(hC*tailT); B_C<-sum(hT*tailC)
  W2<-mT$a*mC$Theta+B_T; L2<-mC$a*mT$Theta+B_C
  Pw<-winL1+W2; Pl<-lossL1+L2; c(Pwin=Pw,Ploss=Pl,WR=Pw/Pl)
}

cat("=== brute-force truth (Nbf=1.5e6/arm) ===\n")
bfT<-simArm(1.5e6,trt,11)$bf; bfC<-simArm(1.5e6,ctrl,12)$bf
TR<-bfWinRatio(bfT,bfC)
cat(sprintf("  TRUTH: P(win)=%.4f P(loss)=%.4f WR=%.3f\n", TR["Pwin"],TR["Ploss"],TR["WR"]))
cat("=== estimator (n=4000/arm) ===\n")
mT<-armMarginal(simArm(4000,trt,1)$est); mC<-armMarginal(simArm(4000,ctrl,2)$est)
ES<-assembleWR(mT,mC)
cat(sprintf("  EST:   P(win)=%.4f P(loss)=%.4f WR=%.3f\n", ES["Pwin"],ES["Ploss"],ES["WR"]))
cat(sprintf("  rel.error WR = %.1f%%\n", 100*(ES["WR"]/TR["WR"]-1)))
