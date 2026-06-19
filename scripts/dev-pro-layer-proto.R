## Standalone prototype: a continuous PRO (KCCQ) bottom layer for the hierarchical
## win statistics, with REACH-WEIGHTED distribution standardization.
##
## A PRO layer is reached only by pairs tied on the higher (death) layer, i.e.
## both alive at the landmark. The marker is undefined off that set, so we do NOT
## use the plain marginal plug-in n^{-1} sum_i Qhat(y | W_i, A=a). We standardize
## the conditional CDF by the arm-specific reach probability:
##
##   G_a^R(y) = E[ p_{R,a}(W) Q_a(y | W, R=1) ] / E[ p_{R,a}(W) ],
##
## where p_{R,a}(W) = P(alive at landmark | W, A=a) and Q_a is the conditional CDF
## of the PRO among reachable subjects (estimated by IPCW-weighted binary-threshold
## regression to correct landmark-observation missingness -- separate from death).
## Then  w_k = int [1 - G_1^R(y+delta)] dG_0^R(y),  W^(k) = (pair reach) * w_k.
##
## Validated against a brute-force pairwise truth on complete data (both alive,
## margin delta). Models are glm here for speed; the binary-threshold SuperLearner
## is a drop-in (same 1{Y<=c} ~ W, binomial API). Sequential (MC=1).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages(library(data.table)))
set.seed(1)
tauL <- 1.0; delta <- 5; B <- 30L; n <- 2000L
muD <- function(W,A) 0.55*exp(0.35*W - 0.7*A)            # death hazard (treatment protects)
muY <- function(W,A) 58 + 6*W + 7*A                       # KCCQ mean (treatment raises; higher=better)
sdY <- 14
piObs <- function(W) plogis(0.7 + 0.5*W)                  # MAR observation prob among the alive

## ---- brute-force truth (complete data, no censoring/missingness) ----
truthWL <- local({
  Ng <- 3e5L
  draw <- function(A){ W<-rnorm(Ng); tD<-rexp(Ng, muD(W,A)); Y<-muY(W,A)+rnorm(Ng,0,sdY)
    list(alive = tD>tauL, Y=Y, pR=mean(tD>tauL)) }
  X<-draw(1); Y<-draw(0)
  rk <- X$pR * Y$pR                                       # pair reach = P(both alive)
  ax <- which(X$alive); ay <- which(Y$alive); M<-2e6L
  i<-sample(ax,M,TRUE); j<-sample(ay,M,TRUE)
  w <- mean(X$Y[i] > Y$Y[j] + delta); l <- mean(X$Y[i] < Y$Y[j] - delta)
  c(W=rk*w, L=rk*l, rk=rk, w=w, l=l)
}); invisible(gc(FALSE))

## ---- estimator for one simulated trial ----
simOne <- function(seed){
  set.seed(seed); N<-2L*n; A<-rep(0:1,each=n); W<-rnorm(N); W2<-rnorm(N)
  tD<-rexp(N, muD(W,A)); alive<-tD>tauL                  # reach (assume follow-up to landmark)
  Y<-ifelse(alive, muY(W,A)+rnorm(N,0,sdY), NA_real_)
  obs<-alive & (runif(N) < piObs(W))                     # MAR landmark observation among the alive
  dat<-data.table(A=A, W=W, W2=W2, alive=as.integer(alive), Y=Y, obs=obs)

  ## p_{R,a}(W): P(alive at landmark | W, A=a), per-arm binary glm, predict for all W
  pR <- function(a){ d<-dat[A==a]; f<-glm(alive~W+W2, binomial, d)
    as.numeric(predict(f, dat, type="response")) }
  pR1<-pR(1); pR0<-pR(0)
  rk_hat <- mean(pR1)*mean(pR0)

  ## landmark-observation IPCW among the alive (MAR): 1 / P(obs | W, alive)
  do<-dat[alive==1]; fo<-glm(obs~W+W2, binomial, do)
  ipcw_all <- 1/pmax(predict(fo, dat, type="response"), 0.05)

  ## Q_a(y | W, R=1): IPCW-weighted binary-threshold regressions among observed-alive
  cgrid <- seq(quantile(dat$Y,0.02,na.rm=TRUE), quantile(dat$Y,0.98,na.rm=TRUE), length.out=24)
  Qa <- function(a){
    sub <- dat[A==a & obs==TRUE]; wsub <- ipcw_all[dat$A==a & dat$obs==TRUE]
    Qmat <- sapply(cgrid, function(c){
      yb<-as.integer(sub$Y<=c)
      if(length(unique(yb))<2) return(rep(mean(yb), nrow(dat)))
      f<-suppressWarnings(glm(yb~W+W2, binomial, data=sub, weights=wsub))
      as.numeric(predict(f, dat, type="response")) })
    t(apply(Qmat,1,cummax))                              # enforce monotone CDF across cutpoints
  }
  Q1<-Qa(1); Q0<-Qa(0)
  ## reach-weighted standardized CDFs G_a^R(cgrid)
  G1 <- colSums(pR1*Q1)/sum(pR1); G0 <- colSums(pR0*Q0)/sum(pR0)
  G1<-cummax(pmin(pmax(G1,0),1)); G0<-cummax(pmin(pmax(G0,0),1))
  Gfun<-function(G,y) approx(cgrid,G,xout=y,rule=2)$y
  dG0<-diff(c(0,G0))                                      # step increments of G0
  w_hat <- sum((1-Gfun(G1,cgrid+delta))*dG0)             # P(Y1 > Y0 + delta)
  l_hat <- sum(Gfun(G1,cgrid-delta)*dG0)                 # P(Y1 < Y0 - delta)
  ## NAIVE: plain marginal plug-in (ignores reach-weighting) -> biased
  G1n<-cummax(pmin(pmax(colMeans(Q1),0),1)); G0n<-cummax(pmin(pmax(colMeans(Q0),0),1))
  dG0n<-diff(c(0,G0n))
  w_naive <- sum((1-Gfun(G1n,cgrid+delta))*dG0n)
  c(W=rk_hat*w_hat, L=rk_hat*l_hat, rk=rk_hat, w=w_hat, l=l_hat, w_naive=w_naive)
}

R<-do.call(rbind, lapply(seq_len(B), function(b) tryCatch(simOne(100+b), error=function(e) rep(NA,6))))
R<-R[stats::complete.cases(R),,drop=FALSE]
cat(sprintf("\n===== PRO layer (reach-weighted) prototype: %d reps, n=%d/arm, delta=%g =====\n", nrow(R), n, delta))
cat(sprintf("  pair reach r_k:  truth %.4f  est %.4f (bias %+.4f)\n", truthWL["rk"], mean(R[,"rk"]), mean(R[,"rk"])-truthWL["rk"]))
cat(sprintf("  stage win w_k:   truth %.4f  est %.4f (bias %+.4f)\n", truthWL["w"],  mean(R[,"w"]),  mean(R[,"w"]) -truthWL["w"]))
cat(sprintf("  stage loss l_k:  truth %.4f  est %.4f (bias %+.4f)\n", truthWL["l"],  mean(R[,"l"]),  mean(R[,"l"]) -truthWL["l"]))
cat(sprintf("  tier win  W^(k): truth %.4f  est %.4f (bias %+.4f)\n", truthWL["W"],  mean(R[,"W"]),  mean(R[,"W"]) -truthWL["W"]))
cat(sprintf("  tier loss L^(k): truth %.4f  est %.4f (bias %+.4f)\n", truthWL["L"],  mean(R[,"L"]),  mean(R[,"L"]) -truthWL["L"]))
## contrast: naive (non-reach-weighted) standardization is biased for w_k
cat(sprintf("\n  reach-weighted w_k:  est %.4f (bias %+.4f)\n", mean(R[,"w"]), mean(R[,"w"])-truthWL["w"]))
cat(sprintf("  NAIVE marginal w_k:  est %.4f (bias %+.4f)  <- ignores reach-weighting\n",
            mean(R[,"w_naive"]), mean(R[,"w_naive"])-truthWL["w"]))
## ---- ORACLE check: true reach + true CDF (no estimation) isolates whether the
## residual bias is the standardization/integral (structural) or learner misspec ----
oracle <- local({
  Ng<-3e5L; W<-rnorm(Ng)
  pR1<-exp(-muD(W,1)*tauL); pR0<-exp(-muD(W,0)*tauL)
  cg<-seq(quantile(muY(W,1)+rnorm(Ng,0,sdY),0.01), quantile(muY(W,0)+rnorm(Ng,0,sdY),0.99), length.out=200)
  Q1<-sapply(cg, function(c) pnorm((c-muY(W,1))/sdY)); Q0<-sapply(cg, function(c) pnorm((c-muY(W,0))/sdY))
  G1<-colSums(pR1*Q1)/sum(pR1); G0<-colSums(pR0*Q0)/sum(pR0)
  Gf<-function(G,y) approx(cg,G,xout=y,rule=2)$y; dG0<-diff(c(0,G0))
  c(w=sum((1-Gf(G1,cg+delta))*dG0), l=sum(Gf(G1,cg-delta)*dG0))
})
cat(sprintf("\n  [oracle, fine grid] w_k %.4f (truth %.4f, gap %+.4f) | l_k %.4f (truth %.4f, gap %+.4f)\n",
            oracle["w"], truthWL["w"], oracle["w"]-truthWL["w"], oracle["l"], truthWL["l"], oracle["l"]-truthWL["l"]))
cat("PRO-PROTO-DONE\n")
