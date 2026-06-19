## PRO bottom-layer: ANALYTIC EIF validation (step 1 of wiring PRO tiers in).
##
## The tier win component is a bilinear form in the un-normalized reach-weighted
## marker measures  H_a(y) = E[ rho_a(W) Q_a(y | W) ]:
##
##     W^(k) = int [ E rho_T - H_T(y + delta) ] dH_C(y),     (T wins at the PRO tier)
##     L^(k) = int [ E rho_T - H_T(y - delta) ] dH_C(y),     mirror (T loses) -- via swap
##
## with rho_a(W) = P(reach | W, a) and Q_a(y|W) = P(marker <= y | W, reach, a).
## Each H_a(y) has a three-term influence function that collapses cleanly when the
## reach indicator is observed (no censoring before the landmark):
##
##   ifH_a,i(c) = Q(c|W_i) * alive_i  -  H_a(c)
##              + rho(W_i) * (alive_i obs_i / pi(W_i)) * ( 1{Y_i<=c} - Q(c|W_i) )
##   ifTot_a,i  = ifH_a,i(Inf) = alive_i - E rho_a .
##
## (main + reach-residual telescope to Q(c|W_i)*alive_i; the third term is the MAR
## IPCW correction for estimating the marker CDF among observed reachers.) In the
## multistate ENGINE the reach is event-free-alive at the horizon and may be
## censored, so the first two terms are replaced by the occupancy adjoint IF
## (survKilled/survIFset, already validated for RMT-IF) carrying terminal reward
## Q(c|W); the third term is unchanged. This script validates the H-functional EIF
## and the bilinear win-component assembly in the directly-observed-reach case.
##
## Reports: (A) one-step W^(k) bias + analytic-SE coverage vs brute-force truth
## with a fast glm threshold-CDF; (B) SL-vs-glm-vs-oracle bias on the marker CDF.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ library(data.table) }))
set.seed(1)
tauL <- 1.0; delta <- 5; B <- 80L; n <- 1500L; NCUT <- 120L
CGRID <- seq(15, 100, length.out = NCUT)        # FIXED marker grid (not data-dependent)
muD <- function(W,A) 0.55*exp(0.35*W - 0.7*A)
muY <- function(W,A) 58 + 6*W + 7*A
sdY <- 14
piObs <- function(W) plogis(0.7 + 0.5*W)

## ---- brute-force truth (complete data) ----
truth <- local({
  Ng <- 4e5L
  draw <- function(A){ W<-rnorm(Ng); tD<-rexp(Ng, muD(W,A)); Y<-muY(W,A)+rnorm(Ng,0,sdY)
    list(alive = tD>tauL, Y=Y, pR=mean(tD>tauL)) }
  X<-draw(1); Y<-draw(0)
  rk <- X$pR * Y$pR
  ax<-which(X$alive); ay<-which(Y$alive); M<-3e6L
  i<-sample(ax,M,TRUE); j<-sample(ay,M,TRUE)
  w<-mean(X$Y[i] > Y$Y[j] + delta); l<-mean(X$Y[i] < Y$Y[j] - delta)
  c(Wk=rk*w, Lk=rk*l, rk=rk, w=w, l=l)
}); invisible(gc(FALSE))

## ---- per-arm reach-weighted measure H_a(cgrid) + per-subject IF rows ----
## cdfFit: "glm" | "sl" | "oracle".  Returns Q (n x G), H (G), Erho, and the
## per-subject ifH (n x G) and ifTot (n) influence rows for this arm.
armMeasure <- function(d, a, cgrid, cdfFit, SL.library = NULL) {
  G <- length(cgrid); sub <- d[A==a]; nA <- nrow(sub)
  ## reach rho(W) = P(alive at landmark | W): per-arm logistic (oracle: true survival)
  rho <- if (cdfFit=="oracle") exp(-muD(sub$W,a)*tauL) else
    as.numeric(predict(glm(alive~W+W2, binomial, sub), sub, type="response"))
  Erho <- mean(rho)
  ## landmark-observation IPCW among the alive (MAR)
  oa <- sub$alive==1
  pio <- if (cdfFit=="oracle") piObs(sub$W) else
    as.numeric(predict(glm(obs~W+W2, binomial, sub[oa]), sub, type="response"))
  ipcw <- ifelse(oa & sub$obs==1, 1/pmax(pio,0.05), 0)
  ## conditional CDF Q(c|W) over cutpoints
  if (cdfFit=="oracle") {
    Q <- sapply(cgrid, function(c) pnorm((c - muY(sub$W,a))/sdY))
  } else {
    obA <- which(oa & sub$obs==1); wsub <- ipcw[obA]
    Q <- sapply(cgrid, function(c){
      yb <- as.integer(sub$Y[obA] <= c)
      if (length(unique(yb))<2) return(rep(mean(yb), nA))
      if (cdfFit=="sl") {
        Xtr <- as.data.frame(sub[obA, .(W,W2)]); Xnew <- as.data.frame(sub[, .(W,W2)])
        pr <- tryCatch(suppressWarnings({
          f <- SuperLearner::SuperLearner(Y=yb, X=Xtr, newX=Xnew, family=binomial(),
                 SL.library=SL.library, obsWeights=wsub, cvControl=list(V=3L))
          as.numeric(f$SL.predict) }), error=function(e) NULL)
        if (is.null(pr) || any(!is.finite(pr)))                       # fall back to glm at this cutpoint
          pr <- as.numeric(predict(suppressWarnings(glm(yb~W+W2, binomial, data=sub[obA], weights=wsub)),
                                   sub, type="response"))
        pr
      } else {
        f <- suppressWarnings(glm(yb~W+W2, binomial, data=sub[obA], weights=wsub))
        as.numeric(predict(f, sub, type="response"))
      }
    })
  }
  Q <- t(apply(Q,1,cummax)); Q <- pmin(pmax(Q,0),1)        # monotone, in [0,1]
  H <- as.numeric(colSums(rho*Q)/nA)                        # H_a(c) = mean rho*Q
  ## per-subject IF rows of H_a(c):  Q(c)*alive - H(c) + rho*ipcw*(1{Y<=c}-Q(c))
  aliveV <- as.numeric(sub$alive)
  Yle <- if (cdfFit=="oracle") matrix(0,nA,G) else outer(ifelse(is.na(sub$Y),Inf,sub$Y), cgrid, `<=`)*1
  ifH <- Q*aliveV - matrix(H,nA,G,byrow=TRUE) + (rho*ipcw)*(Yle - Q)
  ifTot <- aliveV - Erho
  list(Q=Q, H=H, Erho=Erho, ifH=ifH, ifTot=ifTot, nA=nA)
}

## ---- one-step tier win component W^(k) (T wins) + analytic SE ----
## Midpoint rule for  W^(k) = int [E rho_T - H_T(y+delta)] dH_C(y): evaluate the
## (smooth, decreasing) coefficient at interval midpoints against the dH_C mass in
## each interval (second-order accurate; a right-Riemann sum biased it ~+0.006).
tierWin <- function(MT, MC, cgrid, delta) {
  G <- length(cgrid)
  cmid <- (cgrid[-1] + cgrid[-G]) / 2                        # interval midpoints (G-1)
  dHc  <- diff(MC$H)                                         # dH_C mass per interval (G-1)
  HTmid <- approx(cgrid, MT$H, xout=cmid+delta, rule=2)$y    # H_T at midpoint + delta
  coefT <- MT$Erho - HTmid                                   # [E rho_T - H_T(mid+delta)]
  Wk <- sum(coefT * dHc)
  ## IF over arm T:  sum_int dH_C * ( ifTot_T - ifH_T(mid+delta) )
  ifHT_mid <- t(apply(MT$ifH, 1, function(r) approx(cgrid, r, xout=cmid+delta, rule=2)$y))
  IFwinT <- as.numeric((MT$ifTot %o% rep(1,G-1) - ifHT_mid) %*% dHc)
  ## IF over arm C:  sum_int coefT * ( ifH_C(c) - ifH_C(c_prev) )
  difC <- MC$ifH[,-1,drop=FALSE] - MC$ifH[,-G,drop=FALSE]
  IFlosC <- as.numeric(difC %*% coefT)
  Wk_os <- Wk + mean(IFwinT) + mean(IFlosC)
  seWk <- sqrt(var(IFwinT)/MT$nA + var(IFlosC)/MC$nA)
  c(Wk=Wk_os, se=seWk)
}

simOne <- function(seed, cdfFit="glm", SL.library=NULL) {
  set.seed(seed); N<-2L*n; A<-rep(0:1,each=n); W<-rnorm(N); W2<-rnorm(N)
  tD<-rexp(N, muD(W,A)); alive<-tD>tauL
  Y<-ifelse(alive, muY(W,A)+rnorm(N,0,sdY), NA_real_)
  obs<-alive & (runif(N) < piObs(W))
  d<-data.table(A=A, W=W, W2=W2, alive=as.integer(alive), Y=Y, obs=as.integer(obs))
  cgrid <- CGRID                                  # fixed estimand grid
  MT<-armMeasure(d,1,cgrid,cdfFit,SL.library); MC<-armMeasure(d,0,cgrid,cdfFit,SL.library)
  win<-tierWin(MT,MC,cgrid,delta)
  c(Wk=win["Wk"], se=win["se"],
    cov=as.integer(win["Wk"]-1.96*win["se"] <= truth["Wk"] & truth["Wk"] <= win["Wk"]+1.96*win["se"]))
}

## ---- (A) analytic-EIF coverage of W^(k): oracle nuisances vs glm ----
covBlock <- function(fit, lab) {
  R <- do.call(rbind, lapply(seq_len(B), function(b)
    tryCatch(simOne(200+b, fit), error=function(e) rep(NA,3))))
  R <- R[stats::complete.cases(R),,drop=FALSE]
  cat(sprintf("  %-11s est %.4f bias %+.4f | empSD %.4f meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
    lab, mean(R[,1]), mean(R[,1])-truth["Wk"], sd(R[,1]), mean(R[,2]),
    sd(R[,1])/mean(R[,2]), mean(R[,3])))
}
cat(sprintf("\n===== (A) PRO tier W^(k) analytic-EIF coverage (%d reps, n=%d/arm) =====\n", B, n))
cat(sprintf("  truth W^(k) %.4f (reach %.4f, w %.4f)\n", truth["Wk"], truth["rk"], truth["w"]))
covBlock("oracle", "oracle")    # isolates the bilinear win-IF assembly (no nuisance estimation)
covBlock("glm",    "glm")       # standalone: glm reach + glm CDF + glm IPCW

## ---- (B) marker-CDF spec: SL vs glm vs oracle (bias on W^(k)), few reps ----
SLlib <- c("SL.mean","SL.glm","SL.glm.interaction")
have_sl <- requireNamespace("SuperLearner", quietly=TRUE)
nb <- 12L
bias <- function(fit, lib=NULL) {
  v <- sapply(seq_len(nb), function(b) tryCatch(simOne(500+b, fit, lib)["Wk.Wk"], error=function(e) NA))
  mean(v, na.rm=TRUE) - truth["Wk"] }
cat(sprintf("\n===== (B) marker-CDF spec, W^(k) bias (%d reps) =====\n", nb))
cat(sprintf("  oracle CDF : bias %+.4f\n", bias("oracle")))
cat(sprintf("  glm    CDF : bias %+.4f\n", bias("glm")))
if (have_sl) cat(sprintf("  SL     CDF : bias %+.4f  (library: %s)\n", bias("sl", SLlib), paste(SLlib,collapse=", ")))
cat("PRO-EIF-DONE\n")
