## Blocks C and D of the RMT-IF validation (re-run lean: the combined suite ran
## out of memory accumulating brute-force truth vectors + mclapply forks).
## Memory hygiene: small truth Nb, rm()+gc() immediately after each scalar truth.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 2L; tau <- 5; targT <- 4
fl <- function(...) { cat(sprintf(...)); flush(stdout()) }
ts <- seq(0, targT, length.out = 400); dtg <- ts[2] - ts[1]
trapz <- function(z) sum(z) * dtg - 0.5 * dtg * (z[1] + z[length(z)])

## ---------- Block C: getRMTIF first-event K=2 coverage ----------
l1f <- function(W, A) 0.14 * exp(0.3*W - 0.5*A); l2f <- function(W, A) 0.16 * exp(0.2*W - 0.3*A)
truthC <- local({
  set.seed(3); Nb <- 1.5e6
  W <- rnorm(Nb); lt <- l1f(W,1)+l2f(W,1); Tt <- rexp(Nb,lt); Jt <- 1L+rbinom(Nb,1,l2f(W,1)/lt)
  W0 <- rnorm(Nb); lt0 <- l1f(W0,0)+l2f(W0,0); Tc <- rexp(Nb,lt0); Jc <- 1L+rbinom(Nb,1,l2f(W0,0)/lt0)
  fav <- c(`1`=2,`2`=1)
  rk <- function(T,J,t) ifelse(T<=t, fav[as.character(J)], 0L)
  z <- vapply(ts, function(t){rx<-rk(Tt,Jt,t);ry<-rk(Tc,Jc,t);mean(rx<ry)-mean(rx>ry)}, numeric(1))
  trapz(z)
}); gc(FALSE)
simC <- function(seed, n) {
  set.seed(seed); A <- rep(0:1, each=n); W <- rnorm(2*n)
  lt <- l1f(W,A)+l2f(W,A); Te <- rexp(2*n,lt); cause <- 1L+rbinom(2*n,1,l2f(W,A)/lt); C <- rexp(2*n,0.06)
  dat <- data.table(id=1:(2*n), time=pmin(Te,C,tau), status=ifelse(Te<=pmin(C,tau),cause,0L), trt=A, W=W, W2=rnorm(2*n))
  a <- suppressMessages(formatArguments(DataTable=dat, EventTime="time", EventType="status", Treatment="trt",
       ID="id", Intervention=0:1, TargetTime=c(1,2,3,4), TargetEvent=c(1,2), CVArg=list(V=2),
       MaxUpdateIter=8, Verbose=FALSE,
       Model=list(trt="SL.mean","0"=list(Cox=survival::Surv(time,status==0)~.),
         "1"=list(Cox=survival::Surv(time,status==1)~.),"2"=list(Cox=survival::Surv(time,status==2)~.))))
  o <- suppressMessages(getRMTIF(suppressMessages(suppressWarnings(doConcrete(a))),
        Horizon=4, Intervention=c(1,2), TargetEvent=c(1,2)))
  r <- o[Estimand=="RMT-IF"]; c(est=r[["Pt Est"]], se=r[["se"]],
    cov=as.integer(r[["CI Low"]]<=truthC & truthC<=r[["CI Hi"]]))
}
BC <- 120L; nC <- 600L
RC <- do.call(rbind, parallel::mclapply(seq_len(BC), function(b)
  tryCatch({on.exit(gc(FALSE)); simC(700+b, nC)}, error=function(e) rep(NA,3)), mc.cores=MC))
RC <- RC[stats::complete.cases(RC), , drop=FALSE]
fl("\n=== BLOCK C: getRMTIF first-event K=2, %d reps, n=%d/arm (truth %.4f) ===\n", nrow(RC), nC, truthC)
fl("  mean est %.4f (bias %+.4f) | empSD %.4f | meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
   mean(RC[,"est"]), mean(RC[,"est"])-truthC, sd(RC[,"est"]), mean(RC[,"se"]),
   sd(RC[,"est"])/mean(RC[,"se"]), mean(RC[,"cov"]))
fl("BLOCK-C-DONE\n"); rm(RC); gc(FALSE)

## ---------- Block D: clinicalRMTIF K=2 cross-fitted (n.folds=5) ----------
bD2 <- function(W,A) 0.10*exp(0.3*W-0.5*A); bH2 <- function(W,A) 0.20*exp(0.2*W-0.3*A)
truthD <- local({
  set.seed(2); Nb <- 1.5e6
  Wx<-rnorm(Nb); tDx<-rexp(Nb,bD2(Wx,1)); tHx<-rexp(Nb,bH2(Wx,1))
  Wy<-rnorm(Nb); tDy<-rexp(Nb,bD2(Wy,0)); tHy<-rexp(Nb,bH2(Wy,0))
  rk <- function(tD,tH,t) ifelse(tD<=t, 2L, ifelse(tH<=t & tH<tD, 1L, 0L))
  z <- vapply(ts, function(t){rx<-rk(tDx,tHx,t);ry<-rk(tDy,tHy,t);mean(rx<ry)-mean(rx>ry)}, numeric(1))
  trapz(z)
}); gc(FALSE)
BD <- 50L; nD <- 700L
RD <- do.call(rbind, parallel::mclapply(seq_len(BD), function(b) tryCatch({on.exit(gc(FALSE))
  set.seed(900+b); A <- rep(0:1, each=nD); W <- rnorm(2*nD)
  tD <- rexp(2*nD,bD2(W,A)); tH <- rexp(2*nD,bH2(W,A)); C <- rexp(2*nD,0.05); obs <- pmin(tD,C,tau)
  dat <- data.frame(arm=A, t_hosp=ifelse(tH<obs,tH,NA), t_term=obs, died=as.integer(tD<=pmin(C,tau)), W=W, W2=rnorm(2*nD))
  o <- suppressMessages(suppressWarnings(clinicalRMTIF(dat, arm="arm", illness.time="t_hosp",
       terminal.time="t_term", terminal.status="died", covariates=c("W","W2"),
       horizon=targT, n.grid=40, n.folds=5)))
  r <- o[Estimand=="RMT-IF"]; c(est=r[["Pt Est"]], se=r[["se"]],
    cov=as.integer(r[["CI Low"]]<=truthD & truthD<=r[["CI Hi"]]))},
  error=function(e) rep(NA,3)), mc.cores=MC))
RD <- RD[stats::complete.cases(RD), , drop=FALSE]
fl("\n=== BLOCK D: clinicalRMTIF K=2 CROSS-FITTED (n.folds=5), %d reps, n=%d/arm (truth %.4f) ===\n", nrow(RD), nD, truthD)
fl("  bias %+.4f | SD/SE %.3f | cover %.3f\n", mean(RD[,"est"])-truthD,
   sd(RD[,"est"])/mean(RD[,"se"]), mean(RD[,"cov"]))
fl("CD-DONE\n")
