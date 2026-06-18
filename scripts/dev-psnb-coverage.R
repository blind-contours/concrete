## Coverage validation of clinicalPSNB's PSNB and PSWR influence-function
## inference, against a brute-force charter-weighted pairwise truth (complete
## data, the estimand under CAR). K=2 (death > hospitalization), charter (0.7,0.3).
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
data.table::setDTthreads(1L)
MC <- 1L; B <- 60L; n <- 700L; tau <- 5; targT <- 4
alpha <- c(0.7, 0.3)                                   # charter: death 0.7, hosp 0.3
bD <- function(W,A) 0.08*exp(0.3*W - 0.6*A); bH <- function(W,A) 0.25*exp(0.2*W - 0.2*A)

## ---- brute-force truth: pairwise stage-conditional effects on complete data ----
truth <- local({
  set.seed(1); Ng <- 8000L; M <- 2e6L
  drw <- function(a){ W<-rnorm(Ng); list(tD=rexp(Ng,bD(W,a)), tH=rexp(Ng,bH(W,a))) }
  X <- drw(1); Y <- drw(0)
  i <- sample.int(Ng, M, TRUE); j <- sample.int(Ng, M, TRUE)
  tDx<-X$tD[i]; tHx<-X$tH[i]; tDy<-Y$tD[j]; tHy<-Y$tH[j]
  Xa <- tDx>targT; Ya <- tDy>targT                    # alive at horizon
  c1 <- ifelse(Xa & !Ya, 1L, ifelse(!Xa & Ya, -1L, ifelse(!Xa & !Ya, sign(tDx-tDy), 0L)))
  R2 <- (c1==0L)                                       # tied on death -> reach hosp tier
  Xh <- tHx<=targT; Yh <- tHy<=targT
  c2 <- ifelse(!Xh & Yh, 1L, ifelse(Xh & !Yh, -1L, ifelse(Xh & Yh, sign(tHx-tHy), 0L)))
  w1<-mean(c1==1); l1<-mean(c1==-1); D1<-w1-l1
  r2<-mean(R2); w2<-mean(R2 & c2==1)/r2; l2<-mean(R2 & c2==-1)/r2; D2<-w2-l2
  psnb <- alpha[1]*D1 + alpha[2]*D2
  pswr <- (alpha[1]*w1 + alpha[2]*w2)/(alpha[1]*l1 + alpha[2]*l2)
  c(PSNB=psnb, PSWR=pswr, D1=D1, D2=D2, r2=r2)
}); gc(FALSE)

simOne <- function(seed) {
  set.seed(seed); A<-rep(0:1,each=n); W<-rnorm(2*n)
  tD<-rexp(2*n,bD(W,A)); tH<-rexp(2*n,bH(W,A)); C<-rexp(2*n,0.05); obs<-pmin(tD,C,tau)
  dat<-data.frame(arm=A, t_hosp=ifelse(tH<obs,tH,NA), t_term=obs,
                  died=as.integer(tD<=pmin(C,tau)), W=W, W2=rnorm(2*n))
  o<-suppressMessages(suppressWarnings(clinicalPSNB(dat, arm="arm", illness.time="t_hosp",
       terminal.time="t_term", terminal.status="died", covariates=c("W","W2"),
       charter=alpha, horizon=targT, n.grid=40, n.folds=1)))
  pn<-o[Estimand=="PSNB"]; wr<-o[Estimand=="PSWR"]
  c(pn=pn[["Pt Est"]], pnse=pn[["se"]], pncov=as.integer(pn[["CI Low"]]<=truth["PSNB"] & truth["PSNB"]<=pn[["CI Hi"]]),
    wr=wr[["Pt Est"]], wrse=wr[["se"]], wrcov=as.integer(wr[["CI Low"]]<=truth["PSWR"] & truth["PSWR"]<=wr[["CI Hi"]]))
}
R<-do.call(rbind, parallel::mclapply(seq_len(B), function(b)
  tryCatch({on.exit(gc(FALSE)); simOne(300+b)}, error=function(e) rep(NA,6)), mc.cores=MC))
R<-R[stats::complete.cases(R),,drop=FALSE]
cat(sprintf("\n===== clinicalPSNB coverage (%d reps, n=%d/arm, charter 0.7/0.3) =====\n", nrow(R), n))
cat(sprintf("  truth: PSNB %.4f  PSWR %.4f  (Delta_death %.3f, Delta_hosp %.3f, reach_hosp %.3f)\n",
    truth["PSNB"], truth["PSWR"], truth["D1"], truth["D2"], truth["r2"]))
cat(sprintf("  PSNB: mean %.4f bias %+.4f | empSD %.4f meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
    mean(R[,"pn"]), mean(R[,"pn"])-truth["PSNB"], sd(R[,"pn"]), mean(R[,"pnse"]),
    sd(R[,"pn"])/mean(R[,"pnse"]), mean(R[,"pncov"])))
cat(sprintf("  PSWR: mean %.4f bias %+.4f | empSD %.4f meanSE %.4f (SD/SE %.3f) | cover %.3f\n",
    mean(R[,"wr"]), mean(R[,"wr"])-truth["PSWR"], sd(R[,"wr"]), mean(R[,"wrse"]),
    sd(R[,"wr"])/mean(R[,"wrse"]), mean(R[,"wrcov"])))
cat("PSNB-COV-DONE\n")
