## One-step estimator of Theta = P(alive at tau, with a prior HFH) = E_W[p1(tau|W)]
## in a Markov illness-death model, validating the path-probability EIF (incl. xi12).
## v1: no random censoring (admin at tau), so G == 1; isolates the multistate EIF.
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
tau <- 1500; grid <- seq(0, tau, length.out = 61); M <- length(grid) - 1L
ends <- grid[-1L]; starts <- grid[-length(grid)]
b01 <- 7e-4; b02 <- 4e-4; b12 <- 1.0e-3
g01 <- 0.3;  g02 <- 0.2;  g12 <- 0.4

simIllnessDeath <- function(n, seed) {
  set.seed(seed); W <- rnorm(n)
  a01 <- b01*exp(g01*W); a02 <- b02*exp(g02*W); a12 <- b12*exp(g12*W)
  T01 <- rexp(n, a01); T02 <- rexp(n, a02)
  u0 <- pmin(T01, T02); hfh <- T01 < T02
  s <- ifelse(hfh, u0, NA)                       # HFH time (state-1 entry)
  dgap <- rexp(n, a12); dpost <- s + dgap        # post-HFH death (calendar)
  # observed multistate quantities, admin censor at tau
  d01 <- as.integer(hfh & u0 < tau)              # 0->1 occurred (HFH) by tau
  d02 <- as.integer(!hfh & u0 < tau)             # 0->2 (death no HFH) by tau
  t01 <- ifelse(d01==1, u0, Inf)
  t02 <- ifelse(d02==1, u0, Inf)
  d12 <- as.integer(d01==1 & dpost < tau)        # 1->2 (death after HFH) by tau
  t12 <- ifelse(d12==1, dpost, Inf)
  exit0 <- pmin(u0, tau)                          # state-0 exit (HFH/death0/admin)
  entry1 <- ifelse(d01==1, s, Inf)               # state-1 entry
  exit1  <- ifelse(d01==1, pmin(dpost, tau), Inf) # state-1 exit (post-HFH death/admin)
  data.table(W=W, d01=d01, t01=t01, d02=d02, t02=t02, d12=d12, t12=t12,
             exit0=exit0, entry1=entry1, exit1=exit1)
}
truthTheta <- function(N = 2e6) {                # MC: P(in state 1 at tau)
  d <- simIllnessDeath(N, 99)
  mean(d$d01==1 & d$exit1 >= tau & d$d12==0)      # entered state1, still there at tau
}

oneStep <- function(d) {
  n <- nrow(d); W <- d$W
  ## fit the three transition hazards (alpha01, alpha02 from state 0; alpha12 left-truncated)
  Cov <- data.frame(W = W)
  f01 <- fitTransitionSL(rep(0,n), d$exit0, d$d01, Cov, grid, SL.library=c("SL.mean","SL.glm"), V=5)
  f02 <- fitTransitionSL(rep(0,n), d$exit0, d$d02, Cov, grid, SL.library=c("SL.mean","SL.glm"), V=5)
  sub <- d$d01==1                                 # post-HFH subjects for alpha12
  f12 <- fitTransitionSL(d$entry1[sub], d$exit1[sub], d$d12[sub], data.frame(W=W[sub]), grid,
                         SL.library=c("SL.mean","SL.glm"), V=5)
  ## per-subject increments at their own W (M x n)
  I01 <- predictTransitionSL(f01, Cov); I02 <- predictTransitionSL(f02, Cov)
  I12 <- predictTransitionSL(f12, Cov)
  ## curves per subject (work column-wise)
  cumExit <- apply(I01 + I02, 2, cumsum)          # Lambda0 at ends (M x n)
  S0 <- rbind(1, exp(-cumExit))                   # S0 at grid 0..M ((M+1) x n)
  cum12 <- apply(I12, 2, cumsum)                  # Lambda12 at ends
  L12full <- rbind(0, cum12)                      # at grid 0..M
  S12_to_tau <- function() exp(-(matrix(L12full[M+1,], M+1, n, byrow=TRUE) - L12full))  # S12(t_k,tau)
  S12tau <- S12_to_tau()                          # (M+1) x n, row k = S12(t_{k-1}... index 0..M)
  ## pi[k] (enter state1 in interval k, survive to tau) = S0(start)*inc01[k]*S12(end,tau)
  S0start <- S0[1:M, , drop=FALSE]                # S0 at interval starts t_{k-1}
  S12end  <- S12tau[2:(M+1), , drop=FALSE]        # S12(t_k, tau)
  pimat <- S0start * I01 * S12end                 # (M x n)
  Pi <- colSums(pimat)                            # Pi(tau|W_i) plug-in
  Pi_gt <- apply(pimat, 2, function(col) rev(cumsum(rev(col))) - col)  # Pi_{>t_k}: sum_{m>k}
  ## clever covariates phi at u=t_k (G==1); use S0,S12 at t_k (index k -> row k+1 of grid-curves)
  S0_k  <- S0[2:(M+1), , drop=FALSE]              # S0(t_k)
  S12_k <- S12tau[2:(M+1), , drop=FALSE]          # S12(t_k,tau)
  phi01 <- (S0_k*S12_k - Pi_gt) / S0_k
  phi02 <- -Pi_gt / S0_k
  phi12 <- -S12_k
  ## martingale increments per subject per interval
  Y0 <- outer(starts, d$exit0, "<")               # in state 0 at interval start (M x n): starts < exit0
  dN01 <- outer(1:M, 1:M)*0                        # placeholder
  inIv <- function(tt) {                           # (M x n) indicator: event time in interval k
    o <- matrix(0L, M, n); k <- findInterval(tt, grid); ok <- is.finite(tt) & k>=1 & k<=M
    o[cbind(k[ok], which(ok))] <- 1L; o
  }
  dN01 <- inIv(d$t01); dN02 <- inIv(d$t02); dN12 <- inIv(d$t12)
  Y1 <- outer(starts, d$entry1, ">=") & outer(starts, d$exit1, "<")  # in state1 at interval start
  dM01 <- dN01 - Y0*I01; dM02 <- dN02 - Y0*I02; dM12 <- dN12 - Y1*I12
  xi01 <- colSums(phi01*dM01); xi02 <- colSums(phi02*dM02); xi12 <- colSums(phi12*dM12)
  theta_plug <- mean(Pi)
  theta_os <- theta_plug + mean(xi01+xi02+xi12)
  eif <- (Pi - theta_os) + xi01 + xi02 + xi12
  list(plug=theta_plug, os=theta_os, se=sqrt(var(eif)/n),
       Pnxi=c(xi01=mean(xi01), xi02=mean(xi02), xi12=mean(xi12)))
}

TR <- truthTheta(); cat(sprintf("TRUTH  Theta = %.4f\n", TR))
r <- oneStep(simIllnessDeath(4000, 1))
cat(sprintf("plug-in = %.4f | one-step = %.4f (se %.4f) | truth = %.4f\n", r$plug, r$os, r$se, TR))
cat(sprintf("  one-step 95%% CI: [%.4f, %.4f]  covers truth: %s\n",
            r$os-1.96*r$se, r$os+1.96*r$se, (r$os-1.96*r$se<=TR)&(TR<=r$os+1.96*r$se)))
cat(sprintf("  Pn xi: 01=%.4f 02=%.4f 12=%.4f\n", r$Pnxi[1], r$Pnxi[2], r$Pnxi[3]))
