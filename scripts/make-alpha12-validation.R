suppressWarnings(suppressMessages({
  devtools::load_all(".", quiet = TRUE); library(data.table); library(survival)
}))
tau <- 1500; grid <- seq(0, tau, length.out = 60)     # regular grid, 59 intervals
SLlib <- c("SL.mean", "SL.glm")                        # robust small library for the test

# delayed-entry illness-death post-HFH sub-sample: entry = HFH time s, then a
# death gap (hazard a*exp(g*W)) and a dropout gap; left-truncated at s.
simPostHFH <- function(n, a, g, cc, seed) {
  set.seed(seed); W <- rnorm(n); s <- runif(n, 0, 0.5 * tau)
  gapD <- rexp(n, a * exp(g * W)); gapC <- rexp(n, cc)
  gap <- pmin(gapD, gapC); exitcal <- s + gap
  admin <- exitcal > tau
  exit <- ifelse(admin, tau, exitcal)
  event <- ifelse(admin, 0L, as.integer(gapD < gapC))
  data.table(entry = s, exit = exit, event = event, W = W)
}

cat("===== (A) constant hazard, no covariate =====\n")
a <- 1.2e-3
d <- simPostHFH(4000, a = a, g = 0, cc = 4e-4, seed = 1)
cat(sprintf("  n=%d, post-HFH deaths=%d, censored=%d\n", nrow(d), sum(d$event), sum(d$event==0)))
fit <- fitTransitionSL(d$entry, d$exit, d$event, d[, .(W)], grid, SL.library = SLlib, V = 5)
inc <- predictTransitionSL(fit, data.frame(W = 0))     # W has no effect here
for (s in c(200, 500)) for (t in c(900, 1400)) {
  shat <- s12FromIncrements(inc, grid, s, t)[1]
  strue <- exp(-a * (t - s))
  cat(sprintf("  S12(%4d,%4d): SL=%.4f  truth=%.4f  (Lambda SL=%.3f truth=%.3f)\n",
              s, t, shat, strue, -log(shat), a * (t - s)))
}

cat("\n===== (B) covariate-dependent hazard a*exp(g*W) =====\n")
g <- 0.5
d <- simPostHFH(5000, a = a, g = g, cc = 4e-4, seed = 2)
fit <- fitTransitionSL(d$entry, d$exit, d$event, d[, .(W)], grid, SL.library = SLlib, V = 5)
for (w in c(-1, 0, 1)) {
  inc <- predictTransitionSL(fit, data.frame(W = w))
  shat <- s12FromIncrements(inc, grid, 300, 1200)[1]
  strue <- exp(-a * exp(g * w) * (1200 - 300))
  cat(sprintf("  W=%+d: S12(300,1200) SL=%.4f  truth=%.4f\n", w, shat, strue))
}

cat("\n===== (C) left-truncation: SL vs delayed-entry Nelson-Aalen =====\n")
d <- simPostHFH(4000, a = a, g = 0, cc = 4e-4, seed = 3)
fit <- fitTransitionSL(d$entry, d$exit, d$event, d[, .(W)], grid, SL.library = SLlib, V = 5)
inc <- predictTransitionSL(fit, data.frame(W = 0))
# delayed-entry NA via survival::survfit on (entry, exit, event)
naf <- survfit(Surv(entry, exit, event) ~ 1, data = d)
na_cumhaz <- function(t) { i <- findInterval(t, naf$time); if (i < 1) 0 else cumsum(naf$n.event / naf$n.risk)[i] }
for (s in c(200, 500)) for (t in c(900, 1300)) {
  L_sl <- -log(s12FromIncrements(inc, grid, s, t)[1])
  L_na <- na_cumhaz(t) - na_cumhaz(s)
  cat(sprintf("  Lambda12(%4d,%4d): SL=%.3f  delayed-entry NA=%.3f  truth=%.3f\n",
              s, t, L_sl, L_na, a * (t - s)))
}
# control: ignoring left truncation (entry=0) should bias the early hazard
fit0 <- fitTransitionSL(rep(0, nrow(d)), d$exit, d$event, d[, .(W)], grid, SL.library = SLlib, V = 5)
inc0 <- predictTransitionSL(fit0, data.frame(W = 0))
cat(sprintf("  [naive entry=0] Lambda12(200,900)=%.3f  vs truth %.3f (expected biased)\n",
            -log(s12FromIncrements(inc0, grid, 200, 900)[1]), a * 700))
cat("\nDone.\n")
