## Small-sample behavior of clinicalWinRatio(): the win ratio is a ratio, so it
## is mildly biased / anti-conservative at small n and becomes nominal as n grows.
## NULL DGP (both arms identical -> true WR=1). We sweep n/arm and plot coverage.
## Run in-sample (n.folds=1) for speed: cross-fitting does NOT change this curve
## (the n-dependence is a property of the win-ratio functional, not the nuisances).
suppressWarnings(suppressMessages({ devtools::load_all(".", quiet = TRUE); library(data.table) }))
cli <- commandArgs(trailingOnly = TRUE)
B   <- if (length(cli) >= 1) as.integer(cli[1]) else 160L
tau <- 1500; ng <- 40
nper_grid <- c(400L, 800L, 1600L)

simNull <- function(ntot, seed) { set.seed(seed); W <- rnorm(ntot); arm <- rbinom(ntot, 1, 0.5)
  a01 <- 6e-4 * exp(0.3 * W); a02 <- 4e-4 * exp(0.2 * W); a12 <- 1e-3 * exp(0.4 * W)   # no arm term
  T01 <- rexp(ntot, a01); T02 <- rexp(ntot, a02); u0 <- pmin(T01, T02); hfh <- T01 < T02
  s <- ifelse(hfh, u0, Inf); dpost <- ifelse(hfh, u0 + rexp(ntot, a12), Inf); D <- ifelse(hfh, dpost, u0)
  C <- rexp(ntot, 3e-4); term <- pmin(D, C, tau); died <- as.integer(D < C & D < tau)
  t_hfh <- ifelse(hfh & s < term, s, NA)
  data.table(arm = arm, t_hfh = t_hfh, t_term = term, died = died, W = W) }

one <- function(seed, nper) tryCatch({
  d <- simNull(2 * nper, seed)
  r <- clinicalWinRatio(d, arm = "arm", illness.time = "t_hfh", terminal.time = "t_term",
                        terminal.status = "died", covariates = "W", horizon = tau,
                        n.grid = ng, n.folds = 1L)   # in-sample (see header)
  wr <- r[r$Estimand == "Win Ratio", ]
  data.table(cov = as.integer(wr$`CI Low` <= 1 & 1 <= wr$`CI Hi`),
             rej = as.integer(wr$pValue < 0.05), est = wr$`Pt Est`)
}, error = function(e) NULL)

mc <- min(5L, parallel::detectCores() - 1L)
summ <- rbindlist(lapply(nper_grid, function(nper) {
  res <- rbindlist(Filter(Negate(is.null),
    parallel::mclapply(seq_len(B) + 4242L, one, nper = nper, mc.cores = mc)))
  cat(sprintf("n=%4d/arm  reps=%3d  meanWR=%.3f  coverage(1)=%.3f  type-I=%.3f\n",
              nper, nrow(res), mean(res$est), mean(res$cov), mean(res$rej)))
  data.table(nper = nper, reps = nrow(res), meanWR = mean(res$est),
             cov = mean(res$cov), rej = mean(res$rej),
             cov_se = sqrt(mean(res$cov) * (1 - mean(res$cov)) / nrow(res)))
}))

png("vignettes/figures/clinical-winratio-smalln.png", width = 1500, height = 620, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.4, 1))
with(summ, {
  plot(nper, cov, type = "b", pch = 19, ylim = c(0.88, 0.99), log = "x",
       xlab = "sample size per arm", ylab = "95% CI coverage of WR = 1",
       main = "Coverage -> nominal as n grows", xaxt = "n")
  axis(1, at = nper_grid, labels = nper_grid)
  arrows(nper, cov - 1.96 * cov_se, nper, cov + 1.96 * cov_se, angle = 90, code = 3, length = 0.04, col = "grey50")
  abline(h = 0.95, lty = 2, col = "red")
  plot(nper, meanWR, type = "b", pch = 19, ylim = c(0.97, 1.01), log = "x",
       xlab = "sample size per arm", ylab = "mean win-ratio estimate",
       main = "Ratio bias shrinks at ~1/n", xaxt = "n")
  axis(1, at = nper_grid, labels = nper_grid)
  abline(h = 1, lty = 2, col = "red")
})
par(op); dev.off()
cat("wrote vignettes/figures/clinical-winratio-smalln.png\n")
print(summ)
