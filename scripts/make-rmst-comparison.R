#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make-rmst-comparison.R
#
# Simulation comparing two ways of estimating the restricted mean survival time
# from a `concrete` fit:
#   - pointwise: getRMST()   -- integrate pointwise-targeted absolute risks
#   - direct:    targetRMST() -- fluctuate the hazards for the RMST estimating
#                                equation (integrated clever covariate)
#
# The data-generating mechanism uses exponential cause-specific hazards, so the
# true marginal RMST and life-years-lost have closed forms. The fit uses a
# deliberately sparse two-time target grid, which is where the direct method's
# full-grid integration and direct targeting help most.
#
#   Rscript scripts/make-rmst-comparison.R [B] [n]
#
# Writes vignettes/figures/rmst-comparison-{coverage,bias}.png and a summary CSV.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(concrete); library(data.table); library(ggplot2)
}))

cli <- commandArgs(trailingOnly = TRUE)
B <- if (length(cli) >= 1) as.integer(cli[1]) else 120L
n <- if (length(cli) >= 2) as.integer(cli[2]) else 400L
tau <- 1500          # RMST horizon
tau_max <- 2500      # administrative censoring
fig_dir <- "vignettes/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

params <- list(l1_0 = 5e-4, b1 = 0.3, g1 = 0.4, g1b = -0.2,
               l2_0 = 3e-4, b2 = -0.2, g2 = 0.3,
               lc = 4e-4)

haz <- function(A, W1, W2) {
  l1 <- params$l1_0 * exp(params$b1 * A + params$g1 * W1 + params$g1b * W2)
  l2 <- params$l2_0 * exp(params$b2 * A + params$g2 * W1)
  list(l1 = l1, l2 = l2, ltot = l1 + l2)
}

## --- truth: closed-form marginal RMST / LYL under each arm ------------------
truthFor <- function(a, Wn = 4e5) {
  W1 <- stats::rnorm(Wn); W2 <- stats::rnorm(Wn)
  h <- haz(a, W1, W2)
  rmst <- mean((1 - exp(-h$ltot * tau)) / h$ltot)
  lyl1 <- mean((h$l1 / h$ltot) * (tau - (1 - exp(-h$ltot * tau)) / h$ltot))
  c(RMST = rmst, LYL1 = lyl1)
}
set.seed(20260604)
truth1 <- truthFor(1); truth0 <- truthFor(0)
TRUTH <- list(
  RMST_A0 = truth0["RMST"], RMST_A1 = truth1["RMST"],
  RMSTdiff = truth1["RMST"] - truth0["RMST"],
  LYL1_A0 = truth0["LYL1"], LYL1_A1 = truth1["LYL1"])
cat(sprintf("True RMST: A0=%.1f A1=%.1f diff=%.1f\n",
            TRUTH$RMST_A0, TRUTH$RMST_A1, TRUTH$RMSTdiff))

simData <- function(seed) {
  set.seed(seed)
  W1 <- stats::rnorm(n); W2 <- stats::rnorm(n); A <- stats::rbinom(n, 1, 0.5)
  h <- haz(A, W1, W2)
  T1 <- stats::rexp(n, h$l1); T2 <- stats::rexp(n, h$l2); C <- stats::rexp(n, params$lc)
  Tobs <- pmin(T1, T2, C, tau_max)
  event <- ifelse(Tobs >= tau_max, 0L,
            ifelse(T1 <= pmin(T2, C), 1L, ifelse(T2 <= C, 2L, 0L)))
  data.table(id = seq_len(n), time = Tobs, event = event, arm = A, W1 = W1, W2 = W2)
}

cover <- function(pt, se, truth, z = 1.96) {
  is.finite(se) & (pt - z * se <= truth) & (truth <= pt + z * se)
}

oneRep <- function(seed) {
  d <- simData(seed)
  out <- tryCatch({
    args <- formatArguments(
      DataTable = d, EventTime = "time", EventType = "event", Treatment = "arm",
      ID = "id", Intervention = makeITT(), TargetTime = c(tau / 2, tau),
      TargetEvent = c(1, 2), CVArg = list(V = 2), UpdateMethod = "adaptive",
      EICStopRule = "absolute", MaxUpdateIter = 25, Verbose = FALSE)
    est <- doConcrete(args)
    pw <- getRMST(est, Horizon = tau, Intervention = c(1, 2))
    dr <- suppressWarnings(targetRMST(est, Horizon = tau, Intervention = c(1, 2),
                                      MaxUpdateIter = 60))
    grab <- function(tab, estd, ev, intv) {
      r <- tab[tab$Estimand == estd & tab$Event == ev & tab$Intervention == intv, ]
      if (!nrow(r)) return(c(pt = NA, se = NA))
      c(pt = r$`Pt Est`[1], se = r$se[1])
    }
    rbind(
      data.table(method = "pointwise", estimand = "RMST A1",
                 pt = grab(pw, "RMST", -1, "A=1")["pt"], se = grab(pw, "RMST", -1, "A=1")["se"],
                 truth = TRUTH$RMST_A1),
      data.table(method = "direct", estimand = "RMST A1",
                 pt = grab(dr, "RMST", -1, "A=1")["pt"], se = grab(dr, "RMST", -1, "A=1")["se"],
                 truth = TRUTH$RMST_A1),
      data.table(method = "pointwise", estimand = "RMST A0",
                 pt = grab(pw, "RMST", -1, "A=0")["pt"], se = grab(pw, "RMST", -1, "A=0")["se"],
                 truth = TRUTH$RMST_A0),
      data.table(method = "direct", estimand = "RMST A0",
                 pt = grab(dr, "RMST", -1, "A=0")["pt"], se = grab(dr, "RMST", -1, "A=0")["se"],
                 truth = TRUTH$RMST_A0),
      data.table(method = "pointwise", estimand = "LYL1 A1",
                 pt = grab(pw, "Life Years Lost", 1, "A=1")["pt"], se = grab(pw, "Life Years Lost", 1, "A=1")["se"],
                 truth = TRUTH$LYL1_A1),
      data.table(method = "direct", estimand = "LYL1 A1",
                 pt = grab(dr, "Life Years Lost", 1, "A=1")["pt"], se = grab(dr, "Life Years Lost", 1, "A=1")["se"],
                 truth = TRUTH$LYL1_A1)
    )
  }, error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out[, seed := seed]
  out
}

seeds <- 20260604L + seq_len(B)
useCores <- if (.Platform$OS.type == "unix") max(1L, parallel::detectCores() - 1L) else 1L
cat(sprintf("Running %d replicates (n=%d) on %d core(s)...\n", B, n, useCores))
reps <- if (useCores > 1L) {
  parallel::mclapply(seeds, oneRep, mc.cores = useCores)
} else {
  lapply(seeds, oneRep)
}
res <- data.table::rbindlist(Filter(Negate(is.null), reps))

summary <- res[, .(
  Reps = .N,
  Bias = mean(pt - truth, na.rm = TRUE),
  Coverage = mean(cover(pt, se, truth), na.rm = TRUE)
), by = .(estimand, method)]
data.table::fwrite(summary, file.path(fig_dir, "rmst-comparison-summary.csv"))
print(summary)

theme_set(theme_minimal(base_size = 12) + theme(strip.text = element_text(face = "bold")))

g_cov <- ggplot(summary, aes(estimand, Coverage, fill = method)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey30") +
  scale_fill_manual(values = c(direct = "#1b9e77", pointwise = "#d95f02")) +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "RMST confidence-interval coverage: direct vs pointwise targeting",
       subtitle = sprintf("Sparse 2-time target grid, horizon = %d, n = %d, %d reps. Dashed = nominal 0.95.",
                          tau, n, B),
       x = NULL, y = "Empirical coverage", fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "rmst-comparison-coverage.png"), g_cov, width = 8, height = 4,
       dpi = 150, bg = "white")

g_bias <- ggplot(summary, aes(estimand, Bias, fill = method)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c(direct = "#1b9e77", pointwise = "#d95f02")) +
  labs(title = "RMST bias: direct vs pointwise targeting",
       subtitle = sprintf("Sparse 2-time target grid, horizon = %d, n = %d, %d reps.", tau, n, B),
       x = NULL, y = "Bias (days)", fill = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "rmst-comparison-bias.png"), g_bias, width = 8, height = 4,
       dpi = 150, bg = "white")

cat("Done. Wrote rmst-comparison-{coverage,bias}.png and -summary.csv to ", fig_dir, "\n")
