#!/usr/bin/env Rscript
# Coverage / validation figures for the win-ratio family, for the trialist
# vignette. Values are the validated simulation results from the committed
# validation scripts (re-runnable):
#   make-winratio-coverage.R            single-event win ratio
#   make-clinical-wr-eif-validation.R   clinical WR, no random censoring (B=60)
#   make-clinical-wr-censoring-validation.R  clinical WR, ~36% censoring (B=120)
#   make-pathprob-censoring-validation.R     path-probability Theta, censoring (B=150)
# All against closed-form / brute-force pairwise ground truth.
suppressWarnings(suppressMessages({ library(data.table); library(ggplot2) }))
fig_dir <- "vignettes/figures"; dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ---- validated coverage of the 95% influence-function CIs ----
cov <- data.table(
  scenario = c("Single-event WR\n(shipped getWinRatio)",
               "Single-event WR\n(shipped getWinRatio)",
               "Single-event WR\n(shipped getWinRatio)",
               "Clinical WR\n(no censoring)",
               "Clinical WR\n(~36% censoring)",
               "Path-prob Θ\n(~36% censoring)"),
  estimand = c("Win Ratio", "Win Odds", "Net Benefit",
               "Win Ratio", "Win Ratio", "P(alive, prior HFH)"),
  coverage = c(0.967, 0.917, 0.917, 0.950, 0.967, 0.960),
  B        = c(60, 60, 60, 60, 120, 150))
cov[, mc := 1.96 * sqrt(0.95 * 0.05 / B)]                 # Monte-Carlo error band
cov[, lab := factor(scenario, levels = unique(scenario))]

gcov <- ggplot(cov, aes(x = coverage, y = lab, colour = estimand)) +
  geom_vline(xintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = coverage - mc, xmax = coverage + mc), height = 0.18,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.5)) +
  scale_x_continuous(limits = c(0.85, 1.0), breaks = seq(0.85, 1, 0.05)) +
  labs(title = "Win-ratio family: 95% CI coverage vs. ground truth",
       subtitle = "Dashed line = nominal 0.95; bars = Monte-Carlo error. Validated against\nclosed-form / brute-force pairwise truth.",
       x = "Empirical coverage of the 95% CI", y = NULL, colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "winratio-coverage.png"), gcov, width = 7.5, height = 4.2, dpi = 150)

## ---- clinical WR: point estimate vs brute-force truth ----
pt <- data.table(
  setting = c("No censoring", "~36% censoring"),
  truth   = c(1.585, 1.585),
  est     = c(1.596, 1.585),                              # single draw (uncens) / B=120 mean (cens)
  lo      = c(1.490, 1.520), hi = c(1.700, 1.652))        # representative 95% CI
pt[, setting := factor(setting, levels = setting)]
gpt <- ggplot(pt, aes(x = setting)) +
  geom_hline(aes(yintercept = truth[1]), linetype = 2, colour = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, colour = "#3366aa") +
  geom_point(aes(y = est), size = 3, colour = "#3366aa") +
  annotate("text", x = 0.6, y = 1.585, label = "truth = 1.585", vjust = -0.6, size = 3.4, colour = "grey30") +
  labs(title = "Clinical (death-priority) win ratio vs. ground truth",
       subtitle = "Estimate ± 95% CI recovers the brute-force pairwise truth, with and without censoring.",
       x = NULL, y = "Win ratio (treated vs control)") +
  theme_minimal(base_size = 12)
ggsave(file.path(fig_dir, "clinical-winratio-truth.png"), gpt, width = 6.5, height = 4, dpi = 150)
cat("Wrote winratio-coverage.png and clinical-winratio-truth.png to", fig_dir, "\n")
