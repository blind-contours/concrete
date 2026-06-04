#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# make-sim-evidence-figures.R
#
# Turns the committed referee-simulation `metrics.csv` summaries into the
# figures embedded in the "Simulation evidence" vignette. No simulation is run
# here; this only plots the precomputed metrics so the figures are reproducible
# from the committed results.
#
#   Rscript scripts/make-sim-evidence-figures.R
#
# Primary source: a 100-replicate run over three scenarios (n = 500).
# Secondary source: a 20-replicate run that additionally includes a
# positivity-stress scenario.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(data.table)
  library(ggplot2)
}))

base <- "scripts/sim-data/referee-sims/output"
primary <- file.path(base, "event1_b100_n500/metrics.csv")
allscn  <- file.path(base, "survsl_nohal_stabilized_b20_n500_all_events/metrics.csv")

fig_dir <- "vignettes/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

scenario_labels <- c(
  ph_correct = "Proportional hazards",
  nonph      = "Non-proportional hazards",
  rare_early = "Rare early events",
  positivity = "Positivity / informative censoring"
)
relabel <- function(x) factor(scenario_labels[x], levels = scenario_labels)

theme_set(theme_minimal(base_size = 12) +
            theme(panel.grid.minor = element_blank(),
                  strip.text = element_text(face = "bold")))

# Keep the two concrete TMLE libraries and the g-computation plug-in, focusing
# on the absolute-risk estimand for the headline bias/coverage story.
load_concrete <- function(path) {
  m <- as.data.table(read.csv(path, check.names = FALSE))
  m <- m[Package == "concrete"]
  m[, ScenarioLab := relabel(Scenario)]
  m[, EstLabel := fifelse(Estimator == "gcomp", "g-computation (plug-in)",
                   fifelse(ModelLibrary == "rich", "TMLE (rich library)",
                                                    "TMLE (Cox library)"))]
  m[]
}

mp <- load_concrete(primary)

# --- 1. Bias vs follow-up time: TMLE removes plug-in bias -------------------
bias_dat <- mp[Estimand == "Abs Risk" & Intervention %in% c("A=0", "A=1")]
bias_dat <- bias_dat[, .(Bias = mean(Bias)),
                     by = .(ScenarioLab, Time, EstLabel)]
g_bias <- ggplot(bias_dat, aes(Time, Bias, colour = EstLabel, shape = EstLabel)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(alpha = 0.8) + geom_point(size = 2) +
  facet_wrap(~ScenarioLab, nrow = 1) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "Bias of the absolute-risk estimate vs follow-up time",
       subtitle = "100 replicates, n = 500. TMLE targets away the plug-in bias.",
       x = "Target time", y = "Bias (estimate - truth)", colour = NULL, shape = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "sim-bias-vs-time.png"), g_bias, width = 9, height = 4,
       dpi = 150, bg = "white")
message("wrote sim-bias-vs-time.png")

# --- 2. CI coverage vs nominal 95% -----------------------------------------
cov_dat <- mp[Estimator == "tmle" & Estimand == "Abs Risk" &
                Intervention %in% c("A=0", "A=1") & !is.na(Coverage)]
cov_dat <- cov_dat[, .(Coverage = mean(Coverage)),
                   by = .(ScenarioLab, Time, EstLabel)]
g_cov <- ggplot(cov_dat, aes(Time, Coverage, colour = EstLabel, shape = EstLabel)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey30") +
  geom_line(alpha = 0.8) + geom_point(size = 2) +
  facet_wrap(~ScenarioLab, nrow = 1) +
  scale_colour_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(title = "95% confidence-interval coverage of the TMLE estimate",
       subtitle = "Dashed line = nominal 0.95. 100 replicates, n = 500. Coverage degrades at long horizons in the non-PH scenario.",
       x = "Target time", y = "Empirical coverage", colour = NULL, shape = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "sim-coverage.png"), g_cov, width = 9, height = 4,
       dpi = 150, bg = "white")
message("wrote sim-coverage.png")

# --- 3. Standard-error calibration: model SE vs empirical SD ----------------
se_dat <- mp[Estimator == "tmle" & Estimand == "Abs Risk" &
               Intervention %in% c("A=0", "A=1") & !is.na(MeanSE)]
lim <- range(c(se_dat$MeanSE, se_dat$EmpiricalSD), na.rm = TRUE)
g_se <- ggplot(se_dat, aes(EmpiricalSD, MeanSE, colour = ScenarioLab)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey40") +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_brewer(palette = "Set2") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(title = "Influence-function SE vs empirical sampling SD",
       subtitle = "Points on the line indicate well-calibrated standard errors (TMLE).",
       x = "Empirical SD across replicates", y = "Mean estimated SE", colour = NULL) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "sim-se-calibration.png"), g_se, width = 6.5, height = 5.2,
       dpi = 150, bg = "white")
message("wrote sim-se-calibration.png")

# --- 4. Four-scenario coverage including positivity stress (20 reps) --------
if (file.exists(allscn)) {
  ma <- load_concrete(allscn)
  cov4 <- ma[Estimator == "tmle" & Estimand == "Abs Risk" &
               Intervention %in% c("A=0", "A=1") & !is.na(Coverage)]
  cov4 <- cov4[, .(Coverage = mean(Coverage)), by = .(ScenarioLab, Time)]
  g_cov4 <- ggplot(cov4, aes(Time, Coverage)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey30") +
    geom_line(colour = "#1b9e77") + geom_point(size = 2, colour = "#1b9e77") +
    facet_wrap(~ScenarioLab, nrow = 1) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(title = "TMLE coverage across four scenarios, including positivity stress",
         subtitle = "Dashed line = nominal 0.95. 20 replicates, n = 500 (coverage noisier at 20 reps).",
         x = "Target time", y = "Empirical coverage")
  ggsave(file.path(fig_dir, "sim-coverage-4scenario.png"), g_cov4, width = 10, height = 3.6,
         dpi = 150, bg = "white")
  message("wrote sim-coverage-4scenario.png")
}

message("Done. Simulation-evidence figures written to ", fig_dir, ".")
