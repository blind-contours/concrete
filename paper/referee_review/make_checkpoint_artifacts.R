suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(knitr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

root <- normalizePath(getwd(), mustWork = TRUE)
if (!dir.exists(file.path(root, "scripts"))) {
  stop("Run this script from the repository root.")
}

review_dir <- file.path(root, "paper", "referee_review")
artifact_dir <- file.path(review_dir, "checkpoint_artifacts")
table_dir <- file.path(artifact_dir, "tables")
figure_dir <- file.path(artifact_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

event1_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "event1_b100_n500")
survsl_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "survsl_nohal_b20_n500_all_events")
hal_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "hal_sensitivity_nonph_b3_n300")
failed_rep_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "failed_seed_update_sensitivity_representative")
failed_best_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "failed_seed_update_sensitivity_best_all")
rare_stop_rep_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "rare_event_stop_sensitivity_representative")
rare_stop_candidate_dir <- file.path(root, "scripts", "sim-data", "referee-sims", "output", "rare_event_stop_sensitivity_candidate_all")

scenario_label <- function(x) {
  map <- c(
    nonph = "Non-PH",
    ph_correct = "PH correct",
    rare_early = "Rare early",
    positivity = "Positivity"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

estimator_label <- function(x) {
  map <- c(
    concrete_standard_minimal = "Concrete std minimal",
    concrete_standard_rich = "Concrete std rich",
    concrete_adaptive_rich = "Concrete adaptive rich",
    concrete_adaptive_rich_stabilized = "Concrete adaptive rich stabilized",
    concrete_standard_survsl_nohal = "Concrete std survSL no-HAL",
    concrete_adaptive_survsl_nohal = "Concrete adaptive survSL no-HAL",
    concrete_adaptive_survsl_nohal_stabilized = "Concrete adaptive survSL no-HAL stabilized",
    concrete_standard_survsl = "Concrete std survSL + HAL",
    aalen_johansen = "Aalen-Johansen",
    survtmle_6mo = "survtmle 6 mo"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

config_label <- function(x) {
  map <- c(
    standard_rich_e0.10_iter200 = "Std rich, eps 0.10, 200",
    adaptive_rich_e0.10_iter200 = "Adapt rich, eps 0.10, 200",
    adaptive_rich_e0.01_iter300 = "Adapt rich, eps 0.01, 300",
    adaptive_rich_e0.01_min0.05_iter300 = "Adapt rich, eps 0.01, min 0.05, 300"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

stop_config_label <- function(x) {
  map <- c(
    relative_best = "Relative",
    `hybrid_abs2.5e-4` = "Hybrid abs 2.5e-4",
    `hybrid_abs5e-4` = "Hybrid abs 5e-4",
    `hybrid_abs1e-3` = "Hybrid abs 1e-3",
    `absolute_abs1e-3` = "Absolute abs 1e-3"
  )
  unname(ifelse(x %in% names(map), map[x], x))
}

case_label <- function(x) {
  gsub("_", " ", x, fixed = TRUE)
}

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
}

write_kable <- function(dt, file, digits = 3) {
  tex <- knitr::kable(
    as.data.frame(dt),
    format = "latex",
    booktabs = TRUE,
    escape = TRUE,
    linesep = "",
    digits = digits
  )
  writeLines(tex, file.path(table_dir, file))
}

theme_checkpoint <- function() {
  theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "grey92", color = "grey60")
    )
}

save_plot <- function(plot, file, width = 7.2, height = 4.3) {
  ggsave(
    filename = file.path(figure_dir, file),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

event1_diag <- fread(file.path(event1_dir, "diagnostic_metrics.csv"))
event1_metrics <- fread(file.path(event1_dir, "metrics.csv"))
survsl_diag <- fread(file.path(survsl_dir, "diagnostic_metrics.csv"))
survsl_long <- fread(file.path(survsl_dir, "diagnostics_long.csv"))
survsl_metrics <- fread(file.path(survsl_dir, "metrics.csv"))
hal_diag <- fread(file.path(hal_dir, "diagnostic_metrics.csv"))
hal_long <- fread(file.path(hal_dir, "diagnostics_long.csv"))
hal_metrics <- fread(file.path(hal_dir, "metrics.csv"))

event1_conv <- event1_diag[grepl("^concrete_", EstimatorID)]
event1_conv[, `:=`(
  Scenario = scenario_label(Scenario),
  Estimator = estimator_label(EstimatorID),
  `Conv.` = fmt(ConvergenceRate, 2),
  `Median step` = fmt(MedianStep, 1),
  `Mean runtime (s)` = fmt(MeanRuntimeSec, 1)
)]
write_kable(
  event1_conv[, .(Scenario, Estimator, Reps, `Conv.`, `Median step`, `Mean runtime (s)`)],
  "tab_event1_convergence.tex"
)

event1_perf <- event1_metrics[
  Estimator %in% c("tmle", "survtmle_6mo", "Aalen-Johansen") &
    ((Estimand == "Abs Risk" & Intervention == "A=1") | Estimand == "Risk Diff"),
  .(
    `Mean abs. bias` = mean(abs(Bias), na.rm = TRUE),
    `Mean RMSE` = mean(RMSE, na.rm = TRUE),
    `Mean coverage` = mean(Coverage, na.rm = TRUE)
  ),
  by = .(Scenario, Estimand, EstimatorID)
]
event1_perf[, `:=`(
  Scenario = scenario_label(Scenario),
  Estimator = estimator_label(EstimatorID),
  `Mean abs. bias` = fmt(`Mean abs. bias`, 3),
  `Mean RMSE` = fmt(`Mean RMSE`, 3),
  `Mean coverage` = fmt(`Mean coverage`, 2)
)]
setorder(event1_perf, Scenario, Estimand, Estimator)
write_kable(
  event1_perf[, .(Scenario, Estimand, Estimator, `Mean abs. bias`, `Mean RMSE`, `Mean coverage`)],
  "tab_event1_performance.tex"
)

survsl_conv <- survsl_diag[grepl("^concrete_", EstimatorID)]
survsl_conv[, `:=`(
  Scenario = scenario_label(Scenario),
  Estimator = estimator_label(EstimatorID),
  `Conv.` = fmt(ConvergenceRate, 2),
  `Median step` = fmt(MedianStep, 1),
  `Mean runtime (s)` = fmt(MeanRuntimeSec, 1)
)]
write_kable(
  survsl_conv[, .(Scenario, Estimator, Reps, `Conv.`, `Median step`, `Mean runtime (s)`)],
  "tab_survsl_convergence.tex"
)

survsl_abs <- survsl_metrics[
  Estimator %in% c("tmle", "Aalen-Johansen") & Estimand == "Abs Risk",
  .(
    `Mean abs. bias` = mean(abs(Bias), na.rm = TRUE),
    `Mean RMSE` = mean(RMSE, na.rm = TRUE),
    `Mean coverage` = mean(Coverage, na.rm = TRUE)
  ),
  by = .(Scenario, EstimatorID)
]
survsl_abs[, `:=`(
  Scenario = scenario_label(Scenario),
  Estimator = estimator_label(EstimatorID)
)]
best_survsl_abs <- survsl_abs[grepl("^Concrete", Estimator), .SD[which.min(`Mean RMSE`)], by = Scenario]
aj_survsl_abs <- survsl_abs[Estimator == "Aalen-Johansen", .(Scenario, `AJ mean RMSE` = `Mean RMSE`)]
best_survsl_abs <- merge(best_survsl_abs, aj_survsl_abs, by = "Scenario")
best_survsl_abs[, `:=`(
  `Mean abs. bias` = fmt(`Mean abs. bias`, 3),
  `Mean RMSE` = fmt(`Mean RMSE`, 3),
  `Mean coverage` = fmt(`Mean coverage`, 2),
  `AJ mean RMSE` = fmt(`AJ mean RMSE`, 3)
)]
write_kable(
  best_survsl_abs[, .(Scenario, `Best concrete TMLE` = Estimator, `Mean abs. bias`, `Mean RMSE`, `Mean coverage`, `AJ mean RMSE`)],
  "tab_survsl_absrisk_best.tex"
)

survsl_rd <- survsl_metrics[
  Estimator == "tmle" & Estimand == "Risk Diff",
  .(
    `Mean abs. bias` = mean(abs(Bias), na.rm = TRUE),
    `Mean RMSE` = mean(RMSE, na.rm = TRUE),
    `Mean coverage` = mean(Coverage, na.rm = TRUE)
  ),
  by = .(Scenario, EstimatorID)
]
survsl_rd[, `:=`(
  Scenario = scenario_label(Scenario),
  Estimator = estimator_label(EstimatorID),
  `Mean abs. bias` = fmt(`Mean abs. bias`, 3),
  `Mean RMSE` = fmt(`Mean RMSE`, 3),
  `Mean coverage` = fmt(`Mean coverage`, 2)
)]
setorder(survsl_rd, Scenario, Estimator)
write_kable(
  survsl_rd[, .(Scenario, Estimator, `Mean abs. bias`, `Mean RMSE`, `Mean coverage`)],
  "tab_survsl_riskdiff.tex"
)

nonconv <- survsl_long[EstimatorID != "aalen_johansen" & Converged == FALSE,
  .(Estimators = paste(estimator_label(EstimatorID), collapse = "; "), `Number failed` = .N),
  by = .(Scenario, Rep, Seed)
]
nonconv[, Scenario := scenario_label(Scenario)]
setorder(nonconv, Scenario, Rep)
write_kable(
  nonconv[, .(Scenario, Rep, Seed, `Number failed`, Estimators)],
  "tab_survsl_nonconvergence_seeds.tex"
)

hal_table <- hal_diag[EstimatorID != "aalen_johansen"]
hal_table[, `:=`(
  Estimator = estimator_label(EstimatorID),
  `Conv.` = fmt(ConvergenceRate, 2),
  `Median step` = fmt(MedianStep, 1),
  `Mean runtime (s)` = fmt(MeanRuntimeSec, 2)
)]
hal_max <- hal_long[EstimatorID != "aalen_johansen", .(`Max step` = max(Step, na.rm = TRUE)), by = EstimatorID]
hal_table <- merge(hal_table, hal_max, by = "EstimatorID", all.x = TRUE)
write_kable(
  hal_table[, .(Estimator, Reps, `Conv.`, `Median step`, `Max step`, `Mean runtime (s)`)],
  "tab_hal_diagnostics.tex"
)

hal_abs <- hal_metrics[
  Estimator %in% c("tmle", "Aalen-Johansen") & Estimand == "Abs Risk",
  .(
    `Mean abs. bias` = mean(abs(Bias), na.rm = TRUE),
    `Mean RMSE` = mean(RMSE, na.rm = TRUE),
    `Mean coverage` = mean(Coverage, na.rm = TRUE)
  ),
  by = EstimatorID
]
hal_abs[, `:=`(
  Estimator = estimator_label(EstimatorID),
  `Mean abs. bias` = fmt(`Mean abs. bias`, 3),
  `Mean RMSE` = fmt(`Mean RMSE`, 3),
  `Mean coverage` = fmt(`Mean coverage`, 2)
)]
write_kable(
  hal_abs[, .(Estimator, `Mean abs. bias`, `Mean RMSE`, `Mean coverage`)],
  "tab_hal_absrisk.tex"
)

event1_plot <- copy(event1_diag[grepl("^concrete_", EstimatorID)])
event1_plot[, `:=`(Scenario = scenario_label(Scenario), Estimator = estimator_label(EstimatorID))]
save_plot(
  ggplot(event1_plot, aes(x = Scenario, y = ConvergenceRate, fill = Estimator)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    coord_cartesian(ylim = c(0.9, 1.01)) +
    labs(x = NULL, y = "Convergence rate", title = "Event-1 B100 convergence") +
    theme_checkpoint(),
  "fig_event1_convergence.pdf",
  width = 6.8,
  height = 3.8
)

survsl_plot <- copy(survsl_diag[grepl("^concrete_", EstimatorID)])
survsl_plot[, `:=`(Scenario = scenario_label(Scenario), Estimator = estimator_label(EstimatorID))]
save_plot(
  ggplot(survsl_plot, aes(x = Scenario, y = ConvergenceRate, fill = Estimator)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.65) +
    coord_cartesian(ylim = c(0.7, 1.02)) +
    labs(x = NULL, y = "Convergence rate", title = "Survival-SL B20 convergence") +
    theme_checkpoint(),
  "fig_survsl_convergence.pdf",
  width = 7.2,
  height = 4.2
)

survsl_rmse_plot <- copy(survsl_abs)
save_plot(
  ggplot(survsl_rmse_plot, aes(x = Estimator, y = `Mean RMSE`, fill = Estimator)) +
    geom_col(width = 0.68) +
    facet_wrap(~Scenario, scales = "free_x") +
    labs(x = NULL, y = "Mean RMSE", title = "Absolute-risk RMSE in survival-SL grid") +
    theme_checkpoint() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7)),
  "fig_survsl_absrisk_rmse.pdf",
  width = 8.0,
  height = 4.6
)

survsl_runtime_plot <- copy(survsl_diag[grepl("^concrete_", EstimatorID)])
survsl_runtime_plot[, `:=`(Scenario = scenario_label(Scenario), Estimator = estimator_label(EstimatorID))]
save_plot(
  ggplot(survsl_runtime_plot, aes(x = Estimator, y = MeanRuntimeSec, fill = Estimator)) +
    geom_col(width = 0.68) +
    facet_wrap(~Scenario) +
    labs(x = NULL, y = "Mean runtime (s)", title = "Runtime in survival-SL grid") +
    theme_checkpoint() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7)),
  "fig_survsl_runtime.pdf",
  width = 8.0,
  height = 4.6
)

event1_perf_plot <- event1_metrics[
  Estimator %in% c("tmle", "survtmle_6mo", "Aalen-Johansen") &
    Estimand == "Risk Diff",
  .(`Mean RMSE` = mean(RMSE, na.rm = TRUE)),
  by = .(Scenario, EstimatorID)
]
event1_perf_plot[, `:=`(Scenario = scenario_label(Scenario), Estimator = estimator_label(EstimatorID))]
save_plot(
  ggplot(event1_perf_plot, aes(x = Estimator, y = `Mean RMSE`, fill = Estimator)) +
    geom_col(width = 0.68) +
    facet_wrap(~Scenario) +
    labs(x = NULL, y = "Mean RMSE", title = "Event-1 risk-difference RMSE") +
    theme_checkpoint() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7)),
  "fig_event1_riskdiff_rmse.pdf",
  width = 7.6,
  height = 4.2
)

if (file.exists(file.path(failed_rep_dir, "diagnostics.csv"))) {
  failed_rep_diag <- fread(file.path(failed_rep_dir, "diagnostics.csv"))
  failed_rep_diag[, `:=`(
    Case = case_label(CaseID),
    Config = config_label(ConfigID),
    Converged = ifelse(Converged, "yes", "no"),
    Step = fmt(Step, 0),
    `Max ratio` = fmt(MaxRatio, 2),
    `Final norm` = fmt(FinalNormPnEIC, 4),
    `Runtime (s)` = fmt(RuntimeSec, 1)
  )]
  setorder(failed_rep_diag, CaseID, ConfigID)
  write_kable(
    failed_rep_diag[, .(Case, Config, Converged, Step, `Final norm`, `Max ratio`, `Runtime (s)`)],
    "tab_failed_seed_representative.tex"
  )
}

if (file.exists(file.path(failed_best_dir, "diagnostics.csv"))) {
  failed_best_diag <- fread(file.path(failed_best_dir, "diagnostics.csv"))
  failed_best_components <- fread(file.path(failed_best_dir, "final_components.csv"))
  failed_best_trace <- fread(file.path(failed_best_dir, "trace_summary.csv"))

  failed_best_diag[, `:=`(
    Case = case_label(CaseID),
    Scenario = scenario_label(Scenario),
    Converged = ifelse(Converged, "yes", "no"),
    Step = fmt(Step, 0),
    `Final norm` = fmt(FinalNormPnEIC, 4),
    `Max ratio` = fmt(MaxRatio, 2),
    `Runtime (s)` = fmt(RuntimeSec, 1)
  )]
  setorder(failed_best_diag, Scenario, Rep)
  write_kable(
    failed_best_diag[, .(Case, Scenario, Converged, Step, `Final norm`, `Max ratio`, `Runtime (s)`)],
    "tab_failed_seed_best_all.tex"
  )

  top_components <- failed_best_components[, .SD[which.max(Ratio)], by = .(CaseID, ConfigID)]
  top_components[, `:=`(
    Case = case_label(CaseID),
    `Worst component` = paste0(Intervention, ", event ", Event, ", day ", Time),
    `PnEIC` = fmt(PnEIC, 6),
    `Threshold` = fmt(`seEIC/(sqrt(n)log(n))`, 6),
    Ratio = fmt(Ratio, 2),
    seEIC = fmt(seEIC, 4)
  )]
  setorder(top_components, CaseID)
  write_kable(
    top_components[, .(Case, `Worst component`, PnEIC, Threshold, Ratio, seEIC)],
    "tab_failed_seed_worst_components.tex"
  )

  failed_best_trace[, `:=`(
    Case = case_label(CaseID),
    `Accepted steps` = AcceptedSteps,
    `Rejected steps` = RejectedSteps,
    `Initial norm` = fmt(InitialNorm, 4),
    `Final norm` = fmt(FinalNorm, 4),
    `Final max ratio` = fmt(FinalMaxRatio, 2),
    `Median alpha` = fmt(MedianAcceptedAlpha, 4)
  )]
  setorder(failed_best_trace, Scenario, Rep)
  write_kable(
    failed_best_trace[, .(Case, `Accepted steps`, `Rejected steps`, `Initial norm`, `Final norm`, `Final max ratio`, `Median alpha`)],
    "tab_failed_seed_trace_summary.tex"
  )

  failed_best_plot <- copy(failed_best_diag)
  failed_best_plot[, Case := factor(Case, levels = Case)]
  save_plot(
    ggplot(failed_best_plot, aes(x = Case, y = as.numeric(MaxRatio), fill = Scenario)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 1, linetype = 2) +
      scale_y_log10() +
      labs(x = NULL, y = "Final max ratio (log scale)", title = "Best candidate update on all failed seeds") +
      theme_checkpoint(),
    "fig_failed_seed_best_all_maxratio.pdf",
    width = 7.4,
    height = 4.1
  )
}

if (file.exists(file.path(rare_stop_rep_dir, "diagnostics.csv"))) {
  rare_rep_diag <- fread(file.path(rare_stop_rep_dir, "diagnostics.csv"))
  rare_rep_diag[, `:=`(
    Case = case_label(CaseID),
    Config = stop_config_label(ConfigID),
    Converged = ifelse(Converged, "yes", "no"),
    Step = fmt(Step, 0),
    `Final norm` = fmt(FinalNormPnEIC, 4),
    `Stop ratio` = fmt(MaxRatio, 2),
    `Relative ratio` = fmt(MaxRelativeRatio, 2),
    `Max |PnEIC|` = fmt(MaxAbsPnEIC, 6),
    `Runtime (s)` = fmt(RuntimeSec, 1)
  )]
  setorder(rare_rep_diag, EICStopAbsTol, ConfigID)
  write_kable(
    rare_rep_diag[, .(Case, Config, Converged, Step, `Final norm`, `Stop ratio`, `Relative ratio`, `Max |PnEIC|`, `Runtime (s)`)],
    "tab_rare_event_stop_representative.tex"
  )

  rare_rep_plot <- copy(rare_rep_diag)
  rare_rep_plot[, Config := factor(Config, levels = unique(Config))]
  save_plot(
    ggplot(rare_rep_plot, aes(x = Config, y = as.numeric(MaxRatio), fill = Converged)) +
      geom_col(width = 0.68) +
      geom_hline(yintercept = 1, linetype = 2) +
      labs(x = NULL, y = "Final stopping-rule ratio", title = "Representative rare-early stopping-rule sensitivity") +
      theme_checkpoint(),
    "fig_rare_event_stop_representative_ratio.pdf",
    width = 7.0,
    height = 4.0
  )
}

if (file.exists(file.path(rare_stop_candidate_dir, "diagnostics.csv"))) {
  rare_candidate_diag <- fread(file.path(rare_stop_candidate_dir, "diagnostics.csv"))
  rare_candidate_components <- fread(file.path(rare_stop_candidate_dir, "final_components.csv"))

  rare_candidate_diag[, `:=`(
    Case = case_label(CaseID),
    Converged = ifelse(Converged, "yes", "no"),
    Step = fmt(Step, 0),
    `Final norm` = fmt(FinalNormPnEIC, 4),
    `Stop ratio` = fmt(MaxRatio, 2),
    `Relative ratio` = fmt(MaxRelativeRatio, 2),
    `Max |PnEIC|` = fmt(MaxAbsPnEIC, 6),
    `Runtime (s)` = fmt(RuntimeSec, 1)
  )]
  setorder(rare_candidate_diag, Rep)
  write_kable(
    rare_candidate_diag[, .(Case, Converged, Step, `Final norm`, `Stop ratio`, `Relative ratio`, `Max |PnEIC|`, `Runtime (s)`)],
    "tab_rare_event_stop_candidate_all.tex"
  )

  rare_top_components <- rare_candidate_components[, .SD[which.max(RelativeRatio)], by = .(CaseID, ConfigID)]
  rare_top_components[, `:=`(
    Case = case_label(CaseID),
    `Worst relative component` = paste0(Intervention, ", event ", Event, ", day ", Time),
    `PnEIC` = fmt(PnEIC, 6),
    `Stop criterion` = fmt(StopCriteria, 6),
    `Stop ratio` = fmt(Ratio, 2),
    `Relative ratio` = fmt(RelativeRatio, 2),
    `Relative threshold` = fmt(`seEIC/(sqrt(n)log(n))`, 6)
  )]
  setorder(rare_top_components, CaseID)
  write_kable(
    rare_top_components[, .(Case, `Worst relative component`, PnEIC, `Stop criterion`, `Stop ratio`, `Relative threshold`, `Relative ratio`)],
    "tab_rare_event_stop_candidate_components.tex"
  )

  rare_candidate_plot <- copy(rare_candidate_diag)
  rare_candidate_plot[, Case := factor(Case, levels = Case)]
  save_plot(
    ggplot(rare_candidate_plot, aes(x = Case, y = MaxAbsPnEIC)) +
      geom_col(fill = "#4C78A8", width = 0.68) +
      geom_hline(yintercept = 1e-3, linetype = 2) +
      labs(x = NULL, y = "Final max |PnEIC|", title = "Hybrid 1e-3 stopping rule on rare-early failed seeds") +
      theme_checkpoint(),
    "fig_rare_event_stop_candidate_abs_pneic.pdf",
    width = 7.2,
    height = 4.0
  )
}

cat("Wrote checkpoint artifacts to", artifact_dir, "\n")
