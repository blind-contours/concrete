source("scripts/sim-data/referee-sims/referee_sim_functions.R")

args <- parse_cli_args(list(
  B = "6",
  n = "500",
  scenarios = "ph_correct,nonph,rare_early",
  target_times = "180,365,730,1095,1460",
  target_events = "1,2",
  cv_v = "2",
  truth_n = "100000",
  max_update_iter = "200",
  include_survsl = "true",
  include_hal_sl = "false",
  include_stabilized = "false",
  stabilized_min_nuisance = "0.05",
  stabilized_one_step_eps = "0.01",
  stabilized_stop_rule = "hybrid",
  stabilized_stop_abs_tol = "0.001",
  include_survtmle = "false",
  survtmle_months = "6",
  cores = "1",
  out = "scripts/sim-data/referee-sims/output/pilot"
))

B <- as.integer(args$B)
n <- as.integer(args$n)
scenarios <- split_arg(args$scenarios)
target_times <- as.numeric(split_arg(args$target_times))
target_events <- as.numeric(split_arg(args$target_events))
cv_v <- as.integer(args$cv_v)
truth_n <- as.integer(args$truth_n)
max_update_iter <- as.integer(args$max_update_iter)
include_survsl <- as_bool(args$include_survsl)
include_hal_sl <- as_bool(args$include_hal_sl)
include_stabilized <- as_bool(args$include_stabilized)
stabilized_min_nuisance <- as.numeric(args$stabilized_min_nuisance)
stabilized_one_step_eps <- as.numeric(args$stabilized_one_step_eps)
stabilized_stop_rule <- args$stabilized_stop_rule
stabilized_stop_abs_tol <- as.numeric(args$stabilized_stop_abs_tol)
include_survtmle <- as_bool(args$include_survtmle)
survtmle_months <- as.numeric(split_arg(args$survtmle_months))
cores <- max(1L, as.integer(args$cores))
out_dir <- args$out

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "replicates"), recursive = TRUE, showWarnings = FALSE)

cat("Referee simulation pilot\n")
cat("Output:", out_dir, "\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("B:", B, " n:", n, " CV folds:", cv_v, "\n")
cat("Target times:", paste(target_times, collapse = ", "), "\n")
cat("Target events:", paste(target_events, collapse = ", "), "\n")
cat("survtmle:", include_survtmle, "\n\n")
cat("survival SL:", include_survsl, "\n")
cat("HAL in survival SL:", include_hal_sl, "\n\n")
cat("stabilized adaptive update:", include_stabilized, "\n")
if (include_stabilized) {
  cat("  min nuisance:", stabilized_min_nuisance,
      " eps:", stabilized_one_step_eps,
      " stop rule:", stabilized_stop_rule,
      " abs tol:", stabilized_stop_abs_tol, "\n\n")
}
cat("cores:", cores, "\n\n")
cat("max update iterations:", max_update_iter, "\n\n")

truth_list <- list()
for (scenario in scenarios) {
  truth_path <- file.path(out_dir, paste0("truth_", scenario, ".rds"))
  if (file.exists(truth_path)) {
    truth <- readRDS(truth_path)
  } else {
    cat("Estimating truth for", scenario, "with n_truth =", truth_n, "\n")
    truth <- estimate_truth(
      scenario = scenario,
      target_times = target_times,
      target_events = target_events,
      n_truth = truth_n
    )
    saveRDS(truth, truth_path)
    fwrite(truth, sub("\\.rds$", ".csv", truth_path))
  }
  truth_list[[scenario]] <- truth
}
truth_all <- rbindlist(truth_list, fill = TRUE)
saveRDS(truth_all, file.path(out_dir, "truth_all.rds"))
fwrite(truth_all, file.path(out_dir, "truth_all.csv"))

seed_grid_path <- file.path(out_dir, "seed_grid.csv")
if (file.exists(seed_grid_path)) {
  seed_grid <- fread(seed_grid_path)
  if (!all(c("Scenario", "Rep", "Seed") %in% names(seed_grid))) {
    stop("Existing seed grid must contain Scenario, Rep, and Seed columns.")
  }
  seed_grid[, Scenario := as.character(Scenario)]
  seed_grid[, Rep := as.integer(Rep)]
  seed_grid[, Seed := as.integer(Seed)]
  expected_grid <- CJ(Scenario = scenarios, Rep = seq_len(B))
  observed_grid <- seed_grid[, .(Scenario, Rep)][order(Scenario, Rep)]
  expected_grid <- expected_grid[order(Scenario, Rep)]
  if (!isTRUE(all.equal(observed_grid, expected_grid, check.attributes = FALSE))) {
    stop("Existing seed grid does not match requested scenarios and B.")
  }
  cat("Using existing seed grid:", seed_grid_path, "\n")
} else {
  set.seed(20260528)
  seed_grid <- CJ(Scenario = scenarios, Rep = seq_len(B))
  seed_grid[, Seed := sample.int(.Machine$integer.max, .N)]
  fwrite(seed_grid, seed_grid_path)
}

run_seed_row <- function(row) {
  scenario <- seed_grid[row, Scenario]
  rep <- seed_grid[row, Rep]
  seed <- seed_grid[row, Seed]
  rep_path <- file.path(
    out_dir,
    "replicates",
    paste0("scenario=", scenario, "_n=", n, "_rep=", rep, ".rds")
  )
  if (file.exists(rep_path)) {
    cat("Skipping existing", basename(rep_path), "\n")
    return(data.table(Row = row, Scenario = scenario, Rep = rep, Status = "skipped"))
  }

  cat(
    sprintf("[%s] scenario=%s n=%s rep=%s seed=%s\n",
            format(Sys.time(), "%H:%M:%S"), scenario, n, rep, seed)
  )
  result <- run_one_replicate(
    scenario = scenario,
    n = n,
    rep = rep,
    seed = seed,
    target_times = target_times,
    target_events = target_events,
    cv_v = cv_v,
    max_update_iter = max_update_iter,
    include_survsl = include_survsl,
    include_hal_sl = include_hal_sl,
    include_stabilized = include_stabilized,
    stabilized_min_nuisance = stabilized_min_nuisance,
    stabilized_one_step_eps = stabilized_one_step_eps,
    stabilized_stop_rule = stabilized_stop_rule,
    stabilized_stop_abs_tol = stabilized_stop_abs_tol,
    include_survtmle = include_survtmle,
    survtmle_months = survtmle_months
  )
  saveRDS(result, rep_path)
  gc()
  data.table(Row = row, Scenario = scenario, Rep = rep, Status = "done")
}

rows <- seq_len(nrow(seed_grid))
if (cores == 1L) {
  row_status <- rbindlist(lapply(rows, run_seed_row), fill = TRUE)
} else {
  row_status <- rbindlist(parallel::mclapply(
    rows,
    function(row) {
      tryCatch(
        run_seed_row(row),
        error = function(e) {
          data.table(Row = row,
                     Scenario = seed_grid[row, Scenario],
                     Rep = seed_grid[row, Rep],
                     Status = "error",
                     Error = conditionMessage(e))
        }
      )
    },
    mc.cores = cores,
    mc.preschedule = FALSE
  ), fill = TRUE)
}
fwrite(row_status, file.path(out_dir, "row_status.csv"))

summary <- summarize_simulation(output_dir = file.path(out_dir, "replicates"),
                                truth_file = file.path(out_dir, "truth_all.rds"))
saveRDS(summary, file.path(out_dir, "summary.rds"))
fwrite(summary$estimates, file.path(out_dir, "estimates_long.csv"))
fwrite(summary$diagnostics, file.path(out_dir, "diagnostics_long.csv"))
fwrite(summary$metrics, file.path(out_dir, "metrics.csv"))
fwrite(summary$diagnostic_metrics, file.path(out_dir, "diagnostic_metrics.csv"))

cat("\nDiagnostic metrics:\n")
print(summary$diagnostic_metrics)

cat("\nTop-line absolute-risk metrics:\n")
print(summary$metrics[
  Estimand == "Abs Risk",
  .(Scenario, N, EstimatorID, Estimator, Event, Time, Bias, RMSE, Coverage)
][order(Scenario, EstimatorID, Event, Time)][1:min(.N, 30)])

cat("\nDone.\n")
