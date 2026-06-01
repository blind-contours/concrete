source("scripts/sim-data/referee-sims/referee_sim_functions.R")

args <- parse_cli_args(list(
  out = "scripts/sim-data/referee-sims/output/pilot"
))

out_dir <- args$out
summary <- summarize_simulation(
  output_dir = file.path(out_dir, "replicates"),
  truth_file = file.path(out_dir, "truth_all.rds")
)

saveRDS(summary, file.path(out_dir, "summary.rds"))
fwrite(summary$estimates, file.path(out_dir, "estimates_long.csv"))
fwrite(summary$diagnostics, file.path(out_dir, "diagnostics_long.csv"))
fwrite(summary$metrics, file.path(out_dir, "metrics.csv"))
fwrite(summary$diagnostic_metrics, file.path(out_dir, "diagnostic_metrics.csv"))

cat("Wrote summary files to", out_dir, "\n\n")
print(summary$diagnostic_metrics)
