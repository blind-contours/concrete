suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/sim-data/referee-sims/referee_sim_functions.R")

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    out = "scripts/sim-data/referee-sims/output/failed_seed_update_sensitivity",
    cases = "representative",
    configs = "all",
    cores = 1L
  )
  for (arg in args) {
    if (grepl("^--out=", arg)) out$out <- sub("^--out=", "", arg)
    if (grepl("^--cases=", arg)) out$cases <- sub("^--cases=", "", arg)
    if (grepl("^--configs=", arg)) out$configs <- sub("^--configs=", "", arg)
    if (grepl("^--cores=", arg)) out$cores <- as.integer(sub("^--cores=", "", arg))
  }
  out
}

failed_seed_grid <- function(cases = c("representative", "all")) {
  cases <- match.arg(cases)
  grid <- data.table(
    CaseID = c(
      "nonph_rep4",
      "positivity_rep15",
      "positivity_rep18",
      "rare_early_rep8",
      "rare_early_rep10",
      "rare_early_rep11",
      "rare_early_rep13",
      "rare_early_rep17"
    ),
    Scenario = c(
      "nonph",
      "positivity",
      "positivity",
      "rare_early",
      "rare_early",
      "rare_early",
      "rare_early",
      "rare_early"
    ),
    Rep = c(4L, 15L, 18L, 8L, 10L, 11L, 13L, 17L),
    Seed = c(
      120884616L,
      331531145L,
      1559150800L,
      2085160275L,
      51175594L,
      25777032L,
      1186231484L,
      669013693L
    )
  )
  if (identical(cases, "representative")) {
    grid <- grid[CaseID %in% c("nonph_rep4", "positivity_rep15", "rare_early_rep8")]
  }
  grid[]
}

config_grid <- function(configs = c("all", "best")) {
  configs <- match.arg(configs)
  grid <- data.table(
    ConfigID = c(
      "standard_rich_e0.10_iter200",
      "adaptive_rich_e0.10_iter200",
      "adaptive_rich_e0.01_iter300",
      "adaptive_rich_e0.01_min0.05_iter300"
    ),
    UpdateMethod = c("standard", "adaptive", "adaptive", "adaptive"),
    ModelLibrary = c("rich", "rich", "rich", "rich"),
    MaxUpdateIter = c(200L, 200L, 300L, 300L),
    OneStepEps = c(0.1, 0.1, 0.01, 0.01),
    MinNuisance = c(NA_real_, NA_real_, NA_real_, 0.05)
  )
  if (identical(configs, "best")) {
    grid <- grid[ConfigID == "adaptive_rich_e0.01_min0.05_iter300"]
  }
  grid[]
}

extract_component_state <- function(fit) {
  Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL
  stop_rule <- attr(fit, "EICStopRule")
  stop_abs_tol <- attr(fit, "EICStopAbsTol")
  if (is.null(stop_rule)) stop_rule <- "relative"
  if (is.null(stop_abs_tol)) stop_abs_tol <- 0

  summ <- rbindlist(lapply(seq_along(fit), function(a) {
    cbind(Intervention = names(fit)[a], fit[[a]]$SummEIC)
  }), fill = TRUE)
  summ <- summ[Time %in% attr(fit, "TargetTime") & Event %in% attr(fit, "TargetEvent")]
  stop_dt <- data.table::copy(summ)
  data.table::setnames(stop_dt, "Intervention", "Trt")
  stop_dt <- concrete:::makeOneStepStop(stop_dt, stop_rule, stop_abs_tol)
  data.table::setnames(stop_dt, "Trt", "Intervention")
  summ <- merge(
    summ,
    stop_dt[, .(
      Intervention, Time, Event, Ratio = ratio, RelativeRatio, AbsoluteRatio,
      AbsPnEIC, StopCriteria, StopRule, StopAbsTol, StopCheck = check
    )],
    by = c("Intervention", "Time", "Event")
  )
  setorder(summ, -Ratio)
  summ[]
}

extract_nuisance_denoms <- function(fit) {
  rbindlist(lapply(seq_along(fit), function(a) {
    nw <- fit[[a]]$NuisanceWeight
    denom <- as.numeric(1 / nw)
    data.table(
      Intervention = names(fit)[a],
      MinDenom = min(denom, na.rm = TRUE),
      Q01Denom = as.numeric(stats::quantile(denom, probs = 0.01, na.rm = TRUE)),
      MedianDenom = stats::median(denom, na.rm = TRUE),
      MaxWeight = max(as.numeric(nw), na.rm = TRUE)
    )
  }), fill = TRUE)
}

run_failed_config <- function(row) {
  row <- as.list(row)
  target_times <- c(180, 365, 730, 1095, 1460)
  target_events <- c(1, 2)

  dat <- simulate_referee_data(
    n = 500,
    scenario = row$Scenario,
    seed = row$Seed,
    censoring = TRUE,
    max_time = 1826,
    dt = 14
  )

  run <- capture_run({
    model <- make_concrete_model(
      library = row$ModelLibrary,
      include_events = sort(unique(dat$EVENT)),
      ps_library = c("SL.glm", "SL.glmnet")
    )

    format_args <- list(
      DataTable = data.table::copy(dat),
      EventTime = "TIME",
      EventType = "EVENT",
      Treatment = "ARM",
      ID = "id",
      Intervention = concrete::makeITT(),
      TargetTime = target_times,
      TargetEvent = target_events,
      CVArg = list(V = 2),
      Model = model,
      MaxUpdateIter = row$MaxUpdateIter,
      OneStepEps = row$OneStepEps,
      Verbose = FALSE,
      GComp = TRUE,
      ReturnModels = FALSE,
      RenameCovs = FALSE,
      UpdateMethod = row$UpdateMethod
    )
    if (!is.na(row$MinNuisance)) {
      format_args$MinNuisance <- row$MinNuisance
    }

    args <- do.call(concrete::formatArguments, format_args)
    fit <- concrete::doConcrete(args)
    fit
  })

  base <- data.table(
    CaseID = row$CaseID,
    Scenario = row$Scenario,
    Rep = row$Rep,
    Seed = row$Seed,
    ConfigID = row$ConfigID,
    UpdateMethod = row$UpdateMethod,
    ModelLibrary = row$ModelLibrary,
    MaxUpdateIter = row$MaxUpdateIter,
    OneStepEps = row$OneStepEps,
    MinNuisance = row$MinNuisance
  )

  if (inherits(run$value, "error")) {
    diagnostics <- cbind(
      base,
      data.table(
        Converged = FALSE,
        Step = NA_real_,
        FinalNormPnEIC = NA_real_,
        MaxRatio = NA_real_,
        MaxRelativeRatio = NA_real_,
        MaxAbsPnEIC = NA_real_,
        FailingComponents = NA_integer_,
        RelativeFailingComponents = NA_integer_,
        MinNuisanceDenom = NA_real_,
        EICStopRule = NA_character_,
        EICStopAbsTol = NA_real_,
        RuntimeSec = run$elapsed_sec,
        Warnings = paste(run$warnings, collapse = " | "),
        Error = conditionMessage(run$value)
      )
    )
    return(list(
      diagnostics = diagnostics,
      components = data.table(),
      trace = data.table(),
      nuisance = data.table()
    ))
  }

  fit <- run$value
  conv <- attr(fit, "TmleConverged")
  components <- extract_component_state(fit)
  nuisance <- extract_nuisance_denoms(fit)
  trace <- data.table::copy(attr(fit, "TmleUpdateTrace"))
  if (is.null(trace)) trace <- data.table()

  diagnostics <- cbind(
    base,
      data.table(
        Converged = isTRUE(conv$converged),
        Step = suppressWarnings(as.numeric(conv$step)),
        FinalNormPnEIC = tail(attr(fit, "NormPnEICs"), 1),
        MaxRatio = max(components$Ratio, na.rm = TRUE),
        MaxRelativeRatio = max(components$RelativeRatio, na.rm = TRUE),
        MaxAbsPnEIC = max(components$AbsPnEIC, na.rm = TRUE),
        FailingComponents = sum(!components$StopCheck, na.rm = TRUE),
        RelativeFailingComponents = sum(components$RelativeRatio > 1, na.rm = TRUE),
        MinNuisanceDenom = min(nuisance$MinDenom, na.rm = TRUE),
        EICStopRule = attr(fit, "EICStopRule"),
        EICStopAbsTol = attr(fit, "EICStopAbsTol"),
        RuntimeSec = run$elapsed_sec,
      Warnings = paste(run$warnings, collapse = " | "),
      Error = NA_character_
    )
  )

  for (dt in list(components, trace, nuisance)) {
    dt[, `:=`(
      CaseID = row$CaseID,
      Scenario = row$Scenario,
      Rep = row$Rep,
      Seed = row$Seed,
      ConfigID = row$ConfigID
    )]
  }

  list(
    diagnostics = diagnostics,
    components = components,
    trace = trace,
    nuisance = nuisance
  )
}

summarize_trace <- function(trace) {
  if (!nrow(trace)) return(data.table())
  trace[, .(
    TraceRows = .N,
    AcceptedSteps = sum(Status == "accepted", na.rm = TRUE),
    RejectedSteps = sum(grepl("^rejected", Status), na.rm = TRUE),
    NoAcceptedSteps = sum(Status == "no_accepted_step", na.rm = TRUE),
    InitialNorm = NormAfter[Status == "initial"][1],
    FinalNorm = tail(na.omit(NormAfter), 1),
    BestNorm = min(NormAfter, na.rm = TRUE),
    FinalMaxRatio = tail(na.omit(MaxRatio), 1),
    FinalMaxRelativeRatio = tail(na.omit(MaxRelativeRatio), 1),
    FinalMaxAbsPnEIC = tail(na.omit(MaxAbsPnEIC), 1),
    BestMaxRatio = min(MaxRatio, na.rm = TRUE),
    MinAcceptedAlpha = suppressWarnings(min(Alpha[Status == "accepted"], na.rm = TRUE)),
    MedianAcceptedAlpha = suppressWarnings(stats::median(Alpha[Status == "accepted"], na.rm = TRUE)),
    FinalWorstTrt = tail(na.omit(WorstTrt), 1),
    FinalWorstTime = tail(na.omit(WorstTime), 1),
    FinalWorstEvent = tail(na.omit(WorstEvent), 1)
  ), by = .(CaseID, Scenario, Rep, Seed, ConfigID)]
}

main <- function() {
  cli <- parse_cli()
  dir.create(cli$out, recursive = TRUE, showWarnings = FALSE)

  cases <- failed_seed_grid(cli$cases)
  configs <- config_grid(cli$configs)
  cases[, CrossJoinID := 1L]
  configs[, CrossJoinID := 1L]
  jobs <- merge(cases, configs, by = "CrossJoinID", allow.cartesian = TRUE)
  jobs[, CrossJoinID := NULL]
  cases[, CrossJoinID := NULL]
  configs[, CrossJoinID := NULL]
  setorder(jobs, Scenario, Rep, ConfigID)
  fwrite(jobs, file.path(cli$out, "job_grid.csv"))

  cat("Failed-seed update sensitivity\n")
  cat("Output:", cli$out, "\n")
  cat("Cases:", cli$cases, " (", nrow(cases), " seeds )\n", sep = "")
  cat("Configs:", cli$configs, " (", nrow(configs), " )\n", sep = "")
  cat("Jobs:", nrow(jobs), "\n")
  cat("Cores:", cli$cores, "\n\n")

  run_job <- function(i) {
    job <- jobs[i]
    cat(sprintf("[%s] %s / %s\n", format(Sys.time(), "%H:%M:%S"), job$CaseID, job$ConfigID))
    run_failed_config(job)
  }

  if (cli$cores > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = cli$cores, mc.preschedule = FALSE)
  } else {
    results <- lapply(seq_len(nrow(jobs)), run_job)
  }

  diagnostics <- rbindlist(lapply(results, `[[`, "diagnostics"), fill = TRUE)
  components <- rbindlist(lapply(results, `[[`, "components"), fill = TRUE)
  trace <- rbindlist(lapply(results, `[[`, "trace"), fill = TRUE)
  nuisance <- rbindlist(lapply(results, `[[`, "nuisance"), fill = TRUE)
  trace_summary <- summarize_trace(trace)

  fwrite(diagnostics, file.path(cli$out, "diagnostics.csv"))
  fwrite(components, file.path(cli$out, "final_components.csv"))
  fwrite(trace, file.path(cli$out, "trace_long.csv"))
  fwrite(trace_summary, file.path(cli$out, "trace_summary.csv"))
  fwrite(nuisance, file.path(cli$out, "nuisance_denoms.csv"))
  saveRDS(
    list(
      jobs = jobs,
      diagnostics = diagnostics,
      components = components,
      trace = trace,
      trace_summary = trace_summary,
      nuisance = nuisance
    ),
    file.path(cli$out, "failed_seed_update_sensitivity.rds")
  )

  cat("\nDiagnostics:\n")
  print(diagnostics[, .(
    CaseID, ConfigID, Converged, Step, FinalNormPnEIC,
    MaxRatio, MaxRelativeRatio, MaxAbsPnEIC, FailingComponents,
    RelativeFailingComponents, MinNuisanceDenom, RuntimeSec
  )])
  cat("\nTrace summary:\n")
  print(trace_summary)
  cat("\nDone.\n")
}

main()
