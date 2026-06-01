suppressPackageStartupMessages({
  library(data.table)
})

parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    out = "scripts/sim-data/referee-sims/output/alt_convergence_methods",
    cases = "representative",
    configs = "candidate",
    cores = 1L,
    n = 500L,
    max_iter = 300L,
    cv_v = 2L,
    model_library = "rich",
    target_times = "180,365,730,1095,1460",
    target_events = "1,2",
    ps_library = "SL.glm,SL.glmnet"
  )
  for (arg in args) {
    if (grepl("^--out=", arg)) out$out <- sub("^--out=", "", arg)
    if (grepl("^--cases=", arg)) out$cases <- sub("^--cases=", "", arg)
    if (grepl("^--configs=", arg)) out$configs <- sub("^--configs=", "", arg)
    if (grepl("^--cores=", arg)) out$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--n=", arg)) out$n <- as.integer(sub("^--n=", "", arg))
    if (grepl("^--max-iter=", arg)) out$max_iter <- as.integer(sub("^--max-iter=", "", arg))
    if (grepl("^--cv-v=", arg)) out$cv_v <- as.integer(sub("^--cv-v=", "", arg))
    if (grepl("^--model-library=", arg)) out$model_library <- sub("^--model-library=", "", arg)
    if (grepl("^--target-times=", arg)) out$target_times <- sub("^--target-times=", "", arg)
    if (grepl("^--target-events=", arg)) out$target_events <- sub("^--target-events=", "", arg)
    if (grepl("^--ps-library=", arg)) out$ps_library <- sub("^--ps-library=", "", arg)
  }
  out
}

split_int_arg <- function(x) {
  as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}

split_char_arg <- function(x) {
  out <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

load_concrete_for_sim <- function(lib_dir) {
  local_version <- unname(read.dcf("DESCRIPTION")[1, "Version"])

  if (requireNamespace("pkgload", quietly = TRUE) &&
      requireNamespace("pkgbuild", quietly = TRUE)) {
    pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
    return(invisible(TRUE))
  }

  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  lib_dir <- normalizePath(lib_dir, mustWork = TRUE)
  .libPaths(c(lib_dir, .libPaths()))

  installed_version <- tryCatch(
    as.character(utils::packageVersion("concrete")),
    error = function(e) NA_character_
  )
  if (identical(installed_version, local_version)) {
    suppressPackageStartupMessages(library(concrete))
    return(invisible(TRUE))
  }

  cat(
    "Installing concrete ", local_version,
    " into simulation library: ", lib_dir, "\n",
    sep = ""
  )

  pkg_dir <- normalizePath(".", mustWork = TRUE)
  build_dir <- tempfile("concrete_pkg_build_")
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(build_dir)
  on.exit(setwd(old_wd), add = TRUE)
  build_output <- system2(
    file.path(R.home("bin"), "R"),
    c(
      "CMD", "build",
      "--no-build-vignettes", "--no-manual",
      pkg_dir
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  build_status <- attr(build_output, "status")
  if (!is.null(build_status) && build_status != 0) {
    stop(
      "Failed to build local concrete package for simulations:\n",
      paste(build_output, collapse = "\n")
    )
  }
  tarballs <- list.files(build_dir, pattern = "^concrete_.*[.]tar[.]gz$",
                         full.names = TRUE)
  if (!length(tarballs)) {
    stop(
      "Package build did not produce a concrete source tarball:\n",
      paste(build_output, collapse = "\n")
    )
  }
  tarball <- tarballs[which.max(file.info(tarballs)$mtime)]

  install_output <- system2(
    file.path(R.home("bin"), "R"),
    c(
      "CMD", "INSTALL",
      "--no-docs", "--no-help", "--no-html",
      paste0("--library=", lib_dir),
      tarball
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  install_status <- attr(install_output, "status")
  if (!is.null(install_status) && install_status != 0) {
    stop(
      "Failed to install local concrete package for simulations:\n",
      paste(install_output, collapse = "\n")
    )
  }

  suppressPackageStartupMessages(library(concrete))
  loaded_version <- as.character(utils::packageVersion("concrete"))
  if (!identical(loaded_version, local_version)) {
    stop(
      "Loaded concrete version ", loaded_version,
      " but expected local version ", local_version, "."
    )
  }
  invisible(TRUE)
}

enable_internal_coordinate_update <- function() {
  ns <- asNamespace("concrete")
  patched <- function(UpdateMethod) {
    choices <- c("standard", "adaptive", "coordinated")
    disabled <- c("jacobi", "rootsolve")

    if (is.null(UpdateMethod) || length(UpdateMethod) == 0) {
      return("standard")
    }
    if (length(UpdateMethod) > 1) {
      UpdateMethod <- UpdateMethod[1]
    }
    if (!is.character(UpdateMethod) || is.na(UpdateMethod)) {
      stop("UpdateMethod must be one of: ", paste(choices, collapse = ", "))
    }
    UpdateMethod <- trimws(UpdateMethod)
    if (identical(UpdateMethod, "accelerated")) {
      message("UpdateMethod = 'accelerated' has been renamed to 'adaptive'.")
      UpdateMethod <- "adaptive"
    }
    UpdateMethod <- tolower(UpdateMethod)
    if (UpdateMethod %in% disabled) {
      stop("UpdateMethod = '", UpdateMethod, "' is disabled for this experiment.")
    }
    if (!(UpdateMethod %in% choices)) {
      stop("UpdateMethod must be one of: ", paste(choices, collapse = ", "))
    }
    UpdateMethod
  }

  if (bindingIsLocked("getUpdateMethod", ns)) {
    unlockBinding("getUpdateMethod", ns)
    on.exit(lockBinding("getUpdateMethod", ns), add = TRUE)
  }
  assign("getUpdateMethod", patched, envir = ns)
  invisible(TRUE)
}

hard_seed_grid <- function(cases = c("representative", "all", "rare", "smoke")) {
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
  if (identical(cases, "smoke")) {
    grid <- grid[CaseID == "rare_early_rep8"]
  } else if (identical(cases, "representative")) {
    grid <- grid[CaseID %in% c("nonph_rep4", "positivity_rep15", "rare_early_rep8")]
  } else if (identical(cases, "rare")) {
    grid <- grid[Scenario == "rare_early"]
  }
  grid[]
}

method_config_grid <- function(configs = c("candidate", "all", "smoke", "finalist", "primary")) {
  configs <- match.arg(configs)
  grid <- data.table(
    ConfigID = c(
      "relative_adaptive_e0.01_min0.05",
      "hybrid_fixed_5e-4",
      "hybrid_fixed_1e-3",
      "hybrid_nsqrt_0.02",
      "hybrid_nsqrt_0.05",
      "hybrid_nsqrt_0.10",
      "absolute_nsqrt_0.02",
      "coordinated_relative_e0.01",
      "coordinated_hybrid_fixed_1e-3"
    ),
    UpdateMethod = c(
      "adaptive",
      "adaptive",
      "adaptive",
      "adaptive",
      "adaptive",
      "adaptive",
      "adaptive",
      "coordinated",
      "coordinated"
    ),
    EICStopRule = c(
      "relative",
      "hybrid",
      "hybrid",
      "hybrid",
      "hybrid",
      "hybrid",
      "absolute",
      "relative",
      "hybrid"
    ),
    TolScale = c(
      "fixed",
      "fixed",
      "fixed",
      "nsqrt",
      "nsqrt",
      "nsqrt",
      "nsqrt",
      "fixed",
      "fixed"
    ),
    TolConst = c(0, 5e-4, 1e-3, 0.02, 0.05, 0.10, 0.02, 0, 1e-3),
    OneStepEps = c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
    MinNuisance = 0.05,
    ModelLibrary = "rich",
    Experimental = c(rep(FALSE, 7), TRUE, TRUE)
  )
  if (identical(configs, "candidate")) {
    grid <- grid[ConfigID %in% c(
      "relative_adaptive_e0.01_min0.05",
      "hybrid_fixed_1e-3",
      "hybrid_nsqrt_0.02",
      "hybrid_nsqrt_0.05",
      "coordinated_hybrid_fixed_1e-3"
    )]
  } else if (identical(configs, "finalist")) {
    grid <- grid[ConfigID %in% c(
      "relative_adaptive_e0.01_min0.05",
      "hybrid_nsqrt_0.05",
      "absolute_nsqrt_0.02",
      "coordinated_hybrid_fixed_1e-3"
    )]
  } else if (identical(configs, "primary")) {
    grid <- grid[ConfigID %in% c(
      "relative_adaptive_e0.01_min0.05",
      "absolute_nsqrt_0.02"
    )]
  } else if (identical(configs, "smoke")) {
    grid <- grid[ConfigID %in% c(
      "relative_adaptive_e0.01_min0.05",
      "hybrid_fixed_1e-3",
      "hybrid_nsqrt_0.02"
    )]
  }
  grid[]
}

resolve_tol <- function(scale, constant, n) {
  if (identical(scale, "fixed")) return(constant)
  if (identical(scale, "nsqrt")) return(constant / sqrt(n))
  if (identical(scale, "nlogn")) return(constant / (sqrt(n) * log(n)))
  stop("Unknown tolerance scale: ", scale)
}

extract_components <- function(fit) {
  out <- concrete::getTmleDiagnostics(fit, type = "components")
  data.table::setnames(out, "Intervention", "Trt", skip_absent = TRUE)
  out[]
}

extract_trace <- function(fit) {
  trace <- concrete::getTmleDiagnostics(fit, type = "trace")
  if (is.null(trace)) trace <- data.table()
  trace[]
}

extract_nuisance_denoms <- function(fit) {
  rbindlist(lapply(seq_along(fit), function(a) {
    nw <- fit[[a]]$NuisanceWeight
    denom <- as.numeric(1 / nw)
    data.table(
      Trt = names(fit)[a],
      MinDenom = min(denom, na.rm = TRUE),
      Q01Denom = as.numeric(stats::quantile(denom, probs = 0.01, na.rm = TRUE)),
      MedianDenom = stats::median(denom, na.rm = TRUE),
      MaxWeight = max(as.numeric(nw), na.rm = TRUE)
    )
  }), fill = TRUE)
}

run_config <- function(row) {
  row <- as.list(row)
  target_times <- split_int_arg(row$TargetTimes)
  target_events <- split_int_arg(row$TargetEvents)
  abs_tol <- resolve_tol(row$TolScale, row$TolConst, row$N)

  dat <- simulate_referee_data(
    n = row$N,
    scenario = row$Scenario,
    seed = row$Seed,
    censoring = TRUE,
    max_time = 1826,
    dt = 14
  )

  base <- data.table(
    CaseID = row$CaseID,
    Scenario = row$Scenario,
    Rep = row$Rep,
    Seed = row$Seed,
    N = row$N,
    ConfigID = row$ConfigID,
    UpdateMethod = row$UpdateMethod,
    ModelLibrary = row$ModelLibrary,
    OneStepEps = row$OneStepEps,
    MinNuisance = row$MinNuisance,
    EICStopRule = row$EICStopRule,
    TolScale = row$TolScale,
    TolConst = row$TolConst,
    EICStopAbsTol = abs_tol,
    Experimental = row$Experimental,
    TargetTimes = row$TargetTimes,
    TargetEvents = row$TargetEvents,
    PSLibrary = row$PSLibrary
  )

  run <- capture_run({
    model <- make_concrete_model(
      library = row$ModelLibrary,
      include_events = sort(unique(dat$EVENT)),
      ps_library = split_char_arg(row$PSLibrary)
    )
    args <- concrete::formatArguments(
      DataTable = data.table::copy(dat),
      EventTime = "TIME",
      EventType = "EVENT",
      Treatment = "ARM",
      ID = "id",
      Intervention = concrete::makeITT(),
      TargetTime = target_times,
      TargetEvent = target_events,
      CVArg = list(V = row$CVV),
      Model = model,
      MaxUpdateIter = row$MaxUpdateIter,
      OneStepEps = row$OneStepEps,
      MinNuisance = row$MinNuisance,
      Verbose = FALSE,
      GComp = TRUE,
      ReturnModels = FALSE,
      RenameCovs = FALSE,
      UpdateMethod = row$UpdateMethod,
      EICStopRule = row$EICStopRule,
      EICStopAbsTol = abs_tol
    )
    fit <- concrete::doConcrete(args)
    output <- as.data.table(concrete::getOutput(
      fit,
      Estimand = c("Risk", "RD", "RR"),
      GComp = TRUE,
      Simultaneous = FALSE
    ))
    list(fit = fit, output = output)
  })

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
        RuntimeSec = run$elapsed_sec,
        Warnings = paste(run$warnings, collapse = " | "),
        Error = conditionMessage(run$value)
      )
    )
    return(list(
      diagnostics = diagnostics,
      components = data.table(),
      trace = data.table(),
      nuisance = data.table(),
      estimates = data.table()
    ))
  }

  fit <- run$value$fit
  conv <- attr(fit, "TmleConverged")
  components <- extract_components(fit)
  trace <- extract_trace(fit)
  nuisance <- extract_nuisance_denoms(fit)
  estimates <- data.table::copy(run$value$output)

  diagnostics <- cbind(
    base,
    data.table(
      Converged = isTRUE(conv$converged),
      Step = suppressWarnings(as.numeric(conv$step)),
      FinalNormPnEIC = tail(attr(fit, "NormPnEICs"), 1),
      MaxRatio = max(components$ratio, na.rm = TRUE),
      MaxRelativeRatio = max(components$RelativeRatio, na.rm = TRUE),
      MaxAbsPnEIC = max(components$AbsPnEIC, na.rm = TRUE),
      FailingComponents = sum(!components$check, na.rm = TRUE),
      RelativeFailingComponents = sum(components$RelativeRatio > 1, na.rm = TRUE),
      MinNuisanceDenom = min(nuisance$MinDenom, na.rm = TRUE),
      RuntimeSec = run$elapsed_sec,
      Warnings = paste(run$warnings, collapse = " | "),
      Error = NA_character_
    )
  )

  for (dt in list(components, trace, nuisance, estimates)) {
    if (!nrow(dt)) next
    for (nm in names(base)) {
      dt[, (nm) := base[[nm]]]
    }
  }

  list(
    diagnostics = diagnostics,
    components = components,
    trace = trace,
    nuisance = nuisance,
    estimates = estimates
  )
}

summarize_methods <- function(diagnostics, estimates) {
  method_summary <- diagnostics[, .(
    Jobs = .N,
    Errors = sum(!is.na(Error)),
    ConvRate = mean(Converged, na.rm = TRUE),
    MedianStep = stats::median(Step, na.rm = TRUE),
    MedianMaxRatio = stats::median(MaxRatio, na.rm = TRUE),
    MedianMaxRelativeRatio = stats::median(MaxRelativeRatio, na.rm = TRUE),
    MedianMaxAbsPnEIC = stats::median(MaxAbsPnEIC, na.rm = TRUE),
    MeanFailingComponents = mean(FailingComponents, na.rm = TRUE),
    MedianRuntimeSec = stats::median(RuntimeSec, na.rm = TRUE)
  ), by = .(ConfigID, UpdateMethod, EICStopRule, TolScale, TolConst, EICStopAbsTol, Experimental)]
  setorder(method_summary, -ConvRate, MedianMaxRatio, MedianRuntimeSec)

  estimate_delta <- data.table()
  if (nrow(estimates)) {
    tmle <- estimates[Estimator == "tmle"]
    ref <- tmle[ConfigID == "relative_adaptive_e0.01_min0.05",
                .(CaseID, Scenario, Time, Event, Estimand, Intervention,
                  RefPtEst = `Pt Est`)]
    estimate_delta <- merge(
      tmle,
      ref,
      by = c("CaseID", "Scenario", "Time", "Event", "Estimand", "Intervention"),
      all.x = TRUE
    )
    estimate_delta[, AbsDiffFromRelative := abs(`Pt Est` - RefPtEst)]
    estimate_delta <- estimate_delta[, .(
      MaxAbsDiffFromRelative = max(AbsDiffFromRelative, na.rm = TRUE),
      MedianAbsDiffFromRelative = stats::median(AbsDiffFromRelative, na.rm = TRUE)
    ), by = .(ConfigID)]
  }

  list(method_summary = method_summary, estimate_delta = estimate_delta)
}

main <- function() {
  cli <- parse_cli()
  dir.create(cli$out, recursive = TRUE, showWarnings = FALSE)
  load_concrete_for_sim(file.path(cli$out, ".r-lib"))
  source("scripts/sim-data/referee-sims/referee_sim_functions.R")
  enable_internal_coordinate_update()

  cases <- hard_seed_grid(cli$cases)
  configs <- method_config_grid(cli$configs)
  cases[, CrossJoinID := 1L]
  configs[, CrossJoinID := 1L]
  jobs <- merge(cases, configs, by = "CrossJoinID", allow.cartesian = TRUE)
  jobs[, `:=`(
    CrossJoinID = NULL,
    N = cli$n,
    MaxUpdateIter = cli$max_iter,
    CVV = cli$cv_v,
    ModelLibrary = cli$model_library,
    TargetTimes = cli$target_times,
    TargetEvents = cli$target_events,
    PSLibrary = cli$ps_library
  )]
  setorder(jobs, Rep, ConfigID)
  fwrite(jobs, file.path(cli$out, "job_grid.csv"))

  cat("Alternative convergence method comparison\n")
  cat("Output:", cli$out, "\n")
  cat("Cases:", cli$cases, " (", nrow(cases), " seeds )\n", sep = "")
  cat("Configs:", cli$configs, " (", nrow(configs), " )\n", sep = "")
  cat("N:", cli$n, " MaxUpdateIter:", cli$max_iter, " CV V:", cli$cv_v, "\n")
  cat("Model library:", cli$model_library, "\n")
  cat("Target times:", cli$target_times, " Target events:", cli$target_events, "\n")
  cat("PS library:", cli$ps_library, "\n")
  cat("Jobs:", nrow(jobs), " Cores:", cli$cores, "\n\n")

  run_job <- function(i) {
    job <- jobs[i]
    cat(sprintf("[%s] %s / %s\n", format(Sys.time(), "%H:%M:%S"), job$CaseID, job$ConfigID))
    run_config(job)
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
  estimates <- rbindlist(lapply(results, `[[`, "estimates"), fill = TRUE)
  summaries <- summarize_methods(diagnostics, estimates)

  fwrite(diagnostics, file.path(cli$out, "diagnostics.csv"))
  fwrite(components, file.path(cli$out, "final_components.csv"))
  fwrite(trace, file.path(cli$out, "trace_long.csv"))
  fwrite(nuisance, file.path(cli$out, "nuisance_denoms.csv"))
  fwrite(estimates, file.path(cli$out, "estimates_long.csv"))
  fwrite(summaries$method_summary, file.path(cli$out, "method_summary.csv"))
  fwrite(summaries$estimate_delta, file.path(cli$out, "estimate_delta_vs_relative.csv"))
  saveRDS(
    list(
      jobs = jobs,
      diagnostics = diagnostics,
      components = components,
      trace = trace,
      nuisance = nuisance,
      estimates = estimates,
      method_summary = summaries$method_summary,
      estimate_delta = summaries$estimate_delta
    ),
    file.path(cli$out, "alt_convergence_methods.rds")
  )

  cat("\nMethod summary:\n")
  print(summaries$method_summary)
  cat("\nEstimate delta vs relative adaptive baseline:\n")
  print(summaries$estimate_delta)
  cat("\nDiagnostics:\n")
  print(diagnostics[, .(
    CaseID, ConfigID, Converged, Step, MaxRatio, MaxRelativeRatio,
    MaxAbsPnEIC, FailingComponents, RuntimeSec, Error
  )])
  cat("\nDone.\n")
}

main()
