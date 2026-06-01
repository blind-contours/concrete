# Trialist smoke test for the concrete package.
#
# After installing concrete, run:
# source(system.file("examples", "trialist-smoke-test.R", package = "concrete"))
#
# To also try optional hazard learners that are installed locally, run:
# Sys.setenv(CONCRETE_RUN_OPTIONAL_LEARNERS = "true")
# source(system.file("examples", "trialist-smoke-test.R", package = "concrete"))

library(concrete)
library(data.table)

set.seed(20260601)

flag_is_true <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

make_trialist_smoke_data <- function(n = 160) {
  trial <- data.table::as.data.table(survival::pbc)
  trial <- trial[!is.na(trt), .(id, time, status, trt, age, sex, albumin, bili)]
  trial[, arm := as.integer(trt == 2)]
  trial[, event := data.table::fifelse(
    status == 2, 1L,
    data.table::fifelse(status == 1, 2L, 0L)
  )]
  trial <- trial[stats::complete.cases(trial)]
  trial[, sex := as.integer(sex == "f")]
  trial <- trial[sample(.N)]
  trial[seq_len(min(.N, n)), .(id, time, event, arm, age, sex, albumin, bili)]
}

make_hazard_model <- function(extra = character()) {
  event_covariates <- "arm + age + albumin"
  cox_formula <- function(j) {
    stats::as.formula(
      paste0("survival::Surv(time, event == ", j, ") ~ ", event_covariates)
    )
  }
  hazard_candidates <- function(j) {
    out <- list(Cox = cox_formula(j))
    for (learner in extra) {
      out[[learner]] <- learner
    }
    out
  }
  list(
    arm = c("SL.mean", "SL.glm"),
    "0" = hazard_candidates(0),
    "1" = hazard_candidates(1),
    "2" = hazard_candidates(2)
  )
}

run_concrete_smoke <- function(label, model, trial, max_update_iter = 25) {
  message("\n--- Running ", label, " ---")
  start <- proc.time()[["elapsed"]]
  tryCatch({
    args <- formatArguments(
      DataTable = trial,
      EventTime = "time",
      EventType = "event",
      Treatment = "arm",
      ID = "id",
      Intervention = makeITT(),
      TargetTime = c(1000, 2000),
      TargetEvent = 1,
      CVArg = list(V = 2),
      Model = model,
      MaxUpdateIter = max_update_iter,
      UpdateMethod = "adaptive",
      EICStopRule = "absolute",
      EICStopAbsTol = 0.02 / sqrt(nrow(trial)),
      Verbose = FALSE,
      ReturnModels = TRUE
    )
    fit <- doConcrete(args)
    output <- getOutput(
      fit,
      Estimand = c("Risk", "RD", "RR"),
      Intervention = c(1, 2),
      GComp = TRUE,
      Simultaneous = FALSE
    )
    diagnostics <- getTmleDiagnostics(fit, type = "components")
    elapsed <- proc.time()[["elapsed"]] - start
    list(
      label = label,
      status = "ok",
      elapsed_sec = elapsed,
      output = output,
      diagnostics = diagnostics,
      convergence = attr(fit, "TmleConverged")
    )
  }, error = function(e) {
    elapsed <- proc.time()[["elapsed"]] - start
    list(
      label = label,
      status = "error",
      elapsed_sec = elapsed,
      error = conditionMessage(e)
    )
  })
}

summarise_smoke_result <- function(result) {
  if (identical(result$status, "error")) {
    return(data.table::data.table(
      analysis = result$label,
      status = result$status,
      elapsed_sec = round(result$elapsed_sec, 1),
      converged = NA,
      step = NA_integer_,
      max_ratio = NA_real_,
      failing_components = NA_integer_,
      message = result$error
    ))
  }
  diagnostics <- result$diagnostics
  data.table::data.table(
    analysis = result$label,
    status = result$status,
    elapsed_sec = round(result$elapsed_sec, 1),
    converged = isTRUE(result$convergence$converged),
    step = result$convergence$step,
    max_ratio = round(max(diagnostics$ratio, na.rm = TRUE), 3),
    failing_components = sum(!diagnostics$check, na.rm = TRUE),
    message = ""
  )
}

trial_smoke_data <- make_trialist_smoke_data()

cat("\nEvent counts by treatment arm:\n")
print(trial_smoke_data[, .N, by = .(arm, event)][order(arm, event)])

smoke_results <- list()
smoke_results[["cox_only"]] <- run_concrete_smoke(
  label = "cox_only",
  model = make_hazard_model(),
  trial = trial_smoke_data
)

run_optional <- flag_is_true(Sys.getenv("CONCRETE_RUN_OPTIONAL_LEARNERS", "false"))

if (run_optional) {
  smoke_results[["additive_hazards"]] <- run_concrete_smoke(
    label = "additive_hazards",
    model = make_hazard_model("aareg"),
    trial = trial_smoke_data
  )

  if (requireNamespace("glmnet", quietly = TRUE)) {
    smoke_results[["coxnet"]] <- run_concrete_smoke(
      label = "coxnet",
      model = make_hazard_model("coxnet"),
      trial = trial_smoke_data
    )
  } else {
    message("Skipping coxnet: package 'glmnet' is not installed.")
  }

  if (requireNamespace("randomForestSRC", quietly = TRUE)) {
    smoke_results[["rsf"]] <- run_concrete_smoke(
      label = "rsf",
      model = make_hazard_model("rsf"),
      trial = trial_smoke_data
    )
  } else {
    message("Skipping rsf: package 'randomForestSRC' is not installed.")
  }

  if (requireNamespace("hal9001", quietly = TRUE)) {
    smoke_results[["hal"]] <- run_concrete_smoke(
      label = "hal",
      model = make_hazard_model("hal"),
      trial = trial_smoke_data
    )
  } else {
    message("Skipping HAL: package 'hal9001' is not installed.")
  }
} else {
  message(
    "\nOptional learner checks skipped. Set ",
    "CONCRETE_RUN_OPTIONAL_LEARNERS=true to run installed optional learners."
  )
}

smoke_summary <- data.table::rbindlist(lapply(smoke_results, summarise_smoke_result), fill = TRUE)

cat("\nSmoke-test summary:\n")
print(smoke_summary)

cat("\nPrimary TMLE output from cox_only:\n")
if (identical(smoke_results[["cox_only"]]$status, "ok")) {
  print(smoke_results[["cox_only"]]$output)
}

if (!identical(smoke_results[["cox_only"]]$status, "ok")) {
  stop("The required cox_only smoke test failed. See smoke_summary$message.")
}

invisible(smoke_results)
