library(data.table)

expit <- function(x) {
  1 / (1 + exp(-x))
}

parse_cli_args <- function(defaults = list()) {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  for (arg in raw) {
    if (!grepl("^--[^=]+=", arg)) next
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    val <- sub("^--[^=]+=", "", arg)
    out[[key]] <- val
  }
  out
}

as_bool <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

split_arg <- function(x) {
  x <- as.character(x)
  x <- x[nzchar(x)]
  unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
}

scenario_catalog <- function() {
  data.table(
    scenario = c("ph_correct", "nonph", "rare_early", "positivity"),
    description = c(
      "Proportional hazards with a Cox-compatible linear predictor.",
      "Non-proportional treatment and covariate effects.",
      "Low early event rate to stress early target-time convergence.",
      "Strong treatment imbalance and informative censoring."
    )
  )
}

simulate_baseline <- function(n) {
  data.table(
    W1 = rbinom(n, 1, 0.45),
    W2 = rnorm(n),
    W3 = rbinom(n, 1, 0.30),
    W4 = rnorm(n),
    W5 = rbinom(n, 1, 0.55)
  )
}

assign_treatment <- function(W, scenario) {
  lp <- -0.10 + 0.55 * W$W1 - 0.45 * W$W2 + 0.35 * W$W3 - 0.25 * W$W5
  if (identical(scenario, "positivity")) {
    lp <- -0.20 + 1.70 * W$W1 - 1.55 * W$W2 + 1.15 * W$W3 - 0.95 * W$W5
  }
  p <- expit(lp)
  p <- pmin(pmax(p, 0.01), 0.99)
  rbinom(nrow(W), 1, p)
}

hazard_matrix <- function(t, W, A, scenario, censoring = TRUE) {
  years <- pmax(t, 1) / 365.25

  base1 <- (0.095 / 365.25) * 1.15 * years^0.15
  base2 <- (0.060 / 365.25) * 1.05 * years^0.05
  basec <- (0.035 / 365.25) * 1.10 * years^0.10

  lp1 <- log(0.72) * A + 0.42 * W$W1 + 0.34 * W$W2 -
    0.22 * W$W3 + 0.20 * W$W4
  lp2 <- log(1.25) * A + 0.32 * W$W3 - 0.26 * W$W2 +
    0.24 * W$W5
  lpc <- log(0.95) * A + 0.25 * W$W1 - 0.20 * W$W4 + 0.18 * W$W5

  if (identical(scenario, "nonph")) {
    switch_term <- expit((t - 730) / 120)
    lp1 <- 0.42 * W$W1 + 0.34 * W$W2 - 0.22 * W$W3 + 0.20 * W$W4 +
      A * (log(0.50) + log(2.60) * switch_term) +
      0.45 * W$W2 * switch_term - 0.30 * W$W1 * switch_term
    lp2 <- 0.32 * W$W3 - 0.26 * W$W2 + 0.24 * W$W5 +
      A * (log(1.05) + log(1.60) * switch_term)
  } else if (identical(scenario, "rare_early")) {
    base1 <- base1 * 0.35
    lp1 <- lp1 + 0.90 * (years > 2.0)
  } else if (identical(scenario, "positivity")) {
    basec <- basec * 1.65
    lpc <- lpc + 0.70 * W$W1 - 0.60 * W$W2 + 0.45 * A * W$W3
  }

  cens <- if (isTRUE(censoring)) basec * exp(lpc) else rep(0, nrow(W))
  event1 <- base1 * exp(lp1)
  event2 <- base2 * exp(lp2)
  cbind(cens = cens, event1 = event1, event2 = event2)
}

simulate_referee_data <- function(n,
                                  scenario = "ph_correct",
                                  seed = NULL,
                                  intervention = NULL,
                                  censoring = TRUE,
                                  max_time = 1826,
                                  dt = 14) {
  if (!is.null(seed)) set.seed(seed)
  W <- simulate_baseline(n)
  A <- if (is.null(intervention)) {
    assign_treatment(W, scenario)
  } else {
    rep_len(as.integer(intervention), n)
  }

  time <- rep(max_time, n)
  event <- rep(0L, n)
  alive <- rep(TRUE, n)
  starts <- seq(0, max_time - dt, by = dt)
  event_codes <- c(0L, 1L, 2L)

  for (t0 in starts) {
    idx <- which(alive)
    if (!length(idx)) break

    haz <- hazard_matrix(
      t = t0 + dt / 2,
      W = W[idx],
      A = A[idx],
      scenario = scenario,
      censoring = censoring
    )
    total_haz <- rowSums(haz)
    p_event <- 1 - exp(-total_haz * dt)
    hit <- runif(length(idx)) < p_event
    if (!any(hit)) next

    hit_idx <- idx[hit]
    probs <- haz[hit, , drop = FALSE] / total_haz[hit]
    probs[!is.finite(probs)] <- 0
    cum_probs <- t(apply(probs, 1, cumsum))
    draws <- runif(length(hit_idx))
    event_col <- max.col(draws <= cum_probs, ties.method = "first")

    time[hit_idx] <- pmin(t0 + runif(length(hit_idx), 0, dt), max_time)
    event[hit_idx] <- event_codes[event_col]
    alive[hit_idx] <- FALSE
  }

  out <- data.table(
    id = seq_len(n),
    TIME = time,
    EVENT = event,
    ARM = A
  )
  cbind(out, W)
}

estimate_truth <- function(scenario,
                           target_times,
                           target_events = c(1, 2),
                           n_truth = 100000,
                           chunk_size = 50000,
                           seed = 20260528,
                           max_time = 1826,
                           dt = 14) {
  target_times <- sort(unique(target_times))
  target_events <- sort(unique(target_events))
  truth <- data.table()

  for (a in c(1L, 0L)) {
    counts <- data.table(expand.grid(
      Time = target_times,
      Event = target_events
    ))
    counts[, count := 0]

    done <- 0L
    chunk <- 0L
    while (done < n_truth) {
      chunk <- chunk + 1L
      m <- min(chunk_size, n_truth - done)
      dat <- simulate_referee_data(
        n = m,
        scenario = scenario,
        seed = seed + 1000L * a + chunk,
        intervention = a,
        censoring = FALSE,
        max_time = max_time,
        dt = dt
      )
      for (j in target_events) {
        event_times <- dat[EVENT == j, TIME]
        for (tau in target_times) {
          counts[Time == tau & Event == j, count := count + sum(event_times <= tau)]
        }
      }
      done <- done + m
    }

    counts[, `:=`(
      Scenario = scenario,
      Intervention = paste0("A=", a),
      Estimand = "Abs Risk",
      True = count / n_truth
    )]
    truth <- rbind(truth, counts[, .(Scenario, Intervention, Estimand, Event, Time, True)])
  }

  rd <- dcast(truth[Estimand == "Abs Risk"], Scenario + Event + Time ~ Intervention,
              value.var = "True")
  rd[, `:=`(
    Intervention = "[A=1] - [A=0]",
    Estimand = "Risk Diff",
    True = `A=1` - `A=0`
  )]

  rr <- dcast(truth[Estimand == "Abs Risk"], Scenario + Event + Time ~ Intervention,
              value.var = "True")
  rr[, `:=`(
    Intervention = "[A=1] / [A=0]",
    Estimand = "Rel Risk",
    True = `A=1` / `A=0`
  )]

  rbind(
    truth,
    rd[, .(Scenario, Intervention, Estimand, Event, Time, True)],
    rr[, .(Scenario, Intervention, Estimand, Event, Time, True)]
  )
}

make_concrete_model <- function(library = c("minimal", "rich", "coxnet", "survsl", "survsl_nohal"),
                                include_events = c(0, 1, 2),
                                ps_library = c("SL.glm", "SL.glmnet")) {
  library <- match.arg(library)
  covs <- c("W1", "W2", "W3", "W4", "W5")
  main_rhs <- paste(c("ARM", covs), collapse = " + ")
  rich_rhs <- paste(
    c("ARM", covs, "ARM:W1", "ARM:W2", "ARM:W3", "I(W2^2)", "I(W4^2)"),
    collapse = " + "
  )

  model <- list(ARM = ps_library)
  for (j in include_events) {
    lhs <- paste0("Surv(TIME, EVENT == ", j, ") ~ ")
    if (identical(library, "minimal")) {
      model[[as.character(j)]] <- list(
        TrtOnly = as.formula(paste0(lhs, "ARM")),
        MainTerms = as.formula(paste0(lhs, main_rhs))
      )
    } else if (identical(library, "rich")) {
      model[[as.character(j)]] <- list(
        TrtOnly = as.formula(paste0(lhs, "ARM")),
        MainTerms = as.formula(paste0(lhs, main_rhs)),
        RichInteractions = as.formula(paste0(lhs, rich_rhs))
      )
    } else if (identical(library, "coxnet")) {
      model[[as.character(j)]] <- list(
        TrtOnly = as.formula(paste0(lhs, "ARM")),
        MainTerms = as.formula(paste0(lhs, main_rhs)),
        Coxnet = "coxnet"
      )
    } else if (identical(library, "survsl")) {
      model[[as.character(j)]] <- list(
        TrtOnly = as.formula(paste0(lhs, "ARM")),
        MainTerms = as.formula(paste0(lhs, main_rhs)),
        Coxnet = "coxnet",
        RandomForest = "rsf",
        AdditiveHazards = "aareg",
        HAL = "hal"
      )
    } else if (identical(library, "survsl_nohal")) {
      model[[as.character(j)]] <- list(
        TrtOnly = as.formula(paste0(lhs, "ARM")),
        MainTerms = as.formula(paste0(lhs, main_rhs)),
        Coxnet = "coxnet",
        RandomForest = "rsf",
        AdditiveHazards = "aareg"
      )
    }
  }
  model
}

capture_run <- function(expr) {
  warnings <- character()
  start <- Sys.time()
  value <- tryCatch(
    withCallingHandlers(
      suppressMessages(expr),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  list(value = value, warnings = unique(warnings), elapsed_sec = elapsed)
}

extract_concrete_diagnostics <- function(fit, estimator_id, elapsed_sec, warnings) {
  Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL

  conv <- attr(fit, "TmleConverged")
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

  min_denom <- unlist(lapply(fit, function(x) {
    nw <- x$NuisanceWeight
    denom <- 1 / nw
    as.numeric(min(denom, na.rm = TRUE))
  }), use.names = FALSE)

  data.table(
    EstimatorID = estimator_id,
    Converged = isTRUE(conv$converged),
    Step = suppressWarnings(as.numeric(conv$step)),
    FinalNormPnEIC = tail(attr(fit, "NormPnEICs"), 1),
    MaxRatio = max(stop_dt$ratio, na.rm = TRUE),
    MaxRelativeRatio = max(stop_dt$RelativeRatio, na.rm = TRUE),
    MaxAbsPnEIC = max(stop_dt$AbsPnEIC, na.rm = TRUE),
    FailingComponents = sum(!stop_dt$check, na.rm = TRUE),
    RelativeFailingComponents = sum(stop_dt$RelativeRatio > 1, na.rm = TRUE),
    MinNuisanceDenom = min(min_denom, na.rm = TRUE),
    EICStopRule = stop_rule,
    EICStopAbsTol = stop_abs_tol,
    RuntimeSec = elapsed_sec,
    Warnings = paste(warnings, collapse = " | ")
  )
}

run_concrete_estimator <- function(dat,
                                   target_times,
                                   target_events,
                                   estimator_id,
                                   update_method = "standard",
                                   model_library = "rich",
                                   max_update_iter = 200,
                                   one_step_eps = 0.1,
                                   min_nuisance = NULL,
                                   eic_stop_rule = "relative",
                                   eic_stop_abs_tol = 0,
                                   cv_v = 2,
                                   ps_library = c("SL.glm", "SL.glmnet")) {
  run <- capture_run({
    model <- make_concrete_model(
      library = model_library,
      include_events = sort(unique(dat$EVENT)),
      ps_library = ps_library
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
      CVArg = list(V = cv_v),
      Model = model,
      MaxUpdateIter = max_update_iter,
      OneStepEps = one_step_eps,
      Verbose = FALSE,
      GComp = TRUE,
      ReturnModels = FALSE,
      RenameCovs = FALSE,
      UpdateMethod = update_method,
      EICStopRule = eic_stop_rule,
      EICStopAbsTol = eic_stop_abs_tol
    )
    if (!is.null(min_nuisance)) {
      format_args$MinNuisance <- min_nuisance
    }
    args <- do.call(concrete::formatArguments, format_args)
    fit <- concrete::doConcrete(args)
    out <- concrete::getOutput(
      fit,
      Estimand = c("Risk", "RD", "RR"),
      GComp = TRUE,
      Simultaneous = FALSE
    )
    list(fit = fit, output = as.data.table(out))
  })

  if (inherits(run$value, "error")) {
    return(list(
      estimates = data.table(),
      diagnostics = data.table(
        EstimatorID = estimator_id,
        Converged = FALSE,
        Step = NA_real_,
        FinalNormPnEIC = NA_real_,
        MaxRatio = NA_real_,
        MaxRelativeRatio = NA_real_,
        MaxAbsPnEIC = NA_real_,
        FailingComponents = NA_integer_,
        RelativeFailingComponents = NA_integer_,
        MinNuisanceDenom = NA_real_,
        EICStopRule = eic_stop_rule,
        EICStopAbsTol = eic_stop_abs_tol,
        RuntimeSec = run$elapsed_sec,
        Warnings = paste(run$warnings, collapse = " | "),
        Error = conditionMessage(run$value)
      )
    ))
  }

  est <- run$value$output
  est[, `:=`(
    EstimatorID = estimator_id,
    Package = "concrete",
    UpdateMethod = update_method,
    ModelLibrary = model_library,
    MaxUpdateIter = max_update_iter,
    EICStopRule = eic_stop_rule,
    EICStopAbsTol = eic_stop_abs_tol
  )]
  diag <- extract_concrete_diagnostics(
    fit = run$value$fit,
    estimator_id = estimator_id,
    elapsed_sec = run$elapsed_sec,
    warnings = run$warnings
  )
  diag[, Error := NA_character_]
  list(estimates = est, diagnostics = diag)
}

add_contrasts <- function(risk_dt) {
  risk_dt <- copy(risk_dt)
  risk_dt[, Event := as.numeric(Event)]
  risk_dt[, Time := as.numeric(Time)]
  rd <- risk_dt[Estimand == "Abs Risk",
                .(`Pt Est` = `Pt Est`[Intervention == "A=1"] -
                    `Pt Est`[Intervention == "A=0"],
                  se = sqrt(sum(se^2, na.rm = TRUE))),
                by = setdiff(names(risk_dt),
                             c("Intervention", "Estimand", "Pt Est", "se",
                               "CI Low", "CI Hi", "SimCI Low", "SimCI Hi"))]
  rd[, `:=`(
    Intervention = "[A=1] - [A=0]",
    Estimand = "Risk Diff",
    `CI Low` = `Pt Est` - 1.96 * se,
    `CI Hi` = `Pt Est` + 1.96 * se
  )]

  rr <- risk_dt[Estimand == "Abs Risk",
                .(`Pt Est` = `Pt Est`[Intervention == "A=1"] /
                    `Pt Est`[Intervention == "A=0"],
                  se = NA_real_),
                by = setdiff(names(risk_dt),
                             c("Intervention", "Estimand", "Pt Est", "se",
                               "CI Low", "CI Hi", "SimCI Low", "SimCI Hi"))]
  rr[, `:=`(
    Intervention = "[A=1] / [A=0]",
    Estimand = "Rel Risk",
    `CI Low` = NA_real_,
    `CI Hi` = NA_real_
  )]

  rbindlist(list(risk_dt, rd, rr), fill = TRUE)
}

run_aalen_johansen <- function(dat, target_times, target_events) {
  run <- capture_run({
    fit <- survival::survfit(
      survival::Surv(time = TIME, event = as.factor(EVENT)) ~ ARM,
      data = dat,
      ctype = 1
    )
    s <- summary(fit, times = target_times, extend = TRUE)
    pstate <- as.data.table(s$pstate)
    stderr <- as.data.table(s$std.err)
    state_names <- colnames(pstate)
    if (is.null(state_names)) {
      state_names <- as.character(seq_len(ncol(pstate)) - 1L)
    }
    colnames(pstate) <- paste0("risk_", state_names)
    colnames(stderr) <- paste0("se_", state_names)
    out <- cbind(
      data.table(Strata = as.character(s$strata), Time = s$time),
      pstate,
      stderr
    )
    out[, Intervention := paste0("A=", sub("^ARM=", "", Strata))]

    long <- melt(
      out,
      id.vars = c("Intervention", "Time"),
      measure.vars = patterns("^risk_", "^se_"),
      value.name = c("Pt Est", "se"),
      variable.name = "State"
    )
    long[, Event := as.numeric(sub("^risk_", "", names(pstate)[State]))]
    long <- long[Event %in% target_events]
    long[, `:=`(
      Estimand = "Abs Risk",
      Estimator = "Aalen-Johansen",
      EstimatorID = "aalen_johansen",
      Package = "survival",
      UpdateMethod = NA_character_,
      ModelLibrary = NA_character_,
      MaxUpdateIter = NA_real_,
      `CI Low` = `Pt Est` - 1.96 * se,
      `CI Hi` = `Pt Est` + 1.96 * se
    )]
    long[, State := NULL]
    add_contrasts(long)
  })

  if (inherits(run$value, "error")) {
    return(list(
      estimates = data.table(),
      diagnostics = data.table(
        EstimatorID = "aalen_johansen",
        RuntimeSec = run$elapsed_sec,
        Warnings = paste(run$warnings, collapse = " | "),
        Error = conditionMessage(run$value)
      )
    ))
  }

  diag <- data.table(
    EstimatorID = "aalen_johansen",
    RuntimeSec = run$elapsed_sec,
    Warnings = paste(run$warnings, collapse = " | "),
    Error = NA_character_
  )
  list(estimates = run$value, diagnostics = diag)
}

run_survtmle_estimator <- function(dat,
                                   target_times,
                                   target_events,
                                   months = 6,
                                   cv_v = 2,
                                   max_iter = 50) {
  estimator_id <- paste0("survtmle_", months, "mo")
  run <- capture_run({
    target_times_local <- target_times
    target_events_local <- target_events
    library(SuperLearner)
    bins_per_year <- 12 / months
    disc_time <- pmax(1L, ceiling(dat$TIME / 365.25 * bins_per_year))
    disc_targets <- pmax(1L, ceiling(target_times_local / 365.25 * bins_per_year))
    W <- as.data.frame(dat[, .SD, .SDcols = c("W1", "W2", "W3", "W4", "W5")])

    fit <- do.call(survtmle::survtmle, list(
      ftime = disc_time,
      ftype = dat$EVENT,
      trt = dat$ARM,
      adjustVars = W,
      t0 = max(disc_targets),
      SL.ftime = c("SL.glm", "SL.glmnet"),
      SL.ctime = c("SL.glm", "SL.glmnet"),
      SL.trt = c("SL.glm", "SL.glmnet"),
      returnIC = TRUE,
      returnModels = TRUE,
      method = "hazard",
      verbose = FALSE,
      maxIter = max_iter,
      Gcomp = FALSE,
      cvControl = list(V = as.integer(cv_v), stratifyCV = FALSE, shuffle = TRUE),
      ftypeOfInterest = target_events_local,
      trtOfInterest = c(0, 1)
    ))
    tp <- survtmle::timepoints(fit, times = disc_targets)
    invisible(capture.output(tp_print <- print(tp)))
    risk <- as.data.table(tp_print$est, keep.rownames = "Row")
    var <- as.data.table(tp_print$var, keep.rownames = "Row")
    risk_long <- melt(risk, id.vars = "Row", variable.name = "TimeIndex",
                      value.name = "Pt Est")
    var_long <- melt(var, id.vars = "Row", variable.name = "TimeIndex",
                     value.name = "var")
    out <- merge(risk_long, var_long, by = c("Row", "TimeIndex"))
    out[, Time := target_times_local[as.integer(sub("^t", "", TimeIndex))]]
    row_parts <- tstrsplit(trimws(out$Row), "\\s+")
    out[, `:=`(
      Intervention = paste0("A=", row_parts[[1]]),
      Event = as.numeric(row_parts[[2]]),
      se = sqrt(var),
      Row = NULL,
      TimeIndex = NULL,
      var = NULL
    )]
    out <- out[which(out[["Event"]] %in% target_events_local)]
    out[, `:=`(
      Estimand = "Abs Risk",
      Estimator = estimator_id,
      EstimatorID = estimator_id,
      Package = "survtmle",
      UpdateMethod = NA_character_,
      ModelLibrary = paste0(months, "mo-discrete"),
      MaxUpdateIter = max_iter,
      `CI Low` = `Pt Est` - 1.96 * se,
      `CI Hi` = `Pt Est` + 1.96 * se
    )]
    add_contrasts(out)
  })

  if (inherits(run$value, "error")) {
    return(list(
      estimates = data.table(),
      diagnostics = data.table(
        EstimatorID = estimator_id,
        RuntimeSec = run$elapsed_sec,
        Warnings = paste(run$warnings, collapse = " | "),
        Error = conditionMessage(run$value)
      )
    ))
  }

  diag <- data.table(
    EstimatorID = estimator_id,
    RuntimeSec = run$elapsed_sec,
    Warnings = paste(run$warnings, collapse = " | "),
    Error = NA_character_
  )
  list(estimates = run$value, diagnostics = diag)
}

run_one_replicate <- function(scenario,
                              n,
                              rep,
                              seed,
                              target_times,
                              target_events,
                              cv_v = 2,
                              max_update_iter = 200,
                              include_survsl = FALSE,
                              include_hal_sl = FALSE,
                              include_stabilized = FALSE,
                              stabilized_min_nuisance = 0.05,
                              stabilized_one_step_eps = 0.01,
                              stabilized_stop_rule = "hybrid",
                              stabilized_stop_abs_tol = 1e-3,
                              include_survtmle = FALSE,
                              survtmle_months = c(6),
                              max_time = 1826,
                              dt = 14) {
  dat <- simulate_referee_data(
    n = n,
    scenario = scenario,
    seed = seed,
    censoring = TRUE,
    max_time = max_time,
    dt = dt
  )

  concrete_grid <- data.table(
    EstimatorID = c(
      "concrete_standard_minimal",
      "concrete_standard_rich",
      "concrete_adaptive_rich"
    ),
    UpdateMethod = c("standard", "standard", "adaptive"),
    ModelLibrary = c("minimal", "rich", "rich"),
    MaxUpdateIter = rep(max_update_iter, 3),
    OneStepEps = c(0.1, 0.1, 0.1),
    MinNuisance = NA_real_,
    EICStopRule = "relative",
    EICStopAbsTol = 0
  )
  if (isTRUE(include_survsl)) {
    concrete_grid <- rbind(
      concrete_grid,
      data.table(
        EstimatorID = c(
          "concrete_standard_survsl_nohal",
          "concrete_adaptive_survsl_nohal"
        ),
        UpdateMethod = c("standard", "adaptive"),
        ModelLibrary = c("survsl_nohal", "survsl_nohal"),
        MaxUpdateIter = rep(max_update_iter, 2),
        OneStepEps = c(0.1, 0.1),
        MinNuisance = NA_real_,
        EICStopRule = "relative",
        EICStopAbsTol = 0
      ),
      fill = TRUE
    )
  }
  if (isTRUE(include_stabilized)) {
    stabilized_grid <- data.table(
      EstimatorID = "concrete_adaptive_rich_stabilized",
      UpdateMethod = "adaptive",
      ModelLibrary = "rich",
      MaxUpdateIter = max_update_iter,
      OneStepEps = stabilized_one_step_eps,
      MinNuisance = stabilized_min_nuisance,
      EICStopRule = stabilized_stop_rule,
      EICStopAbsTol = stabilized_stop_abs_tol
    )
    if (isTRUE(include_survsl)) {
      stabilized_grid <- rbind(
        stabilized_grid,
        data.table(
          EstimatorID = "concrete_adaptive_survsl_nohal_stabilized",
          UpdateMethod = "adaptive",
          ModelLibrary = "survsl_nohal",
          MaxUpdateIter = max_update_iter,
          OneStepEps = stabilized_one_step_eps,
          MinNuisance = stabilized_min_nuisance,
          EICStopRule = stabilized_stop_rule,
          EICStopAbsTol = stabilized_stop_abs_tol
        ),
        fill = TRUE
      )
    }
    concrete_grid <- rbind(concrete_grid, stabilized_grid, fill = TRUE)
  }
  if (isTRUE(include_hal_sl)) {
    concrete_grid <- rbind(
      concrete_grid,
      data.table(
        EstimatorID = "concrete_standard_survsl",
        UpdateMethod = "standard",
        ModelLibrary = "survsl",
        MaxUpdateIter = max_update_iter,
        OneStepEps = 0.1,
        MinNuisance = NA_real_,
        EICStopRule = "relative",
        EICStopAbsTol = 0
      ),
      fill = TRUE
    )
  }

  results <- list()
  for (i in seq_len(nrow(concrete_grid))) {
    cfg <- concrete_grid[i]
    results[[cfg$EstimatorID]] <- run_concrete_estimator(
      dat = dat,
      target_times = target_times,
      target_events = target_events,
      estimator_id = cfg$EstimatorID,
      update_method = cfg$UpdateMethod,
      model_library = cfg$ModelLibrary,
      max_update_iter = cfg$MaxUpdateIter,
      one_step_eps = cfg$OneStepEps,
      min_nuisance = if (is.na(cfg$MinNuisance)) NULL else cfg$MinNuisance,
      eic_stop_rule = cfg$EICStopRule,
      eic_stop_abs_tol = cfg$EICStopAbsTol,
      cv_v = cv_v
    )
  }

  results[["aalen_johansen"]] <- run_aalen_johansen(
    dat = dat,
    target_times = target_times,
    target_events = target_events
  )

  if (isTRUE(include_survtmle) && requireNamespace("survtmle", quietly = TRUE)) {
    for (months in survtmle_months) {
      key <- paste0("survtmle_", months, "mo")
      results[[key]] <- run_survtmle_estimator(
        dat = dat,
        target_times = target_times,
        target_events = target_events,
        months = months,
        cv_v = cv_v,
        max_iter = 50
      )
    }
  }

  estimates <- rbindlist(lapply(results, `[[`, "estimates"), fill = TRUE)
  diagnostics <- rbindlist(lapply(results, `[[`, "diagnostics"), fill = TRUE)

  estimates[, `:=`(Scenario = scenario, N = n, Rep = rep, Seed = seed)]
  diagnostics[, `:=`(Scenario = scenario, N = n, Rep = rep, Seed = seed)]

  list(
    scenario = scenario,
    n = n,
    rep = rep,
    seed = seed,
    estimates = estimates,
    diagnostics = diagnostics
  )
}

read_replicates <- function(output_dir) {
  files <- list.files(output_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("truth_|summary_|metrics_", basename(files))]
  if (!length(files)) stop("No replicate RDS files found in ", output_dir)
  lapply(files, readRDS)
}

summarize_simulation <- function(output_dir, truth_file = NULL) {
  reps <- read_replicates(output_dir)
  estimates <- rbindlist(lapply(reps, `[[`, "estimates"), fill = TRUE)
  diagnostics <- rbindlist(lapply(reps, `[[`, "diagnostics"), fill = TRUE)

  if (is.null(truth_file)) {
    truth_files <- list.files(output_dir, pattern = "^truth_.*\\.rds$", full.names = TRUE)
    truth <- rbindlist(lapply(truth_files, readRDS), fill = TRUE)
  } else {
    truth <- readRDS(truth_file)
  }

  joined <- merge(
    estimates,
    truth,
    by = c("Scenario", "Intervention", "Estimand", "Event", "Time"),
    all.x = TRUE
  )
  joined[, Error := `Pt Est` - True]
  joined[, Covered := `CI Low` <= True & `CI Hi` >= True]

  metric_cols <- c(
    "Scenario", "N", "EstimatorID", "Package", "Estimator", "UpdateMethod",
    "ModelLibrary", "MaxUpdateIter", "Estimand", "Intervention", "Event", "Time"
  )
  metrics <- joined[, .(
    Reps = uniqueN(Rep),
    MeanEstimate = mean(`Pt Est`, na.rm = TRUE),
    Truth = mean(True, na.rm = TRUE),
    Bias = mean(Error, na.rm = TRUE),
    PercentBias = mean(Error / True, na.rm = TRUE) * 100,
    EmpiricalSD = stats::sd(`Pt Est`, na.rm = TRUE),
    MeanSE = mean(se, na.rm = TRUE),
    RMSE = sqrt(mean(Error^2, na.rm = TRUE)),
    Coverage = mean(Covered, na.rm = TRUE)
  ), by = metric_cols]

  diag_metrics <- diagnostics[, .(
    Reps = uniqueN(Rep),
    ErrorRate = mean(!is.na(Error)),
    ConvergenceRate = mean(Converged, na.rm = TRUE),
    MedianStep = stats::median(Step, na.rm = TRUE),
    MedianMaxRatio = stats::median(MaxRatio, na.rm = TRUE),
    MeanRuntimeSec = mean(RuntimeSec, na.rm = TRUE)
  ), by = .(Scenario, N, EstimatorID)]

  list(
    estimates = joined,
    diagnostics = diagnostics,
    metrics = metrics,
    diagnostic_metrics = diag_metrics
  )
}
