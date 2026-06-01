###############################################################################
# Required libraries (for example):
###############################################################################
# library(data.table)
# library(nleqslv)
# library(foreach)       # for parallel Jacobi approach
# library(doParallel)    # or another parallel backend

###############################################################################
# UTILITY FUNCTIONS
###############################################################################

#' Print diagnostic information about TMLE convergence
#'
#' @param OneStepStop data.table with convergence check results
#' @param NormPnEIC numeric, the current norm of PnEIC
#' @keywords internal
printOneStepDiagnostics <- function(OneStepStop, NormPnEIC) {
  ratio <- NULL
  # Sort by the ratio, pick top 3
  Worst <- OneStepStop[, !"check"][, ratio := round(ratio, 2)][order(ratio, decreasing = TRUE)]
  print(Worst[1:min(nrow(Worst), 3), ])
  cat("Norm PnEIC = ", NormPnEIC, "\n", sep = "")
}

#' Deep-copy an Estimates list so we can do parallel updates safely (Jacobi).
#'
#' @param E The original Estimates object (list of arms).
#'
#' @return A new, independent copy of E.
#' @keywords internal
deepcopyEstimates <- function(E) {
  copy_obj <- function(x) {
    if (inherits(x, "data.table")) {
      return(data.table::copy(x))
    }
    if (is.list(x)) {
      y <- lapply(x, copy_obj)
      attributes(y) <- attributes(x)
      return(y)
    }
    y <- x
    attributes(y) <- attributes(x)
    y
  }

  newE <- lapply(E, copy_obj)
  attributes(newE) <- attributes(E)
  return(newE)
}

summarizeUpdateState <- function(SummEIC, TargetTime, TargetEvent,
                                 EICStopRule = "relative",
                                 EICStopAbsTol = 0) {
  Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- ratio <-
    RelativeRatio <- AbsoluteRatio <- AbsPnEIC <- StopCriteria <- check <- NULL

  target <- data.table::copy(targetSummEIC(SummEIC, TargetTime, TargetEvent))
  if (!nrow(target)) {
    return(list(
      max_ratio = NA_real_,
      failing_components = NA_integer_,
      worst_trt = NA_character_,
      worst_time = NA_real_,
      worst_event = NA_real_,
      worst_ratio = NA_real_,
      worst_pneic = NA_real_,
      worst_threshold = NA_real_,
      max_relative_ratio = NA_real_,
      max_absolute_ratio = NA_real_,
      max_abs_pneic = NA_real_,
      relative_failing_components = NA_integer_,
      worst_relative_ratio = NA_real_,
      worst_absolute_ratio = NA_real_
    ))
  }

  target <- makeOneStepStop(
    SummEIC = target,
    EICStopRule = EICStopRule,
    EICStopAbsTol = EICStopAbsTol
  )
  target <- target[order(ratio, decreasing = TRUE)]
  worst <- target[1]
  list(
    max_ratio = max(target$ratio, na.rm = TRUE),
    failing_components = sum(!target$check, na.rm = TRUE),
    worst_trt = as.character(worst$Trt),
    worst_time = as.numeric(worst$Time),
    worst_event = as.numeric(worst$Event),
    worst_ratio = as.numeric(worst$ratio),
    worst_pneic = as.numeric(worst$PnEIC),
    worst_threshold = as.numeric(worst$StopCriteria),
    max_relative_ratio = max(target$RelativeRatio, na.rm = TRUE),
    max_absolute_ratio = max(target$AbsoluteRatio, na.rm = TRUE),
    max_abs_pneic = max(target$AbsPnEIC, na.rm = TRUE),
    relative_failing_components = sum(target$RelativeRatio > 1, na.rm = TRUE),
    worst_relative_ratio = as.numeric(worst$RelativeRatio),
    worst_absolute_ratio = as.numeric(worst$AbsoluteRatio)
  )
}

appendUpdateTrace <- function(trace, method, step, line_iter, status, alpha,
                              norm_before, norm_after, SummEIC,
                              TargetTime, TargetEvent,
                              EICStopRule = "relative",
                              EICStopAbsTol = 0) {
  state <- if (is.null(SummEIC)) {
    list(
      max_ratio = NA_real_,
      failing_components = NA_integer_,
      worst_trt = NA_character_,
      worst_time = NA_real_,
      worst_event = NA_real_,
      worst_ratio = NA_real_,
      worst_pneic = NA_real_,
      worst_threshold = NA_real_,
      max_relative_ratio = NA_real_,
      max_absolute_ratio = NA_real_,
      max_abs_pneic = NA_real_,
      relative_failing_components = NA_integer_,
      worst_relative_ratio = NA_real_,
      worst_absolute_ratio = NA_real_
    )
  } else {
    summarizeUpdateState(
      SummEIC,
      TargetTime,
      TargetEvent,
      EICStopRule = EICStopRule,
      EICStopAbsTol = EICStopAbsTol
    )
  }

  data.table::rbindlist(
    list(
      trace,
      data.table::data.table(
        Method = method,
        Step = step,
        LineIter = line_iter,
        Status = status,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol,
        Alpha = alpha,
        NormBefore = norm_before,
        NormAfter = norm_after,
        MaxRatio = state$max_ratio,
        MaxRelativeRatio = state$max_relative_ratio,
        MaxAbsoluteRatio = state$max_absolute_ratio,
        MaxAbsPnEIC = state$max_abs_pneic,
        FailingComponents = state$failing_components,
        RelativeFailingComponents = state$relative_failing_components,
        WorstTrt = state$worst_trt,
        WorstTime = state$worst_time,
        WorstEvent = state$worst_event,
        WorstRatio = state$worst_ratio,
        WorstRelativeRatio = state$worst_relative_ratio,
        WorstAbsoluteRatio = state$worst_absolute_ratio,
        WorstPnEIC = state$worst_pneic,
        WorstThreshold = state$worst_threshold
      )
    ),
    fill = TRUE
  )
}

###############################################################################
# TMLE update function
###############################################################################

#' TMLE update with multiple implementation options
#'
#' @param Estimates list of estimates for each arm/treatment, containing Hazards, Surv, SummEIC, etc.
#' @param SummEIC data.table summarizing the current empirical mean EIC (PnEIC), etc.
#' @param Data data.table containing the original data
#' @param TargetEvent numeric vector of event types
#' @param TargetTime numeric vector of target times
#' @param MaxUpdateIter maximum number of steps/iterations
#' @param OneStepEps initial step size for the incremental updates
#' @param NormPnEIC numeric, initial ||PnEIC|| for reference
#' @param Verbose boolean, whether to print debugging output
#' @param Method character - one of
#'   - "standard"
#'   - "adaptive"
#' @param EICStopRule character stopping rule for empirical mean EIC checks.
#'   Supported values are `"relative"`, `"absolute"`, and `"hybrid"`.
#' @param EICStopAbsTol numeric absolute `|PnEIC|` tolerance used by the
#'   `"absolute"` and `"hybrid"` stopping rules.
#'
#' @return Updated \code{Estimates} object, with updated hazards, SummEIC, etc.
#' @importFrom nleqslv nleqslv
#' @import foreach
#' @import doParallel
#' @import parallel
#' @keywords internal
doTmleUpdate <- function(Estimates, SummEIC, Data, TargetEvent, TargetTime,
                         MaxUpdateIter, OneStepEps, NormPnEIC, Verbose,
                         Method = "standard",
                         EICStopRule = "relative",
                         EICStopAbsTol = 0) {

  `.` <- Trt <- Time <- Event <- seEIC <- `seEIC/(sqrt(n)log(n))` <- PnEIC <- NULL
  Method <- getUpdateMethod(Method)
  EICStopRule <- getEICStopRule(EICStopRule)
  EICStopAbsTol <- getEICStopAbsTol(EICStopAbsTol)
  EvalTimes <- attr(Estimates, "Times")
  T.tilde   <- Data[[attr(Data, "EventTime")]]
  Delta     <- Data[[attr(Data, "EventType")]]

  WorkingEps <- OneStepEps
  NormPnEICs <- NormPnEIC
  UpdateTrace <- appendUpdateTrace(
    trace = data.table::data.table(),
    method = Method,
    step = 0L,
    line_iter = NA_integer_,
    status = "initial",
    alpha = NA_real_,
    norm_before = NA_real_,
    norm_after = NormPnEIC,
    SummEIC = SummEIC,
    TargetTime = TargetTime,
    TargetEvent = TargetEvent,
    EICStopRule = EICStopRule,
    EICStopAbsTol = EICStopAbsTol
  )

  # Dimension k = (#Interventions) * (#TargetEvents) * (#TargetTimes)
  k_dim <- length(TargetEvent) * length(TargetTime) * length(Estimates)
  if (Verbose) cat("Problem dimension: k =", k_dim, "\n")

  ###############################################################################
  # 1) "standard" Approach: iterative small steps
  ###############################################################################
  if (Method == "standard") {
    if (Verbose) cat("Using standard universal LFM approach (iterative small steps)\n")
    StepNum <- 1
    IterNum <- 1

    while (StepNum <= MaxUpdateIter & IterNum <= MaxUpdateIter * 2) {
      IterNum <- IterNum + 1
      if (Verbose) {
        cat("Starting step", StepNum, "with update epsilon =", WorkingEps, "\n")
      }

      CurrentObjective <- targetUpdateObjective(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      # Update hazards and EIC for each arm
      newEsts <- lapply(Estimates, function(est.a) {
        NewHazards <- updateHazard(
          GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
          Hazards        = est.a[["Hazards"]],
          TotalSurv      = est.a[["EvntFreeSurv"]],
          NuisanceWeight = est.a[["NuisanceWeight"]],
          EvalTimes      = EvalTimes,
          T.tilde        = T.tilde,
          Delta          = Delta,
          PnEIC          = est.a[["SummEIC"]],
          NormPnEIC      = NormPnEIC,
          OneStepEps     = WorkingEps,
          TargetEvent    = TargetEvent,
          TargetTime     = TargetTime
        )
        # Fix any NA
        NewHazards <- lapply(NewHazards, function(hazmat) {
          if (anyNA(hazmat)) hazmat[is.na(hazmat) | is.nan(hazmat)] <- 0
          hazmat
        })

        # Recompute survival
        NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(hz) exp(-cumsum(hz)))
        NewSurv[NewSurv < 1e-12 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-12

        # Recompute EIC
        NewIC <- getIC(
          GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
          Hazards        = NewHazards,
          TotalSurv      = NewSurv,
          NuisanceWeight = est.a[["NuisanceWeight"]],
          TargetEvent    = TargetEvent,
          TargetTime     = TargetTime,
          T.tilde        = T.tilde,
          Delta          = Delta,
          EvalTimes      = EvalTimes,
          GComp          = FALSE
        )
        list(
          Hazards      = NewHazards,
          EvntFreeSurv = NewSurv,
          SummEIC      = summarizeIC(NewIC),
          IC           = NewIC
        )
      })

      # Check improvement: gather EIC
      NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
        cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])
      }))
      NewNormPnEIC <- getNormPnEIC(targetSummEIC(NewSummEIC, TargetTime, TargetEvent)[, PnEIC])
      NewObjective <- targetUpdateObjective(
        SummEIC = NewSummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (anyNA(NewNormPnEIC) || is.nan(NewNormPnEIC) || is.infinite(NewNormPnEIC) ||
          anyNA(NewObjective) || is.nan(NewObjective) || is.infinite(NewObjective)) {
        UpdateTrace <- appendUpdateTrace(
          trace = UpdateTrace,
          method = Method,
          step = StepNum,
          line_iter = NA_integer_,
          status = "rejected_invalid_norm",
          alpha = WorkingEps,
          norm_before = NormPnEIC,
          norm_after = NewNormPnEIC,
          SummEIC = NewSummEIC,
          TargetTime = TargetTime,
          TargetEvent = TargetEvent,
          EICStopRule = EICStopRule,
          EICStopAbsTol = EICStopAbsTol
        )
        if (Verbose) cat("Update produced an invalid ||PnEIC||, halving OneStepEps\n")
        WorkingEps <- WorkingEps / 2
        next
      }
      if (CurrentObjective < NewObjective) {
        UpdateTrace <- appendUpdateTrace(
          trace = UpdateTrace,
          method = Method,
          step = StepNum,
          line_iter = NA_integer_,
          status = "rejected_increased_objective",
          alpha = WorkingEps,
          norm_before = NormPnEIC,
          norm_after = NewNormPnEIC,
          SummEIC = NewSummEIC,
          TargetTime = TargetTime,
          TargetEvent = TargetEvent,
          EICStopRule = EICStopRule,
          EICStopAbsTol = EICStopAbsTol
        )
        if (Verbose) cat("Update increased the active convergence objective, halving OneStepEps\n")
        WorkingEps <- WorkingEps / 2
        next
      }
      StepNum <- StepNum + 1

      # Commit updates
      for (a in seq_along(Estimates)) {
        Estimates[[a]][["Hazards"]]      <- newEsts[[a]][["Hazards"]]
        Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
        Estimates[[a]][["SummEIC"]]      <- newEsts[[a]][["SummEIC"]]
        Estimates[[a]][["IC"]]           <- newEsts[[a]][["IC"]]
      }

      SummEIC   <- NewSummEIC
      NormPnEIC <- NewNormPnEIC
      NormPnEICs <- c(NormPnEICs, NormPnEIC)
      UpdateTrace <- appendUpdateTrace(
        trace = UpdateTrace,
        method = Method,
        step = StepNum - 1L,
        line_iter = NA_integer_,
        status = "accepted",
        alpha = WorkingEps,
        norm_before = tail(NormPnEICs, 2)[1],
        norm_after = NormPnEIC,
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      OneStepStop <- targetOneStepStop(
        SummEIC = NewSummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)

      # if all converged, stop
      if (all(sapply(OneStepStop[["check"]], isTRUE))) {
        attr(Estimates, "TmleConverged") <- list(converged = TRUE, step = StepNum)
        attr(Estimates, "NormPnEICs")    <- NormPnEICs
        attr(Estimates, "TmleUpdateTrace") <- UpdateTrace
        return(Estimates)
      }
    }
    # if we exit
    warning("TMLE has not converged by step ", MaxUpdateIter,
            " - Estimates may not have the desired asymptotic properties")
    attr(Estimates, "TmleConverged") <- list(converged = FALSE, step = StepNum)
    attr(Estimates, "NormPnEICs")    <- NormPnEICs
    attr(Estimates, "TmleUpdateTrace") <- UpdateTrace
    return(Estimates)

    ###############################################################################
    # 2) "adaptive" Approach: Adaptive Line Search with rollback if norm goes up
    ###############################################################################
  } else if (Method == "adaptive") {
    if (Verbose) cat("Using Adaptive Line Search with rollback\n")

    # -----------------------------
    # Example: bounding parameters
    # -----------------------------
    alphaCap     <- 0.3    # Hard maximum on the step size
    maxHaz       <- 1000    # Example clamp on hazard
    c1           <- 0.1    # Armijo parameter
    tau_reduce   <- 0.5
    tau_increase <- 1.5
    max_line_iterations <- 10

    prev_alpha           <- NULL
    prev_norm            <- NULL
    prev_reduction_ratio <- NULL
    consecutive_gains    <- 0

    StepNum <- 1
    while (StepNum <= MaxUpdateIter) {
      if (Verbose) cat("Starting step", StepNum, "\n")

      # =========== DETERMINE INITIAL STEP SIZE =============
      if (!is.null(prev_alpha)) {
        # If we had multiple consecutive successes, we can try bigger alpha,
        # but never above alphaCap:
        if (consecutive_gains >= 3) {
          if (NormPnEIC < 1) {
            initial_alpha <- min(0.1, prev_alpha * tau_increase)
          } else {
            initial_alpha <- min(1.0, prev_alpha * tau_increase)
          }
        } else if (consecutive_gains > 0) {
          initial_alpha <- min(0.5, prev_alpha * tau_increase)
        } else {
          initial_alpha <- prev_alpha
        }
        if (Verbose) cat("  Using alpha =", initial_alpha, "\n")
      } else {
        # first iteration
        initial_alpha <- min(0.5, 5 * OneStepEps)
        if (Verbose) cat("  Initial alpha =", initial_alpha, "\n")
      }

      # Enforce a hard cap on alpha:
      alpha <- min(initial_alpha, alphaCap)
      if (Verbose && alpha < initial_alpha) {
        cat("  (Capped alpha at", alphaCap, ")\n")
      }

      # =========== LINE SEARCH BACKTRACKING =============
      line_iter <- 0
      improvement_found <- FALSE
      current_norm <- NormPnEIC
      current_objective <- targetUpdateObjective(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      while (line_iter < max_line_iterations && !improvement_found) {
        line_iter <- line_iter + 1
        if (Verbose) cat("  line search iter", line_iter, "with alpha =", alpha, "\n")

        # 1) Keep a deep copy of old state
        oldEstimates <- deepcopyEstimates(Estimates)
        oldSummEIC   <- data.table::copy(SummEIC)
        oldNorm      <- NormPnEIC

        # 2) Attempt an update with alpha
        newEsts <- tryCatch({
          lapply(Estimates, function(est.a) {
            NewHazards <- updateHazard(
              GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
              Hazards        = est.a[["Hazards"]],
              TotalSurv      = est.a[["EvntFreeSurv"]],
              NuisanceWeight = est.a[["NuisanceWeight"]],
              EvalTimes      = EvalTimes,
              T.tilde        = T.tilde,
              Delta          = Delta,
              PnEIC          = est.a[["SummEIC"]],
              NormPnEIC      = NormPnEIC,
              OneStepEps     = alpha,
              TargetEvent    = TargetEvent,
              TargetTime     = TargetTime
            )

            # =========== EXAMPLE: clamp the updated hazards ===========
            # For instance, limit them to [0, maxHaz].
            # If any are NA, set them to 0.
            # You can also clamp *lower* if needed.
            NewHazards <- lapply(NewHazards, function(hzmat) {
              # fix NA
              hzmat[is.na(hzmat) | is.nan(hzmat)] <- 0
              # clamp upper
              hzmat <- pmin(hzmat, maxHaz)
              # optionally clamp lower if you prefer:
              hzmat <- pmax(hzmat, 0)
              hzmat
            })

            # recompute survival
            NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(hz) exp(-cumsum(hz)))
            NewSurv[NewSurv < 1e-8 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-8

            # recompute EIC
            NewIC <- getIC(
              GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
              Hazards        = NewHazards,
              TotalSurv      = NewSurv,
              NuisanceWeight = est.a[["NuisanceWeight"]],
              TargetEvent    = TargetEvent,
              TargetTime     = TargetTime,
              T.tilde        = T.tilde,
              Delta          = Delta,
              EvalTimes      = EvalTimes,
              GComp          = FALSE
            )
            list(
              Hazards      = NewHazards,
              EvntFreeSurv = NewSurv,
              SummEIC      = summarizeIC(NewIC),
              IC           = NewIC
            )
          })
        }, error = function(e) {
          if (Verbose) cat("    Error in trial update:", e$message, "\n")
          NULL
        })

        # If update failed, revert & reduce alpha
        if (is.null(newEsts)) {
          UpdateTrace <- appendUpdateTrace(
            trace = UpdateTrace,
            method = Method,
            step = StepNum,
            line_iter = line_iter,
            status = "rejected_update_error",
            alpha = alpha,
            norm_before = oldNorm,
            norm_after = NA_real_,
            SummEIC = NULL,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )
          if (Verbose) cat("    Update failed, revert & reduce alpha\n")
          Estimates <- oldEstimates
          SummEIC   <- oldSummEIC
          NormPnEIC <- oldNorm
          alpha     <- alpha * tau_reduce
          next
        }

        # Check improvement
        NewSummEIC <- tryCatch({
          do.call(rbind, lapply(seq_along(newEsts), function(a) {
            cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])
          }))
        }, error = function(e) {
          if (Verbose) cat("    EIC gather error:", e$message, "\n")
          NULL
        })

        if (is.null(NewSummEIC)) {
          UpdateTrace <- appendUpdateTrace(
            trace = UpdateTrace,
            method = Method,
            step = StepNum,
            line_iter = line_iter,
            status = "rejected_eic_error",
            alpha = alpha,
            norm_before = oldNorm,
            norm_after = NA_real_,
            SummEIC = NULL,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )
          if (Verbose) cat("    EIC gather failed, revert & reduce alpha\n")
          Estimates <- oldEstimates
          SummEIC   <- oldSummEIC
          NormPnEIC <- oldNorm
          alpha     <- alpha * tau_reduce
          next
        }

        # compute new norm
        NewNormPnEIC <- tryCatch({
          getNormPnEIC(targetSummEIC(NewSummEIC, TargetTime, TargetEvent)[, PnEIC])
        }, error = function(e) {
          if (Verbose) cat("    norm error:", e$message, "\n")
          Inf
        })
        NewObjective <- tryCatch({
          targetUpdateObjective(
            SummEIC = NewSummEIC,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )
        }, error = function(e) {
          if (Verbose) cat("    objective error:", e$message, "\n")
          Inf
        })

        if (anyNA(NewNormPnEIC) || is.nan(NewNormPnEIC) || is.infinite(NewNormPnEIC) ||
            anyNA(NewObjective) || is.nan(NewObjective) || is.infinite(NewObjective)) {
          UpdateTrace <- appendUpdateTrace(
            trace = UpdateTrace,
            method = Method,
            step = StepNum,
            line_iter = line_iter,
            status = "rejected_invalid_norm",
            alpha = alpha,
            norm_before = oldNorm,
            norm_after = NewNormPnEIC,
            SummEIC = NewSummEIC,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )
          if (Verbose) cat("    Invalid new norm, revert & reduce alpha\n")
          Estimates <- oldEstimates
          SummEIC   <- oldSummEIC
          NormPnEIC <- oldNorm
          alpha     <- alpha * tau_reduce
          next
        }

        # Armijo condition on the active convergence objective.
        expected_reduction <- c1 * alpha * current_objective
        actual_reduction   <- current_objective - NewObjective

        if ((NewObjective >= current_objective) || (actual_reduction < expected_reduction)) {
          # NO improvement => revert, reduce alpha
          UpdateTrace <- appendUpdateTrace(
            trace = UpdateTrace,
            method = Method,
            step = StepNum,
            line_iter = line_iter,
            status = "rejected_armijo",
            alpha = alpha,
            norm_before = current_norm,
            norm_after = NewNormPnEIC,
            SummEIC = NewSummEIC,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )
          if (Verbose) {
            cat("    Objective went up or Armijo not satisfied => revert & reduce alpha\n")
            cat("    old objective:", current_objective, " new objective:", NewObjective, "\n")
          }
          Estimates <- oldEstimates
          SummEIC   <- oldSummEIC
          NormPnEIC <- oldNorm
          alpha     <- alpha * tau_reduce

        } else {
          # Improvement => commit
          if (Verbose) cat("    Good improvement => commit\n")
          for (a in seq_along(Estimates)) {
            Estimates[[a]][["Hazards"]]      <- newEsts[[a]][["Hazards"]]
            Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
            Estimates[[a]][["SummEIC"]]      <- newEsts[[a]][["SummEIC"]]
            Estimates[[a]][["IC"]]           <- newEsts[[a]][["IC"]]
          }
          SummEIC   <- NewSummEIC
          NormPnEIC <- NewNormPnEIC
          improvement_found <- TRUE
          UpdateTrace <- appendUpdateTrace(
            trace = UpdateTrace,
            method = Method,
            step = StepNum,
            line_iter = line_iter,
            status = "accepted",
            alpha = alpha,
            norm_before = current_norm,
            norm_after = NormPnEIC,
            SummEIC = SummEIC,
            TargetTime = TargetTime,
            TargetEvent = TargetEvent,
            EICStopRule = EICStopRule,
            EICStopAbsTol = EICStopAbsTol
          )

          # Track momentum
          prev_alpha <- alpha
          if (!is.null(prev_norm) && !is.null(prev_reduction_ratio)) {
            local_ratio <- if (expected_reduction > 0) {
              actual_reduction / expected_reduction
            } else {
              Inf
            }
            if (local_ratio > prev_reduction_ratio) {
              consecutive_gains <- consecutive_gains + 1
            } else {
              consecutive_gains <- 0
            }
            prev_reduction_ratio <- local_ratio
          } else {
            # first time we track a ratio
            prev_reduction_ratio <- if (expected_reduction > 0) {
              actual_reduction / expected_reduction
            } else {
              Inf
            }
          }
          prev_norm <- current_norm
        }
      } # end while line_iter

      # If no improvement found after max_line_iterations => alpha got super small
      if (!improvement_found) {
        UpdateTrace <- appendUpdateTrace(
          trace = UpdateTrace,
          method = Method,
          step = StepNum,
          line_iter = max_line_iterations,
          status = "no_accepted_step",
          alpha = alpha,
          norm_before = current_norm,
          norm_after = NormPnEIC,
          SummEIC = SummEIC,
          TargetTime = TargetTime,
          TargetEvent = TargetEvent,
          EICStopRule = EICStopRule,
          EICStopAbsTol = EICStopAbsTol
        )
        if (Verbose) cat("No suitable alpha found => end iteration with minimal step.\n")
        consecutive_gains <- 0
      }

      # check convergence
      OneStepStop <- targetOneStepStop(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)

      if (all(sapply(OneStepStop[["check"]], isTRUE))) {
        if (Verbose) cat("All equations converged!\n")
        attr(Estimates, "TmleConverged") <- list(converged = TRUE, step = StepNum)
        attr(Estimates, "NormPnEICs")    <- c(NormPnEICs, NormPnEIC)
        attr(Estimates, "TmleUpdateTrace") <- UpdateTrace
        return(Estimates)
      }

      StepNum <- StepNum + 1
      NormPnEICs <- c(NormPnEICs, NormPnEIC)
    } # end while StepNum

    # If we exit the loop without convergence
    warning("TMLE with Adaptive + rollback has not converged by step ", MaxUpdateIter,
            " - Estimates may not have the desired asymptotic properties")
    attr(Estimates, "TmleConverged") <- list(converged = FALSE, step = StepNum)
    attr(Estimates, "NormPnEICs")    <- NormPnEICs
    attr(Estimates, "TmleUpdateTrace") <- UpdateTrace
    return(Estimates)

    ###############################################################################
    # 3) "coordinated" approach (Gauss–Seidel coordinate descent)
    ###############################################################################
  } else if (Method == "coordinated") {
    if (Verbose) cat("Using coordinate-wise targeted update approach\n")
    target_coords <- expand.grid(
      Trt   = names(Estimates),
      Event = TargetEvent,
      Time  = TargetTime
    )
    if (Verbose) cat("Total equations to solve:", nrow(target_coords), "\n")

    StepNum <- 1
    global_converged <- FALSE

    while (StepNum <= MaxUpdateIter && !global_converged) {
      if (Verbose) cat("Starting iteration", StepNum, "\n")

      for (i in seq_len(nrow(target_coords))) {
        coord       <- target_coords[i, ]
        curr_trt    <- coord$Trt
        curr_event  <- coord$Event
        curr_time   <- coord$Time

        if (Verbose) {
          cat("  Targeting eqn: Trt =", curr_trt,
              "Event =", curr_event, "Time =", curr_time, "\n")
        }

        curr_eic <- SummEIC[Trt == curr_trt & Event == curr_event & Time == curr_time, PnEIC]
        curr_stop <- makeOneStepStop(
          SummEIC = SummEIC[Trt == curr_trt & Event == curr_event & Time == curr_time],
          EICStopRule = EICStopRule,
          EICStopAbsTol = EICStopAbsTol
        )
        if (all(sapply(curr_stop[["check"]], isTRUE))) {
          if (Verbose) cat("    Already converged, skip\n")
          next
        }
        # local attempts
        WorkingEps       <- OneStepEps
        coord_converged  <- FALSE
        max_coord_iter   <- 10
        coord_iter       <- 0
        a_idx            <- which(names(Estimates) == curr_trt)

        while (!coord_converged && coord_iter < max_coord_iter) {
          coord_iter <- coord_iter + 1
          est_a      <- Estimates[[a_idx]]

          # SummEIC but only for this single coordinate
          coord_SummEIC <- data.table::copy(est_a[["SummEIC"]])
          coord_SummEIC[, PnEIC := 0]
          coord_SummEIC[Time == curr_time & Event == curr_event, PnEIC := curr_eic]

          NewHazards <- updateHazard(
            GStar          = attr(est_a[["PropScore"]], "g.star.obs"),
            Hazards        = est_a[["Hazards"]],
            TotalSurv      = est_a[["EvntFreeSurv"]],
            NuisanceWeight = est_a[["NuisanceWeight"]],
            EvalTimes      = EvalTimes,
            T.tilde        = T.tilde,
            Delta          = Delta,
            PnEIC          = coord_SummEIC,
            NormPnEIC      = abs(curr_eic),
            OneStepEps     = WorkingEps,
            TargetEvent    = curr_event,
            TargetTime     = curr_time
          )
          NewHazards <- lapply(NewHazards, function(hz) {
            if (anyNA(hz)) hz[is.na(hz) | is.nan(hz)] <- 0
            hz
          })
          NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(hz) exp(-cumsum(hz)))
          NewSurv[NewSurv < 1e-8 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-8

          NewIC <- getIC(
            GStar          = attr(est_a[["PropScore"]], "g.star.obs"),
            Hazards        = NewHazards,
            TotalSurv      = NewSurv,
            NuisanceWeight = est_a[["NuisanceWeight"]],
            TargetEvent    = curr_event,
            TargetTime     = curr_time,
            T.tilde        = T.tilde,
            Delta          = Delta,
            EvalTimes      = EvalTimes,
            GComp          = FALSE
          )
          localSumm <- summarizeIC(NewIC)
          new_eic   <- localSumm[Time == curr_time & Event == curr_event, PnEIC]

          if (abs(new_eic) >= abs(curr_eic)) {
            if (Verbose) cat("    No improvement, halving step size\n")
            WorkingEps <- WorkingEps / 2
            if (WorkingEps < 1e-6 * OneStepEps) {
              if (Verbose) cat("    Step size tiny, stop\n")
              break
            }
          } else {
            if (Verbose) cat("    EIC from", curr_eic, "to", new_eic, "\n")
            # commit
            Estimates[[a_idx]][["Hazards"]]      <- NewHazards
            Estimates[[a_idx]][["EvntFreeSurv"]] <- NewSurv

            # re‐compute EIC across entire dimension for this 'a'
            bigIC <- getIC(
              GStar          = attr(Estimates[[a_idx]][["PropScore"]], "g.star.obs"),
              Hazards        = Estimates[[a_idx]][["Hazards"]],
              TotalSurv      = Estimates[[a_idx]][["EvntFreeSurv"]],
              NuisanceWeight = Estimates[[a_idx]][["NuisanceWeight"]],
              TargetEvent    = TargetEvent,
              TargetTime     = TargetTime,
              T.tilde        = T.tilde,
              Delta          = Delta,
              EvalTimes      = EvalTimes,
              GComp          = FALSE
            )
            Estimates[[a_idx]][["SummEIC"]] <- summarizeIC(bigIC)
            Estimates[[a_idx]][["IC"]]      <- bigIC

            SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(xi) {
              cbind("Trt" = names(Estimates)[xi], Estimates[[xi]][["SummEIC"]])
            }))

            curr_eic <- SummEIC[Trt == curr_trt & Event == curr_event & Time == curr_time, PnEIC]
            curr_stop <- makeOneStepStop(
              SummEIC = SummEIC[Trt == curr_trt & Event == curr_event & Time == curr_time],
              EICStopRule = EICStopRule,
              EICStopAbsTol = EICStopAbsTol
            )
            if (all(sapply(curr_stop[["check"]], isTRUE))) {
              if (Verbose) cat("    coord eqn converged\n")
              coord_converged <- TRUE
            }
          }
        } # while local iteration
      } # end for coords

      # check global
      OneStepStop <- targetOneStepStop(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (Verbose) {
        cat("End of iteration", StepNum, "\n")
        printOneStepDiagnostics(
          OneStepStop,
          getNormPnEIC(targetSummEIC(SummEIC, TargetTime, TargetEvent)[, PnEIC])
        )
      }
      if (all(sapply(OneStepStop[["check"]], isTRUE))) {
        global_converged <- TRUE
        if (Verbose) cat("All coordinates converged!\n")
      }
      StepNum <- StepNum + 1
    }

    if (global_converged) {
      attr(Estimates, "TmleConverged") <- list(converged = TRUE, step = StepNum - 1)
    } else {
      warning("Coordinate-wise TMLE not converged by iteration ", MaxUpdateIter)
      attr(Estimates, "TmleConverged") <- list(converged = FALSE, step = StepNum - 1)
    }
    final_norm <- getNormPnEIC(targetSummEIC(SummEIC, TargetTime, TargetEvent)[, PnEIC])
    attr(Estimates, "NormPnEICs") <- c(NormPnEIC, final_norm)
    return(Estimates)

    ###############################################################################
    # 4) "jacobi" approach: parallel coordinate updates
    ###############################################################################
  } else if (Method == "jacobi") {
    if (Verbose) cat("Using Jacobi-style parallel updates for all coordinates\n")

    cl <- makeCluster(5)
    registerDoParallel(cl)

    target_coords <- expand.grid(
      Trt   = names(Estimates),
      Event = TargetEvent,
      Time  = TargetTime
    )
    if (Verbose) cat("Total coordinates to solve in each iteration:", nrow(target_coords), "\n")

    StepNum         <- 1
    global_converged<- FALSE

    while (StepNum <= MaxUpdateIter && !global_converged) {
      if (Verbose) cat("Starting Jacobi iteration", StepNum, "\n")

      # 1) Copy the current hazards so that parallel updates all read the "old" hazards
      oldEstimates <- deepcopyEstimates(Estimates)

      # 2) Parallel loop over all coordinates
      updates <- foreach(i = seq_len(nrow(target_coords)), .combine = rbind) %dopar% {
        coord       <- target_coords[i, ]
        curr_trt    <- coord$Trt
        curr_event  <- coord$Event
        curr_time   <- coord$Time

        oldSummEIC  <- oldEstimates[[curr_trt]][["SummEIC"]]
        curr_eic    <- oldSummEIC[Time == curr_time & Event == curr_event, PnEIC]
        curr_thresh <- oldSummEIC[Time == curr_time & Event == curr_event, `seEIC/(sqrt(n)log(n))`]

        HazChange <- NULL
        if (abs(curr_eic) > curr_thresh) {
          coord_SummEIC <- data.table::copy(oldSummEIC)
          coord_SummEIC[, PnEIC := 0]
          coord_SummEIC[Time == curr_time & Event == curr_event, PnEIC := curr_eic]

          newHazards <- updateHazard(
            GStar          = attr(oldEstimates[[curr_trt]][["PropScore"]], "g.star.obs"),
            Hazards        = oldEstimates[[curr_trt]][["Hazards"]],
            TotalSurv      = oldEstimates[[curr_trt]][["EvntFreeSurv"]],
            NuisanceWeight = oldEstimates[[curr_trt]][["NuisanceWeight"]],
            EvalTimes      = attr(oldEstimates, "Times"),
            T.tilde        = Data[[attr(Data, "EventTime")]],
            Delta          = Data[[attr(Data, "EventType")]],
            PnEIC          = coord_SummEIC,
            NormPnEIC      = abs(curr_eic),
            OneStepEps     = OneStepEps,
            TargetEvent    = curr_event,
            TargetTime     = curr_time
          )

          HazChange <- list(
            "trt"   = curr_trt,
            "event" = curr_event,
            "time"  = curr_time,
            "update_haz" = newHazards
          )
        }
        if (is.null(HazChange)) {
          HazChange <- list(
            "trt"   = curr_trt,
            "event" = curr_event,
            "time"  = curr_time,
            "update_haz" = NULL
          )
        }
        df <- data.frame(trt = curr_trt, event = curr_event, time = curr_time)
        attr(df, "hazList") <- HazChange[["update_haz"]]
        df
      }

      # 3) Single commit step
      for (row_i in seq_len(nrow(updates))) {
        this_row  <- updates[row_i, ]
        this_trt  <- as.character(this_row[["trt"]])
        hazList   <- attr(this_row, "hazList")
        if (!is.null(hazList)) {
          Estimates[[this_trt]][["Hazards"]] <- hazList
        }
      }

      # 4) Recompute SummEIC across all events/times for each trt
      for (tnm in names(Estimates)) {
        newSurv <- apply(Reduce(`+`, Estimates[[tnm]][["Hazards"]]), 2, function(hz) exp(-cumsum(hz)))
        newSurv[newSurv < 1e-8 | is.na(newSurv) | is.nan(newSurv)] <- 1e-8
        Estimates[[tnm]][["EvntFreeSurv"]] <- newSurv

        bigIC <- getIC(
          GStar          = attr(Estimates[[tnm]][["PropScore"]], "g.star.obs"),
          Hazards        = Estimates[[tnm]][["Hazards"]],
          TotalSurv      = Estimates[[tnm]][["EvntFreeSurv"]],
          NuisanceWeight = Estimates[[tnm]][["NuisanceWeight"]],
          TargetEvent    = TargetEvent,
          TargetTime     = TargetTime,
          T.tilde        = Data[[attr(Data, "EventTime")]],
          Delta          = Data[[attr(Data, "EventType")]],
          EvalTimes      = EvalTimes,
          GComp          = FALSE
        )
        Estimates[[tnm]][["SummEIC"]] <- summarizeIC(bigIC)
        Estimates[[tnm]][["IC"]]      <- bigIC
      }

      SummEIC <- do.call(rbind, lapply(seq_along(Estimates), function(xi) {
        cbind("Trt" = names(Estimates)[xi], Estimates[[xi]][["SummEIC"]])
      }))

      OneStepStop <- targetOneStepStop(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (Verbose) {
        cat("End of Jacobi iteration", StepNum, "\n")
        printOneStepDiagnostics(
          OneStepStop,
          getNormPnEIC(targetSummEIC(SummEIC, TargetTime, TargetEvent)[, PnEIC])
        )
      }
      if (all(sapply(OneStepStop[["check"]], isTRUE))) {
        global_converged <- TRUE
        if (Verbose) cat("All coordinates converged!\n")
      }
      StepNum <- StepNum + 1
    }

    # finalize
    if (!global_converged) {
      warning("Jacobi-style TMLE not converged by iteration ", MaxUpdateIter)
    }
    finalNorm <- getNormPnEIC(targetSummEIC(SummEIC, TargetTime, TargetEvent)[, PnEIC])
    attr(Estimates, "TmleConverged") <- list(converged = global_converged, step = StepNum - 1)
    attr(Estimates, "NormPnEICs")    <- c(NormPnEIC, finalNorm)
    return(Estimates)

  } else if (Method == "rootSolve") {

    if (Verbose) cat("Using hybrid approach: Stage 1 (mini-adaptive) => Stage 2 (nleqslv)\n")

    #---------------------------------------------------------------------------
    # Stage 1: a few small-step updates
    #---------------------------------------------------------------------------
    mini_iter <- min(5, MaxUpdateIter)
    for (st in seq_len(mini_iter)) {
      if (Verbose) cat("  [Stage1] Mini-adaptive step", st,
                       "with OneStepEps =", WorkingEps, "\n")

      # We'll do a single iteration of the "standard" small-step approach:
      newEsts <- lapply(Estimates, function(est.a) {
        NewHazards <- updateHazard(
          GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
          Hazards        = est.a[["Hazards"]],
          TotalSurv      = est.a[["EvntFreeSurv"]],
          NuisanceWeight = est.a[["NuisanceWeight"]],
          EvalTimes      = EvalTimes,
          T.tilde        = T.tilde,
          Delta          = Delta,
          PnEIC          = est.a[["SummEIC"]],
          NormPnEIC      = NormPnEIC,
          OneStepEps     = WorkingEps,
          TargetEvent    = TargetEvent,
          TargetTime     = TargetTime
        )
        # Fix any NA
        NewHazards <- lapply(NewHazards, function(hz) {
          if (anyNA(hz)) hz[is.na(hz) | is.nan(hz)] <- 0
          hz
        })
        # Recompute survival
        NewSurv <- apply(Reduce(`+`, NewHazards), 2, function(hz) exp(-cumsum(hz)))
        NewSurv[NewSurv < 1e-12 | is.na(NewSurv) | is.nan(NewSurv)] <- 1e-12

        # Recompute EIC
        NewIC <- getIC(
          GStar          = attr(est.a[["PropScore"]], "g.star.obs"),
          Hazards        = NewHazards,
          TotalSurv      = NewSurv,
          NuisanceWeight = est.a[["NuisanceWeight"]],
          TargetEvent    = TargetEvent,
          TargetTime     = TargetTime,
          T.tilde        = T.tilde,
          Delta          = Delta,
          EvalTimes      = EvalTimes,
          GComp          = FALSE
        )
        list(
          Hazards      = NewHazards,
          EvntFreeSurv = NewSurv,
          SummEIC      = summarizeIC(NewIC),
          IC           = NewIC
        )
      })

      # Evaluate new SummEIC / Norm
      NewSummEIC <- do.call(rbind, lapply(seq_along(newEsts), function(a) {
        cbind("Trt" = names(newEsts)[a], newEsts[[a]][["SummEIC"]])
      }))
      NewNormPnEIC <- getNormPnEIC(targetSummEIC(NewSummEIC, TargetTime, TargetEvent)[, PnEIC])

      if (NewNormPnEIC > NormPnEIC) {
        if (Verbose) cat("    Norm increased, halving step =>", WorkingEps / 2, "\n")
        WorkingEps <- WorkingEps / 2
        next
      }

      # Commit
      for (a in seq_along(Estimates)) {
        Estimates[[a]][["Hazards"]]      <- newEsts[[a]][["Hazards"]]
        Estimates[[a]][["EvntFreeSurv"]] <- newEsts[[a]][["EvntFreeSurv"]]
        Estimates[[a]][["SummEIC"]]      <- newEsts[[a]][["SummEIC"]]
        Estimates[[a]][["IC"]]           <- newEsts[[a]][["IC"]]
      }
      SummEIC   <- NewSummEIC
      NormPnEIC <- NewNormPnEIC
      NormPnEICs<- c(NormPnEICs, NormPnEIC)

      # Check if converged
      OneStepStop <- targetOneStepStop(
        SummEIC = SummEIC,
        TargetTime = TargetTime,
        TargetEvent = TargetEvent,
        EICStopRule = EICStopRule,
        EICStopAbsTol = EICStopAbsTol
      )

      if (Verbose) printOneStepDiagnostics(OneStepStop, NormPnEIC)
      if (all(sapply(OneStepStop[["check"]], isTRUE))) {
        if (Verbose) cat("**Within mini-adaptive, we have converged!**\n")
        attr(Estimates, "TmleConverged") <- list(converged = TRUE, step = st)
        attr(Estimates, "NormPnEICs")    <- NormPnEICs
        return(Estimates)
      }
    }

    #---------------------------------------------------------------------------
    # Stage 2: run nleqslv / Broyden to solve for eps
    #---------------------------------------------------------------------------
    eps_start <- rep(0, k_dim)
    if (Verbose) cat("  [Stage2] Running nleqslv with method=Broyden\n")

    objectiveFn <- function(epsVec) {
      submodelAndEIC(Estimates, epsVec, Data, TargetEvent, TargetTime, Verbose)
    }

    root_res <- nleqslv::nleqslv(
      x = eps_start,
      fn = objectiveFn,
      method = "Broyden",
      control = list(stepmax = 1, maxit=50, ftol=1e-6, xtol=1e-6)
    )
    if (Verbose) {
      cat("nleqslv convergence code:", root_res$convergence,
          "message:", root_res$message, "\n")
    }

    final_eps <- root_res$x

    # 1) Apply final submodel to real 'Estimates'
    updatedRes <- applyFinalSubmodel(Estimates, final_eps, Data, TargetEvent, TargetTime, Verbose)

    # 2) Now unify SummEIC across arms to ensure consistent columns
    #    e.g. force Time x Event grid, add seEIC = 0 if missing, etc.
    #    Then do a final rbind() to check dimension
    arms <- names(updatedRes)
    AllTimeEvent <- data.table::CJ(
      Time  = TargetTime,
      Event = TargetEvent
    )
    for (arm in arms) {
      seicDT <- updatedRes[[arm]][["SummEIC"]]
      # Force all Time,Event combos
      seicDT <- merge(AllTimeEvent, seicDT, by = c("Time","Event"), all.x = TRUE)
      # Fill missing columns or NAs
      if (!("seEIC" %in% colnames(seicDT))) {
        seicDT[, seEIC := 0]
      }
      seicDT[is.na(PnEIC),  PnEIC  := 0]
      seicDT[is.na(seEIC),  seEIC  := 0]
      updatedRes[[arm]][["SummEIC"]] <- seicDT
      if (Verbose) {
        cat(sprintf("  After unify: arm=%s SummEIC has %d rows, columns= %s\n",
                    arm, nrow(seicDT), paste(colnames(seicDT), collapse=",")))
      }
    }

    # 3) Final check of norm
    finalSummEIC <- do.call(rbind, lapply(seq_along(updatedRes), function(a) {
      cbind("Trt" = names(updatedRes)[a], updatedRes[[a]][["SummEIC"]])
    }))
    finalNorm <- getNormPnEIC(targetSummEIC(finalSummEIC, TargetTime, TargetEvent)[, PnEIC])
    if (Verbose) cat("Final norm of PnEIC after hybrid rootSolve:", finalNorm, "\n")

    # 4) Return with TmleConverged attribute
    attr(updatedRes, "TmleConverged") <- list(
      converged = (root_res$convergence == 1),
      step = paste0("Stage1=", mini_iter, ", Stage2=", root_res$iter)
    )

    attr(updatedRes, "NormPnEICs") <- c(NormPnEICs, finalNorm)
    return(updatedRes)

  } else {
    # if unrecognized method, fallback
    warning("Unrecognized Method '", Method, "'. Using standard universal LFM approach.")
    return(doTmleUpdate(
      Estimates, SummEIC, Data, TargetEvent, TargetTime,
      MaxUpdateIter, OneStepEps, NormPnEIC, Verbose,
      Method = "standard"
    ))
  }
}


###############################################################################
# Hazard update function
###############################################################################

#' Update hazards based on clever covariate and PnEIC
#'
#' @param GStar numeric vector of star probabilities
#' @param Hazards list of hazard matrices
#' @param TotalSurv matrix of survival probabilities
#' @param NuisanceWeight matrix of nuisance weights
#' @param EvalTimes numeric vector of evaluation times
#' @param T.tilde numeric vector of observed event times
#' @param Delta numeric vector of observed event types
#' @param PnEIC data.table with PnEIC values
#' @param NormPnEIC numeric norm of PnEIC
#' @param OneStepEps numeric step size for updates
#' @param TargetEvent numeric vector of target events
#' @param TargetTime numeric vector of target times
#'
#' @return list of updated hazard matrices
#' @keywords internal
updateHazard <- function(GStar, Hazards, TotalSurv, NuisanceWeight, EvalTimes, T.tilde,
                         Delta, PnEIC, NormPnEIC, OneStepEps,
                         TargetEvent, TargetTime) {
  eps <- Time <- Event <- NULL
  Iterative <- FALSE
  GStar <- as.numeric(unlist(GStar))

  if (min(TotalSurv) == 0) {
    stop("Some individual's survival probability = 0, which can make the clever covariate blow up.")
  }
  if (Iterative) {
    warning("Iterative TMLE not yet implemented. Performing one-step TMLE instead.")
  }
  if (!is.numeric(NormPnEIC) || length(NormPnEIC) != 1 ||
      !is.finite(NormPnEIC) || NormPnEIC <= .Machine$double.eps) {
    return(lapply(Hazards, function(haz.al) {
      attr(haz.al, "j") <- attr(haz.al, "j")
      haz.al
    }))
  }

  # We do a single-step universal LFM update
  NewHazards <- lapply(Hazards, function(haz.al) {
    l <- attr(haz.al, "j")
    update.l <- Reduce("+", lapply(TargetEvent, function(j) {
      F.j.t <- apply(Hazards[[as.character(j)]] * TotalSurv, 2, cumsum)
      Reduce("+", lapply(TargetTime, function(tau) {
        ClevCov <- h.FS <- matrix(0, nrow = nrow(F.j.t), ncol = ncol(F.j.t))

        idx_le_tau <- EvalTimes <= tau
        if (any(idx_le_tau)) {
          F_tau <- F.j.t[EvalTimes == tau, ]
          F_tau_mat <- matrix(F_tau, ncol = ncol(F.j.t), nrow = sum(idx_le_tau), byrow = TRUE)

          h.FS[idx_le_tau, ] <-
            (F_tau_mat - F.j.t[idx_le_tau, , drop = FALSE]) /
            TotalSurv[idx_le_tau, , drop = FALSE]

          ClevCov[idx_le_tau, ] <-
            getCleverCovariate(
              GStar          = GStar,
              NuisanceWeight = NuisanceWeight[idx_le_tau, , drop = FALSE],
              hFS            = h.FS[idx_le_tau, , drop = FALSE],
              LeqJ           = as.integer(l == j)
            )
        }
        eic_contrib <- PnEIC[Time == tau & Event == j, PnEIC]
        if (length(eic_contrib) < 1) eic_contrib <- 0
        ClevCov * eic_contrib
      }))
    }))
    newhaz.al <- haz.al * exp(update.l * OneStepEps / NormPnEIC)
    attr(newhaz.al, "j") <- l
    newhaz.al
  })

  return(NewHazards)
}

###############################################################################
# NEW FUNCTIONS FOR THE ROOTSOLVE APPROACH
###############################################################################

#' submodelAndEIC
#'
#' Evaluate the empirical mean EIC (length k) for a guessed epsVec, without
#' modifying your real 'Estimates'.
#'
#' @param Estimates list of current estimates (hazards, surv, etc.)
#' @param epsVec numeric vector (length = k)
#' @param Data data.table containing the original data
#' @param TargetEvent numeric vector
#' @param TargetTime numeric vector
#' @param Verbose boolean
#'
#' @return numeric vector of length k, the Pn(EIC) values
#' @keywords internal
submodelAndEIC <- function(Estimates, epsVec, Data, TargetEvent, TargetTime, Verbose = FALSE) {
  Trt <- Time <- Event <- NULL

  # 1) Make an ephemeral copy
  localEst <- deepcopyEstimates(Estimates)

  # 2) Apply the submodel with the given epsVec
  #    (You can define your own logic for how epsVec modifies hazards)
  localEst <- updateHazardsWithEps(localEst, epsVec, Data, TargetEvent, TargetTime, Verbose)

  # 3) Recompute EIC on this ephemeral submodel
  for (arm in names(localEst)) {
    newSurv <- apply(Reduce(`+`, localEst[[arm]][["Hazards"]]), 2, function(hz) exp(-cumsum(hz)))
    newSurv[newSurv < 1e-12] <- 1e-12

    newIC <- getIC(
      GStar          = attr(localEst[[arm]][["PropScore"]], "g.star.obs"),
      Hazards        = localEst[[arm]][["Hazards"]],
      TotalSurv      = newSurv,
      NuisanceWeight = localEst[[arm]][["NuisanceWeight"]],
      TargetEvent    = TargetEvent,
      TargetTime     = TargetTime,
      T.tilde        = Data[[attr(Data, "EventTime")]],
      Delta          = Data[[attr(Data, "EventType")]],
      EvalTimes      = attr(localEst, "Times"),
      GComp          = FALSE
    )
    localEst[[arm]][["SummEIC"]] <- summarizeIC(newIC)
    localEst[[arm]][["IC"]]      <- newIC
  }

  # 4) Gather EIC in consistent order
  SummAll <- do.call(rbind, lapply(names(localEst), function(arm) {
    cbind("Trt" = arm, localEst[[arm]][["SummEIC"]])
  }))
  SummAll <- SummAll[Time %in% TargetTime & Event %in% TargetEvent]
  setorder(SummAll, Trt, Event, Time)

  eicVector <- SummAll[["PnEIC"]]
  if (Verbose) {
    cat("submodelAndEIC returning vector of length", length(eicVector),
        "; norm=", signif(sqrt(sum(eicVector^2)), 5), "\n")
  }
  return(eicVector)
}

#' updateHazardsWithEps
#'
#' Illustrative helper function that modifies local hazards based on epsVec.
#' You must define how each eps entry maps to a hazard shift or scale.
#'
#' @param localEst ephemeral copy of the Estimates list
#' @param epsVec numeric vector (length = k)
#' @param Data data.table
#' @param TargetEvent numeric
#' @param TargetTime numeric
#' @param Verbose boolean
#'
#' @return The updated localEst (ephemeral)
#' @keywords internal
updateHazardsWithEps <- function(localEst, epsVec, Data, TargetEvent, TargetTime, Verbose=FALSE) {
  # Example approach: assume ordering in epsVec is (arm1, event1, time1), ...
  arms <- names(localEst)
  nA <- length(arms)
  nE <- length(TargetEvent)
  nT <- length(TargetTime)
  k_dim <- nA * nE * nT

  if (length(epsVec) != k_dim) {
    stop("epsVec length does not match number of coordinates = nArms*nEvents*nTimes")
  }

  idx <- 1
  for (arm in arms) {
    for (ev in TargetEvent) {
      for (tt in TargetTime) {
        eps_val <- epsVec[idx]
        idx <- idx + 1
        # For this example, multiply the hazard for event=ev by exp(eps_val) entirely:
        old_mat <- localEst[[arm]][["Hazards"]][[as.character(ev)]]
        localEst[[arm]][["Hazards"]][[as.character(ev)]] <- old_mat * exp(eps_val)

        if (Verbose) {
          cat("  updateHazardsWithEps => arm=", arm,
              "event=", ev, "time=", tt,
              "with eps=", eps_val, "\n")
        }
      }
    }
  }
  return(localEst)
}

#' applyFinalSubmodel
#'
#' Commits the final epsVec solution from nleqslv onto the real 'Estimates'.
#'
#' @param Estimates the real estimates list to be updated in-place
#' @param final_eps numeric vector from root solver
#' @param Data data.table
#' @param TargetEvent numeric
#' @param TargetTime numeric
#' @param Verbose boolean
#'
#' @return Updated \code{Estimates} object with final hazards + EIC
#' @keywords internal
applyFinalSubmodel <- function(Estimates, final_eps, Data, TargetEvent, TargetTime, Verbose=FALSE) {
  if (Verbose) cat("applyFinalSubmodel: applying final eps => updating hazards & EIC in-place\n")

  # ephemeral copy to confirm final hazard updates
  localEst <- deepcopyEstimates(Estimates)

  # same logic as submodelAndEIC, but we store changes in localEst
  localEst <- updateHazardsWithEps(localEst, final_eps, Data, TargetEvent, TargetTime, Verbose)

  # Recompute final EIC
  for (arm in names(localEst)) {
    newSurv <- apply(Reduce(`+`, localEst[[arm]][["Hazards"]]), 2, function(hz) exp(-cumsum(hz)))
    newSurv[newSurv < 1e-12] <- 1e-12

    newIC <- getIC(
      GStar          = attr(localEst[[arm]][["PropScore"]], "g.star.obs"),
      Hazards        = localEst[[arm]][["Hazards"]],
      TotalSurv      = newSurv,
      NuisanceWeight = localEst[[arm]][["NuisanceWeight"]],
      TargetEvent    = TargetEvent,
      TargetTime     = TargetTime,
      T.tilde        = Data[[attr(Data, "EventTime")]],
      Delta          = Data[[attr(Data, "EventType")]],
      EvalTimes      = attr(localEst, "Times"),
      GComp          = FALSE
    )
    localEst[[arm]][["SummEIC"]] <- summarizeIC(newIC)
    localEst[[arm]][["IC"]]      <- newIC
    localEst[[arm]][["EvntFreeSurv"]] <- newSurv
  }

  # commit from localEst back into real 'Estimates'
  for (arm in names(localEst)) {
    Estimates[[arm]] <- localEst[[arm]]
  }

  if (Verbose) cat("applyFinalSubmodel: final hazards + EIC stored in 'Estimates'\n")
  return(Estimates)
}
