#' Directly target the restricted mean survival time
#'
#' @description
#' `targetRMST()` runs an additional one-step TMLE update that targets the
#' restricted mean survival time (RMST) and cause-specific life-years-lost (LYL)
#' estimands **directly**, rather than integrating pointwise-targeted absolute
#' risks (the approach used by [getRMST()]).
#'
#' The estimand \eqn{\mathrm{LYL}_j(\tau)=\int_0^\tau F_j(t)\,dt} is a smooth
#' time-average of the cumulative incidence, so its efficient influence function
#' is the time-integral of the pointwise risk influence functions and its clever
#' covariate is the time-integral of the pointwise clever covariates:
#' \deqn{H_{l,j,\tau}(s) = \frac{\pi^*(A\mid W)}{\pi(A\mid W) S_c(s^-)}
#'   \left[\mathbf 1(l=j)(\tau-s) - \frac{\int_s^\tau F_j(t)\,dt - (\tau-s)F_j(s)}{S(s)}\right].}
#' Targeting this single, well-conditioned functional avoids the dependence on a
#' dense `TargetTime` grid and tends to converge and cover better than the
#' pointwise approach in rare-event, competing-risk, and long-horizon settings.
#'
#' Starting from the hazards in `ConcreteEst` (already fit by [doConcrete()]),
#' `targetRMST()` fluctuates them along \eqn{H} until the empirical mean of the
#' RMST/LYL influence function is small, then reports the directly-targeted
#' estimates with influence-function standard errors, p-values, and (optionally)
#' a non-inferiority assessment.
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Horizon numeric: the restriction horizon \eqn{\tau}. Defaults to the
#'   largest target time; snapped to the nearest target time at or below it.
#' @param Intervention numeric (default `seq_along(ConcreteEst)`): interventions
#'   to summarize; the first two are treated as treatment and control.
#' @param TargetEvent numeric: event types to target. Defaults to the events the
#'   model was fit on.
#' @param MaxUpdateIter integer: maximum fluctuation steps.
#' @param OneStepEps numeric: initial step size for the fluctuation.
#' @param Signif numeric (default 0.05): alpha for confidence intervals and
#'   two-sided Wald p-values.
#' @param NIMargin,NIDirection optional non-inferiority margin and direction,
#'   passed to the contrast estimands (see [getOutput()]).
#' @param Verbose logical: print per-step convergence diagnostics.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the directly-targeted
#'   RMST / life-years-lost estimates. The per-arm convergence status is stored
#'   in `attr(., "RMSTConverged")`.
#'
#' @seealso [getRMST()] for the integrate-pointwise-risks version.
#' @export targetRMST
#' @importFrom stats qnorm pnorm
targetRMST <- function(ConcreteEst, Horizon = NULL, Intervention = seq_along(ConcreteEst),
                       TargetEvent = NULL, MaxUpdateIter = 50L, OneStepEps = 0.1,
                       Signif = 0.05, NIMargin = NULL, NIDirection = c("lower", "upper"),
                       EICStopRule = c("hybrid", "relative", "absolute"),
                       EICStopAbsTol = NULL, Verbose = FALSE) {
  Estimand <- Event <- Intervention.col <- NULL
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  NIDirection <- match.arg(NIDirection)
  EICStopRule <- match.arg(EICStopRule)

  EvalTimes <- attr(ConcreteEst, "Times")
  TargetTime <- attr(ConcreteEst, "TargetTime")
  AllEvents <- setdiff(sort(unique(attr(ConcreteEst, "Delta"))), 0)
  T.tilde <- attr(ConcreteEst, "T.tilde")
  Delta <- attr(ConcreteEst, "Delta")
  if (is.null(TargetEvent)) TargetEvent <- attr(ConcreteEst, "TargetEvent")
  # The fraction-scale estimand (LYL / Horizon) lives in [0, 1] like a risk, so
  # the rare-event absolute tolerance 0.02 / sqrt(n) applies on that scale.
  if (is.null(EICStopAbsTol)) EICStopAbsTol <- 0.02 / sqrt(length(attr(ConcreteEst, "T.tilde")))

  if (is.null(Horizon)) Horizon <- max(TargetTime)
  available <- EvalTimes[EvalTimes <= Horizon & EvalTimes > 0]
  if (!length(available)) stop("Horizon is below the first event time.")
  if (!(Horizon %in% EvalTimes)) {
    Horizon <- max(available)
    message("Horizon snapped to the nearest evaluation time at or below it: ", Horizon)
  }

  ## --- target each arm's hazards for the RMST/LYL functional --------------
  armResults <- lapply(seq_along(ConcreteEst), function(a) {
    rmstTargetArm(ConcreteEst[[a]], TargetEvent = TargetEvent, Horizon = Horizon,
                  T.tilde = T.tilde, Delta = Delta, EvalTimes = EvalTimes,
                  MaxUpdateIter = MaxUpdateIter, OneStepEps = OneStepEps,
                  EICStopRule = EICStopRule, EICStopAbsTol = EICStopAbsTol,
                  Verbose = Verbose, ArmName = names(ConcreteEst)[a])
  })
  names(armResults) <- names(ConcreteEst)
  converged <- vapply(armResults, function(x) x$converged, logical(1))
  if (!all(converged))
    warning("Direct RMST targeting did not fully converge for intervention(s): ",
            paste(names(converged)[!converged], collapse = ", "),
            ". Increase MaxUpdateIter or inspect data support.")

  n <- length(T.tilde)
  haveRMST <- setequal(TargetEvent, AllEvents)

  ## --- assemble per-arm LYL (+ event-free RMST) point estimates and ICs ----
  perArmList <- list()
  icList <- list()
  for (a in names(armResults)) {
    res <- armResults[[a]]
    # fraction-scale -> day-scale: LYL_j = Horizon * (fraction estimand)
    LYLday <- lapply(res$LYL, function(x) x * Horizon)
    ICday <- lapply(res$IC, function(x) x * Horizon)
    for (jj in as.character(TargetEvent)) {
      perArmList[[length(perArmList) + 1L]] <- data.table::data.table(
        Intervention = a, Estimand = "Life Years Lost", Estimator = "tmle",
        Event = as.numeric(jj), Time = Horizon,
        `Pt Est` = LYLday[[jj]], se = sqrt(mean(ICday[[jj]]^2) / n))
      icList[[length(icList) + 1L]] <- data.table::data.table(
        Intervention = a, ID = seq_len(n), Time = Horizon,
        Event = as.numeric(jj), IC = ICday[[jj]])
    }
    if (haveRMST) {
      rmstEst <- Horizon - sum(unlist(LYLday))
      rmstIC <- -Reduce(`+`, ICday)
      perArmList[[length(perArmList) + 1L]] <- data.table::data.table(
        Intervention = a, Estimand = "RMST", Estimator = "tmle",
        Event = -1, Time = Horizon,
        `Pt Est` = rmstEst, se = sqrt(mean(rmstIC^2) / n))
      icList[[length(icList) + 1L]] <- data.table::data.table(
        Intervention = a, ID = seq_len(n), Time = Horizon, Event = -1, IC = rmstIC)
    }
  }
  perArm <- data.table::rbindlist(perArmList)
  data.table::setattr(perArm, "IC", data.table::rbindlist(icList))

  Output <- data.table::copy(perArm)
  if (length(Intervention) >= 2) {
    A1 <- names(ConcreteEst)[Intervention[1]]
    A0 <- names(ConcreteEst)[Intervention[2]]
    contrast <- getRD(Risks = perArm, A1 = A1, A0 = A0, TargetTime = Horizon,
                      TargetEvent = unique(perArm[["Event"]]), GComp = FALSE)
    contrast[, Estimand := data.table::fifelse(Event == -1, "RMST Diff", "LYL Diff")]
    Output <- rbind(Output, contrast, use.names = TRUE, fill = TRUE)
  }

  data.table::setcolorder(
    Output, c("Intervention", "Estimand", "Estimator", "Event", "Time", "Pt Est", "se"))
  Output <- addWaldInference(Output, Signif = Signif, NIMargin = NIMargin,
                             NIDirection = NIDirection)
  data.table::setorderv(Output, c("Event", "Estimand", "Intervention"))
  attr(Output, "Signif") <- Signif
  attr(Output, "Horizon") <- Horizon
  attr(Output, "Estimand") <- "RMST"
  attr(Output, "Targeting") <- "direct"
  attr(Output, "RMSTConverged") <- converged
  attr(Output, "RMSTSteps") <- vapply(armResults, function(x) x$steps, numeric(1))
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  class(Output) <- union("ConcreteOut", class(Output))
  Output
}

#' Integrated (direct) RMST clever covariate and influence function for one arm
#'
#' Computes, for the current hazards, the per-subject life-years-lost influence
#' function for each target event and the integrated clever covariate pieces
#' needed to fluctuate the hazards. See [targetRMST()] for the formula.
#'
#' @keywords internal
rmstArmPieces <- function(GStar, Hazards, TotalSurv, NuisanceWeight, TargetEvent,
                          Horizon, T.tilde, Delta, EvalTimes) {
  GStar <- as.numeric(unlist(GStar))
  idx <- EvalTimes <= Horizon
  s <- EvalTimes[idx]
  M <- length(s)
  n <- ncol(TotalSurv)
  # Work on a rescaled time axis (fraction of the horizon) so the integrated
  # clever covariate is O(1) and the fluctuation is well conditioned, exactly as
  # for the pointwise risk. The fraction-scale estimand is LYL_j / Horizon; the
  # caller multiplies the estimate and its influence function back by Horizon.
  tau_minus_s <- (Horizon - s) / Horizon
  Surv <- TotalSurv[idx, , drop = FALSE]
  Weight <- sweep(NuisanceWeight[idx, , drop = FALSE], 2, GStar, `*`)

  causes <- names(Hazards)
  ## martingale residual (dN_l - Y dLambda_l) restricted to s <= Horizon, per cause
  resid <- lapply(causes, function(l) {
    NLdS <- matrix(0, M, n)
    for (i in which(Delta == as.numeric(l) & T.tilde <= Horizon))
      NLdS[which(s == T.tilde[i]), i] <- 1
    HazLS <- getHazLS(T_Tilde = T.tilde, EvalTimes = s,
                      HazL = Hazards[[l]][idx, , drop = FALSE])
    NLdS - HazLS
  })
  names(resid) <- causes

  IC <- list(); LYL <- list(); baseList <- list(); BList <- list()
  for (jj in as.character(TargetEvent)) {
    F.j <- apply(Hazards[[jj]] * TotalSurv, 2, cumsum)[idx, , drop = FALSE]
    ## G_j(s_m) = int_{s_m}^Horizon F_j(t) dt via reverse-cumsum of interval integrals
    if (M >= 2L) {
      dt <- diff(s) / Horizon
      Iint <- (F.j[-M, , drop = FALSE] + F.j[-1, , drop = FALSE]) * (dt / 2)
      Gj <- rbind(apply(Iint, 2, function(col) rev(cumsum(rev(col)))),
                  matrix(0, 1, n))
    } else {
      Gj <- matrix(0, M, n)
    }
    LYLcond <- Gj[1, ]                                   # int_0^Horizon F_j per subject
    base_j <- Weight * (-(Gj - tau_minus_s * F.j) / Surv)
    B_j <- Weight * tau_minus_s                          # the 1(l=j) part
    ## IC_j = sum_l colSums(H_{l,j} * resid_l) + (LYLcond - mean), H_{l,j}=base_j + 1(l=j)B_j
    mart <- Reduce(`+`, lapply(causes, function(l) {
      H_lj <- if (identical(l, jj)) base_j + B_j else base_j
      colSums(H_lj * resid[[l]])
    }))
    IC[[jj]] <- mart + LYLcond - mean(LYLcond)
    LYL[[jj]] <- mean(LYLcond)
    baseList[[jj]] <- base_j
    BList[[jj]] <- B_j
  }
  list(idx = idx, IC = IC, LYL = LYL, base = baseList, B = BList, causes = causes)
}

#' Fluctuate one arm's hazards to solve the RMST/LYL estimating equations
#' @keywords internal
rmstTargetArm <- function(est, TargetEvent, Horizon, T.tilde, Delta, EvalTimes,
                          MaxUpdateIter, OneStepEps, EICStopRule = "hybrid",
                          EICStopAbsTol = 0, Verbose, ArmName = "") {
  GStar <- attr(est[["PropScore"]], "g.star.obs")
  NuisanceWeight <- est[["NuisanceWeight"]]
  Hazards <- est[["Hazards"]]
  TotalSurv <- est[["EvntFreeSurv"]]
  n <- ncol(TotalSurv)
  thresh_scale <- 1 / (sqrt(n) * log(n))

  pieces <- function(H, S) rmstArmPieces(GStar, H, S, NuisanceWeight, TargetEvent,
                                         Horizon, T.tilde, Delta, EvalTimes)
  summ <- function(p) {
    pn <- vapply(as.character(TargetEvent), function(jj) mean(p$IC[[jj]]), numeric(1))
    se <- vapply(as.character(TargetEvent), function(jj) sqrt(mean(p$IC[[jj]]^2)), numeric(1))
    thresh <- switch(EICStopRule,
                     relative = se * thresh_scale,
                     absolute = rep(EICStopAbsTol, length(pn)),
                     hybrid = pmax(se * thresh_scale, EICStopAbsTol))
    list(PnEIC = pn, se = se, norm = sqrt(sum(pn^2)),
         converged = all(abs(pn) <= thresh))
  }

  p <- pieces(Hazards, TotalSurv)
  st <- summ(p)
  eps <- OneStepEps
  step <- 0L
  while (step < MaxUpdateIter && !st$converged) {
    ## fluctuation direction per cause: dir_l = sum_j H_{l,j} * PnEIC_j (on s<=Horizon)
    idx <- p$idx
    dirFull <- lapply(p$causes, function(l) {
      d <- matrix(0, length(EvalTimes), n)
      block <- Reduce(`+`, lapply(as.character(TargetEvent), function(jj) {
        H_lj <- if (identical(l, jj)) p$base[[jj]] + p$B[[jj]] else p$base[[jj]]
        H_lj * st$PnEIC[[jj]]
      }))
      d[idx, ] <- block
      d
    })
    names(dirFull) <- p$causes
    newHaz <- lapply(p$causes, function(l) {
      h <- Hazards[[l]] * exp(eps * dirFull[[l]] / st$norm)
      h[!is.finite(h)] <- 0
      attr(h, "j") <- attr(Hazards[[l]], "j")
      h
    })
    names(newHaz) <- p$causes
    newSurv <- apply(Reduce(`+`, newHaz), 2, function(hz) exp(-cumsum(hz)))
    newSurv[newSurv < 1e-12 | !is.finite(newSurv)] <- 1e-12

    np <- pieces(newHaz, newSurv)
    nst <- summ(np)
    if (!is.finite(nst$norm) || nst$norm >= st$norm) {
      eps <- eps / 2
      if (eps < 1e-8 * OneStepEps) break
      next
    }
    Hazards <- newHaz; TotalSurv <- newSurv; p <- np; st <- nst
    step <- step + 1L
    if (Verbose)
      cat(sprintf("  [RMST %s] step %d  ||PnEIC||=%.3g  converged=%s\n",
                  ArmName, step, st$norm, st$converged))
  }
  list(Hazards = Hazards, EvntFreeSurv = TotalSurv, IC = p$IC, LYL = p$LYL,
       converged = st$converged, steps = step, norm = st$norm)
}
