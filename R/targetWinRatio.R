#' Directly targeted win ratio, win odds, and net benefit
#'
#' @description
#' A targeted-maximum-likelihood update that fluctuates \strong{both arms'}
#' cause-specific hazards to solve the efficient-influence-function estimating
#' equations of the win probability and loss probability \strong{directly},
#' rather than plugging pointwise-targeted risk curves into the win functional
#' (which is what [getWinRatio()] does).
#'
#' For prioritized events \eqn{k = 1, \ldots, K} (first = highest priority) with
#' treated-arm cumulative incidences \eqn{G_k} and control-arm \eqn{H_k}, the
#' win probability over horizon \eqn{\tau} is
#' \deqn{P(win) = S_G(\tau)\sum_k H_k(\tau)
#'   + \sum_{a \ge 2} G_a(\tau) \sum_{k<a} H_k(\tau)
#'   + \sum_k \int_0^\tau [G_k(\tau) - G_k(t)] \, dH_k(t),}
#' and \eqn{P(loss)} swaps the arms. The clever covariate for each functional is
#' the gradient-weighted combination of the pointwise clever covariates over the
#' \strong{full event-time grid} (the chain rule analogue of [targetRMST()]'s
#' time-integrated clever covariate), so the win integral is evaluated on every
#' event time rather than the [getWinRatio()] target-time grid. Because the
#' gradient coefficients depend on both arms' current curves, they are
#' recomputed at every fluctuation step and the two arms are updated jointly
#' until \eqn{P_n D_{win} \approx P_n D_{loss} \approx 0}.
#'
#' Compared to the plug-in this removes the target-grid discretization of the
#' win integral and solves the functional's own estimating equation, which
#' matters most for sparse target grids and for the win odds / net benefit.
#' Derivation: `notes/target-win-ratio.md` in the source repository.
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Horizon numeric: the comparison horizon \eqn{\tau}. Defaults to the
#'   largest target time; snapped down to the nearest evaluation time.
#' @param Intervention length-2 numeric: `ConcreteEst` list indices of the
#'   treated and control interventions, in that order.
#' @param TargetEvent numeric: event types in priority order (first = highest,
#'   e.g. death). Defaults to the fitted target events in their given order.
#' @param MaxUpdateIter integer: maximum fluctuation steps.
#' @param OneStepEps numeric: initial step size for the fluctuation.
#' @param Signif numeric (default 0.05): alpha for CIs and two-sided p-values.
#' @param EICStopRule one of `"hybrid"` (default), `"relative"`, `"absolute"`:
#'   stopping rule for the (win, loss) estimating equations.
#' @param EICStopAbsTol numeric absolute tolerance for the `"absolute"` and
#'   `"hybrid"` rules. Defaults to `0.02 / sqrt(n)` (the functionals are
#'   probabilities, so the risk-scale default applies).
#' @param Verbose logical: print per-step convergence diagnostics.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
#'   net benefit, and win/loss/tie probabilities, each with influence-function
#'   CIs and (for the comparative statistics) p-values. Convergence status is in
#'   `attr(., "WRConverged")` and the step count in `attr(., "WRSteps")`. If the
#'   fit was built with `Strata` (see [formatArguments()]), the standard errors
#'   are corrected for the stratified / covariate-adaptive randomization design.
#'
#' @seealso [getWinRatio()] for the plug-in version; [targetRMST()] for the same
#'   construction applied to the restricted mean survival time.
#' @export targetWinRatio
targetWinRatio <- function(ConcreteEst, Horizon = NULL, Intervention = c(1, 2),
                           TargetEvent = NULL, MaxUpdateIter = 50L, OneStepEps = 0.1,
                           Signif = 0.05, EICStopRule = c("hybrid", "relative", "absolute"),
                           EICStopAbsTol = NULL, Verbose = FALSE) {
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  if (length(Intervention) < 2)
    stop("Intervention must give the treatment and control indices, e.g. c(1, 2).")
  EICStopRule <- match.arg(EICStopRule)

  EvalTimes <- attr(ConcreteEst, "Times")
  TargetTime <- attr(ConcreteEst, "TargetTime")
  T.tilde <- attr(ConcreteEst, "T.tilde")
  Delta <- attr(ConcreteEst, "Delta")
  if (is.null(TargetEvent)) TargetEvent <- attr(ConcreteEst, "TargetEvent")
  K <- length(TargetEvent)
  n <- length(T.tilde)
  if (is.null(EICStopAbsTol)) EICStopAbsTol <- 0.02 / sqrt(n)

  if (is.null(Horizon)) Horizon <- max(TargetTime)
  available <- EvalTimes[EvalTimes <= Horizon & EvalTimes > 0]
  if (!length(available)) stop("Horizon is below the first event time.")
  if (!(Horizon %in% EvalTimes)) {
    Horizon <- max(available)
    message("Horizon snapped to the nearest evaluation time at or below it: ", Horizon)
  }
  idx <- EvalTimes <= Horizon
  s <- EvalTimes[idx]
  M <- length(s)

  A1 <- names(ConcreteEst)[Intervention[1]]
  A0 <- names(ConcreteEst)[Intervention[2]]
  jks <- as.character(TargetEvent)

  ## --- static per-arm pieces: weights and counting-process jumps -------------
  ## (the compensator part of the martingale residual depends on the CURRENT
  ## hazards, so residuals are rebuilt inside pieces() at every iteration)
  armStatic <- function(est) {
    GStar <- as.numeric(unlist(attr(est[["PropScore"]], "g.star.obs")))
    Weight <- sweep(est[["NuisanceWeight"]][idx, , drop = FALSE], 2, GStar, `*`)
    causes <- names(est[["Hazards"]])
    NLdS <- lapply(causes, function(l) {
      NL <- matrix(0, M, n)
      for (i in which(Delta == as.numeric(l) & T.tilde <= Horizon))
        NL[which(s == T.tilde[i]), i] <- 1
      NL
    })
    names(NLdS) <- causes
    list(Weight = Weight, causes = causes, NLdS = NLdS)
  }
  st1 <- armStatic(ConcreteEst[[Intervention[1]]])
  st0 <- armStatic(ConcreteEst[[Intervention[2]]])
  armResid <- function(stArm, Hazards) {
    resid <- lapply(stArm$causes, function(l) {
      stArm$NLdS[[l]] - getHazLS(T_Tilde = T.tilde, EvalTimes = s,
                                 HazL = Hazards[[l]][idx, , drop = FALSE])
    })
    names(resid) <- stArm$causes
    resid
  }

  ## --- gradient coefficients of P(W beats L) on the grid ---------------------
  ## Wm, Lm: lists (length K) of length-M marginal CIF curves. Returns the win
  ## probability and the M x K coefficient matrices on the winner / loser curves.
  gradCoef <- function(Wm, Lm) {
    Wtau <- vapply(Wm, function(v) v[M], numeric(1))
    Ltau <- vapply(Lm, function(v) v[M], numeric(1))
    SWtau <- 1 - sum(Wtau)
    cW <- matrix(0, M, K); cL <- matrix(0, M, K)
    P <- SWtau * sum(Ltau) +
      (if (K >= 2) sum(vapply(2:K, function(a) Wm[[a]][M] * sum(Ltau[seq_len(a - 1)]),
                              numeric(1))) else 0)
    for (k in seq_len(K)) {
      dL <- diff(c(0, Lm[[k]]))
      P <- P + sum((Wm[[k]][M] - Wm[[k]]) * dL)
      ck <- -dL
      ck[M] <- ck[M] - (if (k < K) sum(Ltau[(k + 1):K]) else 0)
      cW[, k] <- ck
      cl <- numeric(M)
      if (M >= 2) cl[1:(M - 1)] <- Wm[[k]][2:M] - Wm[[k]][1:(M - 1)]
      cl[M] <- SWtau + (if (k < K) sum(Wtau[(k + 1):K]) else 0)
      cL[, k] <- cl
    }
    list(P = P, cW = cW, cL = cL)
  }

  revcs <- function(x) rev(cumsum(rev(x)))

  ## --- everything that changes with the hazards, for one iteration -----------
  ## c1, c0: M x K coefficient matrices on the treated / control curves of one
  ## functional. Returns the per-subject EIC and the per-arm clever-covariate
  ## pieces (base matrix + per-target-cause add-on) used for the fluctuation.
  funcPieces <- function(c1, c0, F1, F0, S1, S0, resid1, resid0) {
    armPiece <- function(cc, Fc, Surv, stt, resid) {
      A <- apply(cc, 2, revcs)                                   # M x K
      B <- matrix(0, M, n)
      for (k in seq_len(K))
        B <- B + apply(cc[, k] * Fc[[k]], 2, revcs) - A[, k] * Fc[[k]]
      base <- stt$Weight * (-B / Surv)
      mart <- Reduce(`+`, lapply(stt$causes, function(l) {
        H_l <- if (l %in% jks) base + stt$Weight * A[, match(l, jks)] else base
        colSums(H_l * resid[[l]])
      }))
      wp <- Reduce(`+`, lapply(seq_len(K), function(k) colSums(cc[, k] * Fc[[k]])))
      list(A = A, base = base, mart = mart, wp = wp)
    }
    p1 <- armPiece(c1, F1, S1, st1, resid1)
    p0 <- armPiece(c0, F0, S0, st0, resid0)
    WP <- p1$wp + p0$wp
    list(EIC = p1$mart + p0$mart + WP - mean(WP), p1 = p1, p0 = p0)
  }

  pieces <- function(H1, S1, H0, S0) {
    F1 <- lapply(jks, function(k) apply(H1[[k]] * S1, 2, cumsum)[idx, , drop = FALSE])
    F0 <- lapply(jks, function(k) apply(H0[[k]] * S0, 2, cumsum)[idx, , drop = FALSE])
    names(F1) <- names(F0) <- jks
    G <- lapply(F1, rowMeans); Hm <- lapply(F0, rowMeans)
    gw <- gradCoef(G, Hm)                       # win:  coeffs cW on G, cL on H
    gl <- gradCoef(Hm, G)                       # loss: coeffs cW on H, cL on G
    S1i <- S1[idx, , drop = FALSE]; S0i <- S0[idx, , drop = FALSE]
    resid1 <- armResid(st1, H1); resid0 <- armResid(st0, H0)
    win <- funcPieces(gw$cW, gw$cL, F1, F0, S1i, S0i, resid1, resid0)
    loss <- funcPieces(gl$cL, gl$cW, F1, F0, S1i, S0i, resid1, resid0)
    list(Pwin = gw$P, Ploss = gl$P, win = win, loss = loss)
  }

  summ <- function(p) {
    pn <- c(win = mean(p$win$EIC), loss = mean(p$loss$EIC))
    se <- c(win = sqrt(mean(p$win$EIC^2)), loss = sqrt(mean(p$loss$EIC^2)))
    thresh_scale <- 1 / (sqrt(n) * log(n))
    thresh <- switch(EICStopRule,
                     relative = se * thresh_scale,
                     absolute = rep(EICStopAbsTol, 2L),
                     hybrid = pmax(se * thresh_scale, EICStopAbsTol))
    list(PnEIC = pn, norm = sqrt(sum(pn^2)), converged = all(abs(pn) <= thresh))
  }

  ## --- joint targeting loop over both arms -----------------------------------
  Haz1 <- ConcreteEst[[Intervention[1]]][["Hazards"]]
  Sur1 <- ConcreteEst[[Intervention[1]]][["EvntFreeSurv"]]
  Haz0 <- ConcreteEst[[Intervention[2]]][["Hazards"]]
  Sur0 <- ConcreteEst[[Intervention[2]]][["EvntFreeSurv"]]

  p <- pieces(Haz1, Sur1, Haz0, Sur0)
  stt <- summ(p)
  eps <- OneStepEps
  step <- 0L
  while (step < MaxUpdateIter && !stt$converged) {
    updArm <- function(Hazards, fp_win, fp_loss, stArm) {
      newHaz <- lapply(stArm$causes, function(l) {
        Hw <- if (l %in% jks) fp_win$base + stArm$Weight * fp_win$A[, match(l, jks)] else fp_win$base
        Hl <- if (l %in% jks) fp_loss$base + stArm$Weight * fp_loss$A[, match(l, jks)] else fp_loss$base
        d <- matrix(0, length(EvalTimes), n)
        d[idx, ] <- Hw * stt$PnEIC[["win"]] + Hl * stt$PnEIC[["loss"]]
        h <- Hazards[[l]] * exp(eps * d / stt$norm)
        h[!is.finite(h)] <- 0
        attr(h, "j") <- attr(Hazards[[l]], "j")
        h
      })
      names(newHaz) <- stArm$causes
      newSurv <- apply(Reduce(`+`, newHaz), 2, function(hz) exp(-cumsum(hz)))
      newSurv[newSurv < 1e-12 | !is.finite(newSurv)] <- 1e-12
      list(Haz = newHaz, Surv = newSurv)
    }
    new1 <- updArm(Haz1, p$win$p1, p$loss$p1, st1)
    new0 <- updArm(Haz0, p$win$p0, p$loss$p0, st0)
    np <- pieces(new1$Haz, new1$Surv, new0$Haz, new0$Surv)
    nst <- summ(np)
    if (!is.finite(nst$norm) || nst$norm >= stt$norm) {
      eps <- eps / 2
      if (eps < 1e-8 * OneStepEps) break
      next
    }
    Haz1 <- new1$Haz; Sur1 <- new1$Surv; Haz0 <- new0$Haz; Sur0 <- new0$Surv
    p <- np; stt <- nst
    step <- step + 1L
    if (Verbose)
      cat(sprintf("  [WR] step %d  ||PnEIC||=%.3g  (win %.3g, loss %.3g)\n",
                  step, stt$norm, stt$PnEIC[["win"]], stt$PnEIC[["loss"]]))
  }
  if (!stt$converged)
    warning("Direct win-ratio targeting did not fully converge ",
            "(||PnEIC|| = ", signif(stt$norm, 3), " after ", step, " steps). ",
            "Increase MaxUpdateIter or inspect data support.")

  ## --- assemble output on the updated curves ---------------------------------
  Pwin <- p$Pwin; Ploss <- p$Ploss
  Dwin <- p$win$EIC; Dloss <- p$loss$EIC
  Ptie <- max(0, 1 - Pwin - Ploss)
  Dtie <- -(Dwin + Dloss)

  ## strata-corrected SEs when available (ICs are in data row order)
  StrataDT <- attr(ConcreteEst, "StrataDT")
  if (!is.null(StrataDT))
    StrataDT <- data.table::data.table(ID = seq_len(n), A = StrataDT[["A"]],
                                       S = StrataDT[["S"]])
  seIF <- function(d) {
    seStrat <- .strataSE(d, seq_len(n), StrataDT)
    if (is.null(seStrat)) sqrt(mean(d^2) / n) else seStrat
  }

  z <- stats::qnorm(1 - Signif / 2)
  ratioRow <- function(label, num, den, Dnum, Dden) {
    est <- num / den
    Dlog <- Dnum / num - Dden / den
    sl <- seIF(Dlog)
    data.table::data.table(Estimand = label, `Pt Est` = est, se = est * sl,
                           `CI Low` = est * exp(-z * sl), `CI Hi` = est * exp(z * sl),
                           pValue = 2 * stats::pnorm(-abs(log(est) / sl)))
  }
  probRow <- function(label, pp, D)
    data.table::data.table(Estimand = label, `Pt Est` = pp, se = seIF(D),
                           `CI Low` = pp - z * seIF(D), `CI Hi` = pp + z * seIF(D),
                           pValue = NA_real_)
  Dnb <- Dwin - Dloss; nb <- Pwin - Ploss
  Output <- data.table::rbindlist(list(
    ratioRow("Win Ratio", Pwin, Ploss, Dwin, Dloss),
    ratioRow("Win Odds", Pwin + Ptie / 2, Ploss + Ptie / 2,
             Dwin + Dtie / 2, Dloss + Dtie / 2),
    data.table::data.table(Estimand = "Net Benefit", `Pt Est` = nb, se = seIF(Dnb),
                           `CI Low` = nb - z * seIF(Dnb), `CI Hi` = nb + z * seIF(Dnb),
                           pValue = 2 * stats::pnorm(-abs(nb / seIF(Dnb)))),
    probRow("P(win)", Pwin, Dwin),
    probRow("P(loss)", Ploss, Dloss),
    probRow("P(tie)", Ptie, Dtie)))
  EventLab <- if (K == 1L) TargetEvent[1] else paste(TargetEvent, collapse = ">")
  Output[, `:=`(Intervention = paste0("[", A1, "] vs [", A0, "]"),
                Estimator = "tmle", Event = EventLab, Time = Horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event",
                                    "Time", "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif
  attr(Output, "Horizon") <- Horizon
  attr(Output, "Estimand") <- "Win Ratio"
  attr(Output, "Priority") <- TargetEvent
  attr(Output, "Targeting") <- "direct"
  attr(Output, "WRConverged") <- stt$converged
  attr(Output, "WRSteps") <- step
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
