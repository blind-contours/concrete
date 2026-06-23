#' Restricted mean time in favor of treatment (first-event competing risks)
#'
#' @description
#' The restricted mean time in favor of treatment (RMT-IF; Mao 2023) is the
#' average amount of time over \eqn{[0,\tau]} that a randomly chosen treated
#' patient spends in a strictly more favorable state than a randomly chosen
#' control. Unlike the (unitless) win ratio, it is reported in time units, which
#' clinicians find easier to interpret. It is the **time integral of the
#' instantaneous net benefit**,
#' \deqn{\mathrm{RMT\text{-}IF}(\tau) = \int_0^\tau \big[w(t) - l(t)\big]\,dt,}
#' where \eqn{w(t)} (resp. \eqn{l(t)}) is the probability that the treated patient
#' is in a better (resp. worse) state than the control at time \eqn{t}.
#'
#' For a **single event** this reduces \emph{exactly} to the restricted mean
#' survival time difference, \eqn{\mathrm{RMT\text{-}IF} = \mathrm{RMST}_1 -
#' \mathrm{RMST}_0}: with states {alive, with-event}, \eqn{w(t)-l(t) = S_1(t) -
#' S_0(t)}. For a **prioritized hierarchy** of competing events the state at
#' \eqn{t} is "event-free" or "first event was type \eqn{k} by \eqn{t}", ranked by
#' priority (highest priority = least favorable), so \eqn{w(t)} and \eqn{l(t)} are
#' bilinear in the two arms' cause-specific cumulative incidences. RMT-IF is
#' therefore a smooth functional of the targeted marginal curves, with the same
#' covariate-adjusted, doubly-robust, censoring-corrected influence-function
#' inference as the other estimands.
#'
#' This is the \strong{first-event} version: a higher-priority event that follows
#' a lower-priority one is not credited (it treats the events as competing risks),
#' exactly as in [getWinRatio()]. For the clinically intended death-priority
#' hierarchy that credits death after a non-fatal event, use the multistate
#' version (see [clinicalWinRatio()] and `?clinicalRMTIF`).
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Horizon numeric: the restriction horizon \eqn{\tau} (default: the
#'   largest target time). RMT-IF is reported in the time units of `Horizon`.
#' @param Intervention length-2 numeric: treatment and control indices; a positive
#'   RMT-IF favors `Intervention[1]`.
#' @param TargetEvent numeric: the event code, or an ordered vector of codes
#'   giving the priority hierarchy from highest to lowest (default: the first
#'   targeted event).
#' @param Signif numeric (default 0.05): alpha for confidence intervals and the
#'   two-sided Wald p-value.
#'
#' @return a `data.table` of class `"ConcreteOut"` with `RMT-IF` (net), the time
#'   in favor (\eqn{\int w}) and against (\eqn{\int l}), each with an
#'   influence-function CI; the net carries a p-value against the null of no
#'   difference. With a `Strata` fit the SEs are design-corrected.
#'
#' @seealso [getRMST()] (the single-event special case), [getWinRatio()],
#'   [getSimultaneousFamily()].
#' @export getRMTIF
#' @importFrom stats qnorm pnorm
getRMTIF <- function(ConcreteEst, Horizon = NULL, Intervention = c(1, 2),
                     TargetEvent = NULL, Signif = 0.05) {
  Estimator <- Event <- Time <- `Pt Est` <- ID <- IC <- NULL
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  if (length(Intervention) < 2)
    stop("Intervention must give the treatment and control indices, e.g. c(1, 2).")

  TargetTime <- attr(ConcreteEst, "TargetTime")
  if (is.null(TargetEvent)) TargetEvent <- attr(ConcreteEst, "TargetEvent")[1]
  TargetEvent <- as.numeric(TargetEvent)              # priority order, highest first
  if (anyDuplicated(TargetEvent))
    stop("TargetEvent must list each event code at most once (its priority order).")
  if (!all(TargetEvent %in% attr(ConcreteEst, "TargetEvent")))
    stop("All TargetEvent codes must have been targeted in doConcrete().")
  K <- length(TargetEvent)
  if (is.null(Horizon)) Horizon <- max(TargetTime)
  grid <- sort(unique(TargetTime[TargetTime <= Horizon]))
  if (!length(grid))
    stop("No target time is at or below Horizon (", Horizon, "); refit doConcrete() ",
         "with TargetTime values up to the horizon.")
  if (length(grid) < 2L)
    warning("RMT-IF is integrated over fewer than two target times; refit with a ",
            "denser TargetTime grid.")
  ## snap the horizon to the integration grid (the largest target time <= Horizon)
  ## and report THAT, like getRMST()/targetRMST() -- the integral cannot reach a
  ## horizon that is not on the TargetTime grid.
  if (length(grid) && Horizon > max(grid) + 1e-9) {
    message("RMT-IF: requested Horizon (", Horizon, ") is not on the TargetTime grid; ",
            "integrating to the last target time below it (", max(grid), "). ",
            "Add ", Horizon, " to TargetTime in doConcrete() to integrate to it exactly.")
    Horizon <- max(grid)
  }
  m <- length(grid)
  A1 <- names(ConcreteEst)[Intervention[1]]
  A0 <- names(ConcreteEst)[Intervention[2]]

  RisksObj <- getRisk(ConcreteEst, TargetTime = grid, TargetEvent = TargetEvent, GComp = FALSE)
  Risks <- data.table::as.data.table(RisksObj)[Estimator == "tmle"]
  ICdt  <- data.table::as.data.table(attr(RisksObj, "IC"))
  n <- length(unique(ICdt[["ID"]]))
  idsSorted <- sort(unique(ICdt[["ID"]]))

  Fcurve <- function(arm, k) {                        # length-m CIF on the grid
    r <- Risks[Risks[["Intervention"]] == arm & Event == k & Time %in% grid, ]
    data.table::setorder(r, Time); r[["Pt Est"]]
  }
  ICmat <- function(arm, k) {                         # m x n IC matrix, grid x sorted-ID
    d <- ICdt[ICdt[["Intervention"]] == arm & Event == k & Time %in% grid, ]
    mm <- data.table::dcast(d, Time ~ ID, value.var = "IC")
    data.table::setorder(mm, Time)
    as.matrix(mm[, setdiff(names(mm), "Time"), with = FALSE])
  }
  ## priority position p = 1..K (1 = highest priority = least favorable state)
  FX <- lapply(TargetEvent, function(k) Fcurve(A1, k))
  FY <- lapply(TargetEvent, function(k) Fcurve(A0, k))
  DX <- lapply(TargetEvent, function(k) ICmat(A1, k))
  DY <- lapply(TargetEvent, function(k) ICmat(A0, k))
  SX <- 1 - Reduce(`+`, FX); SY <- 1 - Reduce(`+`, FY)
  geTail <- function(F, k) if (k <= K) Reduce(`+`, F[k:K]) else numeric(m)   # sum_{m>=k}
  gtTail <- function(F, k) if (k < K)  Reduce(`+`, F[(k + 1):K]) else numeric(m) # sum_{a>k}

  ## instantaneous favor/against probabilities and the net benefit
  w <- SX * (1 - SY); l <- SY * (1 - SX)
  if (K >= 2) for (a in 2:K) for (b in seq_len(a - 1)) {  # a > b: position a less severe
    w <- w + FX[[a]] * FY[[b]]; l <- l + FY[[a]] * FX[[b]]
  }
  nb <- w - l

  wts <- trapezoidWeights(c(0, grid))[-1]               # integrate from 0 (nb(0)=0)
  rmtif <- sum(wts * nb); tFavor <- sum(wts * w); tAgainst <- sum(wts * l)

  ## per-subject influence functions (delta method, integrated over t)
  contract <- function(coef_list_X, coef_list_Y) {
    d <- numeric(n)
    for (k in seq_len(K)) {
      d <- d + colSums((wts * coef_list_X[[k]]) * DX[[k]]) +
               colSums((wts * coef_list_Y[[k]]) * DY[[k]])
    }
    d
  }
  ## gradients of nb (verified to give (-1,+1) at K=1 -> RMST difference IF):
  ##   d nb / d F_{X,k} = -sum_{m>=k} F_{Y,m} - S_Y - sum_{a>k} F_{Y,a}
  ##   d nb / d F_{Y,k} =  S_X + sum_{a>k} F_{X,a} + sum_{m>=k} F_{X,m}
  cX_nb <- lapply(seq_len(K), function(k) -geTail(FY, k) - SY - gtTail(FY, k))
  cY_nb <- lapply(seq_len(K), function(k)  SX + gtTail(FX, k) + geTail(FX, k))
  Dnet <- contract(cX_nb, cY_nb)
  ## time in favor (int w) and against (int l) for the decomposition rows
  cX_w <- lapply(seq_len(K), function(k) -geTail(FY, k))            # dw/dF_{X,k}
  cY_w <- lapply(seq_len(K), function(k)  SX + gtTail(FX, k))       # dw/dF_{Y,k}
  cX_l <- lapply(seq_len(K), function(k)  SY + gtTail(FY, k))       # dl/dF_{X,k}
  cY_l <- lapply(seq_len(K), function(k) -geTail(FX, k))            # dl/dF_{Y,k}
  Dfav <- contract(cX_w, cY_w); Dagn <- contract(cX_l, cY_l)

  StrataDT <- attr(ConcreteEst, "StrataDT")
  se <- function(d) {
    seStrat <- .strataSE(d, idsSorted, StrataDT)
    if (is.null(seStrat)) sqrt(mean(d^2) / n) else seStrat
  }
  z <- stats::qnorm(1 - Signif / 2)
  row <- function(lab, est, D, pv = TRUE) {
    s <- se(D)
    data.table::data.table(Estimand = lab, `Pt Est` = est, se = s,
                           `CI Low` = est - z * s, `CI Hi` = est + z * s,
                           pValue = if (pv) 2 * stats::pnorm(-abs(est / s)) else NA_real_)
  }
  Output <- data.table::rbindlist(list(
    row("RMT-IF", rmtif, Dnet),
    row("Time in favor", tFavor, Dfav, pv = FALSE),
    row("Time against", tAgainst, Dagn, pv = FALSE)))
  EventLab <- if (K == 1L) TargetEvent[1] else paste(TargetEvent, collapse = ">")
  Output[, `:=`(Intervention = paste0("[", A1, "] vs [", A0, "]"),
                Estimator = "tmle", Event = EventLab, Time = Horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event",
                                    "Time", "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif
  attr(Output, "Horizon") <- Horizon
  attr(Output, "Estimand") <- "RMT-IF"
  attr(Output, "Priority") <- TargetEvent
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  Output <- .attachFamily(Output, idsSorted, list(
    list(key = "RMT-IF", Estimand = "RMT-IF", Event = EventLab, Time = Horizon,
         Intervention = paste0("[", A1, "] vs [", A0, "]"),
         est = rmtif, se = se(Dnet), scale = "identity", ic = Dnet)))
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
