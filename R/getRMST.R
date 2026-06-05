#' Restricted mean survival time and cause-specific life-years lost
#'
#' @description
#' `getRMST()` turns a fitted `"ConcreteEst"` object into restricted mean
#' survival time (RMST) and cause-specific life-years-lost (LYL) estimands, which
#' are collapsible, clinically interpretable summaries that regulators
#' increasingly prefer to a hazard ratio. Both are linear functionals of the
#' cumulative-incidence curves that `concrete` already targets, so their
#' influence functions are time-integrals of the per-subject influence functions
#' of the absolute risks. The integral is taken over the **target times the model
#' was fit on**, so request a reasonably dense `TargetTime` grid in
#' `formatArguments()` for an accurate RMST.
#'
#' - **RMST** (event-free): \eqn{\int_0^\tau S(t)\,dt}, the mean amount of
#'   follow-up time spent free of all events up to the horizon \eqn{\tau}. Only
#'   returned when every event type in the data was targeted.
#' - **Life-years lost** to cause \eqn{j}: \eqn{\int_0^\tau F_j(t)\,dt}, the mean
#'   time lost to cause \eqn{j} by \eqn{\tau}.
#'
#' @param ConcreteEst a `"ConcreteEst"` object returned by [doConcrete()].
#' @param Horizon numeric: the restriction horizon \eqn{\tau}. Defaults to the
#'   largest target time. Snapped to the nearest target time at or below it.
#' @param Intervention numeric (default `seq_along(ConcreteEst)`): which
#'   interventions to summarize. For contrasts the first two are treated as
#'   treatment and control.
#' @param Contrasts logical: also return RMST / LYL differences between the first
#'   two interventions.
#' @param Signif numeric (default 0.05): alpha for two-sided confidence intervals
#'   and two-sided Wald p-values.
#' @param NIMargin numeric (optional): a non-inferiority margin for the contrast
#'   estimands. When supplied, a one-sided non-inferiority assessment is added.
#' @param NIDirection one of `"lower"` or `"upper"`: which side of the margin is
#'   "non-inferior". Use `"lower"` when a larger value is better (e.g. an RMST
#'   difference) and `"upper"` when a smaller value is better.
#'
#' @return a `data.table` of class `"ConcreteOut"` with point estimates,
#'   influence-function standard errors, confidence intervals, and p-values.
#'
#' @seealso [targetRMST()] for the directly-targeted version that fluctuates the
#'   hazards for the RMST estimating equation instead of integrating
#'   pointwise-targeted risks; [getOutput()] for absolute risks, risk
#'   differences, and risk ratios.
#' @export getRMST
#' @importFrom stats qnorm pnorm
getRMST <- function(ConcreteEst, Horizon = NULL, Intervention = seq_along(ConcreteEst),
                    Contrasts = TRUE, Signif = 0.05, NIMargin = NULL,
                    NIDirection = c("lower", "upper")) {
  Intervention.col <- Estimand <- Estimator <- Event <- Time <- `Pt Est` <-
    se <- ID <- IC <- w <- NULL
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  NIDirection <- match.arg(NIDirection)

  TargetTime <- attr(ConcreteEst, "TargetTime")
  TargetEvent <- attr(ConcreteEst, "TargetEvent")
  AllEvents <- setdiff(sort(unique(attr(ConcreteEst, "Delta"))), 0)

  if (is.null(Horizon)) Horizon <- max(TargetTime)
  if (!is.numeric(Horizon) || length(Horizon) != 1 || !is.finite(Horizon))
    stop("Horizon must be a single finite numeric value.")
  available <- TargetTime[TargetTime <= Horizon]
  if (!length(available))
    stop("Horizon is below the smallest target time; nothing to integrate.")
  if (!(Horizon %in% TargetTime)) {
    Horizon <- max(available)
    message("Horizon snapped to the nearest target time at or below it: ", Horizon)
  }
  grid <- sort(unique(TargetTime[TargetTime <= Horizon]))
  if (length(grid) < 2L)
    warning("RMST is integrated over fewer than two target times; the integral ",
            "approximation will be crude. Refit with a denser TargetTime grid.")

  ## trapezoidal weights over {0, grid}; the value at t = 0 is 0 for every curve
  wts <- trapezoidWeights(c(0, grid))
  wdt <- data.table::data.table(Time = c(0, grid), w = wts)[Time > 0]

  ## risk curve and per-subject influence functions at the grid times
  RisksObj <- getRisk(ConcreteEst, TargetTime = grid, TargetEvent = TargetEvent,
                      GComp = FALSE)
  Risks <- data.table::as.data.table(RisksObj)[Estimator == "tmle"]
  ICdt <- data.table::as.data.table(attr(RisksObj, "IC"))
  n <- length(unique(ICdt[["ID"]]))

  ## integrate the curve -> life-years lost point estimates
  lyl <- merge(Risks[, list(Intervention, Event, Time, `Pt Est`)], wdt, by = "Time")
  lyl <- lyl[, list("Pt Est" = sum(w * `Pt Est`)), by = c("Intervention", "Event")]

  ## integrate the influence functions -> per-subject LYL influence functions
  ICw <- merge(ICdt, wdt, by = "Time")
  lylIC <- ICw[, list("IC" = sum(w * IC)), by = c("Intervention", "ID", "Event")]

  ## event-free RMST is only well defined when every event was targeted
  haveRMST <- setequal(TargetEvent, AllEvents)
  if (haveRMST) {
    rmst <- lyl[, list("Event" = -1, "Pt Est" = Horizon - sum(`Pt Est`)),
                by = "Intervention"]
    rmstIC <- lylIC[, list("Event" = -1, "IC" = -sum(IC)), by = c("Intervention", "ID")]
    lyl <- rbind(lyl, rmst, use.names = TRUE)
    lylIC <- rbind(lylIC, rmstIC, use.names = TRUE)
  } else {
    message("RMST (event-free survival) is only returned when every event type ",
            "is targeted; returning cause-specific life-years lost only.")
  }

  data.table::setorder(lylIC, Event, Intervention, ID)
  lylIC[, Time := Horizon]

  lylSe <- lylIC[, list("se" = sqrt(mean(IC^2) / n)), by = c("Intervention", "Event")]
  perArm <- merge(lyl, lylSe, by = c("Intervention", "Event"))
  perArm[, `:=`(Estimator = "tmle", Time = Horizon)]
  perArm[, Estimand := data.table::fifelse(Event == -1, "RMST", "Life Years Lost")]
  data.table::setattr(perArm, "IC", lylIC)

  Output <- data.table::copy(perArm)

  if (Contrasts && length(Intervention) >= 2) {
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
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  class(Output) <- union("ConcreteOut", class(Output))
  return(Output)
}

#' Trapezoidal integration weights for a sorted grid
#'
#' @param grid numeric, sorted, length >= 1
#' @return numeric weights `w` such that `sum(w * f(grid))` approximates the
#'   integral of `f` over `[min(grid), max(grid)]` by the trapezoid rule.
#' @keywords internal
trapezoidWeights <- function(grid) {
  m <- length(grid)
  if (m < 2L) return(rep(0, m))
  w <- numeric(m)
  w[1L] <- (grid[2L] - grid[1L]) / 2
  w[m] <- (grid[m] - grid[m - 1L]) / 2
  if (m > 2L) {
    interior <- 2:(m - 1L)
    w[interior] <- (grid[interior + 1L] - grid[interior - 1L]) / 2
  }
  w
}

#' Add Wald confidence intervals, p-values, and optional non-inferiority tests
#'
#' @param dt a data.table with `Estimand`, `Pt Est`, and `se` columns.
#' @param Signif numeric two-sided alpha.
#' @param NIMargin numeric or NULL non-inferiority margin.
#' @param NIDirection `"lower"` (larger is better) or `"upper"` (smaller better).
#' @return `dt` with `CI Low`, `CI Hi`, `pValue`, and (if `NIMargin`) the
#'   non-inferiority columns added.
#' @keywords internal
#' @importFrom stats qnorm pnorm
addWaldInference <- function(dt, Signif = 0.05, NIMargin = NULL,
                             NIDirection = "lower") {
  `Pt Est` <- se <- NULL
  dt <- data.table::as.data.table(dt)
  z <- stats::qnorm(1 - Signif / 2)
  pt <- dt[["Pt Est"]]
  sterr <- dt[["se"]]
  est <- dt[["Estimand"]]

  dt[, "CI Low" := pt - z * sterr]
  dt[, "CI Hi" := pt + z * sterr]

  ratio_est <- grepl("Ratio|Rel Risk", est)
  comparative <- grepl("Diff|Ratio|Rel Risk", est)
  null0 <- ifelse(ratio_est, 1, 0)
  zstat <- (pt - null0) / sterr
  pval <- 2 * stats::pnorm(-abs(zstat))
  pval[!comparative | !is.finite(sterr) | sterr <= 0] <- NA_real_
  dt[, "pValue" := pval]

  if (!is.null(NIMargin)) {
    if (!is.numeric(NIMargin) || length(NIMargin) != 1)
      stop("NIMargin must be a single numeric value.")
    cilo <- dt[["CI Low"]]
    cihi <- dt[["CI Hi"]]
    if (identical(NIDirection, "upper")) {
      # smaller is better: non-inferior if the whole CI sits below the margin
      noninf <- cihi < NIMargin
      nip <- stats::pnorm((pt - NIMargin) / sterr)
    } else {
      # larger is better: non-inferior if the whole CI sits above the margin
      noninf <- cilo > NIMargin
      nip <- stats::pnorm((NIMargin - pt) / sterr)
    }
    noninf[!comparative | !is.finite(sterr) | sterr <= 0] <- NA
    nip[!comparative | !is.finite(sterr) | sterr <= 0] <- NA_real_
    dt[, "NIMargin" := NIMargin]
    dt[, "NIpValue" := nip]
    dt[, "NonInferior" := noninf]
  }
  dt[]
}
