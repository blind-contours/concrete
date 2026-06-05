#' Covariate-adjusted restricted win ratio, win odds, and net benefit
#'
#' @description
#' For a single terminal time-to-event outcome, `getWinRatio()` estimates the
#' **restricted win ratio** and its relatives (win odds, net benefit) from the
#' covariate-adjusted, censoring-corrected counterfactual survival curves that
#' `concrete` targets. Comparing a random treated patient to a random control
#' patient over `[0, Horizon]`:
#'
#' - a treated patient **wins** if the control patient has the event first
#'   (before the treated patient and before the horizon),
#' - **loses** if the treated patient has the event first,
#' - **ties** otherwise (both event-free at the horizon).
#'
#' These probabilities are functionals of the marginal counterfactual survival
#' \eqn{\bar S_a} and cumulative incidence \eqn{\bar F_a = 1 - \bar S_a}:
#' \deqn{P(\text{win}) = \sum_{t_k \le \tau} \bar S_1(t_k)\, d\bar F_0(t_k), \qquad
#'       P(\text{loss}) = \sum_{t_k \le \tau} \bar S_0(t_k)\, d\bar F_1(t_k),}
#' and the win statistics are
#' \deqn{\text{win ratio} = \frac{P(\text{win})}{P(\text{loss})}, \quad
#'       \text{win odds} = \frac{P(\text{win}) + P(\text{tie})/2}{P(\text{loss}) + P(\text{tie})/2}, \quad
#'       \text{net benefit} = P(\text{win}) - P(\text{loss}).}
#' Because the win/loss probabilities are smooth functionals of the targeted
#' curves, their influence functions are weighted combinations of the per-subject
#' curve influence functions; the win ratio, win odds, and net benefit then
#' follow by the delta method, giving doubly-robust, covariate-adjusted inference.
#' Unlike the standard (unadjusted, censoring-sensitive) win ratio, this is
#' restricted to the horizon and corrects for censoring through the same
#' inverse-probability machinery as the rest of the package.
#'
#' This is the single-terminal-event version; a hierarchical / competing-risk
#' win ratio is planned. The integral is taken over the fitted target times, so
#' use a reasonably dense `TargetTime` grid.
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Horizon numeric: the restriction horizon \eqn{\tau} (default: the
#'   largest target time).
#' @param Intervention length-2 numeric: treatment and control indices.
#' @param TargetEvent numeric: the single terminal event code (default: the first
#'   targeted event).
#' @param Signif numeric (default 0.05): alpha for confidence intervals and
#'   p-values. Win ratio and win odds are inferred on the log scale.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
#'   net benefit, and the win/loss/tie probabilities, each with a CI and (for the
#'   comparative statistics) a p-value against the null of no difference.
#'
#' @seealso [getRMST()], [getOutput()]
#' @export getWinRatio
#' @importFrom stats qnorm pnorm
getWinRatio <- function(ConcreteEst, Horizon = NULL, Intervention = c(1, 2),
                        TargetEvent = NULL, Signif = 0.05) {
  Estimator <- Intervention.col <- Event <- Time <- `Pt Est` <- ID <- IC <- NULL
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  if (length(Intervention) < 2)
    stop("Intervention must give the treatment and control indices, e.g. c(1, 2).")

  TargetTime <- attr(ConcreteEst, "TargetTime")
  if (is.null(TargetEvent)) TargetEvent <- attr(ConcreteEst, "TargetEvent")[1]
  if (length(TargetEvent) != 1L)
    stop("getWinRatio() handles a single terminal event; pass one TargetEvent.")
  if (is.null(Horizon)) Horizon <- max(TargetTime)
  grid <- sort(unique(TargetTime[TargetTime <= Horizon]))
  if (length(grid) < 2L)
    warning("The win ratio is integrated over fewer than two target times; ",
            "refit with a denser TargetTime grid.")
  A1 <- names(ConcreteEst)[Intervention[1]]
  A0 <- names(ConcreteEst)[Intervention[2]]
  jj <- as.numeric(TargetEvent)

  ## marginal CIF curves and per-subject influence functions at the grid times
  RisksObj <- getRisk(ConcreteEst, TargetTime = grid, TargetEvent = TargetEvent, GComp = FALSE)
  Risks <- data.table::as.data.table(RisksObj)[Estimator == "tmle" & Event == jj]
  IC <- data.table::as.data.table(attr(RisksObj, "IC"))
  IC <- IC[Event == jj]
  n <- length(unique(IC[["ID"]]))

  Fcurve <- function(arm) {
    sel <- Risks[["Intervention"]] == arm
    r <- Risks[sel, ]
    data.table::setorder(r, Time)
    stats::setNames(c(0, r[["Pt Est"]]), as.character(c(0, grid)))   # F(0) = 0
  }
  ICmat <- function(arm) {
    sel <- IC[["Intervention"]] == arm
    d <- IC[sel, ]
    m <- data.table::dcast(d, Time ~ ID, value.var = "IC")
    data.table::setorder(m, Time)
    as.matrix(m[, setdiff(names(m), "Time"), with = FALSE])           # (#grid) x n, grid order
  }
  F1 <- Fcurve(A1); F0 <- Fcurve(A0)                                  # length m+1 (incl t0=0)
  S1 <- 1 - F1; S0 <- 1 - F0
  D1 <- ICmat(A1); D0 <- ICmat(A0)                                    # m x n (grid, no t0)
  m <- length(grid)
  dF1 <- diff(F1); dF0 <- diff(F0)                                    # increments at grid times, length m
  S1g <- S1[-1]; S0g <- S0[-1]                                        # survival at grid times, length m

  Pwin <- sum(S1g * dF0)
  Ploss <- sum(S0g * dF1)
  Ptie <- max(0, 1 - Pwin - Ploss)

  ## gradient coefficients (functional delta method); index k = 1..m over grid
  nextdF1 <- c(dF1[-1], NA); nextdF0 <- c(dF0[-1], NA)
  coef0 <- ifelse(seq_len(m) < m, nextdF1, S1g[m])   # dP(win)/dF0(t_k)
  coef1 <- ifelse(seq_len(m) < m, nextdF0, S0g[m])   # dP(loss)/dF1(t_k)

  ## per-subject influence functions of P(win) and P(loss)
  Dwin <- colSums((-dF0) * D1) + colSums(coef0 * D0)
  Dloss <- colSums((-dF1) * D0) + colSums(coef1 * D1)
  Dtie <- -(Dwin + Dloss)

  z <- stats::qnorm(1 - Signif / 2)
  se <- function(d) sqrt(mean(d^2) / n)
  ratioRow <- function(label, num, den, Dnum, Dden) {
    est <- num / den
    Dlog <- Dnum / num - Dden / den            # influence function of log(num/den)
    sl <- se(Dlog)
    pv <- 2 * stats::pnorm(-abs(log(est) / sl))
    data.table::data.table(Estimand = label, `Pt Est` = est, se = est * sl,
                           `CI Low` = est * exp(-z * sl), `CI Hi` = est * exp(z * sl),
                           pValue = pv)
  }
  probRow <- function(label, p, D)
    data.table::data.table(Estimand = label, `Pt Est` = p, se = se(D),
                           `CI Low` = p - z * se(D), `CI Hi` = p + z * se(D), pValue = NA_real_)
  Dnb <- Dwin - Dloss; nb <- Pwin - Ploss
  nbRow <- data.table::data.table(Estimand = "Net Benefit", `Pt Est` = nb, se = se(Dnb),
                                  `CI Low` = nb - z * se(Dnb), `CI Hi` = nb + z * se(Dnb),
                                  pValue = 2 * stats::pnorm(-abs(nb / se(Dnb))))

  Output <- data.table::rbindlist(list(
    ratioRow("Win Ratio", Pwin, Ploss, Dwin, Dloss),
    ratioRow("Win Odds", Pwin + Ptie / 2, Ploss + Ptie / 2, Dwin + Dtie / 2, Dloss + Dtie / 2),
    nbRow,
    probRow("P(win)", Pwin, Dwin),
    probRow("P(loss)", Ploss, Dloss),
    probRow("P(tie)", Ptie, Dtie)))
  Output[, `:=`(Intervention = paste0("[", A1, "] vs [", A0, "]"),
                Estimator = "tmle", Event = jj, Time = Horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event", "Time",
                                    "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif
  attr(Output, "Horizon") <- Horizon
  attr(Output, "Estimand") <- "Win Ratio"
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
