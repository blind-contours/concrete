#' Covariate-adjusted restricted win ratio, win odds, and net benefit
#'
#' @description
#' `getWinRatio()` estimates the **restricted win ratio** and its relatives (win
#' odds, net benefit) from the covariate-adjusted, censoring-corrected
#' counterfactual cumulative-incidence curves that `concrete` targets. It supports
#' a single time-to-event outcome or a **prioritized hierarchy** of competing
#' events (e.g.\ death \eqn{>} hospitalization \eqn{>} stroke).
#'
#' \strong{Single event.} Comparing a random treated patient to a random control
#' patient over \eqn{[0,\tau]}: the treated patient \emph{wins} if the control has
#' the event first (and within the horizon), \emph{loses} if the treated has it
#' first, and \emph{ties} if both are event-free at \eqn{\tau}. By independence of
#' the two patients these are functionals of the marginal survival
#' \eqn{\bar S_a} and cumulative incidence \eqn{\bar F_a = 1 - \bar S_a}:
#' \deqn{P(\text{win}) = \int_0^\tau \bar S_1(t)\, d\bar F_0(t), \qquad
#'       P(\text{loss}) = \int_0^\tau \bar S_0(t)\, d\bar F_1(t).}
#'
#' \strong{Hierarchy.} When `TargetEvent` lists several event codes in priority
#' order (highest priority first), the comparison is the prioritized
#' (Pocock / Finkelstein--Schoenfeld) rule applied to the patients' \emph{first}
#' events: a patient who is event-free beats one who had any event; between two
#' patients with events of \emph{different} priority, the one whose event is
#' \emph{lower} priority (less severe) wins; between two with the \emph{same}
#' event, the one whose event is \emph{later} wins. Writing the per-arm
#' cause-specific cumulative incidences \eqn{F_a^{(k)}} for priority \eqn{k}
#' (\eqn{k=1} highest), the win probability is again a smooth functional of those
#' marginal curves,
#' \deqn{P(\text{win}) = S_1(\tau)\bigl(1 - S_0(\tau)\bigr)
#'   + \sum_{a>b} F_1^{(a)}(\tau)F_0^{(b)}(\tau)
#'   + \sum_k \int_0^\tau \bigl[F_1^{(k)}(\tau) - F_1^{(k)}(t)\bigr]\, dF_0^{(k)}(t),}
#' with \eqn{S_a(\tau) = 1 - \sum_k F_a^{(k)}(\tau)}, and \eqn{P(\text{loss})} the
#' mirror image. This reduces exactly to the single-event formula when one event
#' is given.
#'
#' In all cases the win/loss probabilities are smooth functionals of the targeted
#' curves, so their influence functions are weighted combinations of the
#' per-subject curve influence functions and the win ratio, win odds, and net
#' benefit follow by the delta method --- giving doubly-robust, covariate-adjusted,
#' censoring-corrected inference, unlike the standard unadjusted, censoring-
#' sensitive win ratio. The integral is taken over the fitted target times, so use
#' a reasonably dense `TargetTime` grid.
#'
#' \strong{Assumption (hierarchy).} The prioritized comparison uses each patient's
#' \emph{first} observed event, treating the listed events as competing risks ---
#' the structure `concrete` models. Events occurring after a patient's first event
#' (e.g.\ death following a non-fatal hospitalization) are not used; the fully
#' semi-competing version requires the within-patient joint law and is future work.
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Horizon numeric: the restriction horizon \eqn{\tau} (default: the
#'   largest target time).
#' @param Intervention length-2 numeric: treatment and control indices.
#' @param TargetEvent numeric: the event code, or an ordered vector of event codes
#'   giving the priority hierarchy from highest to lowest (default: the first
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
  TargetEvent <- as.numeric(TargetEvent)            # priority order, highest first
  if (anyDuplicated(TargetEvent))
    stop("TargetEvent must list each event code at most once (its priority order).")
  if (!all(TargetEvent %in% attr(ConcreteEst, "TargetEvent")))
    stop("All TargetEvent codes must have been targeted in doConcrete().")
  K <- length(TargetEvent)
  if (is.null(Horizon)) Horizon <- max(TargetTime)
  grid <- sort(unique(TargetTime[TargetTime <= Horizon]))
  if (length(grid) < 2L)
    warning("The win ratio is integrated over fewer than two target times; ",
            "refit with a denser TargetTime grid.")
  A1 <- names(ConcreteEst)[Intervention[1]]
  A0 <- names(ConcreteEst)[Intervention[2]]
  m <- length(grid)

  ## marginal cause-specific CIF curves and per-subject influence functions, at
  ## the grid times, for every event in the hierarchy
  RisksObj <- getRisk(ConcreteEst, TargetTime = grid, TargetEvent = TargetEvent, GComp = FALSE)
  Risks <- data.table::as.data.table(RisksObj)[Estimator == "tmle"]
  ICdt <- data.table::as.data.table(attr(RisksObj, "IC"))
  n <- length(unique(ICdt[["ID"]]))

  Fcurve <- function(arm, k) {                      # length-m CIF on the grid
    r <- Risks[Risks[["Intervention"]] == arm & Event == k & Time %in% grid, ]
    data.table::setorder(r, Time)
    r[["Pt Est"]]
  }
  ICmat <- function(arm, k) {                       # m x n IC matrix, grid order
    d <- ICdt[ICdt[["Intervention"]] == arm & Event == k & Time %in% grid, ]
    mm <- data.table::dcast(d, Time ~ ID, value.var = "IC")
    data.table::setorder(mm, Time)
    as.matrix(mm[, setdiff(names(mm), "Time"), with = FALSE])
  }
  ## index by priority position 1..K (1 = highest priority); G = treated, H = control
  G  <- lapply(TargetEvent, function(k) Fcurve(A1, k))
  H  <- lapply(TargetEvent, function(k) Fcurve(A0, k))
  DG <- lapply(TargetEvent, function(k) ICmat(A1, k))
  DH <- lapply(TargetEvent, function(k) ICmat(A0, k))

  ## win probability and its per-subject influence function for a given
  ## (winner, loser) assignment of the two arms. W/L are lists of length-m CIF
  ## curves (one per priority level); DW/DL the matching m x n IC matrices.
  winProb <- function(W, L, DW, DL) {
    Wtau <- vapply(W, function(v) v[m], numeric(1))   # CIF at tau, per level
    Ltau <- vapply(L, function(v) v[m], numeric(1))
    SWtau <- 1 - sum(Wtau)                            # winner-arm survival at tau
    ## point estimate
    P <- SWtau * sum(Ltau) +
      (if (K >= 2) sum(vapply(2:K, function(a) W[[a]][m] * sum(Ltau[seq_len(a - 1)]),
                              numeric(1))) else 0) +
      sum(vapply(seq_len(K), function(k) {
        dL <- diff(c(0, L[[k]]))
        sum((W[[k]][m] - W[[k]]) * dL)
      }, numeric(1)))
    ## influence function via gradient coefficients, contracted with the IC mats
    D <- numeric(n)
    for (k in seq_len(K)) {
      dL <- diff(c(0, L[[k]]))
      ## d P / d W^{(k)}(t_i)
      cW <- -dL
      cW[m] <- cW[m] - (if (k < K) sum(Ltau[(k + 1):K]) else 0)
      ## d P / d L^{(k)}(t_i)
      cL <- numeric(m)
      if (m >= 2) cL[1:(m - 1)] <- W[[k]][2:m] - W[[k]][1:(m - 1)]
      cL[m] <- SWtau + (if (k < K) sum(Wtau[(k + 1):K]) else 0)
      D <- D + colSums(cW * DW[[k]]) + colSums(cL * DL[[k]])
    }
    list(P = P, D = D)
  }

  win  <- winProb(G, H, DG, DH)                       # treated beats control
  loss <- winProb(H, G, DH, DG)                       # control beats treated
  Pwin <- win$P;  Dwin <- win$D
  Ploss <- loss$P; Dloss <- loss$D
  Ptie <- max(0, 1 - Pwin - Ploss)
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
  EventLab <- if (K == 1L) TargetEvent[1] else paste(TargetEvent, collapse = ">")
  Output[, `:=`(Intervention = paste0("[", A1, "] vs [", A0, "]"),
                Estimator = "tmle", Event = EventLab, Time = Horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event", "Time",
                                    "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif
  attr(Output, "Horizon") <- Horizon
  attr(Output, "Estimand") <- "Win Ratio"
  attr(Output, "Priority") <- TargetEvent
  attr(Output, "Simultaneous") <- FALSE
  attr(Output, "GComp") <- FALSE
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
