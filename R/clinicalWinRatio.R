#' Clinical (death-priority) win ratio for an illness-death outcome [EXPERIMENTAL]
#'
#' @description
#' \strong{Experimental.} Estimates the \emph{clinical} (death-priority) win
#' ratio, win odds, and net benefit for a two-arm trial with a non-fatal
#' intercurrent event (e.g.\ heart-failure hospitalization) and a terminal event
#' (death). Unlike the competing-risks (first-event) win ratio in [getWinRatio()],
#' this estimand counts **death even when it follows the non-fatal event** --- the
#' clinically intended hierarchy "compare on death first; break ties on the
#' non-fatal event." It is built on a Markov illness-death model with three
#' transition intensities (alive\eqn{\to}non-fatal, alive\eqn{\to}death,
#' post-non-fatal\eqn{\to}death), each estimated by a Super Learner (the
#' post-non-fatal death hazard on a left-truncated risk set), and returns
#' influence-function confidence intervals that are doubly-robust, covariate-
#' adjusted, and censoring-corrected (IPCW).
#'
#' The estimator and its inference are validated against ground truth (a
#' brute-force pairwise win ratio on full simulated histories): see the
#' "Win ratio for trialists" article and `scripts/make-clinical-wr-*.R`. It is
#' marked experimental because it currently takes its own multistate data frame
#' (below) rather than the standard [formatArguments()] pipeline, and assumes a
#' single non-fatal event type and conditionally-independent censoring (CAR).
#'
#' @param data a `data.frame`/`data.table`, one row per subject.
#' @param arm character: name of the binary treatment column (1 = active arm).
#' @param illness.time character: name of the non-fatal-event time column; `NA`
#'   (or `Inf`) for subjects who never had the non-fatal event.
#' @param terminal.time character: name of the terminal time column (time of death
#'   or of censoring, whichever came first).
#' @param terminal.status character: name of the terminal status column
#'   (1 = death, 0 = censored).
#' @param covariates character vector: baseline covariate column names.
#' @param horizon numeric: the restriction horizon \eqn{\tau} (default: the
#'   largest terminal time).
#' @param n.grid integer (default 60): number of time intervals for the discrete
#'   hazard / path-probability quadrature.
#' @param SL.library character vector: SuperLearner library for the transition and
#'   censoring hazards (default `c("SL.mean", "SL.glm")`).
#' @param Signif numeric (default 0.05): alpha for confidence intervals.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
#'   net benefit, and the win/loss/tie probabilities, each with an
#'   influence-function standard error, confidence interval, and (for the
#'   comparative statistics) a p-value against the null of no difference.
#'
#' @seealso [getWinRatio()] for the single-event and competing-risks win ratio.
#' @export clinicalWinRatio
#' @importFrom stats qnorm pnorm binomial predict var
#' @examples
#' \dontrun{
#' # data with: arm, t_hfh (NA if none), t_term (death or censoring), died (1/0), age, sex
#' clinicalWinRatio(trial, arm = "arm", illness.time = "t_hfh",
#'                  terminal.time = "t_term", terminal.status = "died",
#'                  covariates = c("age", "sex"), horizon = 1460)
#' }
clinicalWinRatio <- function(data, arm, illness.time, terminal.time, terminal.status,
                             covariates, horizon = NULL, n.grid = 60L,
                             SL.library = c("SL.mean", "SL.glm"), Signif = 0.05) {
  data <- as.data.frame(data)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  if (length(unique(A)) != 2L) stop("arm must contain both 0 and 1.")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  s <- data[[illness.time]]; s[is.na(s)] <- Inf
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  tau <- horizon
  grid <- seq(0, tau, length.out = as.integer(n.grid) + 1L)
  W <- data[, covariates, drop = FALSE]

  ## --- parse to per-subject illness-death observed quantities ---
  parse <- function(idx) {
    si <- s[idx]; ti <- term[idx]; di <- delta[idx]
    hadNF <- is.finite(si) & si < ti
    exit0 <- pmin(ifelse(hadNF, si, ti), tau)
    d01 <- as.integer(hadNF & si < tau)
    d02 <- as.integer(!hadNF & di == 1 & ti < tau)
    d12 <- as.integer(hadNF & di == 1 & ti < tau)
    data.frame(
      W[idx, , drop = FALSE],
      d01 = d01, t01 = ifelse(d01 == 1, si, Inf),
      d02 = d02, t02 = ifelse(d02 == 1, ti, Inf),
      d12 = d12, t12 = ifelse(d12 == 1, ti, Inf),
      exit0 = exit0,
      entry1 = ifelse(d01 == 1, si, Inf),
      exit1  = ifelse(d01 == 1, pmin(ti, tau), Inf),
      obsT = pmin(ti, tau),
      censE = as.integer(di == 0 & ti < tau),
      check.names = FALSE)
  }
  idxT <- which(A == 1); idxC <- which(A == 0)
  AT <- cwrArmFit(parse(idxT), covariates, grid, SL.library)
  AC <- cwrArmFit(parse(idxC), covariates, grid, SL.library)
  out <- cwrAssemble(AT, AC, piT = mean(A == 1), piC = mean(A == 0),
                     Ntot = length(A), Signif = Signif)
  attr(out, "Horizon") <- tau
  attr(out, "Estimand") <- "Clinical Win Ratio"
  attr(out, "Experimental") <- TRUE
  class(out) <- union("ConcreteOut", class(out))
  out[]
}

#' Per-arm transition fits + building-block influence functions (IPCW).
#' @keywords internal
#' @noRd
cwrArmFit <- function(est, covariates, grid, SL.library) {
  M <- length(grid) - 1L; Mp1 <- M + 1L; starts <- grid[-Mp1]; Gmin <- 0.05
  n <- nrow(est); Cov <- est[, covariates, drop = FALSE]
  f01 <- fitTransitionSL(rep(0, n), est$exit0, est$d01, Cov, grid, SL.library = SL.library)
  f02 <- fitTransitionSL(rep(0, n), est$exit0, est$d02, Cov, grid, SL.library = SL.library)
  sub <- est$d01 == 1
  f12 <- fitTransitionSL(est$entry1[sub], est$exit1[sub], est$d12[sub],
                         Cov[sub, , drop = FALSE], grid, SL.library = SL.library)
  fc  <- fitTransitionSL(rep(0, n), est$obsT, est$censE, Cov, grid, SL.library = SL.library)
  I01 <- predictTransitionSL(f01, Cov); I02 <- predictTransitionSL(f02, Cov)
  I12 <- predictTransitionSL(f12, Cov)
  Glag <- pmax(rbind(1, exp(-apply(predictTransitionSL(fc, Cov), 2, cumsum)))[1:M, , drop = FALSE], Gmin)
  cur <- multistateCurves(I01, I02, I12); S0 <- cur$S0; SD <- cur$SD; pimat <- cur$pi
  L12 <- rbind(0, apply(I12, 2, cumsum)); S12tau <- cur$S12toTau
  ## martingale increments (Y already censoring-adjusted via exit times)
  Y0 <- outer(starts, est$exit0, "<")
  inIv <- function(tt) { o <- matrix(0L, M, n); k <- findInterval(tt, grid)
    ok <- is.finite(tt) & k >= 1 & k <= M; o[cbind(k[ok], which(ok))] <- 1L; o }
  dM01 <- inIv(est$t01) - Y0 * I01; dM02 <- inIv(est$t02) - Y0 * I02; dM0 <- dM01 + dM02
  Y1 <- outer(starts, est$entry1, ">=") & outer(starts, est$exit1, "<")
  dM12 <- inIv(est$t12) - Y1 * I12
  S0end <- S0[2:Mp1, , drop = FALSE]
  abar <- mean(S0[Mp1, ]); Da <- (S0[Mp1, ] - abar) - S0[Mp1, ] * colSums(dM0 / (S0end * Glag))
  Pi <- colSums(pimat); Pigt <- apply(pimat, 2, function(c) rev(cumsum(rev(c))) - c)
  phi01 <- ((S0end * S12tau - Pigt) / S0end) / Glag; phi02 <- (-Pigt / S0end) / Glag; phi12 <- (-S12tau) / Glag
  xiTheta <- colSums(phi01 * dM01) + colSums(phi02 * dM02) + colSums(phi12 * dM12)
  Tbar <- mean(Pi); DTheta <- (Pi - Tbar) + xiTheta
  DS0 <- matrix(0, Mp1, n); cumdM0 <- rbind(0, apply(dM0 / (S0end * Glag), 2, cumsum))
  for (j in 1:Mp1) DS0[j, ] <- (S0[j, ] - mean(S0[j, ])) - S0[j, ] * cumdM0[j, ]
  Dp1 <- matrix(0, Mp1, n)
  for (j in 2:Mp1) { H <- j - 1L
    S12toH <- exp(-(matrix(L12[j, ], M, n, byrow = TRUE) - L12[2:Mp1, , drop = FALSE]))
    piH <- cur$S0mid * I01 * S12toH; if (H < M) piH[(H + 1):M, ] <- 0
    PiH <- colSums(piH); PigtH <- apply(piH, 2, function(c) rev(cumsum(rev(c))) - c)
    p01 <- ((S0end * S12toH - PigtH) / S0end) / Glag; p02 <- (-PigtH / S0end) / Glag; p12 <- (-S12toH) / Glag
    if (H < M) { p01[(H + 1):M, ] <- 0; p02[(H + 1):M, ] <- 0; p12[(H + 1):M, ] <- 0 }
    Dp1[j, ] <- (PiH - mean(PiH)) + colSums(p01 * dM01) + colSums(p02 * dM02) + colSums(p12 * dM12) }
  DFD <- -(DS0 + Dp1)
  list(n = n, M = M, Mp1 = Mp1, SDbar = rowMeans(SD), DFD = DFD, a = abar, Da = Da,
       Theta = Tbar, DTheta = DTheta, hbar = rowMeans(pimat), pimat = pimat)
}

#' Assemble the clinical win ratio + influence-function inference.
#' @keywords internal
#' @noRd
cwrAssemble <- function(AT, AC, piT, piC, Ntot, Signif) {
  M <- AT$M; Mp1 <- AT$Mp1
  ## exact single-event win EIF on overall survival S^D (getWinRatio coef logic)
  levelOne <- function(A, B) {
    ST <- A$SDbar; SC <- B$SDbar; S1g <- ST[2:Mp1]
    dF0 <- diff(1 - SC); dF1 <- diff(1 - ST)
    coef0 <- ifelse(seq_len(M) < M, c(dF1[-1L], NA), S1g[M])
    list(point = sum(S1g * dF0),
         Dt = colSums((-dF0) * A$DFD[2:Mp1, , drop = FALSE]),
         Dc = colSums(coef0 * B$DFD[2:Mp1, , drop = FALSE]))
  }
  L1w <- levelOne(AT, AC); L1l <- levelOne(AC, AT)
  hT <- AT$hbar; hC <- AC$hbar
  tailT <- rev(cumsum(rev(hT))) - hT; tailC <- rev(cumsum(rev(hC))) - hC
  W2 <- AT$a * AC$Theta + sum(hC * tailT); L2 <- AC$a * AT$Theta + sum(hT * tailC)
  Pwin <- L1w$point + W2; Ploss <- L1l$point + L2; Ptie <- max(0, 1 - Pwin - Ploss)
  ## per-subject within-arm influence contributions for Pwin and Ploss
  wC_lt <- c(0, cumsum(hC)[seq_len(M - 1L)])
  DPwin_T <- L1w$Dt + AC$Theta * AT$Da + (colSums(wC_lt * AT$pimat) - sum(wC_lt * hT))
  DPwin_C <- L1w$Dc + AT$a * AC$DTheta + (colSums(tailT * AC$pimat) - sum(tailT * hC))
  wT_lt <- c(0, cumsum(hT)[seq_len(M - 1L)])
  DPloss_C <- L1l$Dt + AT$Theta * AC$Da + (colSums(wT_lt * AC$pimat) - sum(wT_lt * hC))
  DPloss_T <- L1l$Dc + AC$a * AT$DTheta + (colSums(tailC * AT$pimat) - sum(tailC * hT))
  z <- stats::qnorm(1 - Signif / 2)
  ## SE for an estimand with gradient (gw, gl) wrt (Pwin, Ploss)
  seGrad <- function(gw, gl) {
    Dt <- (1 / piT) * (gw * DPwin_T + gl * DPloss_T)
    Dc <- (1 / piC) * (gw * DPwin_C + gl * DPloss_C)
    sqrt((sum(Dt^2) + sum(Dc^2)) / Ntot^2)
  }
  ratioRow <- function(label, val, gw, gl) {            # log-scale CI for ratio estimands
    se <- seGrad(gw, gl); sl <- se / val
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val * exp(-z * sl), `CI Hi` = val * exp(z * sl),
      pValue = 2 * stats::pnorm(-abs(log(val) / sl)))
  }
  diffRow <- function(label, val, gw, gl) {             # natural-scale CI
    se <- seGrad(gw, gl)
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val - z * se, `CI Hi` = val + z * se,
      pValue = 2 * stats::pnorm(-abs(val / se)))
  }
  probRow <- function(label, val, gw, gl) {
    se <- seGrad(gw, gl)
    data.table::data.table(Estimand = label, `Pt Est` = val, se = se,
      `CI Low` = val - z * se, `CI Hi` = val + z * se, pValue = NA_real_)
  }
  NB <- Pwin - Ploss; dWO <- 2 / (1 - NB)^2
  out <- data.table::rbindlist(list(
    ratioRow("Win Ratio", Pwin / Ploss, 1 / Ploss, -Pwin / Ploss^2),
    ratioRow("Win Odds", (1 + NB) / (1 - NB), dWO, -dWO),
    diffRow("Net Benefit", NB, 1, -1),
    probRow("P(win)", Pwin, 1, 0),
    probRow("P(loss)", Ploss, 0, 1),
    probRow("P(tie)", Ptie, -1, -1)))
  out[]
}
