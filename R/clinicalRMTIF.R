#' Clinical (death-priority) restricted mean time in favor of treatment
#'
#' @description
#' The death-priority (multistate) version of the restricted mean time in favor of
#' treatment (RMT-IF; Mao 2023): the average time over \eqn{[0,\tau]} that a random
#' treated patient spends in a clinically more favorable state than a random
#' control, under the hierarchy that compares on the most serious event first
#' \emph{whenever} it occurs. Unlike the first-event version [getRMTIF()], this
#' \strong{credits a higher-priority event that follows a lower-priority one}
#' (death after a hospitalization), and unlike the win ratio it is in time units.
#'
#' It is built on the same Markov multistate engine as [clinicalWinRatio()]: each
#' subject's state at time \eqn{t} is the subset of non-fatal events experienced
#' (plus alive/dead), every transition intensity is estimated by a Super Learner
#' with IPCW censoring correction and optional cross-fitting, and the favorability
#' of a state is its most-severe event so far (event-free best, death worst). The
#' point estimate is
#' \deqn{\mathrm{RMT\text{-}IF} = \int_0^\tau \big[w(t) - l(t)\big]\,dt,
#'   \quad w(t)=\sum_{r<r'}\pi^1_r(t)\,\pi^0_{r'}(t),}
#' with \eqn{\pi^a_r(t)} the probability that arm \eqn{a} occupies favorability
#' level \eqn{r} at \eqn{t}.
#'
#' @inheritParams clinicalWinRatio
#' @param nBoot integer (default 0): inference method. With `nBoot = 0` (default)
#'   the **analytic adjoint-value efficient influence function** is used for the
#'   net RMT-IF standard error (fast). If positive, a stratified nonparametric
#'   bootstrap with `nBoot` resamples is used instead, which also gives SEs for
#'   the time-in-favor / time-against decomposition (correct but heavy --- each
#'   resample refits the transition hazards).
#'
#' @details
#' \strong{Inference.} The net RMT-IF is the time integral of a bilinear
#' functional of the two arms' multistate level-occupancy curves. Its efficient
#' influence function is obtained by the chain rule: the per-arm occupancy
#' influence functions (a reward-accumulation adjoint over the multistate process,
#' with the death state absorbing) weighted by the gradient coefficients
#' \eqn{c^1_r = P(\text{control worse than } r) - P(\text{control better than }
#' r)} and the mirror for the control arm. This reuses the same adjoint-value
#' machinery, IPCW correction, and cross-fitting as [clinicalWinRatio()]. The
#' point estimate is exact and validated against a brute-force pairwise ground
#' truth; the analytic SE is validated against the bootstrap and against
#' empirical coverage. SEs for the favor/against decomposition require `nBoot`.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the net `RMT-IF` (with an
#'   influence-function SE, CI and p-value), the time in favor, and the time
#'   against.
#'
#' @seealso [getRMTIF()] (first-event version with closed-form inference),
#'   [clinicalWinRatio()].
#' @export clinicalRMTIF
#' @importFrom stats qnorm
clinicalRMTIF <- function(data, arm, illness.time, terminal.time, terminal.status,
                          covariates, horizon = NULL, n.grid = 60L, n.folds = 5L,
                          SL.library = c("SL.mean", "SL.glm"), Signif = 0.05,
                          id = NULL, censoring.tv = NULL, crossover = NULL,
                          min.cens.surv = 0.05, nBoot = 0L) {
  data <- as.data.frame(data)
  illness.time <- as.character(illness.time)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  if (!is.null(crossover) && !crossover %in% names(data))
    stop("crossover column '", crossover, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  if (length(unique(A)) != 2L) stop("arm must contain both 0 and 1.")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  if (!all(delta %in% c(0, 1))) stop("terminal.status must be coded 0/1 (1 = death).")
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  K <- 1L + length(illness.time)
  grid <- seq(0, horizon, length.out = as.integer(n.grid) + 1L)
  eng <- .msEngine(K, grid)

  sw <- if (is.null(crossover)) rep(Inf, nrow(data)) else { x <- as.numeric(data[[crossover]]); x[is.na(x)] <- Inf; x }
  parseArm <- function(rows) {
    D <- data[rows, covariates, drop = FALSE]
    D$tD <- ifelse(delta[rows] == 1, term[rows], Inf)
    for (ei in seq_along(illness.time)) {
      ti <- data[[illness.time[ei]]][rows]; ti[is.na(ti)] <- Inf; D[[paste0("t", ei)]] <- ti }
    D$C <- pmin(ifelse(delta[rows] == 0, term[rows], Inf), sw[rows])   # re-censor at switch
    D$switch <- sw[rows]
    D
  }
  tvMats <- NULL
  if (!is.null(censoring.tv)) {
    if (is.null(id)) stop("`id` is required when `censoring.tv` is supplied.")
    censoring.tv <- as.data.frame(censoring.tv)
    tvMats <- .tvLOCF(data[[id]], censoring.tv, id, "time", grid[-length(grid)])
  }
  wts <- c(diff(grid) / 2, 0) + c(0, diff(grid) / 2)            # trapezoid weights on M+1 nodes
  ## fit one arm: returns the engine setup + its marginal level occupancy
  fitArm <- function(rows) {
    D <- parseArm(rows)
    tvA <- if (is.null(tvMats)) NULL else lapply(tvMats, function(m) m[rows, , drop = FALSE])
    nu <- .msNuisances(eng, D, covariates, SL.library, n.folds, tvA, xover = D$switch, minG = min.cens.surv)
    arm <- eng$armSetup(D, nu$rmat, nu$Ginv)
    list(arm = arm, L = eng$levelOcc(nu$rmat))
  }
  ## int_0^tau (X better than Y) per node, from two level-occupancy matrices
  wcurve <- function(X, Y) {
    m <- nrow(X); w <- numeric(m)
    for (jr in seq_len(ncol(X)))
      w <- w + X[, jr] * (if (jr < ncol(Y)) rowSums(Y[, (jr + 1L):ncol(Y), drop = FALSE]) else 0)
    w
  }
  T1 <- fitArm(which(A == 1)); T0 <- fitArm(which(A == 0)); LX <- T1$L; LY <- T0$L
  favor <- wcurve(LX, LY); against <- wcurve(LY, LX)
  est <- c(rmtif = sum(wts * (favor - against)),
           tFavor = sum(wts * favor), tAgainst = sum(wts * against))

  ## --- chain-rule gradient coefficients d/d(pi_r) of favor (int w) and against (int l) ---
  better <- function(L, r) if (r > 1) rowSums(L[, seq_len(r - 1), drop = FALSE]) else numeric(nrow(L))
  worse  <- function(L, r) if (r < ncol(L)) rowSums(L[, (r + 1):ncol(L), drop = FALSE]) else numeric(nrow(L))
  cFavX <- sapply(seq_len(K + 1L), function(r) worse(LY, r))    # d(favor)/d(piX_r) = P(Y worse than r)
  cFavY <- sapply(seq_len(K + 1L), function(r) better(LX, r))   # d(favor)/d(piY_r) = P(X better than r)
  cAgnX <- sapply(seq_len(K + 1L), function(r) better(LY, r))   # d(against)/d(piX_r)
  cAgnY <- sapply(seq_len(K + 1L), function(r) worse(LX, r))    # d(against)/d(piY_r)

  Ntot <- T1$arm$n + T0$arm$n; piT <- T1$arm$n / Ntot; piC <- T0$arm$n / Ntot
  z <- stats::qnorm(1 - Signif / 2)
  seVec <- c(rmtif = NA_real_, tFavor = NA_real_, tAgainst = NA_real_)
  if (nBoot > 0L) {
    i1 <- which(A == 1); i0 <- which(A == 0)
    boots <- vapply(seq_len(nBoot), function(b) {
      b1 <- fitArm(sample(i1, replace = TRUE)); b0 <- fitArm(sample(i0, replace = TRUE))
      fv <- wcurve(b1$L, b0$L); ag <- wcurve(b0$L, b1$L)
      c(sum(wts*(fv-ag)), sum(wts*fv), sum(wts*ag))
    }, numeric(3))
    seVec <- c(rmtif = stats::sd(boots[1, ]), tFavor = stats::sd(boots[2, ]), tAgainst = stats::sd(boots[3, ]))
  } else {
    ## per-arm influence functions of favor / against (net = favor - against)
    DfavX <- eng$rmtifArmIF(T1$arm, cFavX, wts); DfavY <- eng$rmtifArmIF(T0$arm, cFavY, wts)
    DagnX <- eng$rmtifArmIF(T1$arm, cAgnX, wts); DagnY <- eng$rmtifArmIF(T0$arm, cAgnY, wts)
    seArm <- function(DXa, DYa) sqrt((sum(((1/piT)*DXa)^2) + sum(((1/piC)*DYa)^2)) / Ntot^2)
    est["tFavor"]   <- est["tFavor"]   + mean(DfavX) + mean(DfavY)        # one-step
    est["tAgainst"] <- est["tAgainst"] + mean(DagnX) + mean(DagnY)
    est["rmtif"]    <- est["tFavor"] - est["tAgainst"]                     # identity preserved
    seVec <- c(rmtif = seArm(DfavX - DagnX, DfavY - DagnY),
               tFavor = seArm(DfavX, DfavY), tAgainst = seArm(DagnX, DagnY))
  }
  mkrow <- function(lab, nm) {
    s <- unname(seVec[nm]); e <- unname(est[nm])
    data.table::data.table(Estimand = lab, `Pt Est` = e, se = s,
      `CI Low` = if (is.na(s)) NA_real_ else e - z * s,
      `CI Hi`  = if (is.na(s)) NA_real_ else e + z * s,
      pValue   = if (is.na(s) || lab != "RMT-IF") NA_real_ else 2 * stats::pnorm(-abs(e / s)))
  }
  Output <- data.table::rbindlist(list(
    mkrow("RMT-IF", "rmtif"), mkrow("Time in favor", "tFavor"), mkrow("Time against", "tAgainst")))
  Output[, `:=`(Intervention = "[arm=1] vs [arm=0]", Estimator = "tmle",
                Event = paste(c("D", eng$NF), collapse = ">"), Time = horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event",
                                    "Time", "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif; attr(Output, "Horizon") <- horizon
  attr(Output, "Estimand") <- "Clinical RMT-IF"; attr(Output, "Tiers") <- K
  attr(Output, "Experimental") <- TRUE; attr(Output, "nBoot") <- nBoot
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
