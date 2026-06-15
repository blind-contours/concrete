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
#' @param nBoot integer (default 0): if positive, a stratified nonparametric
#'   bootstrap with `nBoot` resamples is used for the standard error and
#'   confidence interval (correct but computationally heavy --- each resample
#'   refits the transition hazards). With `nBoot = 0` only the point estimate and
#'   its favor/against decomposition are returned (`se`/CI are `NA`); see Details.
#'
#' @details
#' \strong{Inference.} The analytic adjoint-value efficient influence function for
#' this estimand (the time integral of a bilinear functional of the multistate
#' occupancy) is under development and validation; until it lands, set `nBoot`
#' for resampling-based inference, or use [getRMTIF()] (the first-event version,
#' which has closed-form influence-function inference and reduces exactly to the
#' RMST difference for a single event). The point estimate here is exact and is
#' validated against a brute-force pairwise ground truth.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the net `RMT-IF`, the
#'   time in favor, and the time against; with `nBoot > 0` each carries a
#'   bootstrap SE and CI.
#'
#' @seealso [getRMTIF()] (first-event version with closed-form inference),
#'   [clinicalWinRatio()].
#' @export clinicalRMTIF
#' @importFrom stats qnorm
clinicalRMTIF <- function(data, arm, illness.time, terminal.time, terminal.status,
                          covariates, horizon = NULL, n.grid = 60L, n.folds = 5L,
                          SL.library = c("SL.mean", "SL.glm"), Signif = 0.05,
                          id = NULL, censoring.tv = NULL, nBoot = 0L) {
  data <- as.data.frame(data)
  illness.time <- as.character(illness.time)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  if (length(unique(A)) != 2L) stop("arm must contain both 0 and 1.")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  if (!all(delta %in% c(0, 1))) stop("terminal.status must be coded 0/1 (1 = death).")
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  K <- 1L + length(illness.time)
  grid <- seq(0, horizon, length.out = as.integer(n.grid) + 1L)
  eng <- .msEngine(K, grid)

  parseArm <- function(rows) {
    D <- data[rows, covariates, drop = FALSE]
    D$tD <- ifelse(delta[rows] == 1, term[rows], Inf)
    for (ei in seq_along(illness.time)) {
      ti <- data[[illness.time[ei]]][rows]; ti[is.na(ti)] <- Inf; D[[paste0("t", ei)]] <- ti }
    D$C <- ifelse(delta[rows] == 0, term[rows], Inf)
    D
  }
  tvMats <- NULL
  if (!is.null(censoring.tv)) {
    if (is.null(id)) stop("`id` is required when `censoring.tv` is supplied.")
    censoring.tv <- as.data.frame(censoring.tv)
    tvMats <- .tvLOCF(data[[id]], censoring.tv, id, "time", grid[-length(grid)])
  }
  ## marginal favorability-level occupancy (M+1 x (K+1)) for an arm's rows
  levelOcc <- function(rows) {
    D <- parseArm(rows)
    tvA <- if (is.null(tvMats)) NULL else lapply(tvMats, function(m) m[rows, , drop = FALSE])
    nu <- .msNuisances(eng, D, covariates, SL.library, n.folds, tvA)
    P <- eng$occupancy(nu$rmat)
    L <- matrix(0, eng$M + 1L, K + 1L)
    for (s in eng$ALIVE) L[, eng$stateRank(s) + 1L] <-
      L[, eng$stateRank(s) + 1L] + rowMeans(P[[as.character(s)]])
    L[, K + 1L] <- 1 - rowSums(L[, seq_len(K), drop = FALSE])   # dead level
    L
  }
  ## RMT-IF point estimate from two arms' level-occupancy matrices
  rmtifFrom <- function(Lt, Lc) {
    wfun <- function(X, Y) {                                    # int (X better than Y)
      m <- nrow(X); w <- numeric(m)
      for (jr in seq_len(ncol(X))) {                            # level r (1-indexed)
        above <- if (jr < ncol(Y)) rowSums(Y[, (jr + 1L):ncol(Y), drop = FALSE]) else numeric(m)
        w <- w + X[, jr] * above
      }
      w
    }
    wts <- c(diff(grid) / 2, 0) + c(0, diff(grid) / 2)          # trapezoid weights on nodes
    favor <- wfun(Lt, Lc); against <- wfun(Lc, Lt)
    c(rmtif = sum(wts * (favor - against)),
      tFavor = sum(wts * favor), tAgainst = sum(wts * against))
  }

  est <- rmtifFrom(levelOcc(which(A == 1)), levelOcc(which(A == 0)))

  ## optional stratified bootstrap inference
  seCI <- function(point) c(se = NA_real_, lo = NA_real_, hi = NA_real_)
  if (nBoot > 0L) {
    i1 <- which(A == 1); i0 <- which(A == 0)
    boots <- vapply(seq_len(nBoot), function(b) {
      rmtifFrom(levelOcc(sample(i1, replace = TRUE)),
                levelOcc(sample(i0, replace = TRUE)))
    }, numeric(3))
    z <- stats::qnorm(1 - Signif / 2)
    seCI <- function(nm) { s <- stats::sd(boots[nm, ]); c(se = s, lo = est[nm] - z*s, hi = est[nm] + z*s) }
  }
  z <- stats::qnorm(1 - Signif / 2)
  mkrow <- function(lab, nm) {
    ci <- seCI(nm)
    data.table::data.table(Estimand = lab, `Pt Est` = unname(est[nm]), se = unname(ci["se"]),
      `CI Low` = unname(ci["lo"]), `CI Hi` = unname(ci["hi"]), pValue = NA_real_)
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
