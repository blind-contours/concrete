#' Super Learner for an illness-death transition hazard with left truncation
#'
#' @description
#' Internal building block for the clinical (death-after-non-fatal) win ratio.
#' Estimates a transition intensity on a **left-truncated** risk set --- in
#' particular the post-non-fatal-event death hazard \eqn{\alpha_{12}(u\mid W)}
#' (state 1 = alive after the non-fatal event, state 2 = dead), where subjects
#' enter the risk set only at their non-fatal-event time (delayed entry).
#'
#' The hazard is fit by **pooled (discrete-time) logistic Super Learning**: each
#' subject is expanded into person-time rows over a regular grid, from the
#' interval containing their entry to the interval containing their exit, with a
#' 0/1 event indicator. A [SuperLearner::SuperLearner] binomial ensemble is fit
#' to those rows; the per-interval cumulative-hazard increment is
#' \eqn{-\log(1-p_k)}. Left truncation, right censoring, and (optionally)
#' time-since-entry as a covariate all fall out of the person-time construction.
#'
#' Markov (default): the hazard depends on calendar time \eqn{u} and \eqn{W}.
#' Semi-Markov (`semiMarkov = TRUE`): time-since-entry \eqn{u-s} is added as a
#' covariate (the clinically realistic clock reset).
#'
#' @keywords internal
#' @noRd
makeTransitionLongData <- function(entry, exit, event, Cov, grid, semiMarkov = FALSE) {
  M <- length(grid) - 1L
  starts <- grid[-length(grid)]                       # interval starts grid[1..M]
  widths <- diff(grid)
  scale <- max(grid)
  Cov <- as.data.frame(Cov)
  rows <- lapply(seq_along(exit), function(i) {
    # intervals whose START is at/after entry (left truncation) and before exit
    ks <- which(starts >= entry[i] & starts < exit[i])
    if (!length(ks)) return(NULL)
    endTimes <- grid[ks + 1L]
    di <- Cov[rep(i, length(ks)), , drop = FALSE]
    di[[".hazTime"]] <- endTimes / scale
    di[[".width"]]   <- widths[ks] / scale
    if (semiMarkov) di[[".sinceEntry"]] <- (endTimes - entry[i]) / scale
    yk <- integer(length(ks))
    if (isTRUE(event[i] == 1L)) yk[length(ks)] <- 1L   # event in the last at-risk interval
    di[[".Y"]] <- yk
    di
  })
  data.table::rbindlist(Filter(Negate(is.null), rows))
}

#' @keywords internal
#' @noRd
#' @importFrom SuperLearner SuperLearner
fitTransitionSL <- function(entry, exit, event, Cov, grid,
                            SL.library = c("SL.mean", "SL.glm", "SL.glmnet"),
                            semiMarkov = FALSE, V = 5L) {
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("Transition-hazard Super Learner requires the 'SuperLearner' package.")
  Long <- makeTransitionLongData(entry, exit, event, Cov, grid, semiMarkov)
  if (!nrow(Long) || length(unique(Long[[".Y"]])) < 2L)
    stop("Transition risk set has no events (or no non-events) to learn from.")
  # .width is kept for bookkeeping but excluded from the learner (it is constant on
  # a regular grid -> collinear with the intercept); the per-interval event
  # probability already encodes the interval length.
  covNames <- setdiff(names(Long), c(".Y", ".width"))
  X <- as.data.frame(Long[, .SD, .SDcols = covNames])
  fit <- SuperLearner::SuperLearner(
    Y = Long[[".Y"]], X = X, family = stats::binomial(),
    SL.library = SL.library, cvControl = list(V = as.integer(V)))
  structure(list(SL = fit, covNames = covNames, baseCov = setdiff(covNames,
                 c(".hazTime", ".width", ".sinceEntry")),
                 grid = grid, semiMarkov = semiMarkov, scale = max(grid)),
            class = "ConcreteTransitionHaz")
}

#' Predict the per-interval cumulative-hazard increments of the transition.
#'
#' For Markov fits returns an `M x n` matrix of increments
#' \eqn{\Lambda_{12}(g_{k+1})-\Lambda_{12}(g_k)} for each subject over the full
#' grid. For semi-Markov fits the increment at calendar interval `k` depends on
#' the subject's entry time, so an `entry` vector is required and the increment
#' is evaluated at time-since-entry.
#' @keywords internal
#' @noRd
predictTransitionSL <- function(Fit, NewCov, entry = NULL) {
  grid <- Fit$grid; M <- length(grid) - 1L; scale <- Fit$scale
  ends <- grid[-1L]; widths <- diff(grid)
  NewCov <- as.data.frame(NewCov); n <- nrow(NewCov)
  Long <- data.table::rbindlist(lapply(seq_len(n), function(i) {
    di <- NewCov[rep(i, M), Fit$baseCov, drop = FALSE]
    di[[".hazTime"]] <- ends / scale
    di[[".width"]]   <- widths / scale
    if (Fit$semiMarkov) {
      e <- if (is.null(entry)) 0 else entry[i]
      di[[".sinceEntry"]] <- pmax(ends - e, 0) / scale
    }
    di
  }))
  X <- as.data.frame(Long[, .SD, .SDcols = Fit$covNames])
  p <- as.numeric(stats::predict(Fit$SL, newdata = X, onlySL = TRUE)$pred)
  p <- pmin(pmax(p, 1e-12), 1 - 1e-8)
  inc <- -log1p(-p)                                   # cumulative-hazard increments
  matrix(inc, nrow = M, ncol = n)                     # M x n
}

#' Per-subject illness-death curves from the three transition increments.
#'
#' Assembles, with **midpoint quadrature** in state 0 (so the path-probability
#' integral is accurate on a coarse grid), the per-subject curves needed for the
#' clinical win ratio: state-0 survival \eqn{S_0}, state-1 occupancy \eqn{p_1},
#' overall survival \eqn{S^D = S_0 + p_1}, and the HFH-entry-to-\eqn{\tau} density
#' \eqn{\pi(s) = S_0(s)\alpha_{01}(s)S_{12}(s,\tau)}.
#'
#' @param I01,I02,I12 `M x n` cumulative-hazard increment matrices from
#'   [predictTransitionSL] for transitions 0->1, 0->2, 1->2.
#' @return a list of matrices: `S0`,`p1`,`SD` are `(M+1) x n` at grid times;
#'   `pi` is `M x n` (HFH in interval k, surviving to the horizon).
#' @keywords internal
#' @noRd
multistateCurves <- function(I01, I02, I12) {
  M <- nrow(I01); n <- ncol(I01)
  cumExit <- apply(I01 + I02, 2, cumsum)                       # Lambda0 at g[k+1]
  S0 <- rbind(1, exp(-cumExit))                               # (M+1) x n
  S0mid <- sqrt(S0[1:M, , drop = FALSE] * S0[2:(M + 1), , drop = FALSE])
  L12 <- rbind(0, apply(I12, 2, cumsum))                       # (M+1) x n, Lambda12 at grid
  L12end <- L12[2:(M + 1), , drop = FALSE]
  w <- S0mid * I01 * exp(L12end)
  p1end <- exp(-L12end) * apply(w, 2, cumsum)                  # p1 at g[k+1]
  p1 <- rbind(0, p1end)
  SD <- S0 + p1
  S12toTau <- exp(-(matrix(L12[M + 1, ], M, n, byrow = TRUE) - L12end))  # S12(g[k+1], tau)
  pimat <- S0mid * I01 * S12toTau
  list(S0 = S0, p1 = p1, SD = SD, pi = pimat, S12toTau = S12toTau, S0mid = S0mid)
}

#' Post-entry survival S_12(s, t | W) = exp(-[Lambda_12(t) - Lambda_12(s)]).
#'
#' @param incMat  `M x n` increments from [predictTransitionSL].
#' @param grid    the time grid (length M+1).
#' @param s,t     scalars with s <= t <= max(grid).
#' @keywords internal
#' @noRd
s12FromIncrements <- function(incMat, grid, s, t) {
  ends <- grid[-1L]
  sel <- ends > s & ends <= t                         # intervals in (s, t]
  if (!any(sel)) return(rep(1, ncol(incMat)))
  exp(-colSums(incMat[sel, , drop = FALSE]))
}
