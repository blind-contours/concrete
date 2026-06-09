#' Time-varying covariates in the censoring model (shared internal machinery).
#'
#' Post-randomization measurements (e.g.\ echo / KCCQ / 6-minute walk at follow-up
#' visits) can drive dropout, making censoring informative conditional on baseline
#' covariates alone. Conditioning the censoring hazard on the time-varying
#' covariate history restores conditional independence (CAR) and removes the
#' inverse-probability-of-censoring-weight bias. These covariates are used
#' \strong{only} in the censoring model --- never the outcome hazards --- so the
#' intent-to-treat / marginal estimand is preserved (they are post-treatment
#' mediators).
#'
#' `.tvLOCF()` builds last-observation-carried-forward value and
#' change-from-baseline matrices on the hazard grid; `.tvCensoringInc()` fits the
#' (cross-fitted) discrete-time censoring hazard with those covariates and returns
#' per-interval cumulative-hazard increments. Both the core IPCW (via an override of
#' the lagged censoring survival) and `clinicalWinRatio()` consume these.
#' @keywords internal
#' @name tvCensoring
NULL

#' LOCF value + change-from-baseline matrices (n x M) on interval starts.
#' @param ids subject ids (one per subject, in data order).
#' @param tv long data.frame with columns `idcol`, `timecol`, and value columns.
#' @param starts numeric vector of interval-start times (the hazard grid, length M).
#' @keywords internal
#' @noRd
.tvLOCF <- function(ids, tv, idcol, timecol, starts) {
  M <- length(starts)
  valcols <- setdiff(names(tv), c(idcol, timecol))
  if (!length(valcols)) stop("censoring time-varying data has no value columns.")
  tvById <- split(tv, factor(tv[[idcol]], levels = unique(ids)))
  out <- list()
  for (vc in valcols) {
    val <- matrix(0, length(ids), M); chg <- matrix(0, length(ids), M)
    for (k in seq_along(ids)) {
      sub <- tvById[[as.character(ids[k])]]
      if (is.null(sub) || !nrow(sub)) next
      ord <- order(sub[[timecol]]); st <- sub[[timecol]][ord]; sv <- sub[[vc]][ord]
      base <- sv[1]; idx <- findInterval(starts, st)            # most recent measurement at/before start
      v <- ifelse(idx == 0, base, sv[pmax(idx, 1L)])            # LOCF; baseline before first measurement
      val[k, ] <- v; chg[k, ] <- v - base
    }
    out[[paste0(vc, "_val")]] <- val; out[[paste0(vc, "_chg")]] <- chg
  }
  out
}

#' Cross-fitted discrete-time censoring-hazard increments with time-varying covariates.
#'
#' @param times numeric grid `c(0, t1, ..., tau)` (length M+1).
#' @param obsT,censInd per-subject observed time and censoring indicator (1 = censored).
#' @param baseCov data.frame of baseline covariates (n rows).
#' @param tvMats named list of `n x M` time-varying covariate matrices (from [.tvLOCF()]).
#' @return an `M x n` matrix of cumulative-hazard increments on `times`.
#' @keywords internal
#' @noRd
.tvCensoringInc <- function(times, obsT, censInd, baseCov, tvMats, SL.library, n.folds) {
  M <- length(times) - 1L; starts <- times[-(M + 1L)]; scale <- max(times); n <- length(obsT)
  tvnames <- names(tvMats)
  V <- max(1L, min(as.integer(n.folds), floor(n / 30)))
  fold <- if (V <= 1L) rep(1L, n) else sample(rep(seq_len(V), length.out = n))
  inc <- matrix(1e-10, M, n)
  fitPred <- function(tr, te) {
    lastj <- pmax(pmin(findInterval(obsT[tr], times), M), 1L)   # intervals at risk per training subject
    si <- rep(tr, lastj); jj <- unlist(lapply(lastj, seq_len)); lj <- rep(lastj, lastj)
    X <- baseCov[si, , drop = FALSE]; X[[".t"]] <- starts[jj] / scale
    for (nm in tvnames) X[[nm]] <- tvMats[[nm]][cbind(si, jj)]
    Y <- as.integer(jj == lj & censInd[si] == 1)               # censoring event in the last at-risk interval
    if (length(unique(Y)) < 2L) return(matrix(1e-10, M, length(te)))
    fit <- suppressWarnings(SuperLearner::SuperLearner(Y = Y, X = as.data.frame(X),
              family = stats::binomial(), SL.library = SL.library, cvControl = list(V = 5L)))
    Sj <- rep(te, each = M); Jj <- rep(seq_len(M), times = length(te))
    Xp <- baseCov[Sj, , drop = FALSE]; Xp[[".t"]] <- starts[Jj] / scale
    for (nm in tvnames) Xp[[nm]] <- tvMats[[nm]][cbind(Sj, Jj)]
    p <- as.numeric(stats::predict(fit, newdata = as.data.frame(Xp), onlySL = TRUE)$pred)
    p <- pmin(pmax(p, 1e-12), 1 - 1e-8)
    matrix(-log1p(-p), nrow = M)
  }
  if (V <= 1L) inc[] <- tryCatch(fitPred(seq_len(n), seq_len(n)), error = function(e) inc)
  else for (v in seq_len(V)) { te <- which(fold == v); tr <- which(fold != v)
    inc[, te] <- tryCatch(fitPred(tr, te), error = function(e) inc[, te]) }
  inc
}
