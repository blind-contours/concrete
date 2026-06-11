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

#' Build the lagged censoring-survival matrix with time-varying covariates, in the
#' core's `LaggedCensSurv` convention (matches `getHazSurvPred()`), for overriding
#' the IPCW throughout the pipeline. `CensoringTV` is a long data.frame with the
#' data's id column (`attr(Data, "ID")`), a `time` column, and value columns.
#' @return a `length(times) x n` matrix of lagged censoring survival.
#' @keywords internal
#' @noRd
#' `Crossover` (optional) is a per-subject switch-time vector (aligned to `Data`
#' rows, `Inf`/`NA` if never): when supplied, a \emph{separate} crossover hazard is
#' fit (same covariates as censoring) and combined with the dropout hazard so the
#' IPCW becomes `1/(S_dropout * S_crossover)` -- the hypothetical "no-crossover"
#' estimand. The outcome must already be re-censored at the switch time (the
#' censored rows whose `Crossover` time matches the observed time are the crossover
#' events; the rest of the censored rows are real dropout).
.tvCensLaggedSurv <- function(Data, CensoringTV = NULL, times, Crossover = NULL,
                              SL.library = c("SL.mean", "SL.glm"), n.folds = 5L, nGridCens = 40L) {
  idcol <- attr(Data, "ID"); typecol <- attr(Data, "EventType"); timecol <- attr(Data, "EventTime")
  ids <- Data[[idcol]]; obsT <- Data[[timecol]]; cens <- Data[[typecol]] <= 0
  covcols <- c(attr(Data, "Treatment"), attr(Data, "CovNames")[["ColName"]])
  baseCov <- as.data.frame(Data[, .SD, .SDcols = covcols])
  ## fit the censoring/crossover hazards on a COARSE regular grid (the eval grid is
  ## every unique event time -> O(n^2) long data; a coarse hazard grid + cumulative-
  ## hazard interpolation is fast and accurate for the IPCW).
  G <- max(1L, as.integer(nGridCens)); coarse <- seq(0, max(times), length.out = G + 1L)
  starts <- coarse[-(G + 1L)]
  tvMats <- list()
  if (!is.null(CensoringTV)) {
    CensoringTV <- as.data.frame(CensoringTV)
    if (!idcol %in% names(CensoringTV)) stop("CensoringTV must contain the id column '", idcol, "'.")
    if (!"time" %in% names(CensoringTV)) stop("CensoringTV must contain a 'time' column.")
    tvMats <- .tvLOCF(ids, CensoringTV, idcol, "time", starts)
  }
  if (is.null(Crossover)) {
    incTotal <- .tvCensoringInc(coarse, obsT, as.integer(cens), baseCov, tvMats, SL.library, n.folds)
  } else {                                           # separate dropout + crossover hazards
    xt <- Crossover; xt[is.na(xt)] <- Inf
    xover <- cens & is.finite(xt) & xt <= obsT + 1e-9   # censored at the switch time
    drop  <- cens & !xover
    incD <- .tvCensoringInc(coarse, obsT, as.integer(drop),  baseCov, tvMats, SL.library, n.folds)
    incX <- .tvCensoringInc(coarse, obsT, as.integer(xover), baseCov, tvMats, SL.library, n.folds)
    incTotal <- incD + incX                           # combined cumulative-hazard increments
  }
  cumHcoarse <- rbind(0, apply(incTotal, 2, cumsum))   # Lambda at coarse nodes, (G+1) x n
  LamFine <- vapply(seq_len(ncol(cumHcoarse)),         # interpolate cumulative hazard onto eval grid
                    function(i) stats::approx(coarse, cumHcoarse[, i], xout = times, rule = 2)$y,
                    numeric(length(times)))
  rbind(1, exp(-LamFine))[seq_along(times), , drop = FALSE]   # lagged censoring survival
}
