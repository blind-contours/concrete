#' Continuous / ordinal patient-reported-outcome (PRO) tiers for the hierarchical
#' win statistics (experimental).
#'
#' A PRO tier (e.g.\ KCCQ, NYHA, 6-minute walk) compares a marker measured at a
#' landmark time among pairs that \emph{reach} the tier --- i.e.\ are tied on every
#' higher-priority (hard-event) tier. In this engine two subjects tie on a
#' time-to-event tier iff neither has that event (continuous times never tie), so a
#' pair reaches a bottom PRO tier iff \strong{both are event-free and alive at the
#' horizon} --- the engine's state-0 occupancy. The per-subject reach probability
#' \eqn{\rho_a(W) = P(\text{state 0 at horizon}\mid W, a)} therefore standardizes
#' the marker distribution by \emph{reach weighting}, not a naive marginal:
#' \deqn{G_a^R(y) = E[\rho_a(W)\,Q_a(y\mid W)] / E[\rho_a(W)],}
#' with \eqn{Q_a(y\mid W) = P(Y \le y \mid W, \text{reach}, a)} the conditional CDF
#' of the marker (estimated by IPCW-weighted binary-threshold Super Learning to
#' correct landmark-observation missingness). The tier win/loss components are the
#' bilinear forms in the un-normalized measures \eqn{H_a(y) = E[\rho_a(W)Q_a(y|W)]}:
#' \deqn{W^{(k)} = \int [E\rho_W - H_W(y+\delta)]\,dH_L(y),\quad
#'       (\text{higher better; } W \text{ beats } L \text{ by margin } \delta).}
#'
#' \strong{Working assumption.} The landmark marker is taken conditionally
#' independent of the post-landmark event process given \eqn{(W, \text{arm})}, so
#' the conditional CDF among reachers equals \eqn{P(Y \le y \mid W)} and is
#' estimable from all observed landmark markers. PRO tiers are supported at the
#' \strong{bottom} of the hierarchy (below all hard-event tiers), which is the
#' clinical norm (soft markers rank below hard events); multiple PRO tiers may be
#' stacked in priority order. A PRO ranked \emph{above} a hard event would require
#' restricting that event's win integral to PRO-tie pairs and is not yet supported.
#'
#' @keywords internal
#' @noRd

## ---- fixed cutpoint grid for a PRO marker ----
## continuous: a regular grid over the pooled observed-marker range; ordinal: the
## sorted unique levels. The grid is part of the estimand (not data-adaptive per
## fit): pass the pooled observed values from BOTH arms so it is stable.
.proCutgrid <- function(markerPooled, type, n.grid = 80L) {
  v <- markerPooled[is.finite(markerPooled)]
  if (identical(type, "ordinal")) return(sort(unique(v)))
  rng <- stats::quantile(v, c(0.005, 0.995), names = FALSE)
  seq(rng[1], rng[2], length.out = as.integer(n.grid))
}

## ---- IPCW-weighted binary-threshold conditional CDF Q_a(c | W) ----
## Fit 1{Y <= c} ~ W among observed-at-landmark subjects (weights = obs IPCW) for
## each cutpoint c, predict for all n subjects, enforce monotonicity across c
## (cumulative max) and [0,1]. SuperLearner per cutpoint with a per-cutpoint glm
## fallback (degenerate / non-converged thresholds).
.proMarkerCDF <- function(marker, obs, ipcw, Cov, cgrid, SL.library, V = 5L) {
  n <- nrow(Cov); G <- length(cgrid); ob <- which(obs)
  Xall <- as.data.frame(Cov); Xtr <- Xall[ob, , drop = FALSE]; wtr <- ipcw[ob]
  haveSL <- requireNamespace("SuperLearner", quietly = TRUE) && length(SL.library)
  glmCut <- function(yb) as.numeric(stats::predict(suppressWarnings(
    stats::glm(yb ~ ., binomial(), data = cbind(yb = yb, Xtr), weights = wtr)),
    newdata = Xall, type = "response"))
  Q <- vapply(cgrid, function(c) {
    yb <- as.integer(marker[ob] <= c)
    if (length(unique(yb)) < 2L) return(rep(mean(yb), n))
    if (haveSL) {
      pr <- tryCatch(suppressWarnings({
        f <- SuperLearner::SuperLearner(Y = yb, X = Xtr, newX = Xall, family = binomial(),
               SL.library = SL.library, obsWeights = wtr, cvControl = list(V = as.integer(V)))
        as.numeric(f$SL.predict) }), error = function(e) NULL)
      if (!is.null(pr) && all(is.finite(pr))) return(pr)
    }
    tryCatch(glmCut(yb), error = function(e) rep(mean(yb), n))
  }, numeric(n))
  Q <- t(apply(Q, 1L, cummax)); pmin(pmax(Q, 0), 1)
}

## ---- per-arm reach-weighted marker measure H_a(.) with influence rows ----
## Uses the engine's state-0 occupancy adjoint (eng$survKilled / eng$survIFset) for
## the reach contribution (carrying terminal reward Q(c|W)), plus the IPCW marker
## residual. `arm` is the armSetup object; `marker`,`obs`,`Cov` are per-subject for
## this arm (marker NA when unobserved; obs = reached landmark & marker recorded).
.proArmMeasure <- function(eng, arm, marker, obs, Cov, cgrid, SL.library, n.folds) {
  n <- arm$n; G <- length(cgrid)
  ## landmark-observation IPCW among reach-eligible (state-0 occupancy > 0): model
  ## P(observed | W) on subjects who could have been measured (have finite marker
  ## opportunity = reached landmark). We approximate the eligible set by "obs or a
  ## positive reach"; with no missingness obs==reach and ipcw==1.
  rho0 <- eng$survKilled(arm$rmat, c(0L), b = 1)[["0"]][1L, ]      # rho_a(W) = state-0 surv to horizon
  elig <- obs | (rho0 > 1e-8)
  pio <- rep(1, n)
  if (any(obs) && !all(obs[elig])) {
    fo <- tryCatch(suppressWarnings(stats::glm(o ~ ., binomial(),
            data = cbind(o = as.integer(obs[elig]), as.data.frame(Cov)[elig, , drop = FALSE]))),
          error = function(e) NULL)
    if (!is.null(fo)) pio <- as.numeric(stats::predict(fo, as.data.frame(Cov), type = "response"))
  }
  ipcw <- ifelse(obs, 1 / pmax(pio, 0.025), 0)
  Q <- .proMarkerCDF(marker, obs, ipcw, Cov, cgrid, SL.library, V = max(2L, n.folds))
  ## reach adjoint: H(c) = E[rho(W) Q(c|W)] via survKilled with terminal reward Q(.,c)
  Vbase <- eng$survKilled(arm$rmat, c(0L), b = 1)
  Erho <- mean(Vbase[["0"]][1L, ])
  ifTot <- (Vbase[["0"]][1L, ] - Erho) +
           eng$survIFset(arm$rmat, arm$YN, arm$Ginv, c(0L), Vbase)
  mLE <- ifelse(is.na(marker), Inf, marker)
  H <- numeric(G); ifH <- matrix(0, n, G)
  for (c in seq_len(G)) {
    Vc <- eng$survKilled(arm$rmat, c(0L), b = Q[, c])
    gi <- Vc[["0"]][1L, ]; H[c] <- mean(gi)
    ifReach <- (gi - H[c]) + eng$survIFset(arm$rmat, arm$YN, arm$Ginv, c(0L), Vc)
    markerResid <- rho0 * ipcw * ((mLE <= cgrid[c]) - Q[, c])
    ifH[, c] <- ifReach + markerResid
  }
  list(H = H, Erho = Erho, ifH = ifH, ifTot = ifTot, n = n, cgrid = cgrid)
}

## ---- PRO tier component "MW beats ML at the tier" (P, IFwin over MW, IFlos over ML) ----
## Matches the shape of one tier of eng$tierComponents. `spec$margin` = delta,
## `spec$direction` in {"higher.better","lower.better"}, `spec$type`.
.proTierP <- function(MW, ML, spec) {
  cg <- ML$cgrid; G <- length(cg); d <- spec$margin
  hi <- !identical(spec$direction, "lower.better")
  if (identical(spec$type, "ordinal")) {                      # exact level sum
    mids <- cg; dHc <- diff(c(0, ML$H))                       # mass at each level (left edge 0)
    shift <- if (hi) cg + d else cg - d
    HWs <- approx(cg, MW$H, xout = shift, rule = 2)$y
    coef <- if (hi) MW$Erho - HWs else HWs
    ifW_at <- t(apply(MW$ifH, 1L, function(r) approx(cg, r, xout = shift, rule = 2)$y))
    ifTotW <- if (hi) MW$ifTot else 0
    P <- sum(coef * dHc)
    IFwin <- as.numeric((ifTotW + (if (hi) -ifW_at else ifW_at)) %*% dHc)
    difC  <- ML$ifH - cbind(0, ML$ifH[, -G, drop = FALSE])
    IFlos <- as.numeric(difC %*% coef)
  } else {                                                    # continuous: midpoint rule
    mids <- (cg[-1] + cg[-G]) / 2; dHc <- diff(ML$H)
    shift <- if (hi) mids + d else mids - d
    HWs <- approx(cg, MW$H, xout = shift, rule = 2)$y
    coef <- if (hi) MW$Erho - HWs else HWs
    ifW_at <- t(apply(MW$ifH, 1L, function(r) approx(cg, r, xout = shift, rule = 2)$y))
    P <- sum(coef * dHc)
    if (hi) IFwin <- as.numeric((MW$ifTot %o% rep(1, G - 1) - ifW_at) %*% dHc)
    else    IFwin <- as.numeric(ifW_at %*% dHc)
    difC  <- ML$ifH[, -1, drop = FALSE] - ML$ifH[, -G, drop = FALSE]
    IFlos <- as.numeric(difC %*% coef)
  }
  P <- P + mean(IFwin) + mean(IFlos)                          # one-step
  list(P = P, IFwin = IFwin, IFlos = IFlos)
}

#' Validate and normalize a PRO-tier specification (or list of them).
#' Each spec: list(marker=, landmark=, margin=, direction=, type=, n.grid=, SL.library=).
#' @keywords internal
#' @noRd
.proNormalize <- function(pro, data, horizon) {
  if (is.null(pro)) return(NULL)
  if (!is.null(pro$marker)) pro <- list(pro)                  # single spec -> list
  lapply(pro, function(s) {
    if (is.null(s$marker) || !s$marker %in% names(data))
      stop("each PRO spec needs a 'marker' column present in data.")
    s$landmark  <- if (is.null(s$landmark)) horizon else s$landmark
    s$margin    <- if (is.null(s$margin)) 0 else s$margin
    s$direction <- match.arg(s$direction %||% "higher.better",
                             c("higher.better", "lower.better"))
    s$type      <- match.arg(s$type %||% "continuous", c("continuous", "ordinal"))
    s$n.grid    <- if (is.null(s$n.grid)) 80L else as.integer(s$n.grid)
    s
  })
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Per-subject "reached landmark & marker observed" flags for one arm.
#' reached = event-free, alive and uncensored through the landmark (engine state 0
#' at the landmark); obs = reached AND a finite marker value is recorded.
#' @keywords internal
#' @noRd
.proArmFlags <- function(Darm, NF, spec) {
  tcols <- if (length(NF)) paste0("t", NF) else character(0); L <- spec$landmark
  reached <- Darm$tD > L & Darm$C > L
  for (tc in tcols) reached <- reached & (Darm[[tc]] > L)
  mk <- Darm[[spec$marker]]
  list(reached = reached, obs = reached & is.finite(mk), marker = ifelse(is.finite(mk), mk, NA_real_))
}

#' Build PRO arm measures + tier components, appended after the event tiers.
#' Returns lists keyed to the win direction the callers need:
#'   winT: P/IFwin(treated)/IFlos(control)  for "treated beats control"
#'   losT: P/IFwin(control)/IFlos(treated)  for "control beats treated"
#' plus the tier labels. `Dtrt`,`Dctl` are the parsed per-arm data.frames (with the
#' marker column carried through), `trt`,`ctl` the armSetup objects.
#' @keywords internal
#' @noRd
.proComponents <- function(eng, pros, Dtrt, Dctl, trt, ctl, covariates, SL.library, n.folds) {
  if (is.null(pros)) return(NULL)
  NF <- eng$NF
  winP <- list(); winW <- list(); winL <- list()
  losP <- list(); losW <- list(); losL <- list(); labs <- character(0)
  for (i in seq_along(pros)) {
    s <- pros[[i]]
    fT <- .proArmFlags(Dtrt, NF, s); fC <- .proArmFlags(Dctl, NF, s)
    cg <- .proCutgrid(c(fT$marker, fC$marker), s$type, s$n.grid)
    MT <- .proArmMeasure(eng, trt, fT$marker, fT$obs, Dtrt[, covariates, drop = FALSE], cg, SL.library, n.folds)
    MC <- .proArmMeasure(eng, ctl, fC$marker, fC$obs, Dctl[, covariates, drop = FALSE], cg, SL.library, n.folds)
    w <- .proTierP(MT, MC, s); l <- .proTierP(MC, MT, s)
    winP[[i]] <- w$P; winW[[i]] <- w$IFwin; winL[[i]] <- w$IFlos
    losP[[i]] <- l$P; losW[[i]] <- l$IFwin; losL[[i]] <- l$IFlos
    labs <- c(labs, if (!is.null(s$label)) s$label else paste0("PRO", i))
  }
  list(winP = unlist(winP), winIFwin = winW, winIFlos = winL,
       losP = unlist(losP), losIFwin = losW, losIFlos = losL, labels = labs)
}
