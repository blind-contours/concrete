#' Continuous / ordinal patient-reported-outcome (PRO) tiers for the hierarchical
#' win statistics (experimental).
#'
#' PRO tiers (e.g.\ KCCQ, NYHA, 6-minute walk) compare markers measured at a final
#' landmark visit among pairs that \emph{reach} the tier --- i.e.\ are tied on every
#' higher-priority (hard-event) tier. In this engine two subjects tie on a
#' time-to-event tier iff neither has that event (continuous times never tie), so a
#' pair reaches the PRO block iff \strong{both are event-free, alive and in
#' follow-up at the horizon}. Within that reach the PRO tiers are compared
#' \strong{sequentially} (KCCQ, then NYHA, then 6-minute walk, ...): a pair is
#' decided at the first PRO tier where the two markers differ by more than that
#' tier's win margin \eqn{\delta}; otherwise it descends to the next PRO tier. This
#' is the generalized pairwise comparison (GPC) that the win ratio performs, and it
#' is exactly the TRISCEND II construction (Hahn et al., NEJM 2025).
#'
#' The PRO block is estimated by a \strong{reach-weighted, IPCW-corrected
#' two-sample GPC} on the joint marker vectors: pairs are restricted to reachers,
#' weighted by inverse censoring survival \eqn{1/G(\tau\mid W)} (so an observed
#' reacher stands in for everyone who would have been event-free at the horizon
#' absent censoring) and by inverse landmark-visit-attendance probability (a per-arm
#' logistic model), and the comparison uses the actual paired marker values, so the
#' markers may be arbitrarily correlated (no marker-independence assumption, and the
#' sequential tie-passing among PRO tiers is exact). Influence-function inference is
#' the two-sample U-statistic (Hajek) projection.
#'
#' \strong{Scope.} The PRO block requires the landmark = the horizon (the standard
#' final-visit QoL design, and TRISCEND II); reach is then event-free survival to
#' the horizon. Hard-event tiers above the PRO block keep the engine's
#' covariate-adjusted, doubly-robust influence-function inference. A PRO ranked
#' \emph{above} a hard event is not supported (it would require restricting that
#' event's win integral to PRO-tie pairs).
#'
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Validate and normalize a PRO-tier specification (or list of them).
#' Each spec: list(marker=, landmark=, margin=, direction=, type=, label=).
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
    s
  })
}

#' Reach-weighted, IPCW-corrected sequential GPC for the PRO block of the
#' hierarchy. Returns, for each PRO tier, the win/loss probabilities and their
#' per-subject influence functions over the treated and control arms --- the same
#' (P, IFwin, IFlos) shape one hard-event tier produces, so the PRO tiers append to
#' the engine's tier components and the WR / PSNB assembly is unchanged.
#'
#' `Dtrt`,`Dctl` are the parsed per-arm data.frames (carrying the marker columns,
#' plus `tD`, `t{1..}`, `C`); `trt`,`ctl` the armSetup objects (unused here beyond
#' the arm sizes; the hard-event reach is recovered from the event/censoring times).
#' @keywords internal
#' @noRd
.proComponents <- function(eng, pros, Dtrt, Dctl, trt, ctl, covariates, SL.library, n.folds) {
  if (is.null(pros)) return(NULL)
  H <- eng$tau; NF <- eng$NF; J <- length(pros)
  if (any(abs(vapply(pros, function(s) s$landmark, numeric(1)) - H) > 1e-8))
    stop("PRO tiers require landmark = horizon (final-visit design). ",
         "Set each PRO 'landmark' to the analysis horizon.")
  tcols <- if (length(NF)) paste0("t", NF) else character(0)
  markers <- vapply(pros, function(s) s$marker, character(1))

  armBlock <- function(Darm, armObj) {
    n <- nrow(Darm)
    reached <- Darm$tD > H & Darm$C >= H                      # event-free, alive, in follow-up at horizon
    for (tc in tcols) reached <- reached & (Darm[[tc]] > H)
    Y <- matrix(vapply(markers, function(m) as.numeric(Darm[[m]]), numeric(n)), n, J)
    finite <- rowSums(!is.finite(Y)) == 0L
    blockobs <- reached & finite
    Cov <- as.data.frame(Darm[, covariates, drop = FALSE])
    pio <- rep(1, n)                                          # P(landmark visit observed | reached, W)
    if (any(reached) && !all(blockobs[reached])) {
      fo <- tryCatch(suppressWarnings(stats::glm(o ~ ., binomial(),
              data = cbind(o = as.integer(blockobs[reached]), Cov[reached, , drop = FALSE]))),
            error = function(e) NULL)
      if (!is.null(fo)) pio <- as.numeric(stats::predict(fo, Cov, type = "response"))
    }
    ## inverse-censoring-survival to the horizon, 1/G(H | W): a reacher (uncensored
    ## to H) stands in for everyone who would be event-free at H absent censoring.
    ## The engine carries 1/G as Ginv (n x M); its last column is 1/G(H^-).
    Gci <- armObj$Ginv[, ncol(armObj$Ginv)]
    u <- ifelse(blockobs, Gci / pmax(pio, 0.025), 0)
    list(idx = which(blockobs), Y = Y, u = u, n = n)
  }
  BT <- armBlock(Dtrt, trt); BC <- armBlock(Dctl, ctl)
  nT <- BT$n; nC <- BC$n; iT <- BT$idx; iC <- BC$idx
  uT <- BT$u[iT]; uC <- BC$u[iC]
  YT <- BT$Y[iT, , drop = FALSE]; YC <- BC$Y[iC, , drop = FALSE]
  mT <- length(iT); mC <- length(iC)

  winP <- losP <- numeric(J)
  winIFwin <- winIFlos <- losIFwin <- losIFlos <- vector("list", J)
  undec <- matrix(TRUE, mT, mC)
  scat <- function(vals, idx, nn, w, P) { v <- numeric(nn); v[idx] <- w * vals; v - P }
  for (j in seq_len(J)) {
    s <- pros[[j]]; d <- s$margin; hi <- !identical(s$direction, "lower.better")
    dif <- outer(YT[, j], YC[, j], "-")                      # mT x mC : Y_T - Y_C
    wj <- undec & (if (hi) dif > d else dif < -d)            # treated wins tier j
    lj <- undec & (if (hi) dif < -d else dif > d)            # treated loses tier j
    Wj <- as.numeric(crossprod(uT, wj %*% uC)) / (nT * nC)
    Lj <- as.numeric(crossprod(uT, lj %*% uC)) / (nT * nC)
    winP[j] <- Wj; losP[j] <- Lj
    ## two-sample U-statistic (Hajek) projection IFs
    winIFwin[[j]] <- scat(as.numeric(wj %*% uC) / nC, iT, nT, uT, Wj)        # over treated (W winner)
    winIFlos[[j]] <- scat(as.numeric(crossprod(wj, uT)) / nT, iC, nC, uC, Wj) # over control (W loser)
    losIFwin[[j]] <- scat(as.numeric(crossprod(lj, uT)) / nT, iC, nC, uC, Lj) # over control (L winner)
    losIFlos[[j]] <- scat(as.numeric(lj %*% uC) / nC, iT, nT, uT, Lj)        # over treated (L loser)
    undec <- undec & !(wj | lj)
  }
  labs <- vapply(seq_len(J), function(j)
    if (!is.null(pros[[j]]$label)) pros[[j]]$label else paste0("PRO", j), character(1))
  list(winP = winP, winIFwin = winIFwin, winIFlos = winIFlos,
       losP = losP, losIFwin = losIFwin, losIFlos = losIFlos, labels = labs)
}
