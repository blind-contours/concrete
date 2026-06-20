#' Priority-standardized net benefit and win ratio (charter-weighted)
#'
#' @description
#' Standard hierarchical win statistics aggregate the layer-specific effects using
#' \emph{reach} weights --- the probability that a treated--control pair is tied on
#' all higher-priority layers and is therefore decided at a given layer. Those
#' weights are set by the outcome distribution and the tie/censoring conventions,
#' not by clinical priority, so a frequently-reached (often low-priority) layer can
#' numerically dominate the composite. The priority-standardized net benefit
#' (PSNB) and win ratio (PSWR) replace the implicit reach weights with a
#' prespecified clinical \strong{charter} \eqn{\alpha} (one weight per layer,
#' \eqn{\alpha_k \ge 0}, \eqn{\sum_k \alpha_k = 1}), fixed before unblinding:
#' \deqn{\mathrm{PSNB}(\alpha) = \sum_k \alpha_k \Delta_k,
#'   \qquad \mathrm{PSWR}(\alpha) = \frac{\sum_k \alpha_k w_k}{\sum_k \alpha_k \ell_k},}
#' where \eqn{\Delta_k = w_k - \ell_k} is the stage-conditional net benefit and
#' \eqn{w_k}, \eqn{\ell_k} are the win / loss probabilities among pairs that reach
#' layer \eqn{k}.
#'
#' This is the \strong{covariate-adjusted, doubly-robust, IPCW-censoring-corrected}
#' estimator of PSNB/PSWR for a death-priority time-to-event hierarchy: it reuses
#' the same validated multistate engine as [clinicalWinRatio()], which already
#' produces the per-tier win/loss components \eqn{W^{(k)} = r_k w_k} and
#' \eqn{L^{(k)} = r_k \ell_k}. The standard win ratio sums those tiers; PSNB
#' divides the reach \eqn{r_k = 1 - \sum_{m<k}(W^{(m)}+L^{(m)})} back out and
#' recombines with the charter. Influence-function inference is propagated by the
#' delta method from the engine's per-tier influence functions.
#'
#' @inheritParams clinicalWinRatio
#' @param charter the priority charter: a numeric vector of length \eqn{K} (the
#'   number of layers, layer 1 = death, followed by any non-fatal event tiers and
#'   then any bottom `pro` tiers in order), giving the weight on each layer
#'   (rescaled to sum to 1). The special value `"reach"` uses the realized reach
#'   weights,
#'   which reproduces the standard net benefit and win ratio (useful as a
#'   reference / sanity check). The charter must be \strong{prespecified}: it is
#'   part of the estimand, not a tuning parameter.
#'
#' @return a `data.table` of class `"ConcreteOut"` with `PSNB` and `PSWR` (each
#'   with an influence-function CI; PSWR on the log scale), plus, for each layer,
#'   the reach probability and the stage-conditional net benefit \eqn{\Delta_k} so
#'   the user can see which layer drives the composite.
#'
#' @seealso [clinicalWinRatio()] (the standard reach-weighted win ratio),
#'   [clinicalRMTIF()].
#' @export clinicalPSNB
#' @importFrom stats qnorm pnorm
clinicalPSNB <- function(data, arm, illness.time, terminal.time, terminal.status,
                         covariates, charter, horizon = NULL, n.grid = 60L, n.folds = 5L,
                         SL.library = c("SL.mean", "SL.glm"), Signif = 0.05,
                         id = NULL, censoring.tv = NULL, pro = NULL) {
  data <- as.data.frame(data)
  illness.time <- as.character(illness.time)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  if (!all(delta %in% c(0, 1))) stop("terminal.status must be coded 0/1 (1 = death).")
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  Kev <- 1L + length(illness.time)                     # number of hard-event tiers
  pros <- .proNormalize(pro, data, horizon)            # bottom PRO tiers (or NULL)
  K <- Kev + length(pros)                              # total tiers (charter length)
  useReach <- is.character(charter) && identical(charter[1], "reach")
  if (!useReach) {
    if (!is.numeric(charter) || length(charter) != K)
      stop("charter must be a numeric vector of length K = ", K,
           " (one weight per layer, layer 1 = death), or \"reach\".")
    if (any(charter < 0)) stop("charter weights must be non-negative.")
    if (sum(charter) <= 0) stop("charter weights must not all be zero.")
    alpha <- charter / sum(charter)
  }
  grid <- seq(0, horizon, length.out = as.integer(n.grid) + 1L)
  eng <- .msEngine(Kev, grid)                          # engine spans the hard-event tiers
  proCols <- unique(vapply(pros, function(s) s$marker, character(1)))

  parseArm <- function(rows) {
    D <- data[rows, covariates, drop = FALSE]
    D$tD <- ifelse(delta[rows] == 1, term[rows], Inf)
    for (ei in seq_along(illness.time)) {
      ti <- data[[illness.time[ei]]][rows]; ti[is.na(ti)] <- Inf; D[[paste0("t", ei)]] <- ti }
    D$C <- ifelse(delta[rows] == 0, term[rows], Inf)
    for (mc in proCols) D[[mc]] <- data[[mc]][rows]    # carry PRO markers through
    D
  }
  tvMats <- NULL
  if (!is.null(censoring.tv)) {
    if (is.null(id)) stop("`id` is required when `censoring.tv` is supplied.")
    censoring.tv <- as.data.frame(censoring.tv)
    tvMats <- .tvLOCF(data[[id]], censoring.tv, id, "time", grid[-length(grid)])
  }
  fitArm <- function(rows) {
    D <- parseArm(rows)
    tvA <- if (is.null(tvMats)) NULL else lapply(tvMats, function(m) m[rows, , drop = FALSE])
    nu <- .msNuisances(eng, D, covariates, SL.library, n.folds, tvA)
    list(arm = eng$armSetup(D, nu$rmat, nu$Ginv), D = D)
  }
  fT <- fitArm(which(A == 1)); fC <- fitArm(which(A == 0))
  trt <- fT$arm; ctl <- fC$arm
  nT <- trt$n; nC <- ctl$n; Ntot <- nT + nC; piT <- nT / Ntot; piC <- nC / Ntot
  z <- stats::qnorm(1 - Signif / 2)

  ## per-tier win/loss components + per-arm influence functions (hard-event tiers)
  winT <- eng$tierComponents(trt, ctl)   # W^{(k)}: IFwin over treated, IFlos over control
  losT <- eng$tierComponents(ctl, trt)   # L^{(k)}: IFwin over control, IFlos over treated
  Wk <- winT$P; Lk <- losT$P
  DWk_T <- winT$IFwin; DWk_C <- winT$IFlos        # treated / control IFs of W^{(k)}
  DLk_C <- losT$IFwin; DLk_T <- losT$IFlos        # control / treated IFs of L^{(k)}
  ## append bottom PRO tiers (reach-weighted marker comparison)
  proLab <- character(0)
  if (!is.null(pros)) {
    pc <- .proComponents(eng, pros, fT$D, fC$D, trt, ctl, covariates, SL.library, n.folds)
    ## rescale the PRO block to the hard-tier residual reach (coherent hierarchy)
    s <- max(1e-6, 1 - sum(Wk) - sum(Lk)) / max(pc$reachEmp, 1e-6)
    Wk <- c(Wk, s * pc$winP); Lk <- c(Lk, s * pc$losP)
    DWk_T <- c(DWk_T, lapply(pc$winIFwin, `*`, s)); DWk_C <- c(DWk_C, lapply(pc$winIFlos, `*`, s))
    DLk_C <- c(DLk_C, lapply(pc$losIFwin, `*`, s)); DLk_T <- c(DLk_T, lapply(pc$losIFlos, `*`, s))
    proLab <- pc$labels
  }

  ## reach, stage-conditional w/l, and their influence functions (delta method)
  zT <- numeric(nT); zC <- numeric(nC)
  rk <- numeric(K); wk <- numeric(K); lk <- numeric(K)
  IFw_T <- IFw_C <- IFl_T <- IFl_C <- vector("list", K)
  cumW_T <- zT; cumW_C <- zC; cumWpt <- 0                    # running sum of decided mass below k
  for (k in seq_len(K)) {
    rk[k] <- 1 - cumWpt
    IFr_T <- -cumW_T; IFr_C <- -cumW_C                       # IF of reach r_k
    wk[k] <- Wk[k] / rk[k]; lk[k] <- Lk[k] / rk[k]
    IFw_T[[k]] <- (DWk_T[[k]] - wk[k] * IFr_T) / rk[k]
    IFw_C[[k]] <- (DWk_C[[k]] - wk[k] * IFr_C) / rk[k]
    IFl_T[[k]] <- (DLk_T[[k]] - lk[k] * IFr_T) / rk[k]
    IFl_C[[k]] <- (DLk_C[[k]] - lk[k] * IFr_C) / rk[k]
    cumWpt <- cumWpt + Wk[k] + Lk[k]
    cumW_T <- cumW_T + DWk_T[[k]] + DLk_T[[k]]
    cumW_C <- cumW_C + DWk_C[[k]] + DLk_C[[k]]
  }
  if (useReach) alpha <- rk      # alpha_k = realized reach -> reproduces standard NB / WR

  seGrad <- function(IF_T, IF_C) sqrt((sum(((1/piT)*IF_T)^2) + sum(((1/piC)*IF_C)^2)) / Ntot^2)
  Dk <- wk - lk                                              # stage-conditional net benefit
  ## PSNB = sum alpha_k Delta_k
  psnb <- sum(alpha * Dk)
  IFpsnb_T <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*(IFw_T[[k]] - IFl_T[[k]])))
  IFpsnb_C <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*(IFw_C[[k]] - IFl_C[[k]])))
  sePsnb <- seGrad(IFpsnb_T, IFpsnb_C)
  ## PSWR = wbar / lbar
  wbar <- sum(alpha * wk); lbar <- sum(alpha * lk); pswr <- wbar / lbar
  IFwbar_T <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*IFw_T[[k]]))
  IFwbar_C <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*IFw_C[[k]]))
  IFlbar_T <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*IFl_T[[k]]))
  IFlbar_C <- Reduce(`+`, lapply(seq_len(K), function(k) alpha[k]*IFl_C[[k]]))
  IFlogwr_T <- IFwbar_T/wbar - IFlbar_T/lbar; IFlogwr_C <- IFwbar_C/wbar - IFlbar_C/lbar
  slwr <- seGrad(IFlogwr_T, IFlogwr_C)

  rows <- list(
    data.table::data.table(Estimand = "PSNB", `Pt Est` = psnb, se = sePsnb,
      `CI Low` = psnb - z*sePsnb, `CI Hi` = psnb + z*sePsnb, pValue = 2*stats::pnorm(-abs(psnb/sePsnb))),
    data.table::data.table(Estimand = "PSWR", `Pt Est` = pswr, se = pswr*slwr,
      `CI Low` = pswr*exp(-z*slwr), `CI Hi` = pswr*exp(z*slwr), pValue = 2*stats::pnorm(-abs(log(pswr)/slwr))))
  tierLab <- c("D", eng$NF, proLab)
  for (k in seq_len(K)) {
    se_dk <- seGrad(IFw_T[[k]] - IFl_T[[k]], IFw_C[[k]] - IFl_C[[k]])
    rows[[length(rows)+1L]] <- data.table::data.table(
      Estimand = paste0("Reach[", tierLab[k], "]"), `Pt Est` = rk[k], se = NA_real_,
      `CI Low` = NA_real_, `CI Hi` = NA_real_, pValue = NA_real_)
    rows[[length(rows)+1L]] <- data.table::data.table(
      Estimand = paste0("NetBenefit[", tierLab[k], "]"), `Pt Est` = Dk[k], se = se_dk,
      `CI Low` = Dk[k] - z*se_dk, `CI Hi` = Dk[k] + z*se_dk, pValue = NA_real_)
  }
  Output <- data.table::rbindlist(rows)
  Output[, `:=`(Intervention = "[arm=1] vs [arm=0]", Estimator = "tmle",
                Event = paste(tierLab, collapse = ">"), Time = horizon)]
  data.table::setcolorder(Output, c("Intervention", "Estimand", "Estimator", "Event",
                                    "Time", "Pt Est", "se", "CI Low", "CI Hi", "pValue"))
  attr(Output, "Signif") <- Signif; attr(Output, "Horizon") <- horizon
  attr(Output, "Estimand") <- "PSNB"; attr(Output, "Tiers") <- K
  attr(Output, "charter") <- if (useReach) "reach" else alpha
  attr(Output, "Experimental") <- TRUE
  class(Output) <- union("ConcreteOut", class(Output))
  Output[]
}
