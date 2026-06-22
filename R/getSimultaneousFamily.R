#' Joint simultaneous inference across a family of estimands
#'
#' @description
#' Trial reports rarely contain a single number: a typical analysis presents a
#' risk difference at several horizons, the RMST difference, and a win ratio
#' together. Pointwise 95% intervals do not control the family-wise error across
#' such a set, and a Bonferroni correction is needlessly conservative because the
#' estimands are computed from the *same* targeted curves and are therefore
#' strongly correlated. `getSimultaneousFamily()` builds **simultaneous**
#' (family-wise) confidence bands across an arbitrary collection of `concrete`
#' estimands by exploiting that correlation.
#'
#' Every `concrete` estimate carries its per-subject efficient influence values.
#' Stacking them into one \eqn{n \times q} matrix gives the joint influence
#' function of the whole family; its empirical correlation \eqn{R} is the
#' correlation of the (asymptotically normal) estimators. The simultaneous
#' critical value is the \eqn{1-\alpha} quantile of \eqn{\max_j |Z_j|} for
#' \eqn{Z \sim N(0, R)} (a multiplier/Gaussian-multiplier bootstrap), and each
#' band is \eqn{\hat\psi \pm q\,\widehat{\mathrm{se}}} (or, for ratio estimands,
#' the same on the log scale). When the family is a single estimand at one time
#' this reduces to the usual Wald interval.
#'
#' This composes with everything else in the package: the influence values used
#' are exactly those reported by [getOutput()], [getRMST()], [targetRMST()],
#' [getWinRatio()], and [targetWinRatio()], so the simultaneous bands inherit the
#' covariate adjustment, censoring correction, cross-fitting, and (when present)
#' the stratified-randomization variance correction.
#'
#' @param ... two or more `"ConcreteOut"` objects produced from the **same**
#'   fitted [doConcrete()] object (so that subjects align). Each may be named;
#'   names are carried into the output `family` column.
#' @param Signif numeric (default 0.05): family-wise alpha.
#' @param nSim integer (default 1e4): Monte Carlo draws for the multiplier
#'   bootstrap critical value.
#'
#' @return a `data.table` with one row per scalar estimand in the family:
#'   `family`, `Estimand`, `Event`, `Time`, `Intervention`, `Pt Est`, `se`, the
#'   pointwise `CI Low` / `CI Hi`, and the simultaneous `SimCI Low` / `SimCI Hi`,
#'   plus the shared simultaneous critical value in `attr(., "critValue")`.
#'
#' @seealso [getOutput()], [getRMST()], [targetRMST()], [getWinRatio()],
#'   [targetWinRatio()].
#' @export getSimultaneousFamily
#' @importFrom stats qnorm quantile cor
getSimultaneousFamily <- function(..., Signif = 0.05, nSim = 1e4L) {
  ekey <- scale <- ic <- ukey <- ID <- `Pt Est` <- se <- NULL
  outs <- list(...)
  if (length(outs) < 1L)
    stop("Supply at least one 'ConcreteOut' object.")
  nm <- names(outs)
  if (is.null(nm)) nm <- rep("", length(outs))
  nm[nm == ""] <- paste0("family", seq_along(outs))[nm == ""]

  estParts <- list(); icParts <- list()
  for (i in seq_along(outs)) {
    o <- outs[[i]]
    fe <- attr(o, "famEst"); fi <- attr(o, "famIC")
    if (is.null(fe) || is.null(fi))
      stop("Object ", i, " carries no family influence functions. ",
           "getSimultaneousFamily() needs output from getOutput(), getRMST(), ",
           "targetRMST(), getWinRatio(), or targetWinRatio().")
    fe <- data.table::copy(fe); fi <- data.table::copy(fi)
    fe[, "family" := nm[i]]; fi[, "family" := nm[i]]
    fe[, "ukey" := paste(nm[i], ekey, sep = "::")]
    fi[, "ukey" := paste(nm[i], ekey, sep = "::")]
    estParts[[i]] <- fe; icParts[[i]] <- fi
  }
  famEst <- data.table::rbindlist(estParts, use.names = TRUE)
  famIC  <- data.table::rbindlist(icParts, use.names = TRUE)

  ## stack influence values into an n x q matrix aligned on subject ID
  wide <- data.table::dcast(famIC, ID ~ ukey, value.var = "ic")
  M <- as.matrix(wide[, setdiff(names(wide), "ID"), with = FALSE])
  if (anyNA(M))
    stop("Influence functions do not align across objects; were they produced ",
         "from the same doConcrete() fit?")
  keep <- apply(M, 2, function(col) stats::sd(col) > 1e-12)
  n <- nrow(M)

  z <- stats::qnorm(1 - Signif / 2)
  if (sum(keep) >= 2L) {
    R <- stats::cor(M[, keep, drop = FALSE])
    draws <- MASS::mvrnorm(n = nSim, mu = rep(0, nrow(R)), Sigma = R)
    crit <- as.numeric(stats::quantile(apply(abs(draws), 1, max), 1 - Signif))
  } else {
    crit <- z                                       # single non-degenerate estimand
  }

  famEst[, c("SimCI Low", "SimCI Hi") := {
    lo <- numeric(.N); hi <- numeric(.N)
    for (r in seq_len(.N)) {
      if (identical(scale[r], "log")) {             # se stored on the log scale
        lo[r] <- `Pt Est`[r] * exp(-crit * se[r]); hi[r] <- `Pt Est`[r] * exp(crit * se[r])
      } else {
        lo[r] <- `Pt Est`[r] - crit * se[r];        hi[r] <- `Pt Est`[r] + crit * se[r]
      }
    }
    list(lo, hi)
  }]
  ## report the natural-scale SE for ratio rows (matching the source outputs)
  famEst[, se := data.table::fifelse(scale == "log", `Pt Est` * se, se)]
  out <- famEst[, c("family", "Estimand", "Event", "Time", "Intervention",
                    "Pt Est", "se", "CI Low", "CI Hi", "SimCI Low", "SimCI Hi"),
                with = FALSE]
  data.table::setattr(out, "critValue", crit)
  data.table::setattr(out, "Signif", Signif)
  data.table::setattr(out, "nEstimands", nrow(out))
  class(out) <- union("ConcreteOut", class(out))
  out[]
}

#' Attach per-subject family influence functions to a ConcreteOut object.
#'
#' `parts` is a list of one entry per scalar estimand, each a list with
#' `key`, `Estimand`, `Event`, `Time`, `Intervention`, `est`, `se`, `scale`
#' (\code{"identity"} or \code{"log"}; for \code{"log"}, `se` is on the log
#' scale and `ic` is the log-scale influence function), and `ic` (a numeric
#' vector). `ids` are the subject IDs aligned with every `ic`. Also records the
#' pointwise CI consistent with the source output.
#' @keywords internal
#' @noRd
.attachFamily <- function(out, ids, parts) {
  sig <- attr(out, "Signif"); if (is.null(sig) || !is.numeric(sig) || length(sig) != 1) sig <- 0.05
  z <- stats::qnorm(1 - sig / 2)              # pointwise CI at the source output's alpha
  est <- data.table::rbindlist(lapply(parts, function(p) {
    ciLo <- if (identical(p$scale, "log")) p$est * exp(-z * p$se) else p$est - z * p$se
    ciHi <- if (identical(p$scale, "log")) p$est * exp(z * p$se)  else p$est + z * p$se
    data.table::data.table(ekey = p$key, Estimand = p$Estimand,
                           Event = as.character(p$Event), Time = p$Time,
                           Intervention = p$Intervention, `Pt Est` = p$est,
                           se = p$se, scale = p$scale,
                           `CI Low` = ciLo, `CI Hi` = ciHi)
  }))
  ic <- data.table::rbindlist(lapply(parts, function(p)
    data.table::data.table(ekey = p$key, ID = ids, ic = as.numeric(p$ic))))
  data.table::setattr(out, "famEst", est)
  data.table::setattr(out, "famIC", ic)
  out
}
