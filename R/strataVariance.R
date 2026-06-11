#' Variance under stratified / covariate-adaptive randomization
#'
#' @description
#' Influence-function standard errors assume simple (iid) randomization. Under
#' covariate-adaptive schemes that achieve strong balance -- stratified permuted
#' blocks, stratified biased coin, minimization (under conditions) -- the iid
#' variance is generically conservative: it includes a between-arm-within-stratum
#' component that the randomization scheme removes (Bugni, Canay & Shaikh 2018,
#' JASA; Ye, Yi & Shao 2022, Biometrika; Ye, Shao, Yi & Zhao 2023, JASA).
#'
#' For an asymptotically linear estimator with influence function psi, target
#' allocation pi, strata S with frequencies p_s, and
#' \eqn{\Delta(s) = E(\psi \mid A=1, S=s) - E(\psi \mid A=0, S=s)}:
#' \deqn{\sigma^2_{CAR} = \sigma^2_{iid} - \pi(1-\pi)\sum_s p_s \Delta(s)^2.}
#' The estimator below computes the equivalent residual + between-stratum form
#' \deqn{\hat\sigma^2 = P_n[(\psi - \bar\psi_{A,S})^2] +
#'       \sum_s \hat p_s (\hat m(s) - \bar m)^2,}
#' with \eqn{\hat m(s) = \hat\pi \bar\psi_{1,s} + (1-\hat\pi)\bar\psi_{0,s}},
#' which is nonnegative by construction. Derivation: notes/strata-variance.md.
#'
#' Because every reported estimand (absolute risk, RD, RR, RMST/LYL, win ratio)
#' stores per-subject influence values, supplying `Strata` to [formatArguments()]
#' corrects the standard errors of all of them through this one helper. When the
#' working models adjust for the stratification variables (recommended; the
#' strata columns stay in the data as covariates), Delta(s) is approximately 0
#' and the correction is approximately 0 -- the iid variance is then already
#' asymptotically correct.
#'
#' @param IC numeric: per-subject influence values (mean approximately 0).
#' @param A numeric/integer: per-subject randomized arm (binary, > 0 = treated).
#' @param S per-subject randomization-stratum labels.
#' @return the adjusted variance of the influence function (so
#'   `se = sqrt(. / n)`), or `NULL` when a stratum-arm cell has fewer than 2
#'   subjects (caller falls back to the iid variance).
#' @keywords internal
#' @name strataVariance
.strataAdjSigma2 <- function(IC, A, S) {
  if (length(IC) != length(A) || length(A) != length(S) || anyNA(IC)) return(NULL)
  Atrt <- as.integer(A > 0)
  S <- as.character(S)
  cellN <- table(Atrt, S)
  if (nrow(cellN) < 2L || any(cellN < 2L)) return(NULL)
  pihat <- mean(Atrt)
  cell <- paste0(Atrt, ".", S)
  mAS <- tapply(IC, cell, mean)
  eps <- IC - as.numeric(mAS[cell])
  ps <- prop.table(table(S))
  sLev <- names(ps)
  m1 <- tapply(IC[Atrt == 1L], S[Atrt == 1L], mean)[sLev]
  m0 <- tapply(IC[Atrt == 0L], S[Atrt == 0L], mean)[sLev]
  ms <- pihat * m1 + (1 - pihat) * m0
  as.numeric(mean(eps^2) + sum(ps * (ms - sum(ps * ms))^2))
}

#' Strata-adjusted SE for a per-subject influence vector keyed by ID; NULL when
#' the correction is unavailable (unmatched IDs or degenerate strata), in which
#' case the caller keeps the iid SE.
#' @keywords internal
#' @noRd
.strataSE <- function(IC, IDs, StrataDT) {
  if (is.null(StrataDT)) return(NULL)
  idx <- match(IDs, StrataDT[["ID"]])
  if (anyNA(idx)) return(NULL)
  s2 <- .strataAdjSigma2(IC, StrataDT[["A"]][idx], StrataDT[["S"]][idx])
  if (is.null(s2)) return(NULL)
  sqrt(s2 / length(IC))
}
