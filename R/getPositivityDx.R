#' Positivity / inverse-weight diagnostics for a fitted estimate
#'
#' @description
#' Reports the practical-positivity health of the inverse-probability weights that
#' every `concrete` estimand relies on. The nuisance weight is
#' \eqn{1/(g(A\mid W)\, S_C(t)\, S_X(t))} --- the inverse of the probability of
#' being assigned the regime's treatment \emph{and} remaining uncensored
#' (\eqn{S_C}) \emph{and}, when a crossover model is used, not yet switched
#' (\eqn{S_X}). Because these are multiplied, the denominator can become very small
#' at later times, which (i) inflates the influence-function variance and (ii)
#' triggers truncation that can bias the estimate. This is exactly the regime to
#' watch with informative censoring and crossover.
#'
#' For each intervention it returns the **effective sample size**
#' \eqn{\mathrm{ESS}(t) = (\sum_i w_{it})^2 / \sum_i w_{it}^2} (as a fraction of
#' \eqn{n}), the largest weight, the smallest observation probability (the
#' positivity floor), and the share of weights sitting at the truncation bound ---
#' overall and at the worst time point. Read it alongside any estimate to judge
#' whether the inference is trustworthy or weight-limited.
#'
#' @param ConcreteEst a `"ConcreteEst"` object from [doConcrete()].
#' @param Verbose logical (default TRUE): print a short interpreted summary.
#'
#' @return invisibly, a list with `summary` (one row per intervention x nuisance:
#'   propensity (g), dropout censoring, crossover if modeled, and OVERALL combined
#'   IPCW) and `byTime` (the per-evaluation-time overall ESS fraction and max
#'   weight for each intervention).
#' @export getPositivityDx
#' @examples
#' \dontrun{
#' est <- doConcrete(formatArguments(...))
#' getPositivityDx(est)
#' }
getPositivityDx <- function(ConcreteEst, Verbose = TRUE) {
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("getPositivityDx takes a 'ConcreteEst' object from doConcrete().")
  mn <- attr(ConcreteEst, "MinNuisance")
  arms <- names(ConcreteEst)
  rows <- list(); byTime <- list(); xoverSuspect <- character(0)

  asMat <- function(x) if (is.matrix(x)) x else matrix(x, nrow = 1)
  ## ESS(worst eval time), max weight, min observation probability for a weight matrix
  wstat <- function(w) {
    w <- asMat(w); essT <- rowSums(w)^2 / rowSums(w^2) / ncol(w)
    list(ess = min(essT), maxw = max(w), minp = min(1 / w), n = ncol(w), essT = essT, maxwT = apply(w, 1, max))
  }
  addrow <- function(a, nuis, w, atbound = NA_real_) {
    if (is.null(w)) return(invisible())
    s <- wstat(w)
    rows[[length(rows) + 1L]] <<- data.frame(
      Intervention = a, Nuisance = nuis, n = s$n,
      ESS_overall = round(s$ess, 3), max_weight = round(s$maxw, 2),
      min_obs_prob = signif(s$minp, 3), pct_at_bound = atbound, stringsAsFactors = FALSE)
  }

  for (a in arms) {
    E <- ConcreteEst[[a]]
    g <- E[["PropScore"]]; cS <- E[["CensSurv"]]; xS <- E[["XoverSurv"]]; w <- E[["NuisanceWeight"]]
    nT <- if (!is.null(w) && is.matrix(w)) nrow(w) else 1L
    ## --- propensity (g): the treatment-assignment weight ---
    if (!is.null(g)) addrow(a, "propensity (g)", 1 / as.numeric(g))
    ## --- censoring: dropout (+ crossover) survival weights ---
    if (!is.null(cS)) {
      addrow(a, "dropout censoring", 1 / cS)
    } else if (!is.null(w) && !is.null(g)) {
      ## components not stored: derive the COMBINED censoring weight = (1/w)/g_inv = w * g
      gmat <- matrix(as.numeric(g), nrow = nT, ncol = length(g), byrow = TRUE)
      addrow(a, "censoring (combined)", asMat(w) * gmat)
    }
    if (!is.null(xS)) {
      addrow(a, "crossover", 1 / xS)
      if (max(1 / xS) < 1.01) xoverSuspect <- union(xoverSuspect, a)   # present but ~no reweighting
    }
    ## --- OVERALL combined IPCW (what the inference actually uses) ---
    if (!is.null(w)) {
      atb <- if (is.numeric(mn) && length(mn) == 1L)
               round(100 * mean((1 / asMat(w)) <= mn * (1 + 1e-8)), 1) else 0
      addrow(a, "OVERALL", w, atb)
      s <- wstat(w)
      byTime[[a]] <- data.frame(time_index = seq_len(nT), ESS_frac = round(s$essT, 3),
                                max_weight = round(s$maxwT, 1))
    }
  }
  summ <- do.call(rbind, rows); rownames(summ) <- NULL
  if (isTRUE(Verbose)) {
    cat("Positivity / inverse-weight diagnostics by nuisance\n")
    cat("  ESS = effective sample size (fraction of n) at the worst eval time; lower = more weight-limited.\n")
    cat("  Components: propensity (g), dropout censoring, crossover (if modeled). OVERALL = combined IPCW used for inference.\n\n")
    print(summ, row.names = FALSE)
    ov <- summ[summ$Nuisance == "OVERALL", , drop = FALSE]
    flag <- ov[ov$ESS_overall < 0.5 | (is.finite(ov$pct_at_bound) & ov$pct_at_bound > 5) | ov$max_weight > 20, , drop = FALSE]
    if (nrow(flag))
      cat("\n  CAUTION: low overall ESS / heavy truncation / large weights for: ",
          paste(flag$Intervention, collapse = ", "),
          ". Inference there is weight-limited (near-positivity violation) -- interpret with caution.\n", sep = "")
    if (length(xoverSuspect))
      cat("\n  WARNING: a crossover model is present for [", paste(xoverSuspect, collapse = ", "),
          "] but its weights are ~1 (no reweighting). The crossover/censoring hazard may have failed ",
          "or be intercept-only -- check the time-varying covariates (e.g. unimputed missingness) and learner library.\n", sep = "")
  }
  invisible(list(summary = summ, byTime = byTime))
}
