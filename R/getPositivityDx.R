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
#' @return invisibly, a list with `summary` (one row per intervention) and
#'   `byTime` (the per-evaluation-time ESS fraction, max weight, and minimum
#'   observation probability for each intervention).
#' @export getPositivityDx
#' @examples
#' \dontrun{
#' est <- doConcrete(formatArguments(...))
#' getPositivityDx(est)
#' }
getPositivityDx <- function(ConcreteEst, Verbose = TRUE) {
  if (!inherits(ConcreteEst, "ConcreteEst"))
    stop("getPositivityDx takes a 'ConcreteEst' object from doConcrete().")
  arms <- names(ConcreteEst)
  summ <- list(); byTime <- list()
  for (a in arms) {
    w <- ConcreteEst[[a]][["NuisanceWeight"]]
    if (is.null(w) || (!is.matrix(w) && length(w) == 1L)) {     # no censoring -> weights ~ 1
      summ[[a]] <- data.frame(Intervention = a, n = NA_integer_, ESS_overall = 1,
        ESS_worst = 1, max_weight = 1, min_obs_prob = 1, pct_at_bound = 0)
      next
    }
    if (!is.matrix(w)) w <- matrix(w, nrow = 1)
    n <- ncol(w)
    essT  <- rowSums(w)^2 / rowSums(w^2) / n                    # ESS fraction at each eval time
    maxwT <- apply(w, 1, max)
    denom <- 1 / w                                             # truncated observation probability
    minpT <- apply(denom, 1, min)
    floorVal <- min(denom)                                     # truncation bound (MinNuisance)
    atBound <- denom <= floorVal * (1 + 1e-8)
    byTime[[a]] <- data.frame(time_index = seq_len(nrow(w)), ESS_frac = round(essT, 3),
                              max_weight = round(maxwT, 1), min_obs_prob = signif(minpT, 3))
    summ[[a]] <- data.frame(Intervention = a, n = n,
      ESS_overall = round(min(essT), 3),                        # whole-window ESS = its worst point
      ESS_worst = round(min(essT), 3),
      max_weight = round(max(w), 1),
      min_obs_prob = signif(floorVal, 3),
      pct_at_bound = round(100 * mean(atBound), 1))
  }
  summ <- do.call(rbind, summ); rownames(summ) <- NULL
  if (isTRUE(Verbose)) {
    cat("Positivity / inverse-weight diagnostics\n")
    cat("  (ESS = effective sample size as a fraction of n; lower = more weight-limited)\n\n")
    print(summ, row.names = FALSE)
    flag <- summ[summ$ESS_worst < 0.5 | summ$pct_at_bound > 5 | summ$max_weight > 20, , drop = FALSE]
    if (nrow(flag)) {
      cat("\n  CAUTION: low ESS / heavy truncation / large weights for: ",
          paste(flag$Intervention, collapse = ", "), ".\n", sep = "")
      cat("  Inference there is weight-limited (often near-positivity violation at later\n",
          "  times). Consider a shorter horizon, fewer/Stabler censoring-or-crossover\n",
          "  covariates, or interpret with caution.\n", sep = "")
    }
  }
  invisible(list(summary = summ, byTime = byTime))
}
