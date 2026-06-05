#' Relative efficiency of a covariate-adjusted vs unadjusted analysis
#'
#' @description
#' In a randomized trial, covariate adjustment does not change the target
#' estimand but it can sharpen it. `getRelativeEfficiency()` quantifies that
#' precision gain by comparing the influence-function standard errors of an
#' adjusted analysis to those of an unadjusted analysis of the **same** estimand,
#' which is the quantity the FDA's 2023 covariate-adjustment guidance is about.
#'
#' Supply two `"ConcreteOut"` tables (from [getOutput()] or [getRMST()]): one
#' from a covariate-adjusted fit and one from an unadjusted fit. The unadjusted
#' fit is the same workflow with treatment-only nuisance models, e.g. a
#' marginal propensity (`"SL.mean"`) and hazard formulas of the form
#' `Surv(time, event == j) ~ arm`.
#'
#' For each matched estimand the function reports the relative efficiency
#' \eqn{\mathrm{RE} = \mathrm{Var}_{\text{unadj}} / \mathrm{Var}_{\text{adj}}}
#' (values above 1 favor adjustment), the implied percentage variance reduction,
#' and the effective sample-size multiplier: an adjusted analysis on \eqn{n}
#' subjects has the precision of an unadjusted analysis on \eqn{\mathrm{RE}\,n}.
#'
#' @param Adjusted a `"ConcreteOut"` table from a covariate-adjusted fit.
#' @param Unadjusted a `"ConcreteOut"` table from an unadjusted (treatment-only)
#'   fit of the same estimands, interventions, events, and times.
#'
#' @return a `data.table` keyed by `Intervention`, `Estimand`, `Event`, and
#'   `Time` with columns `seAdjusted`, `seUnadjusted`, `RelEfficiency`,
#'   `VarReductionPct`, and `EffSampleSizeMult`.
#'
#' @seealso [getOutput()], [getRMST()]
#' @export getRelativeEfficiency
getRelativeEfficiency <- function(Adjusted, Unadjusted) {
  Estimator <- se <- seAdjusted <- seUnadjusted <- RelEfficiency <- NULL
  if (!inherits(Adjusted, "ConcreteOut") || !inherits(Unadjusted, "ConcreteOut"))
    stop("Adjusted and Unadjusted must both be 'ConcreteOut' tables returned by ",
         "getOutput() or getRMST().")

  keys <- c("Intervention", "Estimand", "Event", "Time")
  prep <- function(x, nm) {
    x <- data.table::as.data.table(x)
    if ("Estimator" %in% names(x)) x <- x[Estimator == "tmle"]
    miss <- setdiff(c(keys, "se"), names(x))
    if (length(miss))
      stop("The ", nm, " table is missing required column(s): ",
           paste(miss, collapse = ", "))
    x <- x[is.finite(se) & se > 0, c(keys, "se"), with = FALSE]
    data.table::setnames(x, "se", if (nm == "Adjusted") "seAdjusted" else "seUnadjusted")
    x
  }
  adj <- prep(Adjusted, "Adjusted")
  unadj <- prep(Unadjusted, "Unadjusted")

  out <- merge(adj, unadj, by = keys)
  if (!nrow(out))
    stop("No matching estimand rows (Intervention/Estimand/Event/Time) with ",
         "usable standard errors were found between the two tables.")
  if (nrow(out) < min(nrow(adj), nrow(unadj)))
    warning("Some estimand rows did not match between the adjusted and unadjusted ",
            "tables and were dropped; check that both used the same estimands, ",
            "interventions, events, and times.")

  out[, RelEfficiency := (seUnadjusted^2) / (seAdjusted^2)]
  out[, "VarReductionPct" := 100 * (1 - 1 / RelEfficiency)]
  out[, "EffSampleSizeMult" := RelEfficiency]
  data.table::setorderv(out, c("Estimand", "Event", "Time", "Intervention"))
  out[]
}
