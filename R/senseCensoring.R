#' Tipping-point sensitivity analysis for informative censoring
#'
#' @description
#' `concrete`'s primary analysis assumes censoring is independent of the event
#' given the measured covariates. `senseCensoring()` probes robustness to
#' departures from that assumption with a transparent tipping-point (bounds)
#' analysis: a fraction `delta` of the subjects who are censored before the
#' target time are assumed to have actually experienced the event of interest,
#' and the analysis is re-fit for each `delta`. `delta = 0` is the optimistic
#' bound (censored subjects never have the event), `delta = 1` the pessimistic
#' bound (all do), and the primary inverse-probability-of-censoring analysis sits
#' between them. The **tipping point** is the smallest `delta` at which the
#' conclusion changes -- a sensitivity artifact recommended by ICH E9(R1).
#'
#' Unlike scaling the censoring weight (which leaves a doubly-robust estimator's
#' target unchanged), imputing the event status of the censored changes the
#' estimand and so produces a genuine, interpretable sensitivity curve. Because
#' it re-fits the estimator for every `delta`, it is computationally heavier than
#' the primary analysis.
#'
#' When the fit carries a treatment-switching (crossover) model (i.e.
#' [formatArguments()] was called with `Crossover`), the censored subjects are a
#' mix of two intercurrent events with two *different* untestable assumptions:
#' ordinary **dropout** (conditionally-independent censoring) and **crossover**
#' (the no-switching counterfactual). `mechanism` selects which pool is imputed,
#' so the two assumptions can be probed individually or jointly:
#' \itemize{
#'   \item `"dropout"` -- tip only the genuinely-censored (drop-out) subjects,
#'     holding the switching handled by the crossover hazard;
#'   \item `"crossover"` -- tip only the subjects re-censored at their switch
#'     time, i.e. probe \dQuote{what if switchers would have had the event had
#'     they not switched}, holding dropout handled by the censoring weight;
#'   \item `"all"` (default) -- tip both pools jointly (the original behaviour).
#' }
#' With no crossover model, all censored subjects are dropout and the three modes
#' coincide.
#'
#' @param ConcreteArgs a `"ConcreteArgs"` object from [formatArguments()].
#' @param deltas numeric in `[0, 1]`: fractions of pre-target-time censored
#'   subjects imputed as having the event of interest. Should include 0.
#' @param Estimand one of `"RD"` (default), `"RR"`, `"Risk"`.
#' @param Intervention length-2 numeric: treatment and control indices.
#' @param mechanism one of `"all"` (default), `"dropout"`, `"crossover"`: which
#'   pool of censored subjects to impute (see Details). `"crossover"` requires a
#'   fit built with `Crossover`.
#' @param Signif numeric (default 0.05): two-sided alpha.
#' @param Verbose logical.
#'
#' @return a `data.table` of estimate / CI / p-value by `delta` x event x time,
#'   with a leading `mechanism` column, the tipping point in
#'   `attr(., "tippingPoint")`, and the mechanism in `attr(., "mechanism")`.
#' @export senseCensoring
senseCensoring <- function(ConcreteArgs, deltas = seq(0, 1, by = 0.25),
                           Estimand = c("RD", "RR", "Risk"), Intervention = c(1, 2),
                           mechanism = c("all", "dropout", "crossover"),
                           Signif = 0.05, Verbose = FALSE) {
  Estimand <- match.arg(Estimand)
  mechanism <- match.arg(mechanism)
  if (!inherits(ConcreteArgs, "ConcreteArgs"))
    stop("ConcreteArgs must be a 'ConcreteArgs' object returned by formatArguments().")
  if (any(deltas < 0 | deltas > 1)) stop("deltas must lie in [0, 1].")
  deltas <- sort(unique(c(0, deltas)))

  Data <- ConcreteArgs[["DataTable"]]
  TimeCol <- attr(Data, "EventTime")
  TypeCol <- attr(Data, "EventType")
  TargetEvent <- ConcreteArgs[["TargetEvent"]]
  TargetTime <- ConcreteArgs[["TargetTime"]]
  tau <- max(TargetTime)
  ev1 <- TargetEvent[1]
  EstLabel <- c(RD = "Risk Diff", RR = "Rel Risk", Risk = "Abs Risk")[[Estimand]]

  ## label each censored row as crossover (re-censored at its switch time) vs
  ## dropout, using the per-subject switch times stashed by formatArguments().
  xt <- attr(Data, "CrossoverTime")
  if (is.null(xt)) {
    if (mechanism == "crossover")
      stop("mechanism = \"crossover\" requires a fit built with a crossover model ",
           "(pass `Crossover` to formatArguments()).")
    isXover <- rep(FALSE, nrow(Data))
  } else {
    isXover <- is.finite(xt) & abs(Data[[TimeCol]] - xt) < 1e-9
  }

  ## subjects censored before the horizon, restricted to the chosen pool and
  ## ordered for a reproducible imputation
  censored <- Data[[TypeCol]] == 0 & Data[[TimeCol]] < tau
  pool <- switch(mechanism,
                 all       = censored,
                 dropout   = censored & !isXover,
                 crossover = censored &  isXover)
  censIdx <- which(pool)
  censIdx <- censIdx[order(Data[[TimeCol]][censIdx])]
  m <- length(censIdx)
  if (m == 0) warning("No ", mechanism, "-censored subjects before the target time; ",
                      "the sensitivity analysis is degenerate for this mechanism.")

  runDelta <- function(delta) {
    # copy the arg environment and the data so the imputation never mutates the
    # caller's ConcreteArgs (data.table assigns by reference); reuse the same
    # folds across delta for comparability.
    A2 <- as.environment(as.list(ConcreteArgs, all.names = TRUE))
    class(A2) <- class(ConcreteArgs)
    D2 <- data.table::copy(Data)
    k <- floor(delta * m)
    if (k > 0) D2[censIdx[seq_len(k)], (TypeCol) := ev1]   # impute as the event of interest
    A2[["DataTable"]] <- D2
    A2[["Verbose"]] <- FALSE
    est <- suppressMessages(suppressWarnings(doConcrete(A2)))
    out <- getOutput(est, Estimand = Estimand, Intervention = Intervention,
                     GComp = FALSE, Simultaneous = FALSE, Signif = Signif)
    out <- data.table::as.data.table(out)[Estimator == "tmle" & Estimand == EstLabel]
    out[, list(delta = delta, n_imputed = k, Event, Time, `Pt Est`, se, `CI Low`, `CI Hi`, pValue)]
  }

  if (Verbose) message("Tipping-point ", mechanism, " sensitivity over ", length(deltas),
                       " values of delta (re-fit each) ...")
  res <- data.table::rbindlist(lapply(deltas, runDelta))
  res[, "mechanism" := mechanism]
  data.table::setcolorder(res, c("mechanism", setdiff(names(res), "mechanism")))

  null0 <- if (Estimand == "RR") 1 else 0
  res[, crosses := (`CI Low` <= null0 & `CI Hi` >= null0)]
  ## reference significance at delta = 0, per (Event, Time) stratum
  sig0 <- res[abs(delta) < 1e-12 & crosses == FALSE, list(Event, Time)]
  tip <- if (nrow(sig0)) {
    ## smallest delta that overturns the conclusion, only within strata that were
    ## significant at delta = 0 (do not borrow a flip from a different stratum)
    flip <- res[crosses == TRUE & delta > 0][sig0, on = c("Event", "Time"), nomatch = 0L]
    if (nrow(flip)) {
      flip[order(Event, Time, delta), list(delta = delta[1L]), by = list(Event, Time)]
    } else {
      data.table::data.table(note = "significant conclusion(s) robust over the delta grid")
    }
  } else {
    data.table::data.table(note = "primary analysis is already non-significant at delta = 0")
  }
  data.table::setattr(res, "tippingPoint", tip)
  data.table::setattr(res, "Estimand", Estimand)
  data.table::setattr(res, "mechanism", mechanism)
  res[]
}
