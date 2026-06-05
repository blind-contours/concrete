#' Specify an ICH E9(R1) estimand and apply an intercurrent-event strategy
#'
#' @description
#' `makeEstimand()` records the analysis target using the five attributes of the
#' ICH E9(R1) estimand framework (treatment condition, population, variable /
#' endpoint, intercurrent-event handling, and population-level summary). It does
#' not change the estimation; it documents the target for a statistical analysis
#' plan and travels with the results for reproducibility and regulatory clarity.
#'
#' `applyIntercurrentEvent()` implements the data handling for the
#' intercurrent-event (ICE) strategy. An ICE is coded as one of the event-type
#' values in the data (the observed time is the earliest of the outcome,
#' competing, censoring, and ICE times). Supported strategies:
#'
#' \describe{
#'   \item{`"treatment policy"`}{The ICE is ignored and follow-up continues to
#'     the outcome -- the default intent-to-treat target. The data are returned
#'     unchanged (the ICE code, if any, is treated as a competing/censoring event
#'     exactly as supplied).}
#'   \item{`"hypothetical"`}{The ICE is recoded as censoring (event type `0`).
#'     `concrete`'s inverse-probability-of-censoring weighting then targets the
#'     risk that would be observed if the ICE had not occurred, **under the
#'     assumption that the ICE acts as conditionally independent censoring given
#'     the measured covariates** -- a strong, untestable assumption that should be
#'     accompanied by a censoring sensitivity analysis.}
#'   \item{`"composite"`}{The ICE is recoded as an occurrence of the event of
#'     interest `TargetEvent`, defining a composite endpoint.}
#' }
#'
#' The `"while on treatment"` and `"principal stratum"` strategies are not
#' supported: the former is a different (descriptive, non-extrapolated) estimand
#' that an IPCW estimator does not target, and the latter requires latent
#' stratification beyond the current scope.
#'
#' @param treatment character: the treatment condition / intervention contrast.
#' @param population character: the target population.
#' @param endpoint character: the variable / endpoint (e.g. cause-specific
#'   absolute risk, RMST) and its horizon.
#' @param summary character: the population-level summary (e.g. risk difference,
#'   risk ratio, RMST difference).
#' @param strategy one of `"treatment policy"`, `"hypothetical"`, `"composite"`.
#' @param intercurrent character or NULL: a label for the intercurrent event.
#'
#' @return `makeEstimand()` returns a `"ConcreteEstimand"` object.
#' @export makeEstimand
makeEstimand <- function(treatment = "as randomized (intent-to-treat)",
                         population = "as enrolled",
                         endpoint = "cause-specific absolute risk",
                         summary = "risk difference",
                         strategy = c("treatment policy", "hypothetical", "composite"),
                         intercurrent = NULL) {
  strategy <- match.arg(strategy)
  out <- list(treatment = treatment, population = population, endpoint = endpoint,
              `intercurrent event strategy` = strategy, intercurrent = intercurrent,
              `population-level summary` = summary)
  class(out) <- "ConcreteEstimand"
  out
}

#' @rdname makeEstimand
#' @param x a `"ConcreteEstimand"` object.
#' @param ... ignored.
#' @exportS3Method print ConcreteEstimand
print.ConcreteEstimand <- function(x, ...) {
  cat("ICH E9(R1) estimand\n")
  cat("  Treatment condition       :", x$treatment, "\n")
  cat("  Population                 :", x$population, "\n")
  cat("  Variable / endpoint        :", x$endpoint, "\n")
  cat("  Intercurrent-event strategy:", x$`intercurrent event strategy`,
      if (!is.null(x$intercurrent)) paste0(" (", x$intercurrent, ")") else "", "\n")
  cat("  Population-level summary    :", x$`population-level summary`, "\n")
  invisible(x)
}

#' @rdname makeEstimand
#' @param Data a `data.table`/`data.frame` of one row per subject.
#' @param EventTime character: name of the observed-time column.
#' @param EventType character: name of the event-type column (`0` = censoring,
#'   positive integers = events).
#' @param Intercurrent numeric: the event-type value in `EventType` that codes
#'   the intercurrent event.
#' @param TargetEvent numeric: for the `"composite"` strategy, the event-of-
#'   interest code that the intercurrent event is merged into.
#' @param Verbose logical: report what was recoded.
#'
#' @return `applyIntercurrentEvent()` returns a copy of `Data` with `EventType`
#'   recoded for the chosen `strategy`, carrying the estimand in
#'   `attr(., "Estimand")`.
#' @export applyIntercurrentEvent
#' @importFrom data.table as.data.table copy
applyIntercurrentEvent <- function(Data, EventTime, EventType, Intercurrent,
                                   strategy = c("treatment policy", "hypothetical", "composite"),
                                   TargetEvent = 1, intercurrent = "intercurrent event",
                                   Verbose = TRUE) {
  strategy <- match.arg(strategy)
  if (identical(strategy, "composite") && length(TargetEvent) != 1L)
    stop("For the composite strategy, TargetEvent must be a single event code ",
         "into which the intercurrent event is merged.")
  DT <- data.table::copy(data.table::as.data.table(Data))
  if (!EventType %in% names(DT)) stop("EventType column '", EventType, "' not found.")
  ev <- DT[[EventType]]
  n_ice <- sum(ev == Intercurrent, na.rm = TRUE)
  if (n_ice == 0 && !identical(strategy, "treatment policy"))
    warning("No rows have EventType == ", Intercurrent,
            "; the intercurrent-event recoding changed nothing.")

  if (identical(strategy, "hypothetical")) {
    DT[get(EventType) == Intercurrent, (EventType) := 0]
    if (Verbose) message(n_ice, " intercurrent-event row(s) recoded as censoring ",
                         "(hypothetical strategy via IPCW).")
  } else if (identical(strategy, "composite")) {
    DT[get(EventType) == Intercurrent, (EventType) := TargetEvent]
    if (Verbose) message(n_ice, " intercurrent-event row(s) merged into event ",
                         TargetEvent, " (composite endpoint).")
  } else if (Verbose) {
    message("Treatment-policy strategy: data returned unchanged.")
  }

  attr(DT, "Estimand") <- makeEstimand(strategy = strategy, intercurrent = intercurrent)
  DT
}
