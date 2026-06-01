getEICStopRule <- function(EICStopRule) {
  choices <- c("relative", "absolute", "hybrid")

  if (is.null(EICStopRule) || length(EICStopRule) == 0) {
    return("relative")
  }
  if (length(EICStopRule) > 1) {
    EICStopRule <- EICStopRule[1]
  }
  if (!is.character(EICStopRule) || is.na(EICStopRule)) {
    stop("EICStopRule must be one of: ", paste(choices, collapse = ", "))
  }
  EICStopRule <- tolower(trimws(EICStopRule))
  if (!(EICStopRule %in% choices)) {
    stop("EICStopRule must be one of: ", paste(choices, collapse = ", "))
  }
  EICStopRule
}

getEICStopAbsTol <- function(EICStopAbsTol) {
  EICStopAbsTolOK <- try(all(
    is.numeric(EICStopAbsTol),
    length(EICStopAbsTol) == 1,
    EICStopAbsTol >= 0,
    is.finite(EICStopAbsTol)
  ))
  if (any(inherits(EICStopAbsTolOK, "try-error"), !EICStopAbsTolOK, is.null(EICStopAbsTolOK))) {
    message("EICStopAbsTol must be a non-negative finite scalar, so has been set to 0\n")
    EICStopAbsTol <- 0
  }
  EICStopAbsTol
}

safeEICRatio <- function(abs_pneic, threshold) {
  out <- rep(Inf, length(abs_pneic))
  positive <- is.finite(threshold) & threshold > 0
  out[positive] <- abs_pneic[positive] / threshold[positive]
  out[!positive & abs_pneic == 0] <- 0
  out
}

makeOneStepStop <- function(SummEIC,
                            EICStopRule = "relative",
                            EICStopAbsTol = 0) {
  Trt <- Time <- Event <- PnEIC <- `seEIC/(sqrt(n)log(n))` <- NULL

  EICStopRule <- getEICStopRule(EICStopRule)
  EICStopAbsTol <- getEICStopAbsTol(EICStopAbsTol)
  out <- SummEIC[, .(
    Trt = Trt,
    Time = Time,
    Event = Event,
    PnEIC = PnEIC,
    RelativeCriteria = `seEIC/(sqrt(n)log(n))`
  )]
  out[, AbsPnEIC := abs(PnEIC)]
  out[, AbsoluteCriteria := EICStopAbsTol]

  if (identical(EICStopRule, "relative")) {
    out[, StopCriteria := RelativeCriteria]
  } else if (identical(EICStopRule, "absolute")) {
    out[, StopCriteria := AbsoluteCriteria]
  } else {
    out[, StopCriteria := pmax(RelativeCriteria, AbsoluteCriteria, na.rm = TRUE)]
  }

  out[, RelativeRatio := safeEICRatio(AbsPnEIC, RelativeCriteria)]
  out[, AbsoluteRatio := safeEICRatio(AbsPnEIC, AbsoluteCriteria)]
  out[, ratio := safeEICRatio(AbsPnEIC, StopCriteria)]
  out[, check := ratio <= 1]
  out[, StopRule := EICStopRule]
  out[, StopAbsTol := EICStopAbsTol]
  out[]
}
