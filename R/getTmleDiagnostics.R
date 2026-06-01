#' Extract TMLE convergence diagnostics
#'
#' @description
#' `getTmleDiagnostics()` returns the empirical efficient influence curve (EIC)
#' convergence diagnostics stored on a fitted `"ConcreteEst"` object. Use
#' `type = "components"` to inspect the final component-wise stopping criteria,
#' `type = "trace"` to inspect the update history, or `type = "norm"` to inspect
#' the norm of the empirical mean EIC across update steps.
#'
#' @param ConcreteEst `"ConcreteEst"` object returned by [doConcrete()].
#' @param type character; one of `"components"`, `"trace"`, or `"norm"`.
#'
#' @return A `data.table` containing the requested diagnostics.
#'
#' @examples
#' library(data.table)
#' library(concrete)
#'
#' data <- as.data.table(survival::pbc)
#' data <- data[1:200, .SD, .SDcols = c("id", "time", "status", "trt", "age", "sex")]
#' data[, trt := sample(0:1, nrow(data), TRUE)]
#'
#' concrete.args <- formatArguments(DataTable = data,
#'                                  EventTime = "time",
#'                                  EventType = "status",
#'                                  Treatment = "trt",
#'                                  ID = "id",
#'                                  TargetTime = 2500,
#'                                  TargetEvent = c(1, 2),
#'                                  Intervention = makeITT(),
#'                                  CVArg = list(V = 2),
#'                                  MaxUpdateIter = 2,
#'                                  Verbose = FALSE)
#'
#' \donttest{
#' concrete.est <- doConcrete(concrete.args)
#' getTmleDiagnostics(concrete.est, type = "components")
#' }
#'
#' @export
getTmleDiagnostics <- function(ConcreteEst,
                               type = c("components", "trace", "norm")) {
  Trt <- Time <- Event <- check <- NULL

  if (!inherits(ConcreteEst, "ConcreteEst")) {
    stop("ConcreteEst must be a 'ConcreteEst' object returned by doConcrete().")
  }
  type <- match.arg(type)

  if (identical(type, "trace")) {
    trace <- attr(ConcreteEst, "TmleUpdateTrace")
    if (is.null(trace)) {
      return(data.table::data.table())
    }
    return(data.table::copy(trace))
  }

  if (identical(type, "norm")) {
    norm <- attr(ConcreteEst, "NormPnEICs")
    if (is.null(norm)) {
      return(data.table::data.table())
    }
    return(data.table::data.table(
      Step = seq_along(norm) - 1L,
      NormPnEIC = as.numeric(norm)
    ))
  }

  TargetTime <- attr(ConcreteEst, "TargetTime")
  TargetEvent <- attr(ConcreteEst, "TargetEvent")
  EICStopRule <- attr(ConcreteEst, "EICStopRule")
  EICStopAbsTol <- attr(ConcreteEst, "EICStopAbsTol")
  if (is.null(EICStopRule)) EICStopRule <- "relative"
  if (is.null(EICStopAbsTol)) EICStopAbsTol <- 0

  SummEIC <- data.table::rbindlist(
    lapply(seq_along(ConcreteEst), function(a) {
      out <- data.table::copy(ConcreteEst[[a]][["SummEIC"]])
      out[, Trt := names(ConcreteEst)[a]]
      out
    }),
    fill = TRUE
  )
  data.table::setcolorder(SummEIC, c("Trt", setdiff(names(SummEIC), "Trt")))

  diagnostics <- targetOneStepStop(
    SummEIC = SummEIC,
    TargetTime = TargetTime,
    TargetEvent = TargetEvent,
    EICStopRule = EICStopRule,
    EICStopAbsTol = EICStopAbsTol
  )
  data.table::setnames(diagnostics, "Trt", "Intervention")
  diagnostics[, Converged := isTRUE(attr(ConcreteEst, "TmleConverged")$converged)]
  diagnostics[, ConvergenceStep := attr(ConcreteEst, "TmleConverged")$step]
  diagnostics[]
}
