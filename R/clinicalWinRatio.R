#' Hierarchical (death-priority) clinical win ratio (experimental)
#'
#' @description
#' \strong{The recommended win ratio for most trials.} Estimates the
#' \emph{clinical}, death-priority win ratio, win odds, and net benefit for a
#' two-arm trial with an ordered hierarchy of a terminal event (death) and one or
#' more non-fatal events (e.g.\ heart-failure hospitalization, stroke, valve
#' intervention). Unlike the first-event / competing-risks win ratio in
#' [getWinRatio()], this estimand counts **a higher-priority event even when it
#' follows a lower-priority one** --- death after a non-fatal event, or a stroke
#' after a hospitalization. That is the clinically intended hierarchy ("compare on
#' the most serious event first; break ties on the next"), and it is the win ratio
#' the first-event version cannot produce.
#'
#' It is built on a Markov multistate model whose states are the subsets of
#' non-fatal events a subject has experienced; every transition intensity (each
#' non-fatal event out of each reachable state, and death out of every state) is
#' estimated by a Super Learner, with doubly-robust, covariate-adjusted,
#' censoring-corrected (IPCW) influence-function inference and optional
#' cross-fitting. The estimator and its inference are validated against ground
#' truth (a brute-force pairwise win ratio on full simulated histories) for
#' hierarchies up to four time-to-event tiers: see the "Win ratios for trialists"
#' article and `scripts/genwr-*.R`.
#'
#' It is marked experimental because it currently takes its own per-subject event
#' columns (below) rather than the standard [formatArguments()] pipeline, and
#' assumes non-recurrent events, conditionally-independent censoring (CAR), and a
#' Markov model. Recurrent-event tiers (repeated hospitalizations) and
#' continuous/ordinal tiers (e.g.\ KCCQ) are not yet supported.
#'
#' @param data a `data.frame`/`data.table`, one row per subject.
#' @param arm character: name of the binary treatment column (1 = active arm).
#' @param illness.time character vector: the non-fatal-event time columns, **ordered
#'   highest priority first** (e.g.\ `c("t_stroke", "t_hosp")` for stroke > hosp).
#'   Each entry is the time of that subject's first such event, `NA` (or `Inf`) if
#'   it never occurred. A single column reproduces the two-tier illness-death case.
#'   Death is always the top-priority tier.
#' @param terminal.time character: name of the terminal time column (time of death
#'   or of censoring, whichever came first).
#' @param terminal.status character: name of the terminal status column
#'   (1 = death, 0 = censored).
#' @param covariates character vector: baseline covariate column names.
#' @param horizon numeric: the restriction horizon \eqn{\tau} (default: the
#'   largest terminal time).
#' @param n.grid integer (default 60): number of time intervals for the discrete
#'   hazard / path-probability quadrature.
#' @param n.folds integer (default 5): number of cross-fitting folds. The
#'   transition and censoring hazards are fit out-of-fold, which gives honest
#'   inference when the `SL.library` contains flexible learners that could over-fit
#'   in sample; with simple parametric learners it makes little difference. Set to
#'   1 to disable cross-fitting (faster). \strong{Note:} cross-fitting does
#'   \emph{not} fix the mild small-sample anti-conservatism described below --- that
#'   is a finite-sample property of the win ratio itself.
#' @param SL.library character vector: SuperLearner library for the transition and
#'   censoring hazards (default `c("SL.mean", "SL.glm")`).
#' @param Signif numeric (default 0.05): alpha for confidence intervals.
#' @param id character (optional): name of a subject id column, required only when
#'   `censoring.tv` is supplied (to link the longitudinal measurements to subjects).
#' @param censoring.tv optional `data.frame` of \strong{time-varying covariates for
#'   the censoring model} (e.g.\ post-randomization echo / KCCQ / 6-minute-walk
#'   measured at follow-up visits), in long form with the id column (named as `id`),
#'   a `time` column, and one or more value columns. When supplied, the censoring
#'   hazard is conditioned on the last-observation-carried-forward value and
#'   change-from-baseline of each, which corrects inverse-probability-of-censoring
#'   bias when dropout is driven by these measurements. They enter \strong{only} the
#'   censoring model (never the outcome hazards), so the marginal/ITT estimand is
#'   preserved (they are post-treatment mediators). No effect on the result when
#'   omitted.
#'
#' @return a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
#'   net benefit, and the win/loss/tie probabilities, each with an
#'   influence-function standard error, confidence interval, and (for the
#'   comparative statistics) a p-value against the null of no difference.
#'
#' @section Small-sample behavior:
#' Like the win ratio in general (including the unadjusted Pocock win ratio), the
#' point estimate is a \emph{ratio} and is therefore mildly biased and
#' anti-conservative in small samples. In a null simulation (true win ratio 1,
#' both arms identical) the estimator is biased downward by \eqn{\approx}1\% at
#' \eqn{\sim}400/arm, with Wald coverage \eqn{\approx}0.93--0.94 and type-I error
#' \eqn{\approx}0.06--0.07; this is a finite-sample property of the win-ratio
#' functional, not of the nuisance estimation (cross-fitting does not change it).
#' The bias and under-coverage shrink at the usual \eqn{O(1/n)} rate, and inference
#' is nominal (coverage 0.95--0.97) by \eqn{\sim}800/arm. For small trials,
#' interpret the interval as mildly optimistic, or use a resampling interval.
#'
#' @seealso [getWinRatio()] for the first-event / competing-risks win ratio (the
#'   special case where events are mutually exclusive and a higher-priority event
#'   can never follow a lower-priority one).
#' @export clinicalWinRatio
#' @examples
#' \dontrun{
#' # Two-tier (death > hospitalization):
#' clinicalWinRatio(trial, arm = "arm", illness.time = "t_hosp",
#'                  terminal.time = "t_term", terminal.status = "died",
#'                  covariates = c("age", "sex"), horizon = 1460)
#' # Three-tier hierarchy (death > stroke > hospitalization):
#' clinicalWinRatio(trial, arm = "arm", illness.time = c("t_stroke", "t_hosp"),
#'                  terminal.time = "t_term", terminal.status = "died",
#'                  covariates = c("age", "sex"), horizon = 1460)
#' }
clinicalWinRatio <- function(data, arm, illness.time, terminal.time, terminal.status,
                             covariates, horizon = NULL, n.grid = 60L, n.folds = 5L,
                             SL.library = c("SL.mean", "SL.glm"), Signif = 0.05,
                             id = NULL, censoring.tv = NULL) {
  data <- as.data.frame(data)
  illness.time <- as.character(illness.time)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  if (length(unique(A)) != 2L) stop("arm must contain both 0 and 1.")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  if (!all(delta %in% c(0, 1))) stop("terminal.status must be coded 0/1 (1 = death).")
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  K <- 1L + length(illness.time)
  grid <- seq(0, horizon, length.out = as.integer(n.grid) + 1L)

  ## --- parse to per-subject observed quantities: tD, t1..t{K-1}, C + covariates ---
  D <- data[, covariates, drop = FALSE]
  D$tD <- ifelse(delta == 1, term, Inf)
  for (ei in seq_along(illness.time)) {
    ti <- data[[illness.time[ei]]]; ti[is.na(ti)] <- Inf; D[[paste0("t", ei)]] <- ti
  }
  D$C <- ifelse(delta == 0, term, Inf)

  ## --- optional time-varying censoring covariates (LOCF value + change) ---
  tvMats <- NULL
  if (!is.null(censoring.tv)) {
    if (is.null(id)) stop("`id` (subject id column) is required when `censoring.tv` is supplied.")
    if (!id %in% names(data)) stop("id column '", id, "' not found in data.")
    censoring.tv <- as.data.frame(censoring.tv)
    if (!"time" %in% names(censoring.tv)) stop("`censoring.tv` must have a 'time' column.")
    if (!id %in% names(censoring.tv)) stop("`censoring.tv` must have the id column '", id, "'.")
    tvMats <- .tvLOCF(data[[id]], censoring.tv, id, "time", grid[-length(grid)])
  }

  eng <- .msEngine(K, grid)
  buildArm <- function(av) {
    sel <- A == av; Da <- D[sel, , drop = FALSE]
    tvA <- if (is.null(tvMats)) NULL else lapply(tvMats, function(m) m[sel, , drop = FALSE])
    nu <- .msNuisances(eng, Da, covariates, SL.library, n.folds, tvA)
    eng$armSetup(Da, nu$rmat, nu$Ginv)
  }
  out <- .msWinRatioOut(eng, buildArm(1), buildArm(0), Signif)
  attr(out, "Horizon") <- horizon
  attr(out, "Estimand") <- "Clinical Win Ratio"
  attr(out, "Tiers") <- K
  attr(out, "Experimental") <- TRUE
  class(out) <- union("ConcreteOut", class(out))
  out[]
}
