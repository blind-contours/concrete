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
#' Markov model. Recurrent-event tiers (repeated hospitalizations) are not yet
#' supported. \strong{Continuous / ordinal patient-reported-outcome (PRO) tiers}
#' (e.g.\ KCCQ, NYHA, 6-minute walk) measured at a landmark \emph{are} supported as
#' bottom tiers via the `pro` argument --- see Details.
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
#' @param crossover optional character: name of a column giving each subject's
#'   \strong{treatment-switch time} (e.g.\ days from randomization to crossover),
#'   `NA`/`Inf` for those who never switched. When supplied, switchers are
#'   re-censored at their switch time and a \strong{separate covariate-adjusted
#'   crossover hazard} is fit and combined with the censoring hazard, so the IPCW
#'   becomes \eqn{1/(S_{\mathrm{drop}} S_{\mathrm{cross}})} --- the doubly-robust
#'   \strong{hypothetical no-switching} win ratio (what the contrast would be had
#'   no one crossed over), rather than the ITT treatment-policy estimand. Validated
#'   to recover the no-switching truth under informative switching
#'   (`scripts/dev-crossover-winratio.R`). Requires conditional independence of the
#'   switch and the outcome given the covariates; with heavy late crossover at
#'   small n the estimate can be high-variance.
#' @param min.cens.surv numeric (default 0.05): lower bound (truncation floor) on
#'   the combined censoring (and crossover) survival used in the IPCW, i.e.\ the
#'   inverse-probability weight is capped at `1/min.cens.surv`. Stabilizes the
#'   weights when follow-up / non-switching becomes rare; raise it for more
#'   stability (more bias) or lower it for less truncation. Matters most with heavy
#'   crossover, where the no-switching weights can otherwise blow up.
#' @param pro optional continuous / ordinal patient-reported-outcome (PRO) tier(s)
#'   appended at the \strong{bottom} of the hierarchy (below all hard-event tiers),
#'   the clinical norm for soft markers. A single spec (a named `list`) or a `list`
#'   of specs in priority order, each with: `marker` (column of the final-visit
#'   value, `NA` if not measured), `landmark` (measurement time; must equal the
#'   horizon --- a final-visit design), `margin` (the win margin \eqn{\delta};
#'   default 0), `direction` (`"higher.better"` (default) or `"lower.better"`),
#'   `type` (`"continuous"` (default) or `"ordinal"`), and optional `label`. A pair
#'   reaches the PRO block iff tied on all higher tiers (both event-free, alive and
#'   in follow-up at the horizon); within reach the PRO tiers are compared
#'   \strong{sequentially}, each with its margin \eqn{\delta}. This is the TRISCEND
#'   II construction (death > RV-assist/transplant > tricuspid reintervention > HF
#'   hospitalization > KCCQ \eqn{\ge} 10 > NYHA \eqn{\ge} 1 class > 6-min walk
#'   \eqn{\ge} 30 m). See Details and [clinicalPSNB()].
#'
#' @details
#' \strong{PRO tiers (experimental).} Markers measured at the final visit are
#' compared among pairs that reach the PRO block (tied on every hard-event tier,
#' i.e.\ both event-free and alive at the horizon), \emph{sequentially}: a pair is
#' decided at the first PRO tier whose two markers differ by more than that tier's
#' margin \eqn{\delta}, else it descends to the next PRO tier. The block is
#' estimated by a \strong{reach-weighted, IPCW-corrected two-sample generalized
#' pairwise comparison} on the joint marker vectors --- pairs restricted to
#' reachers and inverse-probability weighted for both pre-horizon censoring
#' (\eqn{1/G(\tau\mid W)}) and missing landmark visits (a per-arm attendance model)
#' --- so the markers may be arbitrarily correlated and the sequential
#' tie-passing among PRO tiers is exact. Inference is the two-sample U-statistic
#' (Hajek) influence function; hard-event tiers above keep the engine's
#' covariate-adjusted, doubly-robust influence-function inference. \strong{Scope}:
#' the PRO landmark must equal the horizon (final-visit design); a PRO ranked above
#' a hard event is not supported.
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
                             id = NULL, censoring.tv = NULL, crossover = NULL, pro = NULL,
                             min.cens.surv = 0.05) {
  data <- as.data.frame(data)
  illness.time <- as.character(illness.time)
  for (col in c(arm, illness.time, terminal.time, terminal.status, covariates))
    if (!col %in% names(data)) stop("column '", col, "' not found in data.")
  if (!is.null(crossover) && !crossover %in% names(data))
    stop("crossover column '", crossover, "' not found in data.")
  A <- data[[arm]]
  if (!all(A %in% c(0, 1))) stop("arm must be coded 0/1 (1 = active arm).")
  if (length(unique(A)) != 2L) stop("arm must contain both 0 and 1.")
  term <- data[[terminal.time]]; delta <- data[[terminal.status]]
  if (!all(delta %in% c(0, 1))) stop("terminal.status must be coded 0/1 (1 = death).")
  if (is.null(horizon)) horizon <- max(term[is.finite(term)])
  K <- 1L + length(illness.time)
  pros <- .proNormalize(pro, data, horizon)
  proCols <- unique(vapply(pros, function(s) s$marker, character(1)))
  grid <- seq(0, horizon, length.out = as.integer(n.grid) + 1L)

  ## --- parse to per-subject observed quantities: tD, t1..t{K-1}, C + covariates ---
  D <- data[, covariates, drop = FALSE]
  D$tD <- ifelse(delta == 1, term, Inf)
  for (ei in seq_along(illness.time)) {
    ti <- data[[illness.time[ei]]]; ti[is.na(ti)] <- Inf; D[[paste0("t", ei)]] <- ti
  }
  D$C <- ifelse(delta == 0, term, Inf)
  ## crossover (treatment-switching): re-censor at the switch time and carry it so
  ## a separate crossover hazard reweights the IPCW (hypothetical no-switching).
  sw <- if (is.null(crossover)) rep(Inf, nrow(data)) else { x <- as.numeric(data[[crossover]]); x[is.na(x)] <- Inf; x }
  D$switch <- sw; D$C <- pmin(D$C, sw)
  for (mc in proCols) D[[mc]] <- data[[mc]]                 # carry PRO markers through

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
    nu <- .msNuisances(eng, Da, covariates, SL.library, n.folds, tvA, xover = Da$switch, minG = min.cens.surv)
    list(arm = eng$armSetup(Da, nu$rmat, nu$Ginv), D = Da)
  }
  bT <- buildArm(1); bC <- buildArm(0)
  proC <- if (is.null(pros)) NULL else
    .proComponents(eng, pros, bT$D, bC$D, bT$arm, bC$arm, covariates, SL.library, n.folds)
  out <- .msWinRatioOut(eng, bT$arm, bC$arm, Signif, pro = proC)
  attr(out, "Horizon") <- horizon
  attr(out, "Estimand") <- "Clinical Win Ratio"
  attr(out, "Tiers") <- K + length(pros)
  attr(out, "Experimental") <- TRUE
  class(out) <- union("ConcreteOut", class(out))
  out[]
}
