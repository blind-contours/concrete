# Hierarchical (death-priority) clinical win ratio (experimental)

**The recommended win ratio for most trials.** Estimates the *clinical*,
death-priority win ratio, win odds, and net benefit for a two-arm trial
with an ordered hierarchy of a terminal event (death) and one or more
non-fatal events (e.g.\\ heart-failure hospitalization, stroke, valve
intervention). Unlike the first-event / competing-risks win ratio in
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md),
this estimand counts **a higher-priority event even when it follows a
lower-priority one** — death after a non-fatal event, or a stroke after
a hospitalization. That is the clinically intended hierarchy ("compare
on the most serious event first; break ties on the next"), and it is the
win ratio the first-event version cannot produce.

It is built on a Markov multistate model whose states are the subsets of
non-fatal events a subject has experienced; every transition intensity
(each non-fatal event out of each reachable state, and death out of
every state) is estimated by a Super Learner, with doubly-robust,
covariate-adjusted, censoring-corrected (IPCW) influence-function
inference and optional cross-fitting. The estimator and its inference
are validated against ground truth (a brute-force pairwise win ratio on
full simulated histories) for hierarchies up to four time-to-event
tiers: see the "Win ratios for trialists" article and
`scripts/genwr-*.R`.

It is marked experimental because it currently takes its own per-subject
event columns (below) rather than the standard
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
pipeline, and assumes non-recurrent events, conditionally-independent
censoring (CAR), and a Markov model. Recurrent-event tiers (repeated
hospitalizations) and continuous/ordinal tiers (e.g.\\ KCCQ) are not yet
supported.

## Usage

``` r
clinicalWinRatio(
  data,
  arm,
  illness.time,
  terminal.time,
  terminal.status,
  covariates,
  horizon = NULL,
  n.grid = 60L,
  n.folds = 5L,
  SL.library = c("SL.mean", "SL.glm"),
  Signif = 0.05
)
```

## Arguments

- data:

  a `data.frame`/`data.table`, one row per subject.

- arm:

  character: name of the binary treatment column (1 = active arm).

- illness.time:

  character vector: the non-fatal-event time columns, **ordered highest
  priority first** (e.g.\\ `c("t_stroke", "t_hosp")` for stroke \>
  hosp). Each entry is the time of that subject's first such event, `NA`
  (or `Inf`) if it never occurred. A single column reproduces the
  two-tier illness-death case. Death is always the top-priority tier.

- terminal.time:

  character: name of the terminal time column (time of death or of
  censoring, whichever came first).

- terminal.status:

  character: name of the terminal status column (1 = death, 0 =
  censored).

- covariates:

  character vector: baseline covariate column names.

- horizon:

  numeric: the restriction horizon \\\tau\\ (default: the largest
  terminal time).

- n.grid:

  integer (default 60): number of time intervals for the discrete hazard
  / path-probability quadrature.

- n.folds:

  integer (default 5): number of cross-fitting folds. The transition and
  censoring hazards are fit out-of-fold, which gives honest inference
  when the `SL.library` contains flexible learners that could over-fit
  in sample; with simple parametric learners it makes little difference.
  Set to 1 to disable cross-fitting (faster). **Note:** cross-fitting
  does *not* fix the mild small-sample anti-conservatism described below
  — that is a finite-sample property of the win ratio itself.

- SL.library:

  character vector: SuperLearner library for the transition and
  censoring hazards (default `c("SL.mean", "SL.glm")`).

- Signif:

  numeric (default 0.05): alpha for confidence intervals.

## Value

a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
net benefit, and the win/loss/tie probabilities, each with an
influence-function standard error, confidence interval, and (for the
comparative statistics) a p-value against the null of no difference.

## Small-sample behavior

Like the win ratio in general (including the unadjusted Pocock win
ratio), the point estimate is a *ratio* and is therefore mildly biased
and anti-conservative in small samples. In a null simulation (true win
ratio 1, both arms identical) the estimator is biased downward by
\\\approx\\1\\ \\\sim\\400/arm, with Wald coverage \\\approx\\0.93–0.94
and type-I error \\\approx\\0.06–0.07; this is a finite-sample property
of the win-ratio functional, not of the nuisance estimation
(cross-fitting does not change it). The bias and under-coverage shrink
at the usual \\O(1/n)\\ rate, and inference is nominal (coverage
0.95–0.97) by \\\sim\\800/arm. For small trials, interpret the interval
as mildly optimistic, or use a resampling interval.

## See also

[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
for the first-event / competing-risks win ratio (the special case where
events are mutually exclusive and a higher-priority event can never
follow a lower-priority one).

## Examples

``` r
if (FALSE) { # \dontrun{
# Two-tier (death > hospitalization):
clinicalWinRatio(trial, arm = "arm", illness.time = "t_hosp",
                 terminal.time = "t_term", terminal.status = "died",
                 covariates = c("age", "sex"), horizon = 1460)
# Three-tier hierarchy (death > stroke > hospitalization):
clinicalWinRatio(trial, arm = "arm", illness.time = c("t_stroke", "t_hosp"),
                 terminal.time = "t_term", terminal.status = "died",
                 covariates = c("age", "sex"), horizon = 1460)
} # }
```
