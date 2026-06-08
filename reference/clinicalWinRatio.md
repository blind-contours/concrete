# Clinical (death-priority) win ratio for an illness-death outcome (experimental)

**Experimental.** Estimates the *clinical* (death-priority) win ratio,
win odds, and net benefit for a two-arm trial with a non-fatal
intercurrent event (e.g.\\ heart-failure hospitalization) and a terminal
event (death). Unlike the competing-risks (first-event) win ratio in
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md),
this estimand counts **death even when it follows the non-fatal event**
— the clinically intended hierarchy "compare on death first; break ties
on the non-fatal event." It is built on a Markov illness-death model
with three transition intensities (alive\\\to\\non-fatal,
alive\\\to\\death, post-non-fatal\\\to\\death), each estimated by a
Super Learner (the post-non-fatal death hazard on a left-truncated risk
set), and returns influence-function confidence intervals that are
doubly-robust, covariate- adjusted, and censoring-corrected (IPCW).

The estimator and its inference are validated against ground truth (a
brute-force pairwise win ratio on full simulated histories): see the
"Win ratio for trialists" article and `scripts/make-clinical-wr-*.R`. It
is marked experimental because it currently takes its own multistate
data frame (below) rather than the standard
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
pipeline, and assumes a single non-fatal event type and
conditionally-independent censoring (CAR).

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

  character: name of the non-fatal-event time column; `NA` (or `Inf`)
  for subjects who never had the non-fatal event.

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

## See also

[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
for the single-event and competing-risks win ratio.

## Examples

``` r
if (FALSE) { # \dontrun{
# data with: arm, t_hfh (NA if none), t_term (death or censoring), died (1/0), age, sex
clinicalWinRatio(trial, arm = "arm", illness.time = "t_hfh",
                 terminal.time = "t_term", terminal.status = "died",
                 covariates = c("age", "sex"), horizon = 1460)
} # }
```
