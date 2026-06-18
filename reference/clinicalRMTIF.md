# Clinical (death-priority) restricted mean time in favor of treatment

The death-priority (multistate) version of the restricted mean time in
favor of treatment (RMT-IF; Mao 2023): the average time over
\\\[0,\tau\]\\ that a random treated patient spends in a clinically more
favorable state than a random control, under the hierarchy that compares
on the most serious event first *whenever* it occurs. Unlike the
first-event version
[`getRMTIF()`](https://blind-contours.github.io/concrete/reference/getRMTIF.md),
this **credits a higher-priority event that follows a lower-priority
one** (death after a hospitalization), and unlike the win ratio it is in
time units.

It is built on the same Markov multistate engine as
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md):
each subject's state at time \\t\\ is the subset of non-fatal events
experienced (plus alive/dead), every transition intensity is estimated
by a Super Learner with IPCW censoring correction and optional
cross-fitting, and the favorability of a state is its most-severe event
so far (event-free best, death worst). The point estimate is
\$\$\mathrm{RMT\text{-}IF} = \int_0^\tau \big\[w(t) - l(t)\big\]\\dt,
\quad w(t)=\sum\_{r\<r'}\pi^1_r(t)\\\pi^0\_{r'}(t),\$\$ with
\\\pi^a_r(t)\\ the probability that arm \\a\\ occupies favorability
level \\r\\ at \\t\\.

## Usage

``` r
clinicalRMTIF(
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
  Signif = 0.05,
  id = NULL,
  censoring.tv = NULL,
  nBoot = 0L
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

- id:

  character (optional): name of a subject id column, required only when
  `censoring.tv` is supplied (to link the longitudinal measurements to
  subjects).

- censoring.tv:

  optional `data.frame` of **time-varying covariates for the censoring
  model** (e.g.\\ post-randomization echo / KCCQ / 6-minute-walk
  measured at follow-up visits), in long form with the id column (named
  as `id`), a `time` column, and one or more value columns. When
  supplied, the censoring hazard is conditioned on the
  last-observation-carried-forward value and change-from-baseline of
  each, which corrects inverse-probability-of-censoring bias when
  dropout is driven by these measurements. They enter **only** the
  censoring model (never the outcome hazards), so the marginal/ITT
  estimand is preserved (they are post-treatment mediators). No effect
  on the result when omitted.

- nBoot:

  integer (default 0): inference method. With `nBoot = 0` (default) the
  **analytic adjoint-value efficient influence function** is used for
  the net RMT-IF standard error (fast). If positive, a stratified
  nonparametric bootstrap with `nBoot` resamples is used instead, which
  also gives SEs for the time-in-favor / time-against decomposition
  (correct but heavy — each resample refits the transition hazards).

## Value

a `data.table` of class `"ConcreteOut"` with the net `RMT-IF` (with an
influence-function SE, CI and p-value), the time in favor, and the time
against.

## Details

**Inference.** The net RMT-IF is the time integral of a bilinear
functional of the two arms' multistate level-occupancy curves. Its
efficient influence function is obtained by the chain rule: the per-arm
occupancy influence functions (a reward-accumulation adjoint over the
multistate process, with the death state absorbing) weighted by the
gradient coefficients \\c^1_r = P(\text{control worse than } r) -
P(\text{control better than } r)\\ and the mirror for the control arm.
This reuses the same adjoint-value machinery, IPCW correction, and
cross-fitting as
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md).
The point estimate is exact and validated against a brute-force pairwise
ground truth; the analytic SE is validated against the bootstrap and
against empirical coverage. SEs for the favor/against decomposition
require `nBoot`.

## See also

[`getRMTIF()`](https://blind-contours.github.io/concrete/reference/getRMTIF.md)
(first-event version with closed-form inference),
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md).
