# Priority-standardized net benefit and win ratio (charter-weighted)

Standard hierarchical win statistics aggregate the layer-specific
effects using *reach* weights — the probability that a treated–control
pair is tied on all higher-priority layers and is therefore decided at a
given layer. Those weights are set by the outcome distribution and the
tie/censoring conventions, not by clinical priority, so a
frequently-reached (often low-priority) layer can numerically dominate
the composite. The priority-standardized net benefit (PSNB) and win
ratio (PSWR) replace the implicit reach weights with a prespecified
clinical **charter** \\\alpha\\ (one weight per layer, \\\alpha_k \ge
0\\, \\\sum_k \alpha_k = 1\\), fixed before unblinding:
\$\$\mathrm{PSNB}(\alpha) = \sum_k \alpha_k \Delta_k, \qquad
\mathrm{PSWR}(\alpha) = \frac{\sum_k \alpha_k w_k}{\sum_k \alpha_k
\ell_k},\$\$ where \\\Delta_k = w_k - \ell_k\\ is the stage-conditional
net benefit and \\w_k\\, \\\ell_k\\ are the win / loss probabilities
among pairs that reach layer \\k\\.

This is the **covariate-adjusted, doubly-robust,
IPCW-censoring-corrected** estimator of PSNB/PSWR for a death-priority
time-to-event hierarchy: it reuses the same validated multistate engine
as
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md),
which already produces the per-tier win/loss components \\W^{(k)} = r_k
w_k\\ and \\L^{(k)} = r_k \ell_k\\. The standard win ratio sums those
tiers; PSNB divides the reach \\r_k = 1 -
\sum\_{m\<k}(W^{(m)}+L^{(m)})\\ back out and recombines with the
charter. Influence-function inference is propagated by the delta method
from the engine's per-tier influence functions.

## Usage

``` r
clinicalPSNB(
  data,
  arm,
  illness.time,
  terminal.time,
  terminal.status,
  covariates,
  charter,
  horizon = NULL,
  n.grid = 60L,
  n.folds = 5L,
  SL.library = c("SL.mean", "SL.glm"),
  Signif = 0.05,
  id = NULL,
  censoring.tv = NULL,
  pro = NULL
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

- charter:

  the priority charter: a numeric vector of length \\K\\ (the number of
  layers, layer 1 = death, followed by any non-fatal event tiers and
  then any bottom `pro` tiers in order), giving the weight on each layer
  (rescaled to sum to 1). The special value `"reach"` uses the realized
  reach weights, which reproduces the standard net benefit and win ratio
  (useful as a reference / sanity check). The charter must be
  **prespecified**: it is part of the estimand, not a tuning parameter.

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

- pro:

  optional continuous / ordinal patient-reported-outcome (PRO) tier(s)
  appended at the **bottom** of the hierarchy (below all hard-event
  tiers), the clinical norm for soft markers. A single spec (a named
  `list`) or a `list` of specs, each with: `marker` (column of the
  landmark value, `NA` if not measured), `landmark` (measurement time;
  default = horizon), `margin` (the win margin \\\delta\\; default 0),
  `direction` (`"higher.better"` (default) or `"lower.better"`), `type`
  (`"continuous"` (default) or `"ordinal"`), `n.grid` (cutpoint
  resolution for continuous markers; default 80), and optional `label`.
  A pair reaches a PRO tier iff tied on all higher tiers (both
  event-free and alive at the horizon); within reach the markers are
  compared with margin \\\delta\\. The marker distribution is
  **reach-weighted** standardized and landmark-missingness is
  IPCW-corrected; see Details and `clinicalPSNB()`.

## Value

a `data.table` of class `"ConcreteOut"` with `PSNB` and `PSWR` (each
with an influence-function CI; PSWR on the log scale), plus, for each
layer, the reach probability and the stage-conditional net benefit
\\\Delta_k\\ so the user can see which layer drives the composite.

## See also

[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
(the standard reach-weighted win ratio),
[`clinicalRMTIF()`](https://blind-contours.github.io/concrete/reference/clinicalRMTIF.md).
