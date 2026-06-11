# Directly targeted win ratio, win odds, and net benefit

A targeted-maximum-likelihood update that fluctuates **both arms'**
cause-specific hazards to solve the efficient-influence-function
estimating equations of the win probability and loss probability
**directly**, rather than plugging pointwise-targeted risk curves into
the win functional (which is what
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
does).

For prioritized events \\k = 1, \ldots, K\\ (first = highest priority)
with treated-arm cumulative incidences \\G_k\\ and control-arm \\H_k\\,
the win probability over horizon \\\tau\\ is \$\$P(win) =
S_G(\tau)\sum_k H_k(\tau) + \sum\_{a \ge 2} G_a(\tau) \sum\_{k\<a}
H_k(\tau) + \sum_k \int_0^\tau \[G_k(\tau) - G_k(t)\] \\ dH_k(t),\$\$
and \\P(loss)\\ swaps the arms. The clever covariate for each functional
is the gradient-weighted combination of the pointwise clever covariates
over the **full event-time grid** (the chain rule analogue of
[`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md)'s
time-integrated clever covariate), so the win integral is evaluated on
every event time rather than the
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
target-time grid. Because the gradient coefficients depend on both arms'
current curves, they are recomputed at every fluctuation step and the
two arms are updated jointly until \\P_n D\_{win} \approx P_n D\_{loss}
\approx 0\\.

Compared to the plug-in this removes the target-grid discretization of
the win integral and solves the functional's own estimating equation,
which matters most for sparse target grids and for the win odds / net
benefit. Derivation: `notes/target-win-ratio.md` in the source
repository.

## Usage

``` r
targetWinRatio(
  ConcreteEst,
  Horizon = NULL,
  Intervention = c(1, 2),
  TargetEvent = NULL,
  MaxUpdateIter = 50L,
  OneStepEps = 0.1,
  Signif = 0.05,
  EICStopRule = c("hybrid", "relative", "absolute"),
  EICStopAbsTol = NULL,
  Verbose = FALSE
)
```

## Arguments

- ConcreteEst:

  a `"ConcreteEst"` object from
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- Horizon:

  numeric: the comparison horizon \\\tau\\. Defaults to the largest
  target time; snapped down to the nearest evaluation time.

- Intervention:

  length-2 numeric: `ConcreteEst` list indices of the treated and
  control interventions, in that order.

- TargetEvent:

  numeric: event types in priority order (first = highest, e.g. death).
  Defaults to the fitted target events in their given order.

- MaxUpdateIter:

  integer: maximum fluctuation steps.

- OneStepEps:

  numeric: initial step size for the fluctuation.

- Signif:

  numeric (default 0.05): alpha for CIs and two-sided p-values.

- EICStopRule:

  one of `"hybrid"` (default), `"relative"`, `"absolute"`: stopping rule
  for the (win, loss) estimating equations.

- EICStopAbsTol:

  numeric absolute tolerance for the `"absolute"` and `"hybrid"` rules.
  Defaults to `0.02 / sqrt(n)` (the functionals are probabilities, so
  the risk-scale default applies).

- Verbose:

  logical: print per-step convergence diagnostics.

## Value

a `data.table` of class `"ConcreteOut"` with the win ratio, win odds,
net benefit, and win/loss/tie probabilities, each with
influence-function CIs and (for the comparative statistics) p-values.
Convergence status is in `attr(., "WRConverged")` and the step count in
`attr(., "WRSteps")`. If the fit was built with `Strata` (see
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)),
the standard errors are corrected for the stratified /
covariate-adaptive randomization design.

## See also

[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md)
for the plug-in version;
[`targetRMST()`](https://blind-contours.github.io/concrete/reference/targetRMST.md)
for the same construction applied to the restricted mean survival time.
