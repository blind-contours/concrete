# Restricted mean time in favor of treatment (first-event competing risks)

The restricted mean time in favor of treatment (RMT-IF; Mao 2023) is the
average amount of time over \\\[0,\tau\]\\ that a randomly chosen
treated patient spends in a strictly more favorable state than a
randomly chosen control. Unlike the (unitless) win ratio, it is reported
in time units, which clinicians find easier to interpret. It is the
**time integral of the instantaneous net benefit**,
\$\$\mathrm{RMT\text{-}IF}(\tau) = \int_0^\tau \big\[w(t) -
l(t)\big\]\\dt,\$\$ where \\w(t)\\ (resp. \\l(t)\\) is the probability
that the treated patient is in a better (resp. worse) state than the
control at time \\t\\.

For a **single event** this reduces *exactly* to the restricted mean
survival time difference, \\\mathrm{RMT\text{-}IF} = \mathrm{RMST}\_1 -
\mathrm{RMST}\_0\\: with states alive, with-event, \\w(t)-l(t) =
S_1(t) - S_0(t)\\. For a **prioritized hierarchy** of competing events
the state at \\t\\ is "event-free" or "first event was type \\k\\ by
\\t\\", ranked by priority (highest priority = least favorable), so
\\w(t)\\ and \\l(t)\\ are bilinear in the two arms' cause-specific
cumulative incidences. RMT-IF is therefore a smooth functional of the
targeted marginal curves, with the same covariate-adjusted,
doubly-robust, censoring-corrected influence-function inference as the
other estimands.

This is the **first-event** version: a higher-priority event that
follows a lower-priority one is not credited (it treats the events as
competing risks), exactly as in
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md).
For the clinically intended death-priority hierarchy that credits death
after a non-fatal event, use the multistate version (see
[`clinicalWinRatio()`](https://blind-contours.github.io/concrete/reference/clinicalWinRatio.md)
and
[`?clinicalRMTIF`](https://blind-contours.github.io/concrete/reference/clinicalRMTIF.md)).

## Usage

``` r
getRMTIF(
  ConcreteEst,
  Horizon = NULL,
  Intervention = c(1, 2),
  TargetEvent = NULL,
  Signif = 0.05
)
```

## Arguments

- ConcreteEst:

  a `"ConcreteEst"` object from
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- Horizon:

  numeric: the restriction horizon \\\tau\\ (default: the largest target
  time). RMT-IF is reported in the time units of `Horizon`.

- Intervention:

  length-2 numeric: treatment and control indices; a positive RMT-IF
  favors `Intervention[1]`.

- TargetEvent:

  numeric: the event code, or an ordered vector of codes giving the
  priority hierarchy from highest to lowest (default: the first targeted
  event).

- Signif:

  numeric (default 0.05): alpha for confidence intervals and the
  two-sided Wald p-value.

## Value

a `data.table` of class `"ConcreteOut"` with `RMT-IF` (net), the time in
favor (\\\int w\\) and against (\\\int l\\), each with an
influence-function CI; the net carries a p-value against the null of no
difference. With a `Strata` fit the SEs are design-corrected.

## See also

[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
(the single-event special case),
[`getWinRatio()`](https://blind-contours.github.io/concrete/reference/getWinRatio.md),
[`getSimultaneousFamily()`](https://blind-contours.github.io/concrete/reference/getSimultaneousFamily.md).
