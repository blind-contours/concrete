# Restricted mean survival time and cause-specific life-years lost

`getRMST()` turns a fitted `"ConcreteEst"` object into restricted mean
survival time (RMST) and cause-specific life-years-lost (LYL) estimands,
which are collapsible, clinically interpretable summaries that
regulators increasingly prefer to a hazard ratio. Both are linear
functionals of the cumulative-incidence curves that `concrete` already
targets, so their influence functions are time-integrals of the
per-subject influence functions of the absolute risks. The integral is
taken over the **target times the model was fit on**, so request a
reasonably dense `TargetTime` grid in
[`formatArguments()`](https://blind-contours.github.io/concrete/reference/formatArguments.md)
for an accurate RMST.

- **RMST** (event-free): \\\int_0^\tau S(t)\\dt\\, the mean amount of
  follow-up time spent free of all events up to the horizon \\\tau\\.
  Only returned when every event type in the data was targeted.

- **Life-years lost** to cause \\j\\: \\\int_0^\tau F_j(t)\\dt\\, the
  mean time lost to cause \\j\\ by \\\tau\\.

## Usage

``` r
getRMST(
  ConcreteEst,
  Horizon = NULL,
  Intervention = seq_along(ConcreteEst),
  Contrasts = TRUE,
  Signif = 0.05,
  NIMargin = NULL,
  NIDirection = c("lower", "upper")
)
```

## Arguments

- ConcreteEst:

  a `"ConcreteEst"` object returned by
  [`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md).

- Horizon:

  numeric: the restriction horizon \\\tau\\. Defaults to the largest
  target time. Snapped to the nearest target time at or below it.

- Intervention:

  numeric (default `seq_along(ConcreteEst)`): which interventions to
  summarize. For contrasts the first two are treated as treatment and
  control.

- Contrasts:

  logical: also return RMST / LYL differences between the first two
  interventions.

- Signif:

  numeric (default 0.05): alpha for two-sided confidence intervals and
  two-sided Wald p-values.

- NIMargin:

  numeric (optional): a non-inferiority margin for the contrast
  estimands. When supplied, a one-sided non-inferiority assessment is
  added.

- NIDirection:

  one of `"lower"` or `"upper"`: which side of the margin is
  "non-inferior". Use `"lower"` when a larger value is better (e.g. an
  RMST difference) and `"upper"` when a smaller value is better.

## Value

a `data.table` of class `"ConcreteOut"` with point estimates,
influence-function standard errors, confidence intervals, and p-values.

## See also

[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
for absolute risks, risk differences, and risk ratios.
