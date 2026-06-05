# Directly target the restricted mean survival time

`targetRMST()` runs an additional one-step TMLE update that targets the
restricted mean survival time (RMST) and cause-specific life-years-lost
(LYL) estimands **directly**, rather than integrating pointwise-targeted
absolute risks (the approach used by
[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)).

The estimand \\\mathrm{LYL}\_j(\tau)=\int_0^\tau F_j(t)\\dt\\ is a
smooth time-average of the cumulative incidence, so its efficient
influence function is the time-integral of the pointwise risk influence
functions and its clever covariate is the time-integral of the pointwise
clever covariates: \$\$H\_{l,j,\tau}(s) = \frac{\pi^\*(A\mid
W)}{\pi(A\mid W) S_c(s^-)} \left\[\mathbf 1(l=j)(\tau-s) -
\frac{\int_s^\tau F_j(t)\\dt - (\tau-s)F_j(s)}{S(s)}\right\].\$\$
Targeting this single, well-conditioned functional avoids the dependence
on a dense `TargetTime` grid and tends to converge and cover better than
the pointwise approach in rare-event, competing-risk, and long-horizon
settings.

Starting from the hazards in `ConcreteEst` (already fit by
[`doConcrete()`](https://blind-contours.github.io/concrete/reference/doConcrete.md)),
`targetRMST()` fluctuates them along \\H\\ until the empirical mean of
the RMST/LYL influence function is small, then reports the
directly-targeted estimates with influence-function standard errors,
p-values, and (optionally) a non-inferiority assessment.

## Usage

``` r
targetRMST(
  ConcreteEst,
  Horizon = NULL,
  Intervention = seq_along(ConcreteEst),
  TargetEvent = NULL,
  MaxUpdateIter = 50L,
  OneStepEps = 0.1,
  Signif = 0.05,
  NIMargin = NULL,
  NIDirection = c("lower", "upper"),
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

  numeric: the restriction horizon \\\tau\\. Defaults to the largest
  target time; snapped to the nearest target time at or below it.

- Intervention:

  numeric (default `seq_along(ConcreteEst)`): interventions to
  summarize; the first two are treated as treatment and control.

- TargetEvent:

  numeric: event types to target. Defaults to the events the model was
  fit on.

- MaxUpdateIter:

  integer: maximum fluctuation steps.

- OneStepEps:

  numeric: initial step size for the fluctuation.

- Signif:

  numeric (default 0.05): alpha for confidence intervals and two-sided
  Wald p-values.

- NIMargin, NIDirection:

  optional non-inferiority margin and direction, passed to the contrast
  estimands (see
  [`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)).

- EICStopRule:

  one of `"hybrid"` (default), `"relative"`, or `"absolute"`: the
  stopping rule for the RMST/LYL estimating equation, evaluated on the
  rescaled fraction-of-horizon scale.

- EICStopAbsTol:

  numeric absolute tolerance for the `"absolute"` and `"hybrid"` rules
  on the fraction scale. Defaults to `0.02 / sqrt(n)`.

- Verbose:

  logical: print per-step convergence diagnostics.

## Value

a `data.table` of class `"ConcreteOut"` with the directly-targeted RMST
/ life-years-lost estimates. The per-arm convergence status is stored in
`attr(., "RMSTConverged")`.

## See also

[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
for the integrate-pointwise-risks version.
