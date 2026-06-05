# Relative efficiency of a covariate-adjusted vs unadjusted analysis

In a randomized trial, covariate adjustment does not change the target
estimand but it can sharpen it. `getRelativeEfficiency()` quantifies
that precision gain by comparing the influence-function standard errors
of an adjusted analysis to those of an unadjusted analysis of the
**same** estimand, which is the quantity the FDA's 2023
covariate-adjustment guidance is about.

Supply two `"ConcreteOut"` tables (from
[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md)
or
[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)):
one from a covariate-adjusted fit and one from an unadjusted fit. The
unadjusted fit is the same workflow with treatment-only nuisance models,
e.g. a marginal propensity (`"SL.mean"`) and hazard formulas of the form
`Surv(time, event == j) ~ arm`.

For each matched estimand the function reports the relative efficiency
\\\mathrm{RE} = \mathrm{Var}\_{\text{unadj}} /
\mathrm{Var}\_{\text{adj}}\\ (values above 1 favor adjustment), the
implied percentage variance reduction, and the effective sample-size
multiplier: an adjusted analysis on \\n\\ subjects has the precision of
an unadjusted analysis on \\\mathrm{RE}\\n\\.

## Usage

``` r
getRelativeEfficiency(Adjusted, Unadjusted)
```

## Arguments

- Adjusted:

  a `"ConcreteOut"` table from a covariate-adjusted fit.

- Unadjusted:

  a `"ConcreteOut"` table from an unadjusted (treatment-only) fit of the
  same estimands, interventions, events, and times.

## Value

a `data.table` keyed by `Intervention`, `Estimand`, `Event`, and `Time`
with columns `seAdjusted`, `seUnadjusted`, `RelEfficiency`,
`VarReductionPct`, and `EffSampleSizeMult`.

## See also

[`getOutput()`](https://blind-contours.github.io/concrete/reference/getOutput.md),
[`getRMST()`](https://blind-contours.github.io/concrete/reference/getRMST.md)
